#' Calculate total fertility rates
#'
#' Sum age-specific fertility rates times their actual interval widths.
#' If parity-bound columns are present, delegate to [calculate_patfr()] instead
#' and return its `patfr` column rather than `tfr`, with the same other strata.
#' Open-ended age bands are omitted, warning if an omitted female rate is known
#' and nonzero that TFR underestimates fertility over the full age range. Zero or
#' missing omitted rates do not trigger that warning. The remaining finite, contiguous age range
#' defines the TFR reported. Any missing ASFR in those bands makes the TFR missing.
#' When `union_type`, `previous_union` (or a legacy `union_status` column) is
#' present, calculate a separate sum for each stratum
#' within each period. These conditional sums describe fertility under continuous
#' exposure to that status; they are not components to add into an overall TFR.
#' Any other column that is not an age band or a measure, such as a group
#' dimension (`race`, `region`), adds separate strata, with its type and labels
#' preserved; strata are not pooled. Age bands are checked independently within
#' each stratum.
#' @param rates A table returned by [calculate_fertility_rates()], a fertility
#'   process object, or a data frame with `age_start`, `age_end`, and `rate_fertility`
#'   or `rate` in births per woman-year. Bare `rate` tables are interpreted as
#'   fertility rates; mortality processes are not valid inputs.
#' @return With parity columns, the result of [calculate_patfr()] with `patfr`.
#'   Otherwise, a tibble with period columns, `union_type`, `previous_union` or
#'   legacy `union_status` when supplied,
#'   other stratum columns (such as group dimensions) when supplied,
#'   `age_start`, `age_end`, and `tfr`.
#'   If `sex` is present, only female rows are used.
#' @export
calculate_tfr <- function(rates) {
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "fertility") ||
        !identical(rates$rate_unit, "person_year"))
      stop("TFR requires a fertility process with person-year rates.", call. = FALSE)
    rates <- rates$rates
  }
  if (any(c("parity_start", "parity_end") %in% names(rates)))
    return(calculate_patfr(rates))
  check_columns(rates, c("age_start", "age_end"), "rates")
  if (!"rate_fertility" %in% names(rates)) {
    if (!"rate" %in% names(rates)) stop("Supply `rate_fertility` or fertility `rate`.", call. = FALSE)
    rates$rate_fertility <- rates$rate
  }
  check_groups(rates, "rates")
  unknown <- attr(rates, "unassigned_births")
  if ("sex" %in% names(rates)) rates <- rates[rates$sex == "f", , drop = FALSE]
  check_numeric_column(rates$age_start, "rates$age_start")
  if (!is.numeric(rates$age_end) || anyNA(rates$age_end))
    stop("`age_end` must contain numeric age boundaries without missing values.", call. = FALSE)
  check_numeric_column(rates$rate_fertility, "rates$rate_fertility", missing = TRUE)
  if (any(rates$rate_fertility < 0, na.rm = TRUE)) stop("ASFRs must be nonnegative.", call. = FALSE)
  if (any(!is.finite(rates$age_end) & rates$rate_fertility != 0, na.rm = TRUE))
    warning("TFR omits an open age band with non-zero fertility; TFR is underestimated over the full age range.", call. = FALSE)
  rates <- rates[is.finite(rates$age_end), , drop = FALSE]
  groups <- summary_strata(rates, by_union_status = TRUE)
  if (!nrow(rates)) {
    result <- rates[FALSE, c(groups, "age_start", "age_end"), drop = FALSE]
    result$tfr <- numeric()
    return(tibble::as_tibble(result))
  }
  tables <- lapply(summary_rows(rates, by_union_status = TRUE), function(rows) {
    r <- rates[rows, , drop = FALSE]
    check_summary_intervals(r, open = FALSE)
    key <- r[1, groups, drop = FALSE]
    cbind(key, data.frame(age_start = min(r$age_start), age_end = max(r$age_end),
                         tfr = sum(r$rate_fertility * (r$age_end - r$age_start))))
  })
  result <- do.call(rbind, tables)
  result <- tibble::as_tibble(result)
  attr(result, "unassigned_births") <- unknown
  result
}

#' Calculate parity- and age-adjusted total fertility rates
#'
#' Progress a synthetic cohort from parity zero through conditional age/parity
#' fertility schedules using the HFD period fertility-table midpoint convention.
#' Women in the final open parity band remain there and can have repeated births.
#' This is a deterministic fertility table without mortality or explicit birth
#' spacing, not a microsimulation or an exact continuous-time hazard solution.
#'
#' Separate the same period, group and union-type strata as [calculate_tfr()];
#' select female rows when sex is present. Parity itself is integrated over.
#' Each stratum has its own cohort: a union-type-specific index assumes that
#' type throughout, without partnership or group transitions.
#' Every age must have complete, nonoverlapping parity coverage from zero to
#' one open final category. Finite parity bounds must be nonnegative integers.
#' Each parity must cover the same contiguous finite age range. Missing included
#' rates yield NA for that stratum. No missing-rate zero fill is performed here.
#'
#' Upper-open age rows are omitted, warning on known nonzero omitted female
#' rates. The returned finite age range defines the calculation: all women
#' start at parity zero at its lower bound. Split at integer ages and input
#' boundaries into steps of at most one year. On each step of length d, use
#' progression probability `d*m / (1 + d*m/2)` and midpoint entry/exit exposure,
#' following HFD Methods Protocol (2015), sections 8.3-8.4. Probabilities above
#' one are rejected. Fractional steps generalize the annual convention.
#' HFD endpoint mapping and missing-cell choices can change agreement with its
#' published indices. See \url{https://www.humanfertility.org/File/GetDocumentFree/Docs/methods.pdf}.
#'
#' @param rates A fertility process or a table with `age_start`, `age_end`,
#'   `parity_start`, `parity_end` and `rate` or `rate_fertility` in births per
#'   woman-year at that current parity. Optional `time_start`, `time_end`,
#'   `union_type`, legacy `union_status` and any other non-measure column (such
#'   as a group dimension) define separate strata.
#' @return A tibble with the supplied stratification columns, finite `age_start`,
#'   `age_end`, and `patfr`. Empty results retain their column types. A supplied
#'   `unassigned_births` attribute is preserved.
#' @export
calculate_patfr <- function(rates) {
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "fertility") || !identical(rates$rate_unit, "person_year"))
      stop("PATFR requires a fertility process with person-year rates.", call. = FALSE)
    rates <- rates$rates
  }
  check_columns(rates, c("age_start", "age_end", "parity_start", "parity_end"), "rates")
  if (!"rate_fertility" %in% names(rates)) {
    if (!"rate" %in% names(rates)) stop("Supply `rate_fertility` or fertility `rate`.", call. = FALSE)
    rates$rate_fertility <- rates$rate
  }
  check_groups(rates, "rates")
  unknown <- attr(rates, "unassigned_births")
  if ("sex" %in% names(rates)) rates <- rates[rates$sex == "f", , drop = FALSE]
  check_numeric_column(rates$age_start, "rates$age_start")
  if (!is.numeric(rates$age_end) || !is.null(dim(rates$age_end)) || anyNA(rates$age_end) ||
      any(rates$age_start < 0 | rates$age_end <= rates$age_start))
    stop("Age bands need nonnegative starts and greater ends, possibly Inf.", call. = FALSE)
  check_numeric_column(rates$parity_start, "rates$parity_start")
  end <- rates$parity_end
  if (!is.numeric(end) || !is.null(dim(end)) || anyNA(end) ||
      any(rates$parity_start < 0 | rates$parity_start != floor(rates$parity_start)) ||
      any(end <= rates$parity_start) || any(is.finite(end) & end != floor(end)))
    stop("Parity bands need nonnegative integer starts and greater integer ends or Inf.", call. = FALSE)
  check_numeric_column(rates$rate_fertility, "rates$rate_fertility", missing = TRUE)
  if (any(rates$rate_fertility < 0, na.rm = TRUE)) stop("Fertility rates must be nonnegative.", call. = FALSE)
  if (any(!is.finite(rates$age_end) & rates$rate_fertility != 0, na.rm = TRUE))
    warning("PATFR omits an open age band with non-zero fertility; PATFR is underestimated over the full age range.", call. = FALSE)
  rates <- rates[is.finite(rates$age_end), , drop = FALSE]
  groups <- summary_strata(rates, by_union_status = TRUE)
  if (!nrow(rates)) {
    result <- tibble::as_tibble(rates[FALSE, c(groups, "age_start", "age_end"), drop = FALSE])
    result$patfr <- numeric()
    attr(result, "unassigned_births") <- unknown
    return(result)
  }
  tables <- lapply(summary_rows(rates, by_union_status = TRUE), function(rows) {
    r <- rates[rows, , drop = FALSE]
    lower <- min(r$age_start); upper <- max(r$age_end)
    parity <- sort(unique(c(r$parity_start, r$parity_end[is.finite(r$parity_end)])))
    if (!any(is.infinite(r$parity_end)) || min(parity) != 0)
      stop("Parity coverage must start at zero and include an open final category.", call. = FALSE)
    top <- max(parity)
    states <- 0:top
    ages <- sort(unique(c(r$age_start, r$age_end, seq(ceiling(lower), floor(upper)))))
    ages <- ages[ages >= lower & ages <= upper]
    hazard <- matrix(NA_real_, nrow = length(ages)-1L, ncol = length(states))
    for (a in seq_len(nrow(hazard))) for (k in seq_along(states)) {
      cell <- which(r$age_start <= ages[a] & r$age_end >= ages[a+1L] &
                    r$parity_start <= states[k] & r$parity_end > states[k])
      if (length(cell) != 1L)
        stop("Each age/parity cell must be covered exactly once, without gaps or overlaps.", call. = FALSE)
      hazard[a,k] <- r$rate_fertility[cell]
    }
    value <- if (anyNA(hazard)) NA_real_ else patfr_midpoint(hazard, diff(ages))
    key <- r[1, groups, drop = FALSE]
    cbind(key, data.frame(age_start = lower, age_end = upper, patfr = value))
  })
  result <- tibble::as_tibble(do.call(rbind, tables))
  attr(result, "unassigned_births") <- unknown
  result
}

patfr_midpoint <- function(hazard, widths) {
  n <- ncol(hazard)
  state <- c(1, rep(0, n-1L))
  total <- 0
  for (a in seq_len(nrow(hazard))) {
    d <- widths[a]
    exposure <- state * d
    if (n > 1L) {
      outgoing <- seq_len(n-1L); incoming <- outgoing + 1L
      q <- d * hazard[a,outgoing] / (1 + d * hazard[a,outgoing]/2)
      if (any(!is.finite(q)) || any(q > 1))
        stop("Rates are too large for valid HFD midpoint progression probabilities.", call. = FALSE)
      exposure[outgoing] <- state[outgoing] * d * (1-q/2)
      exposure[incoming] <- exposure[incoming] + state[outgoing] * q * d/2
    }
    births <- exposure * hazard[a, ]
    total <- total + sum(births)
    if (n > 1L) {
      state[outgoing] <- state[outgoing] - births[outgoing]
      state[incoming] <- state[incoming] + births[outgoing]
    }
  }
  total
}
