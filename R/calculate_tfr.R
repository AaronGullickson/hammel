# Calculate TFR and parity-adjusted TFR from fertility rate schedules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Calculate total fertility and parity-adjusted total fertility rates
#'
#' @description
#'
#' Constructs total fertility rates or parity-adjusted total fertility rats from
#' input tables that provide age-specific fertility rates. Input data can come
#' from user-provided input rate schedules, fertility process objects, or
#' estimated rates from a simulation via [calculate_fertility_rates], and may be
#' stratified by additional variables.
#'
#' @details
#'
#' When the supplied `rates` table contains columns other than `age_start`,
#' `age_end`, and either `rate_fertility` or `rate`, the function will treat
#' these additional columns as stratifying dimensions and will stratify the
#' output appropriately. It will also stratify output by time when `time_start`
#' and `time_end` are included. The function also ignores additional calculated
#' columns from the output of [calculate_fertility_rates], so the results of
#' that function can be directly piped into this function.
#'
#' When the input `rates` table is stratified by parity with `parity_start` and
#' `parity_end` columns, the function `calculate_patfr` will be automatically
#' used in place of `calculate_tfr` to calculate parity-adjusted total fertility
#' rates. See See
#' \url{https://www.humanfertility.org/File/GetDocumentFree/Docs/methods.pdf}
#' for details on the estimation of parity-adjusted total fertility rates.
#'
#' The function will automatically remove any open-ended age bands because TFR
#' cannot be estimated across open-ended age bands. If this band contains a
#' non-zero rate, the function will produce a warning that TFR may be
#' underestimated.
#'
#' Missing values on rates will result in missing values for the total fertility
#' rate.
#'
#' @param rates A table returned by [calculate_fertility_rates()], a fertility
#'   process object from [fertility()], or a data frame with at least
#'   `age_start`, `age_end`, and `rate_fertility` or `rate` columns. Input may
#'   be further stratified by other variables.
#' @return A tibble with the supplied period and other stratification columns,
#'   finite `age_start` and `age_end`, and either `tfr` or `patfr`.
#'
#' @seealso [calculate_fertility_rates()], [fertility()]
#' @rdname calculate_tfr
#' @order 1
#' @export
#' @examples
#' # TFR from age-specific rates
#' calculate_tfr(usa_2000_fertility)
#'
#' # TFR stratified by union type
#' calculate_tfr(usa_2000_fertility_union)
calculate_tfr <- function(rates) {
  # Errors raised while summarizing each stratum name this function.
  call <- rlang::current_env()
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "fertility")) {
      abort_process_kind("fertility")
    }
    rates <- rates$rates
  }
  if (any(c("parity_start", "parity_end") %in% names(rates))) {
    return(calculate_patfr(rates))
  }
  check_columns(rates, c("age_start", "age_end"), "rates")
  if (!"rate_fertility" %in% names(rates)) {
    if (!"rate" %in% names(rates)) {
      abort_rate_column(c("rate_fertility", "rate"))
    }
    rates$rate_fertility <- rates$rate
  }
  check_groups(rates, "rates")
  if (!nrow(rates)) {
    abort_no_bands()
  }
  unknown <- attr(rates, "unassigned_births")
  if ("sex" %in% names(rates)) {
    # TFR describes women only, so sex is never a stratum or output column.
    rates <- rates[rates$sex == "f", setdiff(names(rates), "sex"), drop = FALSE]
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  if (!is.numeric(rates$age_end) || anyNA(rates$age_end)) {
    abort_hammel(
      c(
        "Column {.field age_end} of {.arg rates} must contain numbers.",
        "x" = "It is {describe_value(rates$age_end)}.",
        "i" = "Use {.code Inf} for an open-ended band."
      ),
      "table",
      argument = "rates$age_end"
    )
  }
  check_numeric_column(
    rates$rate_fertility, "rates$rate_fertility",
    missing = TRUE
  )
  if (any(rates$rate_fertility < 0, na.rm = TRUE)) {
    check_nonnegative(rates$rate_fertility, "rates$rate_fertility")
  }
  if (any(
        !is.finite(rates$age_end) & rates$rate_fertility != 0,
        na.rm = TRUE
      )) {
    warn_open_fertility_band("TFR")
  }
  rates <- rates[is.finite(rates$age_end), , drop = FALSE]
  groups <- summary_strata(rates)
  if (!nrow(rates)) {
    abort_no_finite_bands("TFR")
  }
  tables <- lapply(
    summary_rows(rates), function(rows) {
      r <- rates[rows, , drop = FALSE]
      check_summary_intervals(r, strata = groups, call = call)
      key <- r[1, groups, drop = FALSE]
      cbind(
        key, data.frame(
          age_start = min(r$age_start),
          age_end = max(r$age_end),
          tfr = sum(r$rate_fertility * (r$age_end - r$age_start))
        )
      )
    }
  )
  result <- do.call(rbind, tables)
  result <- tibble::as_tibble(result)
  attr(result, "unassigned_births") <- unknown
  result
}

#' @rdname calculate_tfr
#' @order 2
#' @export
calculate_patfr <- function(rates) {
  # Errors raised while summarizing each stratum name this function.
  call <- rlang::current_env()
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "fertility")) {
      abort_process_kind("fertility")
    }
    rates <- rates$rates
  }
  check_columns(
    rates, c("age_start", "age_end", "parity_start", "parity_end"), "rates"
  )
  if (!"rate_fertility" %in% names(rates)) {
    if (!"rate" %in% names(rates)) {
      abort_rate_column(c("rate_fertility", "rate"))
    }
    rates$rate_fertility <- rates$rate
  }
  check_groups(rates, "rates")
  if (!nrow(rates)) {
    abort_no_bands()
  }
  unknown <- attr(rates, "unassigned_births")
  if ("sex" %in% names(rates)) {
    # TFR describes women only, so sex is never a stratum or output column.
    rates <- rates[rates$sex == "f", setdiff(names(rates), "sex"), drop = FALSE]
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  if (!is.numeric(rates$age_end) || !is.null(dim(rates$age_end)) ||
        anyNA(rates$age_end) ||
        any(rates$age_start < 0 | rates$age_end <= rates$age_start)) {
    abort_hammel(
      c(
        paste(
          "Age bands in {.arg rates} must have nonnegative starts and",
          "greater ends."
        ),
        "i" = "Use {.code Inf} for an open-ended band."
      ),
      "rates",
      argument = "rates"
    )
  }
  check_numeric_column(rates$parity_start, "rates$parity_start")
  end <- rates$parity_end
  if (!is.numeric(end) || !is.null(dim(end)) || anyNA(end) ||
        any(
          rates$parity_start < 0 |
            rates$parity_start != floor(rates$parity_start)
        ) ||
        any(end <= rates$parity_start) ||
        any(is.finite(end) & end != floor(end))) {
    abort_hammel(
      c(
        paste(
          "Parity bands in {.arg rates} must have nonnegative whole-number",
          "starts and greater whole-number ends."
        ),
        "i" = "Use {.code Inf} for the open-ended last parity band."
      ),
      "rates",
      argument = "rates"
    )
  }
  check_numeric_column(
    rates$rate_fertility, "rates$rate_fertility",
    missing = TRUE
  )
  if (any(rates$rate_fertility < 0, na.rm = TRUE)) {
    check_nonnegative(rates$rate_fertility, "rates$rate_fertility")
  }
  if (any(
        !is.finite(rates$age_end) & rates$rate_fertility != 0,
        na.rm = TRUE
      )) {
    warn_open_fertility_band("PATFR")
  }
  rates <- rates[is.finite(rates$age_end), , drop = FALSE]
  groups <- summary_strata(rates)
  if (!nrow(rates)) {
    abort_no_finite_bands("PATFR")
  }
  tables <- lapply(
    summary_rows(rates), function(rows) {
      r <- rates[rows, , drop = FALSE]
      lower <- min(r$age_start)
      upper <- max(r$age_end)
      parity <- sort(
        unique(c(r$parity_start, r$parity_end[is.finite(r$parity_end)]))
      )
      if (!any(is.infinite(r$parity_end)) || min(parity) != 0) {
        abort_hammel(
          c(
            "Parity bands must start at 0 and end with an open band.",
            "x" = "The parity boundaries are {.val {parity}}."
          ),
          "rates",
          argument = "rates",
          call = call
        )
      }
      top <- max(parity)
      states <- 0:top
      ages <- sort(
        unique(c(r$age_start, r$age_end, seq(ceiling(lower), floor(upper))))
      )
      ages <- ages[ages >= lower & ages <= upper]
      hazard <- matrix(
        NA_real_,
        nrow = length(ages) - 1L,
        ncol = length(states)
      )
      for (a in seq_len(nrow(hazard))) {
        for (k in seq_along(states)) {
          cell <- which(
            r$age_start <= ages[a] & r$age_end >= ages[a + 1L] &
              r$parity_start <= states[k] & r$parity_end > states[k]
          )
          if (length(cell) != 1L) {
            abort_hammel(
              c(
                "Each age and parity combination must have exactly one rate.",
                "i" = "Age and parity bands must meet without gaps or overlaps."
              ),
              "rates",
              argument = "rates",
              call = call
            )
          }
          hazard[a, k] <- r$rate_fertility[cell]
        }
      }
      value <- if (anyNA(hazard)) {
        NA_real_
      } else {
        patfr_midpoint(hazard, diff(ages), call = call)
      }
      key <- r[1, groups, drop = FALSE]
      cbind(key, data.frame(age_start = lower, age_end = upper, patfr = value))
    }
  )
  result <- tibble::as_tibble(do.call(rbind, tables))
  attr(result, "unassigned_births") <- unknown
  result
}

# HELPER FUNCTIONS ---------------------------------------------------------

patfr_midpoint <- function(hazard, widths, call = rlang::caller_env()) {
  n <- ncol(hazard)
  state <- c(1, rep(0, n - 1L))
  total <- 0
  for (a in seq_len(nrow(hazard))) {
    d <- widths[a]
    exposure <- state * d
    if (n > 1L) {
      outgoing <- seq_len(n - 1L)
      incoming <- outgoing + 1L
      q <- d * hazard[a, outgoing] / (1 + d * hazard[a, outgoing] / 2)
      if (any(!is.finite(q)) || any(q > 1)) {
        abort_hammel(
          c(
            "Fertility rates are too large for PATFR.",
            "i" = paste(
              "Midpoint parity-progression probabilities would exceed 1;",
              "check the rates and band widths."
            )
          ),
          "rates",
          argument = "rates",
          call = call
        )
      }
      exposure[outgoing] <- state[outgoing] * d * (1 - q / 2)
      exposure[incoming] <- exposure[incoming] + state[outgoing] * q * d / 2
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
