# Calculate period life tables from mortality rate schedules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Construct period life tables from age-specific mortality rates
#'
#' @description
#'
#' Constructs a period life table from input tables that provide age-specific
#' mortality rates. Input data can come from user-provided input rate schedules,
#' mortality process objects, or estimated rates from a simulation via
#' [calculate_mortality_rates()], and may be stratified by additional variables.
#'
#' @details
#'
#' When the supplied `rates` table contains character or factor columns other
#' than `age_start`, `age_end`, `sex`, and either `rate_mortality` or `rate`,
#' the function will treat these additional columns as stratifying dimensions
#' and will attempt to stratify the output appropriately. It will also stratify
#' output by time when `time_start` and `time_end` are included.
#'
#' Constructed life tables make no adjustment for the ax values in the first
#' year of life, which is common in demographic methods. This is to ensure that
#' results are consistent with input rates used by the simulation. Users who
#' want to model changes in mortality within the first year of life can provide
#' fractional year age bands for the first year of life.
#'
#' Missing values on rates will result in missing values in the life table.
#'
#' @param rates A table returned by [calculate_mortality_rates()], a mortality
#'   process object from [mortality()], or a data frame with at least
#'   `age_start`, `age_end`, `sex`, and `rate_mortality` or `rate` columns. Must
#'   contain an open-ended age group specified by an `Inf` value for `age_end`.
#'   Input may be further stratified by other variables.
#' @param radix Initial number of life-table survivors. Defaults to 1.
#' @return A tibble containing standard life table parameters as well as `sex`
#'   and additional stratifying columns (when supplied).
#' @export
#' @examples
#' calculate_life_table(data.frame(
#'   age_start = c(0, 1, 5),
#'   age_end = c(1, 5, Inf), rate = 0.1
#' ))
calculate_life_table <- function(rates, radix = 1) {
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "mortality") ||
      !identical(rates$rate_unit, "person_year")) {
      stop("Life tables require a mortality process with person-year rates.",
        call. = FALSE
      )
    }
    rates <- rates$rates
  }
  check_columns(rates, c("age_start", "age_end"), "rates")
  if ("rate_mortality" %in% names(rates)) {
    rate_column <- "rate_mortality"
  } else if ("rate" %in% names(rates)) {
    rate_column <- "rate"
  } else {
    stop("Supply `rate_mortality` or `rate`.", call. = FALSE)
  }
  check_groups(rates, "rates")
  check_numeric_column(
    rates[[rate_column]], paste0("rates$", rate_column), missing = TRUE
  )
  rates$mx <- rates[[rate_column]]
  if (any(rates$mx < 0, na.rm = TRUE)) {
    stop(
      "Mortality rates must be nonnegative.",
      call. = FALSE
    )
  }
  if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix) || radix <=
    0) {
    stop("`radix` must be positive and finite.", call. = FALSE)
  }
  rates <- rates[c(
    summary_strata(rates, by_sex = TRUE, by_union_status = TRUE),
    "age_start", "age_end", "mx"
  )]
  if (!nrow(rates)) {
    for (column in c("n", "ax", "qx", "lx", "lx_end", "dx", "Lx", "Tx", "ex")) {
      rates[[column]] <- numeric()
    }
    return(tibble::as_tibble(rates))
  }
  tables <- lapply(
    summary_rows(rates, by_sex = TRUE, by_union_status = TRUE),
    function(rows) {
      r <- rates[rows, , drop = FALSE]
      check_summary_intervals(
        r,
        open = TRUE,
        strata = summary_strata(r, by_sex = TRUE, by_union_status = TRUE)
      )
      r <- r[order(r$age_start), , drop = FALSE]
      n <- r$age_end - r$age_start
      m <- r$mx
      size <- nrow(r)
      q <- ax <- L <- lx <- lend <- rep(NA_real_, size)
      survivors <- radix
      for (i in seq_len(size)) {
        lx[i] <- survivors
        if (is.na(m[i])) {
          survivors <- NA_real_
          next
        }
        if (is.infinite(n[i])) {
          q[i] <- if (m[i] > 0) 1 else 0
          survival <- if (m[i] > 0) 0 else 1
          ax[i] <- if (m[i] > 0) 1 / m[i] else Inf
          L[i] <- if (is.na(survivors)) {
            NA_real_
          } else if (survivors == 0) {
            0
          } else {
            survivors * ax[i]
          }
        } else {
          z <- m[i] * n[i]
          q[i] <- -expm1(-z)
          survival <- exp(-z)
          L[i] <- survivors * if (m[i] == 0) n[i] else q[i] / m[i]
          # Stable mean time lived by deaths within an exponential age interval.
          ax[i] <- if (z < 1e-5) {
            n[i] * (0.5 - z / 12 + z^3 / 720)
          } else {
            1 / m[i] - n[i] / expm1(z)
          }
        }
        lend[i] <- survivors * survival
        survivors <- lend[i]
      }
      r$n <- n; r$ax <- ax; r$qx <- q; r$lx <- lx; r$lx_end <- lend
      r$dx <- lx - lend; r$Lx <- L
      r$Tx <- rev(cumsum(rev(L)))
      r$ex <- ifelse(lx > 0, r$Tx / lx, NA_real_)
      r
    }
  )
  result <- do.call(rbind, tables)
  tibble::as_tibble(result)
}
