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
#' When the supplied `rates` table contains columns other than `age_start`,
#' `age_end`, and either `rate_mortality` or `rate`, the function will treat
#' these additional columns as stratifying dimensions and will stratify the
#' output appropriately. It will also stratify output by time when `time_start`
#' and `time_end` are included. The function also ignores additional calculated
#' columns from the output of [calculate_mortality_rates()], so the results of
#' that function can be directly piped into this function. Stratification
#' produces separate life tables by strata, not multistate life tables.
#'
#' Constructed life tables do not make any adjustment to the \eqn{a_x} or
#' \eqn{q_x} values for the 0 to 1 age group. commonly used in demographic
#' methods. This is to ensure that results are consistent with input rates used
#' by the simulation. Users who want to model changes in mortality within the
#' first year of life can provide fractional year age bands for the first year
#' of life.
#'
#' Missing values on rates will result in missing values for lx, Lx, Tx, and ex
#' for any later age bands. A zero mortality rate in the open age band will
#' result in infinite life expectancy.
#'
#' @param rates A table returned by [calculate_mortality_rates()], a mortality
#'   process object from [mortality()], or a data frame with at least
#'   `age_start`, `age_end`, and `rate_mortality` or `rate` columns. Must
#'   contain an open-ended age group specified by an `Inf` value for `age_end`.
#'   Input may be further stratified by other variables.
#' @param radix Initial number of life-table survivors. Defaults to 1.
#' @return A tibble containing standard life table parameters as well as `sex`
#'   and additional stratifying columns (when supplied).
#' @seealso [calculate_life_expectancy()] for life expectancy at selected ages.
#' @export
#' @examples
#' data(usa_2000_mortality)
#'
#' calculate_life_table(usa_2000_mortality)
#'
#' # stratify by groups to get stratified life tables
#' stratify_rates(
#'   usa_2000_mortality,
#'   by = list(ethnicity = c(A = 1.2, B = 0.9))
#' ) |>
#'   calculate_life_table()
calculate_life_table <- function(rates, radix = 1) {
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "mortality")) {
      stop("Life tables require a mortality process.", call. = FALSE)
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
    rates[[rate_column]], paste0("rates$", rate_column),
    missing = TRUE
  )
  rates$mx <- rates[[rate_column]]
  if (any(rates$mx < 0, na.rm = TRUE)) {
    stop("Mortality rates must be nonnegative.", call. = FALSE)
  }
  if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix) ||
        radix <= 0) {
    stop("`radix` must be positive and finite.", call. = FALSE)
  }
  strata <- summary_strata(rates)
  rates <- rates[c(strata, "age_start", "age_end", "mx")]
  if (!nrow(rates)) {
    stop("A rate table needs at least one age band.", call. = FALSE)
  }
  for (column in strata) {
    if (anyNA(rates[[column]])) {
      stop(
        sprintf("`rates$%s` must not contain missing values.", column),
        call. = FALSE
      )
    }
  }
  # split rates by stratifying variables and calculate life tables for each
  tables <- lapply(
    summary_rows(rates),
    function(rows) {
      r <- rates[rows, , drop = FALSE]
      check_summary_intervals(r, strata = strata)
      r <- r[order(r$age_start), , drop = FALSE]
      if (r$age_start[1] != 0 || tail(r$age_end, 1) != Inf) {
        stop(
          paste0(
            "Life tables require ages from zero through an ",
            "explicit open interval", stratum_label(r, strata), "."
          ),
          call. = FALSE
        )
      }
      n <- r$age_end - r$age_start
      m <- r$mx
      size <- nrow(r)
      q <- ax <- L <- lx <- lend <- rep(NA_real_, size)
      survivors <- radix
      # step through each age and calculate life table measures
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
      # put it all together in a table and calculate summary measures
      r$n <- n
      r$ax <- ax
      r$qx <- q
      r$lx <- lx
      r$lx_end <- lend
      r$dx <- lx - lend
      r$Lx <- L
      r$Tx <- rev(cumsum(rev(L)))
      r$ex <- ifelse(lx > 0, r$Tx / lx, NA_real_)
      r
    }
  )
  result <- do.call(rbind, tables)
  tibble::as_tibble(result[c(strata, life_table_columns)])
}

#' Calculate life expectancy at given ages
#'
#' @description
#'
#' A convenience wrapper around [calculate_life_table()] that returns only
#' remaining life expectancy (`ex`) at the requested ages, for each stratum of
#' the input rates.
#'
#' @inheritParams calculate_life_table
#' @param age Ages at which to report life expectancy, in years. Each must be
#'   the start of an age band in `rates`. Defaults to `0`, life expectancy at
#'   birth.
#' @return A tibble with the stratifying columns of the input (when supplied),
#'   `age`, and `ex`, with one row per stratum and requested age.
#' @seealso [calculate_life_table()]
#' @export
#' @examples
#' calculate_life_expectancy(usa_2000_mortality)
#' calculate_life_expectancy(usa_2000_mortality, age = c(0, 30, 65))
calculate_life_expectancy <- function(rates, age = 0) {
  table <- calculate_life_table(rates)
  if (!is.numeric(age) || !length(age) || !all(age %in% table$age_start)) {
    stop(
      "`age` must contain ages that start an age band in `rates`.",
      call. = FALSE
    )
  }
  strata <- setdiff(names(table), life_table_columns)
  result <- table[table$age_start %in% age, c(strata, "age_start", "ex")]
  names(result)[names(result) == "age_start"] <- "age"
  result
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Age-band and measure columns of a life table; any other column identifies a
# stratum.
life_table_columns <- c(
  "age_start", "age_end", "mx", "n", "ax", "qx", "lx", "lx_end", "dx", "Lx",
  "Tx", "ex"
)
