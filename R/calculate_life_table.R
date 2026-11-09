# Calculate period life tables from mortality rate schedules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Construct period life tables from mortality rates
#'
#' Treat mortality rates as piecewise constant hazards. Closed-interval survival
#' is exponential, and person-years are integrated under that same assumption.
#' This does not apply conventional HMD infant adjustments. Missing rates remain
#' unknown. A zero terminal rate gives infinite remaining lifetime when terminal
#' survivorship is positive. Require a final explicit open age interval.
#' Tables are constructed separately for each period, sex, group, and union status
#' present. Status-specific tables hold that status fixed; they do not model
#' transitions between statuses.
#'
#' @param rates A table returned by [calculate_mortality_rates()], a mortality
#'   process object, or a data frame with `age_start`, `age_end`, and `rate_mortality`,
#'   `mx`, or `rate` per person-year.
#' @param radix Initial life-table survivors, a positive finite number.
#' @return A tibble with period, sex, group, and union-status identifiers (when supplied),
#'   `age_start`, `age_end`, `mx`, `n`, `ax`, `qx`, `lx`, `lx_end`, `dx`,
#'   `Lx`, `Tx`, and `ex`. Summary counts and exposure columns are omitted.
#'   Life expectancy is `Tx / lx` in years.
#' @export
#' @examples
#' calculate_life_table(data.frame(age_start = c(0, 1, 5),
#'   age_end = c(1, 5, Inf), mx = 0.1))
calculate_life_table <- function(rates, radix = 1) {
  if (inherits(rates, "hammel_process")) {
    if (!identical(rates$name, "mortality") ||
        !identical(rates$rate_unit, "person_year"))
      stop("Life tables require a mortality process with person-year rates.", call. = FALSE)
    rates <- rates$rates
  }
  check_columns(rates, c("age_start", "age_end"), "rates")
  if (!"mx" %in% names(rates)) {
    if ("rate_mortality" %in% names(rates)) rates$mx <- rates$rate_mortality
    else if ("rate" %in% names(rates)) rates$mx <- rates$rate
    else stop("Supply `rate_mortality`, `mx`, or `rate`.", call. = FALSE)
  }
  check_groups(rates, "rates")
  check_numeric_column(rates$mx, "rates$mx", missing = TRUE)
  if (any(rates$mx < 0, na.rm = TRUE)) stop("Mortality rates must be nonnegative.", call. = FALSE)
  if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix) || radix <= 0)
    stop("`radix` must be positive and finite.", call. = FALSE)
  rates <- rates[c(summary_strata(rates, by_sex = TRUE, by_union_status = TRUE),
                   "age_start", "age_end", "mx")]
  if (!nrow(rates)) {
    for (column in c("n", "ax", "qx", "lx", "lx_end", "dx", "Lx", "Tx", "ex"))
      rates[[column]] <- numeric()
    return(tibble::as_tibble(rates))
  }
  tables <- lapply(summary_rows(rates, by_sex = TRUE, by_union_status = TRUE), function(rows) {
    r <- rates[rows, , drop = FALSE]
    check_summary_intervals(r, open = TRUE)
    r <- r[order(r$age_start), , drop = FALSE]
    n <- r$age_end - r$age_start
    m <- r$mx
    size <- nrow(r)
    q <- ax <- L <- lx <- lend <- rep(NA_real_, size)
    survivors <- radix
    for (i in seq_len(size)) {
      lx[i] <- survivors
      if (is.na(m[i])) { survivors <- NA_real_; next }
      if (is.infinite(n[i])) {
        q[i] <- if (m[i] > 0) 1 else 0
        survival <- if (m[i] > 0) 0 else 1
        ax[i] <- if (m[i] > 0) 1/m[i] else Inf
        L[i] <- if (is.na(survivors)) NA_real_ else if (survivors == 0) 0 else survivors * ax[i]
      } else {
        z <- m[i] * n[i]
        q[i] <- -expm1(-z)
        survival <- exp(-z)
        L[i] <- survivors * if (m[i] == 0) n[i] else q[i]/m[i]
        # Stable mean time lived by deaths within an exponential age interval.
        ax[i] <- if (z < 1e-5) n[i] * (0.5 - z/12 + z^3/720) else
          1/m[i] - n[i]/expm1(z)
      }
      lend[i] <- survivors * survival
      survivors <- lend[i]
    }
    r$n <- n; r$ax <- ax; r$qx <- q; r$lx <- lx; r$lx_end <- lend
    r$dx <- lx - lend; r$Lx <- L
    r$Tx <- rev(cumsum(rev(L)))
    r$ex <- ifelse(lx > 0, r$Tx/lx, NA_real_)
    r
  })
  result <- do.call(rbind, tables)
  tibble::as_tibble(result)
}
