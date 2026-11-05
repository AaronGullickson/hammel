#' Calculate the living population at time points
#'
#' Count persons with `birth_time <= time` and missing `death_time` or
#' `death_time > time`. Counts include births at the requested time and exclude
#' deaths at that time, representing the population after recorded events there.
#' States in `by` are also those after all events at each time.
#'
#' @param population A `hammel_population` or persons data frame representing
#'   one independent history. For bare tables, ensure follow-up is complete
#'   through every requested time. Union variables require a population.
#' @param times A numeric vector of finite time points in years. A population
#'   object cannot be queried beyond its known `end_time`.
#' @param by `NULL` (default) or a character vector of variables to count
#'   separately: `"sex"`, `"group"`, `"union_type"` and `"previous_union"`. See
#'   [calculate_mortality_rates()] for their definitions.
#' @param pool_dissolved Logical, default `FALSE`. With `TRUE` and
#'   `"previous_union"` in `by`, combine separated and widowed as `"dissolved"`.
#' @return A tibble with `time`, the `by` columns and `population`. Rows follow
#'   the requested times (repeats kept), with time varying fastest and earlier
#'   `by` variables varying faster than later ones.
#' @export
#' @examples
#' population <- create_baseline_population(100, end_time = 0)
#' calculate_population_size(population, times = c(-10, -5, 0))
#' calculate_population_size(population, times = 0, by = "sex")
calculate_population_size <- function(population, times, by = NULL, pool_dissolved = FALSE) {
  if (!is.numeric(times) || !is.null(dim(times)) || any(!is.finite(times)))
    stop("`times` must be a numeric vector of finite times in years.", call. = FALSE)
  by <- check_by(by, c("sex", "group", "union_type", "previous_union"), "calculate_population_size")
  check_pool_dissolved(pool_dissolved, by)
  obs <- person_observation(population, by, if (length(times)) max(times) else -Inf)
  previous <- "previous_union" %in% by
  utype <- "union_type" %in% by
  map <- if (pool_dissolved) c(0L, 1L, 1L) else 0:2
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  options <- list(typed = obs$typed, exclusive = obs$typed || previous,
    counts_cohab = obs$counts_cohab, sex = "sex" %in% by, group = "group" %in% by,
    union_type = utype, previous_union = previous, previous_map = map)
  counts <- person_snapshot_cpp(obs$persons, obs$unions, obs$group_indices, ng,
    as.numeric(times), options)
  sizes <- c(t = length(times), a = 1L, sex = if ("sex" %in% by) 2L else 1L, parity = 1L,
    group = if ("group" %in% by) ng else 1L,
    union_type = if (utype) length(union_type_labels(obs$typed)) else 1L,
    previous_union = if (previous) length(unique(map)) else 1L)
  strides <- cumprod(c(1, head(sizes, -1)))
  names(strides) <- names(sizes)
  dims <- c(list(dim_of(data.frame(time = times), stride = 1)),
    person_by_dims(by, obs, list(strides = strides), pool_dissolved))
  rows <- build_rows(dims)
  out <- finish_rows(rows[!grepl("^\\.key", names(rows))])
  out$population <- as.integer(counts[rows$.index])
  out
}

# Shared validation of period and age breaks; a finite final age gains Inf.
person_breaks <- function(time_breaks, age_breaks) {
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(age_breaks, "age_breaks", ages = TRUE)
  if (is.finite(tail(age_breaks, 1))) age_breaks <- c(age_breaks, Inf)
  age_breaks
}

#' Calculate age-specific mortality rates
#'
#' Report deaths per person-year, with exact half-open period and age splitting.
#'
#' @section Stratification with `by`:
#' `by` lists categorical variables to split exposure and events over. Each is
#' the person's state at the time, changing only at recorded events:
#' * `"group"`: current group, following the transition history. Requires a
#'   declared group dimension (see [set_groups()]).
#' * `"union_type"`: current union state, `"unpartnered"`, `"cohab"` or
#'   `"marriage"` for typed histories, or `"unpartnered"` and `"partnered"`
#'   for untyped ones.
#' * `"previous_union"`: how the most recent ended union ended, `"none"`,
#'   `"separated"` or `"widowed"`, including during later unions. Conversion
#'   is not an ending; cohabitation endings follow the population's
#'   cohabitation policy.
#'
#' Union variables require a `hammel_population` with union histories;
#' `"previous_union"` (and typed `"union_type"`) require nonoverlapping
#' histories for each person. Each event counts under the state just before it
#' in the simulation's event order: deaths see no same-time changes, so a death
#' that ends a union counts as partnered; births count in a union starting at
#' the same time and in one ending then by dissolution or conversion, but not
#' one ending by the partner's death, and use the mother's group before any
#' same-time transition.
#'
#' Rows vary by time fastest, then age, then sex, then the `by` variables in
#' the order given, earlier ones varying faster. Strata sum to pooled counts and
#' exposure; recompute rates after pooling.
#'
#' @param population A hammel_population or persons data frame.
#' @param time_breaks Strictly increasing finite period boundaries in years.
#' @param age_breaks Strictly increasing nonnegative age boundaries in years.
#'   If the final boundary is finite, Inf is appended to create an open-ended
#'   final age band starting at that boundary.
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over: `"group"`, `"union_type"` and `"previous_union"`. See
#'   below.
#' @param conf.level NULL (default) to omit confidence intervals, or a finite
#'   numeric scalar strictly between 0 and 1. Adds `conf_low` and `conf_high`
#'   from two-sided exact [stats::poisson.test()] intervals for observed event
#'   counts over exposure. Bounds are NA at zero exposure. For fertility, these
#'   describe the observed rate, not the unrealized-adjusted rate.
#' @param pool_dissolved Logical, default `FALSE`. With `TRUE` and
#'   `"previous_union"` in `by`, combine separated and widowed as
#'   `"dissolved"`, pooling counts and exposure before calculating rates.
#' @return A tibble containing period, sex, age, the `by` columns,
#'   person_years, deaths and rate_mortality.
#' @export
calculate_mortality_rates <- function(population, time_breaks, age_breaks, by = NULL,
                                      conf.level = NULL, pool_dissolved = FALSE) {
  validate_rate_conf_level(conf.level)
  by <- check_by(by, c("group", "union_type", "previous_union"), "calculate_mortality_rates")
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(obs, by, time_breaks, age_breaks, sex = TRUE,
                             pool_dissolved = pool_dissolved)
  out <- person_rate_rows(obs, by, counts, time_breaks, age_breaks, c("m", "f"), pool_dissolved)
  out$deaths <- counts$deaths[attr(out, "index")]
  out$rate_mortality <- as.numeric(ifelse(out$person_years > 0,
    out$deaths / out$person_years, NA_real_))
  attr(out, "index") <- NULL
  add_rate_confidence(out, "deaths", "person_years", conf.level)
}

# Rows for person rates: time, sex, age, optional parity, then `by`.
person_rate_rows <- function(obs, by, counts, time_breaks, age_breaks, sexes,
                             pool_dissolved, parity_breaks = NULL) {
  s <- counts$strides
  sex_dim <- label_dim("sex", sexes, s[["sex"]], code = match(sexes, c("m", "f")) - 1L)
  dims <- c(list(band_dim(time_breaks, "time", 1), band_dim(age_breaks, "age", s[["a"]]), sex_dim),
    if (!is.null(parity_breaks)) list(band_dim(parity_breaks, "parity", s[["parity"]])),
    person_by_dims(by, obs, counts, pool_dissolved))
  rows <- build_rows(dims)
  columns <- c("time_start", "time_end", "sex", "age_start", "age_end",
    if (!is.null(parity_breaks)) c("parity_start", "parity_end"), by)
  out <- tibble::as_tibble(rows[columns])
  out$person_years <- counts$person_years[rows$.index]
  attr(out, "index") <- rows$.index
  out
}

#' Calculate age-specific fertility rates, optionally by current parity
#'
#' Report maternal births per woman-year, including time excluded by minimum
#' birth spacing. Only female rows are returned. Current parity is reconstructed
#' from all recorded maternal births, including pre-window births and deceased
#' children; custom person parity columns are ignored. Each birth counts at the
#' parity immediately before it and advances subsequent exposure. Simultaneous
#' births progress successively in recorded person-row order.
#'
#' `by` works as in [calculate_mortality_rates()]: births count under the
#' mother's state just before them.
#' @inheritParams calculate_mortality_rates
#' @param parity_breaks NULL to pool parities, or strictly increasing nonnegative
#'   integer boundaries starting at zero. Bands are half-open. A finite final
#'   boundary receives Inf, as for age breaks.
#' @param unrealized Logical, default `FALSE`. With `TRUE`, the population must
#'   come from `simulate()`; add
#'   `unrealized_births` (counted births that could not occur because their
#'   fertility cell had members but no eligible woman, from
#'   `simulation$unrealized_births`) and
#'   `rate_fertility_adjusted = (births + unrealized_births) / person_years`.
#'   Values are `NA` where they cannot be attributed exactly: periods before the
#'   most recent simulation's start (restarts keep only the latest run's
#'   record), breaks that cut through a recorded one-year period or fertility
#'   cell, and `by` or parity strata the fertility schedule did not use
#'   (`"previous_union"` always).
#' @return A tibble with period, sex (female), age, optional parity_start and
#'   parity_end, the `by` columns, person_years, births and rate_fertility, and
#'   with `unrealized = TRUE`, unrealized_births and rate_fertility_adjusted.
#'   Unknown-mother births are in the period-level unassigned_births attribute.
#' @export
calculate_fertility_rates <- function(population, time_breaks, age_breaks, by = NULL,
                                      parity_breaks = NULL, unrealized = FALSE,
                                      conf.level = NULL, pool_dissolved = FALSE) {
  validate_rate_conf_level(conf.level)
  if (!is.logical(unrealized) || length(unrealized) != 1L || is.na(unrealized))
    stop("`unrealized` must be TRUE or FALSE.", call. = FALSE)
  if (unrealized) {
    record <- if (inherits(population, "hammel_population")) population$simulation else NULL
    if (is.null(record$unrealized_births))
      stop(paste("`unrealized = TRUE` requires a population returned by simulate()",
                 "with an unrealized-births record."), call. = FALSE)
  }
  by <- check_by(by, c("group", "union_type", "previous_union"), "calculate_fertility_rates")
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  if (!is.null(parity_breaks)) {
    check_demographic_breaks(parity_breaks, "parity_breaks", ages = TRUE)
    finite <- parity_breaks[is.finite(parity_breaks)]
    if (parity_breaks[1] != 0 || any(finite != floor(finite)))
      stop("`parity_breaks` must start at zero and contain integer boundaries or Inf.", call. = FALSE)
    if (is.finite(tail(parity_breaks, 1))) parity_breaks <- c(parity_breaks, Inf)
  }
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(obs, by, time_breaks, age_breaks, sex = TRUE,
    parity_breaks = if (is.null(parity_breaks)) numeric() else parity_breaks,
    pool_dissolved = pool_dissolved)
  out <- person_rate_rows(obs, by, counts, time_breaks, age_breaks, "f", pool_dissolved,
                          parity_breaks)
  out$births <- counts$births[attr(out, "index")]
  out$rate_fertility <- as.numeric(ifelse(out$person_years > 0,
    out$births / out$person_years, NA_real_))
  attr(out, "index") <- NULL
  if (unrealized) out <- add_unrealized_births(out, population$simulation)
  attr(out, "unassigned_births") <- tibble::tibble(time_start = head(time_breaks, -1),
    time_end = tail(time_breaks, -1), births = counts$unassigned_births)
  add_rate_confidence(out, "births", "person_years", conf.level)
}

# Attribute recorded unrealized births (D101) to observation rows. A record is
# one run-anchored annual period (or a fragment split at a fertility-rate
# change) and one fertility cell; it is added to
# the row whose period, age, parity, group and union type contain it. Rows that
# a record only partly overlaps, or that a stratum the schedule did not use
# would have to split, become NA, as do periods before the run start.
add_unrealized_births <- function(out, simulation) {
  u <- simulation$unrealized_births
  n <- nrow(out)
  value <- numeric(n)
  covered <- out$time_start >= simulation$start_time
  unknown <- !covered
  if (nrow(u) && n) {
    # Interval indices of each record's bounds within each observation axis.
    within <- function(lo, hi, starts, ends) {
      breaks <- sort(unique(c(starts, ends)))
      i <- findInterval(lo, breaks)
      j <- findInterval(hi, breaks, left.open = TRUE)
      list(lo = breaks[pmax(i, 1L)], ok = i >= 1L & i == j & i < length(breaks))
    }
    t <- within(u$period_start, u$period_end, out$time_start, out$time_end)
    a <- within(u$age_start, u$age_end, out$age_start, out$age_end)
    keys <- list(t$lo, a$lo)
    row_keys <- list(out$time_start, out$age_start)
    ok <- t$ok & a$ok
    # Records outside every observation row contribute nothing anywhere.
    inside <- u$period_start < max(out$time_end) & u$period_end > min(out$time_start) &
      u$age_start < max(out$age_end) & u$age_end > min(out$age_start)
    if ("parity_start" %in% names(out)) {
      if ("parity_start" %in% names(u)) {
        p <- within(u$parity_start, u$parity_end, out$parity_start, out$parity_end)
        keys <- c(keys, list(p$lo)); row_keys <- c(row_keys, list(out$parity_start))
        ok <- ok & p$ok
      } else ok[] <- FALSE # pooled-parity schedule cannot be split by parity
    }
    if ("group" %in% names(out)) {
      keys <- c(keys, list(as.character(u$group)))
      row_keys <- c(row_keys, list(as.character(out$group)))
    }
    # The fertility schedule never splits by previous union.
    if ("previous_union" %in% names(out)) ok[] <- FALSE
    if ("union_type" %in% names(out)) {
      if ("union_type" %in% names(u)) {
        ok <- ok & u$union_type %in% out$union_type
        keys <- c(keys, list(u$union_type)); row_keys <- c(row_keys, list(out$union_type))
      } else ok[] <- FALSE # overall schedule cannot be split by union type
    }
    key <- do.call(paste, c(keys, sep = "\r"))
    row_key <- do.call(paste, c(row_keys, sep = "\r"))
    row <- match(key, row_key)
    exact <- inside & ok & !is.na(row)
    sums <- tapply(u$unrealized_births[exact], row[exact], sum)
    value[as.integer(names(sums))] <- as.numeric(sums)
    # Records that overlap rows without fitting one row make those rows unknown.
    for (r in which(inside & !exact)) {
      hit <- out$time_start < u$period_end[r] & out$time_end > u$period_start[r] &
        out$age_start < u$age_end[r] & out$age_end > u$age_start[r]
      if ("parity_start" %in% names(out) && "parity_start" %in% names(u))
        hit <- hit & out$parity_start < u$parity_end[r] & out$parity_end > u$parity_start[r]
      if ("group" %in% names(out)) hit <- hit & as.character(out$group) == as.character(u$group[r])
      if ("union_type" %in% names(out) && "union_type" %in% names(u) &&
          u$union_type[r] %in% out$union_type)
        hit <- hit & out$union_type == u$union_type[r]
      unknown <- unknown | hit
    }
  }
  value[unknown] <- NA_real_
  out$unrealized_births <- value
  out$rate_fertility_adjusted <- as.numeric(ifelse(out$person_years > 0,
    (out$births + out$unrealized_births) / out$person_years, NA_real_))
  out
}

# Summary calculations use explicit built-in grouping columns.
summary_rows <- function(table, by_sex = FALSE, by_union_status = FALSE) {
  columns <- intersect(c("time_start", "time_end", "group", if (by_sex) "sex",
                         if (by_union_status) c("union_status", "union_type", "previous_union")),
                     names(table))
  if (!length(columns)) return(list(seq_len(nrow(table))))
  codes <- lapply(table[columns], function(column) match(column, unique(column)))
  keys <- do.call(paste, c(codes, sep = ":"))
  split(seq_len(nrow(table)), factor(keys, levels = unique(keys)))
}

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
#' A `group` column automatically adds separate group strata, with its type and
#' labels preserved; groups are not pooled. Age bands are checked independently
#' within each period/group/status combination.
#' @param rates A table returned by [calculate_fertility_rates()], a fertility
#'   process object, or a data frame with `age_start`, `age_end`, and `rate_fertility`
#'   or `rate` in births per woman-year. Bare `rate` tables are interpreted as
#'   fertility rates; mortality processes are not valid inputs.
#' @return With parity columns, the result of [calculate_patfr()] with `patfr`.
#'   Otherwise, a tibble with period columns, `union_type`, `previous_union` or
#'   legacy `union_status` when supplied,
#'   `group` when supplied, `age_start`, `age_end`, and `tfr`.
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
  groups <- intersect(c("time_start", "time_end", "group", "union_status", "union_type",
                        "previous_union"), names(rates))
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
  rates <- rates[c(intersect(c("time_start", "time_end", "sex", "group", "union_status", "union_type",
                                "previous_union"), names(rates)),
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

check_demographic_breaks <- function(breaks, name, ages = FALSE) {
  if (!is.numeric(breaks) || !is.null(dim(breaks)) || length(breaks) < 2L || anyNA(breaks) ||
      any(!is.finite(head(breaks, -1))) || any(diff(breaks) <= 0) ||
      (!ages && !is.finite(tail(breaks, 1))) || (ages && breaks[1] < 0)) {
    stop(sprintf("`%s` must be strictly increasing %sboundaries.", name,
                 if (ages) {
                   if (name == "dur_breaks") "nonnegative duration " else "nonnegative age "
                 } else "finite time "), call. = FALSE)
  }
}

check_summary_intervals <- function(rates, open) {
  if (!nrow(rates)) stop("A rate table needs at least one age band.", call. = FALSE)
  check_numeric_column(rates$age_start, "rates$age_start")
  ends <- rates$age_end
  if (!is.numeric(ends) || anyNA(ends) || any(rates$age_start < 0) ||
      any(ends <= rates$age_start)) stop("Invalid age intervals.", call. = FALSE)
  ordered <- rates[order(rates$age_start), , drop = FALSE]
  if (nrow(ordered) > 1L && any(head(ordered$age_end, -1) != tail(ordered$age_start, -1)))
    stop("Age intervals must be contiguous and nonoverlapping within each group.", call. = FALSE)
  if (open && (ordered$age_start[1] != 0 || tail(ordered$age_end, 1) != Inf))
    stop("Life tables require ages from zero through an explicit open interval.", call. = FALSE)
  if (!open && any(!is.finite(ends)))
    stop("TFR requires finite age intervals.", call. = FALSE)
}

# Population validation already checks label types and missing values.
observation_groups <- function(persons, by_group) {
  if (!is.logical(by_group) || length(by_group) != 1L ||
      !is.null(dim(by_group)) || is.na(by_group))
    stop("`by_group` must be TRUE or FALSE.", call. = FALSE)
  if (!by_group) return(NULL)
  if (!"group" %in% names(persons))
    stop("`by_group = TRUE` requires a persons `group` column.", call. = FALSE)
  groups <- attr(persons, ".hammel_group_values")
  if (is.null(groups)) unique(persons$group) else groups
}
