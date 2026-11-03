#' Calculate the living population at time points
#'
#' Count persons with `birth_time <= time` and missing `death_time` or
#' `death_time > time`. Counts include births at the requested time and exclude
#' deaths at that time, representing the population after recorded events there.
#' Results preserve the order of the requested times, including repeated times.
#'
#' @param population A `hammel_population` or persons data frame representing
#'   one independent history. For bare tables, ensure follow-up is complete
#'   through every requested time.
#' @param times A numeric vector of finite time points in years. A population
#'   object cannot be queried beyond its known `end_time`.
#' @param by_group Logical, default `FALSE`. With `TRUE`, require a persons
#'   `group` column and count separately by group. Groups retain their input type
#'   and first-occurrence order within each requested time. No pooled group row
#'   is added. Empty histories return a typed zero-row table.
#' @return A tibble with `time` and `population` columns, plus `group` when
#'   `by_group = TRUE`.
#' @export
#' @examples
#' population <- create_baseline_population(100, end_time = 0)
#' calculate_population_size(population, times = c(-10, -5, 0))
calculate_population_size <- function(population, times, by_group = FALSE) {
  if (!is.numeric(times) || !is.null(dim(times)) || any(!is.finite(times)))
    stop("`times` must be a numeric vector of finite times in years.", call. = FALSE)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || any(times > population$end_time))
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    persons <- population$persons
  } else {
    if (!is.data.frame(population))
      stop("`population` must be a hammel_population or persons data frame.", call. = FALSE)
    persons <- population
    hammel_population(persons)
  }
  persons <- prepare_group_observations(persons, population, by_group)
  groups <- observation_groups(persons, by_group)
  if (by_group) {
    if (inherits(population, "hammel_population") && !is.null(population$transitions) && nrow(population$transitions)) {
      # Births enter the initial group, deaths leave the final group, and
      # every transition adds one arrival and one departure. Count sorted
      # event dates instead of rebuilding each person's identity per snapshot.
      tx <- population$transitions
      history <- attr(persons, ".hammel_group_history")
      initial <- match(persons$group, groups)
      ordered <- order(tx$time)
      first <- ordered[!duplicated(history$row[ordered])]
      initial[history$row[first] + 1L] <- history$from[first] + 1L
      final <- match(persons$group, groups)
      counts <- lapply(seq_along(groups), function(g) {
        enter <- sort(c(persons$birth_time[initial == g], tx$time[history$to == g - 1L]))
        exit <- sort(c(persons$death_time[final == g & !is.na(persons$death_time)],
                       tx$time[history$from == g - 1L]))
        findInterval(times, enter) - findInterval(times, exit)
      })
      grid <- expand.grid(g = seq_along(groups), t = seq_along(times), KEEP.OUT.ATTRS = FALSE)
      indices <- grid$t + length(times) * (grid$g - 1L)
      return(tibble::tibble(time = times[grid$t], group = groups[grid$g],
        population = as.integer(unlist(counts, use.names = FALSE)[indices])))
    }
    grid <- expand.grid(g = seq_along(groups), t = seq_along(times),
                        KEEP.OUT.ATTRS = FALSE)
    counts <- lapply(seq_along(groups), function(g) {
      rows <- persons$group == groups[g]
      born <- sort(persons$birth_time[rows])
      died <- sort(persons$death_time[rows & !is.na(persons$death_time)])
      findInterval(times, born) - findInterval(times, died)
    })
    # Each group vector follows requested time order; the grid interleaves groups.
    indices <- grid$t + length(times) * (grid$g - 1L)
    return(tibble::tibble(time = times[grid$t], group = groups[grid$g],
      population = as.integer(unlist(counts, use.names = FALSE)[indices])))
  }
  # Sorted event dates let us count each time without scanning every person.
  born <- sort(persons$birth_time)
  died <- sort(persons$death_time[!is.na(persons$death_time)])
  counts <- findInterval(times, born) - findInterval(times, died)
  tibble::tibble(time = times, population = counts)
}

# Shared validation, grids, and native observation counts for process summaries.
collect_demographic_counts <- function(population, time_breaks, age_breaks,
                                   by_union_type = FALSE, by_group = FALSE) {
  if (!is.logical(by_union_type) || length(by_union_type) != 1L ||
      !is.null(dim(by_union_type)) || is.na(by_union_type))
    stop("`by_union_type` must be TRUE or FALSE.", call. = FALSE)
  if (by_union_type && !inherits(population, "hammel_population"))
    stop("`by_union_type = TRUE` requires a hammel_population with union histories.",
         call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(age_breaks, "age_breaks", ages = TRUE)
  if (is.finite(tail(age_breaks, 1))) age_breaks <- c(age_breaks, Inf)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    persons <- population$persons
  } else {
    if (!is.data.frame(population))
      stop("`population` must be a hammel_population or persons data frame.", call. = FALSE)
    persons <- population
    hammel_population(persons)
  }
  persons <- prepare_group_observations(persons, population, by_group)
  groups <- observation_groups(persons, by_group)
  group_count <- if (by_group) length(groups) else 1L
  group_indices <- if (by_group) match(persons$group, groups) - 1L else
    rep(0L, nrow(persons))
  nt <- length(time_breaks) - 1L
  na <- length(age_breaks) - 1L
  grid <- expand.grid(t = seq_len(nt), a = seq_len(na), sex = c("m", "f"),
                      g = seq_len(group_count),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  table <- tibble::tibble(time_start = time_breaks[grid$t], time_end = time_breaks[grid$t + 1L],
    sex = grid$sex, age_start = age_breaks[grid$a], age_end = age_breaks[grid$a + 1L])
  if (by_group) table$group <- groups[grid$g]
  counts <- demographic_counts_cpp(persons, time_breaks, age_breaks,
                                   group_indices, group_count)
  for (column in c("person_years", "deaths", "births")) table[[column]] <- counts[[column]]
  table$births[table$sex == "m"] <- NA_real_
  table$rate_mortality <- as.numeric(ifelse(table$person_years > 0,
                                table$deaths / table$person_years, NA_real_))
  table$rate_fertility <- as.numeric(ifelse(table$person_years > 0,
                                table$births / table$person_years, NA_real_))
  if (by_union_type) {
    total <- table
    typed <- "union_type" %in% names(population$unions)
    types <- if (typed) unique(population$unions$union_type) else "partnered"
    if (typed) {
      split <- typed_demographic_counts_cpp(persons, population$unions,
        time_breaks, age_breaks, group_indices, group_count, types)
    } else {
      x <- union_counts_cpp(persons, population$unions, time_breaks, age_breaks,
        group_indices, group_count)
      split <- list(person_years = x$person_years_partnered,
        births = x$births_partnered, deaths = x$deaths_partnered)
    }
    parts <- lapply(seq_along(types), function(k) {
      x <- total
      rows <- seq_len(nrow(total)) + (k - 1L) * nrow(total)
      for (column in c("person_years", "births", "deaths"))
        x[[column]] <- split[[column]][rows]
      x$union_type <- types[k]
      x
    })
    unpartnered <- total
    for (column in c("person_years", "births", "deaths")) {
      used <- Reduce(`+`, lapply(parts, `[[`, column), init = rep(0, nrow(total)))
      unpartnered[[column]] <- total[[column]] - used
    }
    unpartnered$person_years <- pmax(0, unpartnered$person_years)
    unpartnered$union_type <- "unpartnered"
    parts <- c(list(unpartnered), parts)
    table <- do.call(rbind, parts)
    # Keep types adjacent within each original demographic cell.
    order <- as.vector(t(matrix(seq_len(nrow(table)), nrow = nrow(total),
      ncol = length(parts))))
    table <- table[order, ]
    table$births[table$sex == "m"] <- NA_real_
    table$rate_mortality <- as.numeric(ifelse(table$person_years > 0,
      table$deaths / table$person_years, NA_real_))
    table$rate_fertility <- as.numeric(ifelse(table$person_years > 0,
      table$births / table$person_years, NA_real_))
    table <- table[c("time_start", "time_end", "sex", "age_start", "age_end",
      if (by_group) "group", "union_type", "person_years", "deaths", "births", "rate_mortality", "rate_fertility")]
  }
  attr(table, "unassigned_births") <- tibble::tibble(time_start = head(time_breaks, -1),
    time_end = tail(time_breaks, -1), births = counts$unassigned_births)
  table
}

#' Calculate age-specific mortality rates
#'
#' Report deaths per person-year, with exact half-open period and age splitting.
#' @param population A hammel_population or persons data frame.
#' @param time_breaks Strictly increasing finite period boundaries in years.
#' @param age_breaks Strictly increasing nonnegative age boundaries in years.
#'   If the final boundary is finite, Inf is appended to create an open-ended
#'   final age band starting at that boundary.
#' @param by_union_type Logical, default FALSE. Require population histories and
#'   return `union_type` rows for unpartnered and each discovered union type,
#'   without overall rows. Untyped histories use partnered/unpartnered. Births
#'   tied to dissolution or conversion retain the ending type; deaths use the
#'   type immediately before death. Typed overlapping histories are rejected;
#'   untyped overlapping coverage is merged. For formation rates use
#'   [calculate_union_rates()] or [calculate_marriage_rates()].
#' @param by_group Logical, default `FALSE`. With `TRUE`, require a persons
#'   `group` column and return the complete period/age/sex grid for each group
#'   appearing in the person history, in first-occurrence order. Labels retain
#'   their type. Exposure and deaths use the person's group; births use the
#'   mother's group, even when the child has a different group.
#'   Combines with `by_union_type`. No pooled
#'   group row is added. Unknown-mother births remain in the single period-level
#'   `unassigned_births` attribute, without a group assignment.
#' @param conf.level NULL (default) to omit confidence intervals, or a finite
#'   numeric scalar strictly between 0 and 1. Adds `conf_low` and `conf_high`
#'   from two-sided exact [stats::poisson.test()] intervals for observed event
#'   counts over exposure. Bounds are NA at zero exposure. For fertility, these
#'   describe the observed rate, not the unrealized-adjusted rate.
#' @return A tibble containing period, age, sex, person_years, deaths and
#'   rate_mortality, with optional group and union_type strata.
#' @export
calculate_mortality_rates <- function(population, time_breaks, age_breaks,
                                     by_union_type = FALSE, by_group = FALSE,
                                     conf.level = NULL) {
  validate_rate_conf_level(conf.level)
  x <- collect_demographic_counts(population, time_breaks, age_breaks,
                            by_union_type, by_group)
  x$births <- x$rate_fertility <- NULL
  attr(x, "unassigned_births") <- NULL
  add_rate_confidence(x, "deaths", "person_years", conf.level)
}

#' Calculate age-specific fertility rates, optionally by current parity
#'
#' Report maternal births per woman-year, including time excluded by minimum
#' birth spacing. Only female rows are returned. Current parity is reconstructed
#' from all recorded maternal births, including pre-window births and deceased
#' children; custom person parity columns are ignored. Each birth counts at the
#' parity immediately before it and advances subsequent exposure. Simultaneous
#' births progress successively in recorded person-row order.
#' @inheritParams calculate_mortality_rates
#' @param parity_breaks NULL to pool parities, or strictly increasing nonnegative
#'   integer boundaries starting at zero. Bands are half-open. A finite final
#'   boundary receives Inf, as for age breaks.
#' @param unrealized Logical, default `FALSE`. With `TRUE`, the population must
#'   come from `simulate()`; add
#'   `unrealized_births` (expected births dropped from fertility cells that had
#'   members but no eligible woman, from `simulation$unrealized_births`) and
#'   `rate_fertility_adjusted = (births + unrealized_births) / person_years`.
#'   Values are `NA` where they cannot be attributed exactly: periods before the
#'   most recent simulation's start (restarts keep only the latest run's
#'   record), breaks that cut through a recorded one-year period or fertility
#'   cell, and `union_type` or parity strata the fertility schedule did not use.
#' @return A tibble with period, age, sex (female), person_years, births and
#'   rate_fertility; optional parity_start/parity_end, group and union_type, and
#'   with `unrealized = TRUE`, unrealized_births and rate_fertility_adjusted.
#'   Unknown-mother births are in the period-level unassigned_births attribute.
#' @export
calculate_fertility_rates <- function(population, time_breaks, age_breaks,
                                     by_union_type = FALSE, by_group = FALSE,
                                     parity_breaks = NULL, unrealized = FALSE,
                                     conf.level = NULL) {
  validate_rate_conf_level(conf.level)
  if (!is.logical(unrealized) || length(unrealized) != 1L || is.na(unrealized))
    stop("`unrealized` must be TRUE or FALSE.", call. = FALSE)
  if (unrealized) {
    record <- if (inherits(population, "hammel_population")) population$simulation else NULL
    if (is.null(record$unrealized_births))
      stop(paste("`unrealized = TRUE` requires a population returned by simulate()",
                 "with an unrealized-births record."), call. = FALSE)
  }
  # Shared validation and pooled schema; retain complete family/union histories.
  pooled <- collect_demographic_counts(population, time_breaks, age_breaks,
                                 by_union_type, by_group)
  unknown <- attr(pooled, "unassigned_births")
  female <- pooled[pooled$sex == "f", ]
  female$deaths <- female$rate_mortality <- NULL
  if (is.null(parity_breaks)) {
    if (unrealized) female <- add_unrealized_births(female, population$simulation)
    return(add_rate_confidence(female, "births", "person_years", conf.level))
  }
  check_demographic_breaks(parity_breaks, "parity_breaks", ages = TRUE)
  finite <- parity_breaks[is.finite(parity_breaks)]
  if (parity_breaks[1] != 0 || any(finite != floor(finite)))
    stop("`parity_breaks` must start at zero and contain integer boundaries or Inf.", call. = FALSE)
  if (is.finite(tail(parity_breaks, 1))) parity_breaks <- c(parity_breaks, Inf)
  if (is.finite(tail(age_breaks, 1))) age_breaks <- c(age_breaks, Inf)
  persons <- if (inherits(population, "hammel_population")) population$persons else population
  persons <- prepare_group_observations(persons, population, by_group)
  groups <- observation_groups(persons, by_group)
  ng <- if (by_group) length(groups) else 1L
  gi <- if (by_group) match(persons$group, groups) - 1L else rep(0L, nrow(persons))
  types <- if (by_union_type) unique(female$union_type) else "pooled"
  # Empty grouped populations still return the full typed schema.
  grid <- expand.grid(t = seq_len(length(time_breaks)-1L),
    a = seq_len(length(age_breaks)-1L), g = seq_len(ng),
    p = seq_len(length(parity_breaks)-1L), k = seq_along(types))
  unions <- if (by_union_type) population$unions else data.frame()
  counts <- fertility_parity_counts_cpp(persons, unions, time_breaks,
    age_breaks, parity_breaks, gi, ng, types, by_union_type)
  out <- tibble::tibble(time_start = time_breaks[grid$t],
    time_end = time_breaks[grid$t+1L], sex = rep("f", nrow(grid)),
    age_start = age_breaks[grid$a], age_end = age_breaks[grid$a+1L])
  if (by_group) out$group <- groups[grid$g]
  if (by_union_type) out$union_type <- types[grid$k]
  out$parity_start <- parity_breaks[grid$p]
  out$parity_end <- parity_breaks[grid$p+1L]
  out$person_years <- counts$person_years
  out$births <- counts$births
  out$rate_fertility <- as.numeric(ifelse(out$person_years > 0,
    out$births / out$person_years, NA_real_))
  if (unrealized) out <- add_unrealized_births(out, population$simulation)
  attr(out, "unassigned_births") <- unknown
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
                         if (by_union_status) c("union_status", "union_type")), names(table))
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
#' When `union_type` (or a legacy `union_status` column) is present, calculate
#' a separate sum for each stratum
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
#'   Otherwise, a tibble with period columns, `union_type` or legacy `union_status`
#'   when supplied,
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
  groups <- intersect(c("time_start", "time_end", "group", "union_status", "union_type"), names(rates))
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
  rates <- rates[c(intersect(c("time_start", "time_end", "sex", "group", "union_status", "union_type"), names(rates)),
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
