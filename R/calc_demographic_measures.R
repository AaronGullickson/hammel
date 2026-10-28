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
  groups <- observation_groups(persons, by_group)
  if (by_group) {
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

#' Summarize demographic exposure and rates
#'
#' Calculate person-years, deaths, births, ASMR, and ASFR for each time period,
#' age band, and sex. Births are classified by the mother's age and use all
#' female exposure, including time excluded by birth spacing. Births and ASFR
#' are NA for males. Empty cells are retained; rates with no exposure are NA.
#'
#' Lifetimes are split directly at time and age boundaries in years, using
#' half-open intervals. Missing death times are censored at the observation end.
#' For bare tables, ensure complete history through that endpoint. Population
#' objects cannot be observed beyond their known end_time. Analyze one independent
#' history per call; person IDs must be unique.
#'
#' Births with unknown mothers cannot be assigned maternal ages. Their counts
#' by period are returned in the unassigned_births attribute.
#'
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
#' @return A tibble with time_start, time_end, sex, age_start, age_end,
#'   person_years, deaths, births, rate_mortality, and rate_fertility. With `by_union_type = TRUE`,
#'   also union_type. With `by_group = TRUE`, also
#'   group; groups form separate blocks of rows.
#' @export
#' @examples
#' population <- create_baseline_population(100, end_time = 0)
#' calculate_demog_rates(population, c(-1, 0), c(0, 20, 40, Inf))
calculate_demog_rates <- function(population, time_breaks, age_breaks,
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
#' Open-ended age bands are omitted. The remaining finite, contiguous age range
#' defines the TFR reported. Any missing ASFR in those bands makes the TFR missing.
#' When `union_type` (or a legacy `union_status` column) is present, calculate
#' a separate sum for each stratum
#' within each period. These conditional sums describe fertility under continuous
#' exposure to that status; they are not components to add into an overall TFR.
#' A `group` column automatically adds separate group strata, with its type and
#' labels preserved; groups are not pooled. Age bands are checked independently
#' within each period/group/status combination.
#' @param rates A table returned by [calculate_demog_rates()], a fertility
#'   process object, or a data frame with `age_start`, `age_end`, and `rate_fertility`
#'   or `rate` in births per woman-year. Bare `rate` tables are interpreted as
#'   fertility rates; mortality processes are not valid inputs.
#' @return A tibble with period columns, `union_type` or legacy `union_status`
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
  rates <- rates[is.finite(rates$age_end), , drop = FALSE]
  groups <- intersect(c("time_start", "time_end", "group", "union_status", "union_type"), names(rates))
  check_numeric_column(rates$rate_fertility, "rates$rate_fertility", missing = TRUE)
  if (any(rates$rate_fertility < 0, na.rm = TRUE)) stop("ASFRs must be nonnegative.", call. = FALSE)
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
#' @param rates A table returned by [calculate_demog_rates()], a mortality
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
  unique(persons$group)
}
