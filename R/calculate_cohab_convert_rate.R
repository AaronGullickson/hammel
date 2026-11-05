#' Calculate cohabitation-to-marriage conversion rates
#'
#' Count conversions over cohabitation union-years, once per episode. Marriage
#' episodes contribute neither exposure nor events. Conversion, dissolution and
#' partner death end cohabitation exposure; only `end_reason = "conversion"`
#' counts as an event. Ongoing episodes are censored at the observation endpoint.
#' Cells are half-open; boundary events enter the cell beginning there, and
#' events at the final time boundary are excluded. Zero-duration episodes can
#' contribute events without exposure. Rates are NA when exposure is zero.
#'
#' @param population A `hammel_population` with typed union histories and a known
#'   endpoint. Observation cannot extend beyond it.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param dur_breaks Optional strictly increasing nonnegative duration boundaries,
#'   in years since cohabitation began. Supply exactly one of `dur_breaks` and
#'   `age_breaks`. A finite final boundary gains an open-ended band through Inf.
#' @param age_breaks Optional strictly increasing nonnegative boundaries for the
#'   woman's attained age, matching the conversion schedule's age convention.
#' @param by `NULL` (default) or a character vector: `"group_w"` (the woman's
#'   current group, matching conversion input schedules) and `"group_h"` (the
#'   man's). Requires person groups. Rows vary by time fastest, then the axis
#'   band, then the `by` variables in the order given.
#' @inheritParams calculate_mortality_rates
#' @return A tibble with `time_start`, `time_end`, `dur_start`, `dur_end` (or
#'   `age_start`, `age_end`), the `by` columns, `union_years`, `conversions`
#'   and `rate_cohab_convert`. Empty cells are retained. Rates are conversions per
#'   cohabitation union-year, not conversion probabilities with competing risks.
#'   With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.
#' @export
calculate_cohab_convert_rate <- function(population, time_breaks,
                                         dur_breaks = NULL, age_breaks = NULL,
                                         by = NULL, conf.level = NULL) {
  validate_rate_conf_level(conf.level)
  by <- check_by(by, c("group_w", "group_h"), "calculate_cohab_convert_rate")
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a hammel_population with typed union histories.", call. = FALSE)
  if (is.null(dur_breaks) == is.null(age_breaks))
    stop("Supply exactly one of `dur_breaks` and `age_breaks`.", call. = FALSE)
  validate_population(population)
  if (!"union_type" %in% names(population$unions))
    stop("Conversion observation requires typed union histories with `union_type`.", call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(breaks, if (duration) "dur_breaks" else "age_breaks", ages = TRUE)
  if (is.finite(tail(breaks, 1))) breaks <- c(breaks, Inf)
  if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
    stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
  unions <- population$unions[population$unions$union_type == "cohab", ]
  wives <- match(unions$wife_id, population$persons$person_id)
  origins <- if (duration) unions$start_time else population$persons$birth_time[wives]
  nt <- length(time_breaks) - 1L
  na <- length(breaks) - 1L
  strata <- union_strata(population, unions, by, stride = nt * na, typed_split = FALSE)
  axis <- if (duration) "dur" else "age"
  rows <- build_rows(c(list(band_dim(time_breaks, "time", 1), band_dim(breaks, axis, nt)),
                       strata$dims))
  counts <- dissolution_counts_cpp(strata$unions, time_breaks, breaks, strata$index, strata$count,
    origins = origins, event_reason = "conversion")
  result <- tibble::as_tibble(rows[c("time_start", "time_end", paste0(axis, c("_start", "_end")), by)])
  result$union_years <- counts$union_years[rows$.index]
  result$conversions <- counts$dissolutions[rows$.index]
  result$rate_cohab_convert <- as.numeric(ifelse(result$union_years > 0,
    result$conversions / result$union_years, NA_real_))
  add_rate_confidence(result, "conversions", "union_years", conf.level)
}
