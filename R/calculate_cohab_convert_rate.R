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
#' @param by_group Logical, default FALSE. Separate by the woman's group, matching
#'   conversion input schedules. Require person groups; retain their input type.
#' @return A tibble with `time_start`, `time_end`, `dur_start`, `dur_end` (or
#'   `age_start`, `age_end`), `union_years`, `conversions`, `rate_cohab_convert`,
#'   and optional `group`. Empty cells are retained. Rates are conversions per
#'   cohabitation union-year, not conversion probabilities with competing risks.
#' @export
calculate_cohab_convert_rate <- function(population, time_breaks,
                                         dur_breaks = NULL, age_breaks = NULL,
                                         by_group = FALSE) {
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
  persons <- population$persons
  persons <- prepare_group_observations(persons, population, by_group)
  groups <- observation_groups(persons, by_group)
  group_count <- if (by_group) length(groups) else 1L
  unions <- population$unions[population$unions$union_type == "cohab", ]
  wives <- match(unions$wife_id, persons$person_id)
  indices <- if (by_group) match(persons$group[wives], groups) - 1L else rep(0L, nrow(unions))
  origins <- if (duration) unions$start_time else persons$birth_time[wives]
  grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
    a = seq_len(length(breaks) - 1L), g = seq_len(group_count), KEEP.OUT.ATTRS = FALSE)
  result <- tibble::tibble(time_start = time_breaks[grid$t], time_end = time_breaks[grid$t + 1L])
  axis <- if (duration) "dur" else "age"
  result[[paste0(axis, "_start")]] <- breaks[grid$a]
  result[[paste0(axis, "_end")]] <- breaks[grid$a + 1L]
  if (by_group) result$group <- groups[grid$g]
  if (by_group) unions <- prepare_union_group_observations(unions, population, groups)
  counts <- dissolution_counts_cpp(unions, time_breaks, breaks, indices, group_count,
    origins = origins, event_reason = "conversion")
  result$union_years <- counts$union_years
  result$conversions <- counts$dissolutions
  result$rate_cohab_convert <- as.numeric(ifelse(result$union_years > 0,
    result$conversions / result$union_years, NA_real_))
  result
}
