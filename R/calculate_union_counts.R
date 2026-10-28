#' Calculate union-years, dissolutions and duration-specific rates
#'
#' Split each union episode across calendar periods and elapsed-duration bands,
#' counting union-years once per episode, not once per spouse. Death and
#' dissolution both stop exposure, but only `end_reason = "dissolution"` counts
#' toward dissolutions. Missing ends are censored at the observation endpoint.
#' Repeated episodes are separate observations; distinct overlapping records
#' are not merged. Analyze one independent history per call, then combine counts.
#'
#' Exposure and event cells are half-open. Dissolutions at internal boundaries
#' belong to the cell beginning there, and events at the last time boundary are
#' excluded. Zero-duration episodes may contribute events with no exposure.
#' Rates are dissolutions per union-year: zero with positive exposure and no
#' events, and NA when exposure is zero. This is a cause-specific rate, not the
#' probability of dissolution in the presence of competing partner deaths.
#'
#' @param population A `hammel_population`, or a union data frame containing
#'   unique positive whole-number `union_id`, `start_time`, `end_time`, and
#'   `end_reason`. Ends and reasons are NA for ongoing episodes; ended reasons
#'   are `"dissolution"`, `"death_husband"`, or `"death_wife"`. Population objects
#'   need a known endpoint; observation cannot extend beyond it. For bare tables,
#'   ensure complete follow-up through the requested endpoint.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param dur_breaks Strictly increasing nonnegative duration boundaries in years.
#'   If the last boundary is finite, append Inf for an open final band.
#' @param by_group Logical, default FALSE. With TRUE, population inputs require
#'   person groups and attribute each episode to the woman's group. Bare union
#'   tables require a `group` column specifying her group. Groups retain input
#'   type and first-occurrence order, with a complete grid for each and no pooled
#'   row. Empty grouped histories return a typed zero-row tibble.
#' @return A tibble with `time_start`, `time_end`, `dur_start`, `dur_end`,
#'   `union_years`, `dissolutions`, and `rate_dissolution`, plus `group` when
#'   requested. Period varies fastest within duration bands and group blocks.
#' @export
#' @examples
#' unions <- data.frame(union_id = 1:3, start_time = c(0, 1, 2),
#'   end_time = c(3, 4, NA),
#'   end_reason = c("dissolution", "death_husband", NA))
#' calculate_union_counts(unions, c(0, 2, 5), c(0, 2, Inf))
calculate_union_counts <- function(population, time_breaks, dur_breaks,
                                   by_group = FALSE) {
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
  if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    unions <- population$unions
    groups <- observation_groups(population$persons, by_group)
    group_indices <- if (by_group)
      match(population$persons$group[match(unions$wife_id, population$persons$person_id)], groups) - 1L else
      rep(0L, nrow(unions))
  } else {
    unions <- population
    check_columns(unions, c("union_id", "start_time", "end_time", "end_reason"), "unions")
    check_numeric_column(unions$union_id, "unions$union_id", id = TRUE)
    if (anyDuplicated(unions$union_id))
      stop("`unions$union_id` must be unique; analyze replications separately.", call. = FALSE)
    check_numeric_column(unions$start_time, "unions$start_time")
    check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
    if (any(unions$end_time < unions$start_time, na.rm = TRUE))
      stop("Union end times must not precede start times.", call. = FALSE)
    reason <- unions$end_reason
    if (!is.character(reason) || !is.null(dim(reason)) ||
        any(!is.na(reason) & !reason %in% c("dissolution", "death_husband", "death_wife")) ||
        any(is.na(reason) != is.na(unions$end_time)))
      stop("Ended unions require a valid end reason; ongoing unions require NA end and reason.", call. = FALSE)
    check_groups(unions, "unions")
    if (isTRUE(by_group) && !"group" %in% names(unions))
      stop("`by_group = TRUE` requires a union `group` column specifying the woman's group.", call. = FALSE)
    groups <- observation_groups(unions, by_group)
    group_indices <- if (by_group) match(unions$group, groups) - 1L else rep(0L, nrow(unions))
  }
  group_count <- if (by_group) length(groups) else 1L
  grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
    d = seq_len(length(dur_breaks) - 1L), g = seq_len(group_count),
    KEEP.OUT.ATTRS = FALSE)
  result <- tibble::tibble(time_start = time_breaks[grid$t],
    time_end = time_breaks[grid$t + 1L], dur_start = dur_breaks[grid$d],
    dur_end = dur_breaks[grid$d + 1L])
  if (by_group) result$group <- groups[grid$g]
  counts <- dissolution_counts_cpp(unions, time_breaks, dur_breaks,
                                   group_indices, group_count)
  result$union_years <- counts$union_years
  result$dissolutions <- counts$dissolutions
  result$rate_dissolution <- as.numeric(ifelse(result$union_years > 0,
    result$dissolutions / result$union_years, NA_real_))
  result
}
