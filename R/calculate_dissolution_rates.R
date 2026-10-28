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
#'   are `"dissolution"`, `"death_husband"`, `"death_wife"`, or `"conversion"`. Population objects
#'   need a known endpoint; observation cannot extend beyond it. For bare tables,
#'   ensure complete follow-up through the requested endpoint.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param dur_breaks Strictly increasing nonnegative duration boundaries in years.
#'   If the last boundary is finite, append Inf for an open final band.
#' @param by_group Logical, default FALSE. With TRUE, use unordered group pairs
#'   and return `group_pair` (e.g. `"A:A"`, `"A:B"`, `"B:B"`). Both directions
#'   of a cross-group union share one stratum. Population inputs require person
#'   groups and retain every unordered pair in the full period/duration grid.
#'   Bare union tables require `group_pair`; reversed labels are normalized.
#'   Pair labels use nonempty colon-free character representations of groups,
#'   ordered alphabetically (radix). Empty grouped inputs return a typed zero-row
#'   tibble; no pooled row is added.
#' @param by_union_type Logical, default FALSE. Separate cohabitation and marriage
#'   episodes using each type's own exposure. Requires a valid `union_type` column.
#'   Conversion stops cohabitation exposure without counting as dissolution;
#'   duration resets at the start of the linked marriage episode.
#' @return A tibble with `time_start`, `time_end`, `dur_start`, `dur_end`,
#'   `union_years`, `dissolutions`, and `rate_dissolution`, plus `group_pair` and
#'   `union_type` when requested. Period varies fastest within duration bands and group blocks.
#' @export
#' @examples
#' unions <- data.frame(union_id = 1:3, start_time = c(0, 1, 2),
#'   end_time = c(3, 4, NA),
#'   end_reason = c("dissolution", "death_husband", NA))
#' calculate_dissolution_rates(unions, c(0, 2, 5), c(0, 2, Inf))
calculate_dissolution_rates <- function(population, time_breaks, dur_breaks,
                                   by_group = FALSE, by_union_type = FALSE) {
  if (!is.logical(by_group) || length(by_group) != 1L ||
      !is.null(dim(by_group)) || is.na(by_group))
    stop("`by_group` must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(by_union_type) || length(by_union_type) != 1L ||
      !is.null(dim(by_union_type)) || is.na(by_union_type))
    stop("`by_union_type` must be TRUE or FALSE.", call. = FALSE)
  groups <- NULL
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
  if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    unions <- population$unions
    if (by_group) {
      if (!"group" %in% names(population$persons))
        stop("`by_group = TRUE` requires population group identities.", call. = FALSE)
      persons <- population$persons
      groups <- possible_union_pairs(unique(persons$group))
      labels <- canonical_union_pairs(
        persons$group[match(unions$wife_id, persons$person_id)],
        persons$group[match(unions$husband_id, persons$person_id)])
      group_indices <- match(labels, groups) - 1L
    } else group_indices <- rep(0L, nrow(unions))
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
        any(!is.na(reason) & !reason %in% c("dissolution", "death_husband", "death_wife", "conversion")) ||
        any(is.na(reason) != is.na(unions$end_time)))
      stop("Ended unions require a valid end reason; ongoing unions require NA end and reason.", call. = FALSE)
    if (by_group) {
      if (!"group_pair" %in% names(unions))
        stop("`by_group = TRUE` requires a union `group_pair` column.", call. = FALSE)
      labels <- canonical_union_pairs(unions$group_pair)
      groups <- sort(unique(labels), method = "radix")
      group_indices <- match(labels, groups) - 1L
    } else group_indices <- rep(0L, nrow(unions))
  }
  group_count <- if (by_group) length(groups) else 1L
  types <- c("cohab", "marriage")
  if (by_union_type) {
    if (!"union_type" %in% names(unions))
      stop("`by_union_type = TRUE` requires a union `union_type` column.", call. = FALSE)
    if (!is.character(unions$union_type) || !is.null(dim(unions$union_type)) ||
        anyNA(unions$union_type) || any(!unions$union_type %in% types))
      stop("`union_type` must contain only cohab or marriage without missing values.", call. = FALSE)
    group_indices <- group_indices + group_count * (match(unions$union_type, types) - 1L)
  }
  grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
    d = seq_len(length(dur_breaks) - 1L), g = seq_len(group_count),
    type = if (by_union_type) types else "pooled",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  result <- tibble::tibble(time_start = time_breaks[grid$t],
    time_end = time_breaks[grid$t + 1L], dur_start = dur_breaks[grid$d],
    dur_end = dur_breaks[grid$d + 1L])
  if (by_group) result$group_pair <- groups[grid$g]
  if (by_union_type) result$union_type <- grid$type
  counts <- dissolution_counts_cpp(unions, time_breaks, dur_breaks,
                                   group_indices, group_count * if (by_union_type) 2L else 1L)
  result$union_years <- counts$union_years
  result$dissolutions <- counts$dissolutions
  result$rate_dissolution <- as.numeric(ifelse(result$union_years > 0,
    result$dissolutions / result$union_years, NA_real_))
  result
}
