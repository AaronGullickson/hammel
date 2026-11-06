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
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over: `"<group>_w"` and `"<group>_h"` for declared group
#'   dimensions, such as `"race_w"` and `"region_h"` (the wife's and
#'   husband's current categories, following transition histories;
#'   directional, so sum A-wife/B-husband and B-wife/A-husband rows for an
#'   unordered pair; any combination of dimensions and partners) and
#'   `"union_type"` (cohab or marriage, each with its own exposure;
#'   conversion stops cohabitation exposure without counting as dissolution,
#'   and duration resets at the start of the linked marriage). Population
#'   inputs require declared groups for the group variables; bare union tables
#'   require `group_w`, `group_h` or `union_type` columns. A dissolution counts
#'   under the groups before any same-time transition.
#' @inheritParams calculate_mortality_rates
#' @return A tibble with `time_start`, `time_end`, `dur_start`, `dur_end`,
#'   the `by` columns, `union_years`, `dissolutions`, and `rate_dissolution`.
#'   Rows vary by time fastest, then duration, then the `by` variables in the
#'   order given.
#'   With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.
#' @export
#' @examples
#' unions <- data.frame(union_id = 1:3, start_time = c(0, 1, 2),
#'   end_time = c(3, 4, NA),
#'   end_reason = c("dissolution", "death_husband", NA))
#' calculate_dissolution_rates(unions, c(0, 2, 5), c(0, 2, Inf))
calculate_dissolution_rates <- function(population, time_breaks, dur_breaks, by = NULL,
                                        conf.level = NULL) {
  validate_rate_conf_level(conf.level)
  alias <- observation_alias(population, by, c("_w", "_h"))
  population <- alias$population
  by <- alias$by
  by <- check_by(by, c("group_w", "group_h", "union_type"), "calculate_dissolution_rates", alias$dims)
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
  if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    unions <- population$unions
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
  }
  nt <- length(time_breaks) - 1L
  nd <- length(dur_breaks) - 1L
  strata <- union_strata(population, unions, by, stride = nt * nd)
  rows <- build_rows(c(list(band_dim(time_breaks, "time", 1), band_dim(dur_breaks, "dur", nt)),
                       strata$dims))
  counts <- dissolution_counts_cpp(strata$unions, time_breaks, dur_breaks,
                                   strata$index, strata$count)
  result <- tibble::as_tibble(rows[c("time_start", "time_end", "dur_start", "dur_end", by)])
  result$union_years <- counts$union_years[rows$.index]
  result$dissolutions <- counts$dissolutions[rows$.index]
  result$rate_dissolution <- as.numeric(ifelse(result$union_years > 0,
    result$dissolutions / result$union_years, NA_real_))
  result <- finish_strata(result, alias, "dissolutions", "union_years", "rate_dissolution")
  add_rate_confidence(result, "dissolutions", "union_years", conf.level)
}
