#' Calculate marriage entry rates including cohabitation conversions
#'
#' Count direct marriages and conversions to marriage over nonmarried
#' person-years: unpartnered and cohabiting time both contribute exposure.
#' Married time contributes none. Each marriage episode counts once per spouse,
#' including a marriage beginning at conversion. Untyped unions are interpreted
#' as marriages. Men describe realized outcomes, not male search-rate inputs.
#'
#' Complete, nonoverlapping histories are required. Cells are half-open; starts
#' on boundaries enter the cell beginning there. Exposure is censored at death
#' and the observation endpoint. Rates are NA when exposure is zero. This
#' age-specific measure does not provide duration strata. With
#' `by = "previous_union"`, cohabiting exposure and conversion retain the
#' category held immediately before that cohabitation began. Ending categories
#' follow the recorded cohabitation policy, as in [calculate_union_rates()].
#' For entry into any partnership over unpartnered exposure, use
#' [calculate_union_rates()].
#'
#' @inheritParams calculate_union_rates
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over: `"group"` and `"previous_union"`, as in
#'   [calculate_union_rates()].
#' @return A tibble with `time_start`, `time_end`, `sex`, `age_start`, `age_end`,
#'   the `by` columns, `person_years`, `marriages`, and `rate_marriage`.
#'   Empty cells are retained. Person-years include cohabiting time.
#'   With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.
#' @export
calculate_marriage_rates <- function(population, time_breaks, age_breaks, by = NULL,
                                     conf.level = NULL, pool_dissolved = FALSE) {
  validate_rate_conf_level(conf.level)
  by <- check_by(by, c("group", "previous_union"), "calculate_marriage_rates")
  calculate_entry_rates(population, time_breaks, age_breaks, by, dur_breaks = NULL,
    marriage_measure = TRUE, conf.level = conf.level, pool_dissolved = pool_dissolved)
}
