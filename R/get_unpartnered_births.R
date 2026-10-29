#' Extract births with known parents to unpartnered mothers
#'
#' Return births where the mother was outside both cohabitation and marriage,
#' with parent ages and optional groups supplied by [get_parent_info()]. Known
#' mother and father IDs are required. Each birth is retained separately,
#' including repeated parental pairs; missing-parent records and founders are
#' omitted. These rows can be used to examine nonpartnership father matching.
#'
#' Histories must be complete for the intended interpretation. Partnership
#' coverage is half-open: births at formation are partnered. Births exactly at
#' dissolution or conversion retain the ending partnership (birth-first ties).
#' Births at death endings are unpartnered unless another episode covers that
#' time. Overlapping episodes do not duplicate births. Untyped partnerships are
#' treated as partnered too. No time window is imposed and no endpoint is needed.
#'
#' @param population A validated `hammel_population` with person and union
#'   histories, such as the result of [simulate()].
#' @return A separate tibble containing selected person rows and `age_mother`,
#'   `age_father`, plus `group_mother` and `group_father` when groups are present.
#'   Input order, column types and custom person columns are preserved. Empty
#'   results retain the enriched schema. The population is not modified.
#' @export
#' @examples
#' population <- create_baseline_population(10)
#' get_unpartnered_births(population)
get_unpartnered_births <- function(population) {
  parents <- get_parent_info(population)
  rows <- unpartnered_birth_rows_cpp(population$persons, population$unions)
  parents[rows, , drop = FALSE]
}
