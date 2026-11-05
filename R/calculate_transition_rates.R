#' Calculate age- or duration-specific group-transition rates
#'
#' Count individual moves and person-years in their origin group in half-open
#' Lexis cells. Each outgoing destination shares its origin's exposure: do not
#' sum denominators across destinations. Duration is years since the latest
#' group entry, or age before the first recorded transition; return moves reset
#' duration and union changes do not. Events use the exiting group's duration.
#' Death and competing moves end exposure. Rates are cause-specific hazards,
#' not transition probabilities. Internal-boundary events enter the cell
#' starting there; the final calendar boundary is excluded.
#'
#' @param population A grouped `hammel_population` with a known endpoint and
#'   complete person, transition and (when requested) union histories.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param age_breaks Optional strictly increasing nonnegative age boundaries.
#' @param dur_breaks Optional strictly increasing nonnegative group-duration
#'   boundaries. Supply exactly one of age_breaks and dur_breaks. A finite final
#'   axis boundary gains Inf.
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over: `"sex"`, `"union_type"` (current union state; typed
#'   histories use unpartnered/cohab/marriage, untyped unpartnered/partnered)
#'   and `"previous_union"`, as in [calculate_mortality_rates()]. Union
#'   variables require exclusive histories. A transition counts under the union
#'   state after same-time union changes, which precede transitions.
#' @param pool_dissolved Logical, default `FALSE`. With `TRUE` and
#'   `"previous_union"` in `by`, combine separated and widowed as `"dissolved"`.
#' @inheritParams calculate_mortality_rates
#' @return A tibble with `time_start`, `time_end`, age or duration bounds,
#'   `group_origin`, `group_dest`, the `by` columns, `person_years`,
#'   `transitions`, and `rate_transition` (NA at zero exposure). All distinct
#'   directed pairs among current and historical identities are retained,
#'   including zero-event cells. Fewer than two known groups returns a typed
#'   empty tibble. Rows vary by time fastest, then the axis band, then
#'   destination, then origin, then the `by` variables in the order given.
#'   With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1, sex = "f", birth_time = 0,
#'   death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_, group = "B")
#' population <- hammel_population(persons, end_time = 5,
#'   transitions = data.frame(person_id = 1, time = 2,
#'     group_origin = "A", group_dest = "B"))
#' calculate_transition_rates(population, c(0, 5), age_breaks = c(0, Inf))
calculate_transition_rates <- function(population, time_breaks, age_breaks = NULL,
                                       dur_breaks = NULL, by = NULL, conf.level = NULL,
                                       pool_dissolved = FALSE) {
  validate_rate_conf_level(conf.level)
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a grouped hammel_population.", call. = FALSE)
  validate_population(population)
  if (!declared_group(population))
    stop("`calculate_transition_rates()` requires a declared group dimension.", call. = FALSE)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(by, c("sex", "union_type", "previous_union"), "calculate_transition_rates")
  check_pool_dissolved(pool_dissolved, by)
  if (is.null(age_breaks) == is.null(dur_breaks))
    stop("Supply exactly one of `age_breaks` and `dur_breaks`.", call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(breaks, if (duration) "dur_breaks" else "age_breaks", ages = TRUE)
  if (is.finite(tail(breaks, 1))) breaks <- c(breaks, Inf)
  obs <- person_observation(population, by, tail(time_breaks, 1), grouped = TRUE,
                            validate = FALSE)
  union_vars <- any(c("union_type", "previous_union") %in% by)
  counts <- run_person_lexis(obs, setdiff(by, "sex"), time_breaks, breaks,
    axis = if (duration) 1L else 0L, sex = "sex" %in% by, group = TRUE,
    pool_dissolved = pool_dissolved, exclusive = union_vars)
  s <- counts$strides
  groups <- obs$groups
  dims <- c(list(band_dim(time_breaks, "time", 1),
    band_dim(breaks, if (duration) "dur" else "age", s[["a"]]),
    label_dim("group_dest", groups, counts$cells),
    label_dim("group_origin", groups, s[["group"]])),
    person_by_dims(by, obs, counts, pool_dissolved))
  rows <- build_rows(dims)
  rows <- rows[rows$group_origin != rows$group_dest, , drop = FALSE]
  dest <- rows$.key3 - 1L
  axis <- if (duration) "dur" else "age"
  out <- tibble::as_tibble(rows[c("time_start", "time_end", paste0(axis, c("_start", "_end")),
                                  "group_origin", "group_dest", by)])
  out$person_years <- counts$person_years[rows$.index - dest * counts$cells]
  out$transitions <- counts$transitions[rows$.index]
  out$rate_transition <- as.numeric(ifelse(out$person_years > 0, out$transitions / out$person_years, NA_real_))
  add_rate_confidence(out, "transitions", "person_years", conf.level)
}
