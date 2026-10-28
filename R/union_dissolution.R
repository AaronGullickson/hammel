#' Specify duration-dependent union dissolution
#'
#' Construct or attach an optional dissolution process. One continuous hazard
#' per union is carried by the female partner, using years since the union
#' began, independently of either partner's age. Rates are per partnered
#' woman-year (equivalently union-year), not annual probabilities. Omitting
#' this process retains death-only union endings.
#'
#' Births preserve the pending dissolution date. Exact same-person ties select
#' death, then birth, then dissolution; a tied birth occurs while partnered.
#' Dissolution ends the episode, restores both partners' availability, and
#' reschedules the woman under her unpartnered fertility and existing search
#' rules, retaining birth spacing. Separate remarriage rates are not introduced.
#' Restarts resample at attained union duration rather than reset duration.
#'
#' @param rates A nonempty data frame with numeric `dur_start`, `dur_end`, and
#'   `rate`. Bands are half-open `[dur_start, dur_end)` in years. Starts and
#'   rates must be finite and nonnegative; ends must exceed starts and may be
#'   `Inf`. Bands cannot overlap within a group. Gaps and durations outside
#'   support have zero hazard. Sex, age, union-status and calendar-time columns
#'   are not supported. Input order, custom columns and provenance are retained.
#'   An optional `group` column selects the woman's group; shared rates apply
#'   to all groups. Group coverage and ungrouped fallback follow [mortality()].
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   rate multipliers by the woman's group. Missing groups use 1; zero suppresses
#'   risk. With no population groups, simulation warns and ignores multipliers.
#' @param model A `hammel_model` object.
#' @return `union_dissolution()` returns an inspectable
#'   `hammel_union_dissolution` process with a tibble of rates;
#'   `set_union_dissolution()` returns an updated model.
#' @export
#' @examples
#' # Illustrative rates, not demographic defaults.
#' rates <- data.frame(dur_start = c(0, 5), dur_end = c(5, Inf),
#'                     rate = c(0.04, 0.02))
#' model <- hammel_model() |> set_union_dissolution(rates)
union_dissolution <- function(rates, group_prop_hazard = NULL) {
  process <- new_rate_process("union_dissolution", rates)
  process$group_prop_hazard <- group_prop_hazard
  validate_union_dissolution(process)
  process
}

#' @rdname union_dissolution
#' @export
set_union_dissolution <- function(model, rates, group_prop_hazard = NULL) {
  set_process(model, union_dissolution(rates, group_prop_hazard))
}

validate_union_dissolution <- function(process) {
  if (!identical(process$name, "union_dissolution") ||
      !identical(process$rate_unit, "person_year"))
    stop("Invalid union-dissolution specification or units.", call. = FALSE)
  validate_group_prop_hazard(process$group_prop_hazard)
  rates <- process$rates
  check_columns(rates, c("dur_start", "dur_end", "rate"), "rates")
  check_numeric_column(rates$dur_start, "rates$dur_start")
  check_numeric_column(rates$rate, "rates$rate")
  ends <- rates$dur_end
  if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) ||
      anyNA(ends) || any(ends <= rates$dur_start) ||
      any(rates$dur_start < 0) || any(rates$rate < 0) ||
      any(c("sex", "age_start", "age_end", "union_status", "time_start", "time_end") %in% names(rates)))
    stop("Use nonnegative dissolution rates in valid duration bands, without sex, age, status or calendar columns.",
         call. = FALSE)
  for (rows in split_rate_groups(rates)) {
    ordered <- rates[rows, , drop = FALSE]
    ordered <- ordered[order(ordered$dur_start), , drop = FALSE]
    if (nrow(ordered) > 1L &&
        any(tail(ordered$dur_start, -1) < head(cummax(ordered$dur_end), -1)))
      stop("Union-dissolution duration bands must not overlap.", call. = FALSE)
  }
  invisible(process)
}

# Use the shared Band reader internally; public schedules retain duration names.
prepare_dissolution_rates <- function(process, persons, group_values) {
  if (is.null(process))
    return(data.frame(age_start = numeric(), age_end = numeric(), rate = numeric()))
  rates <- prepare_group_rates(process$rates, persons, "union dissolution", group_values)
  rates <- apply_group_prop_hazard(rates, persons, process$group_prop_hazard,
                                  "union dissolution", group_values)
  rates$age_start <- rates$dur_start
  rates$age_end <- rates$dur_end
  rates
}
