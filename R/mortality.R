# Construct mortality processes and attach them to model specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify an age-specific mortality process
#'
#' `mortality()` constructs a process specification without running a simulation.
#' `set_mortality()` constructs the same specification and attaches it to a
#' model, equivalent to `set_process(model, mortality(rates, weights))`.
#' The input model is not modified and neither function runs a simulation.
#'
#' Rates are piecewise constant annual hazards over half-open age intervals.
#' Supply separate schedules for both sexes. Overlapping intervals are rejected;
#' adjacent intervals are allowed. Gaps and uncovered ages are not filled or
#' interpreted as zero rates. Simulation requires coverage of [0, Inf).
#' Input order and custom columns are preserved, rate tables become tibbles,
#' and a `source` attribute on the input is copied into the process.
#'
#' Columns named by declared group dimensions (see [set_groups()]), such as
#' `race` or `region`, supply separate schedules for each combination of their
#' categories. Dimensions without a column share the schedule. Each stratum
#' needs both sexes; overlap is checked within strata. Simulation requires
#' schedules for every combination of categories of the dimensions a table
#' uses. A table using an undeclared persons column is an error.
#'
#' @section Calendar rates:
#' Optional `time_start`/`time_end` columns select half-open calendar intervals
#' in absolute years, including fractional years. Both columns are required
#' together; starts are finite and ends exceed starts (possibly Inf). Active
#' slices must satisfy the strata and coverage rules and cover the requested
#' run. Rates are constant within intervals, with no interpolation. Calendar
#' changes precede same-time events and preserve recorded histories.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`, and character `sex` values `"m"` and `"f"`, with both represented.
#'   Starts must be finite and nonnegative; ends must exceed starts and may be
#'   `Inf` for upper-open bands. Rates must be finite and nonnegative, per
#'   person-year. Use [stratify_rates()] to add proportional group differences.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate.
#' @param model A `hammel_model` object, used by `set_mortality()`.
#' @return `mortality()` returns an inspectable process inheriting from
#'   `hammel_mortality` and `hammel_process`, with `name`, `rates`, `rate_unit`,
#'   and `weights`. `set_mortality()` returns an updated `hammel_model`.
#' @seealso [set_process()], [fertility()], [union_formation()],
#'   [union_dissolution()], [group_transition()]
#' @export
#' @examples
#' rates <- data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008))
#' process <- mortality(rates)
#' model <- hammel_model() |> set_process(process)
#' # Equivalent convenience setter:
#' model <- hammel_model() |> set_mortality(rates)
#' @order 1
mortality <- function(rates, weights = NULL) {
  validate_rates(rates, mortality = TRUE)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(weights, weight_strata("mortality", rates), "mortality")
  process <- new_rate_process("mortality", rates)
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname mortality
#' @export
#' @order 2
set_mortality <- function(model, rates, weights = NULL) {
  set_process(model, mortality(rates, weights = weights))
}
