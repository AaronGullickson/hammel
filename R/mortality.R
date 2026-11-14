# Construct mortality processes and attach them to model specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify an age-specific mortality process
#'
#' @description
#'
#' Constructs a mortality process specification for a simulation model. The
#' convenience function `set_mortality` is equivalent to `set_process(model,
#' mortality(...))` and allows for direct piping from a model.
#'
#' @details
#'
#' The primary and only required input to the function is a mortality rate
#' schedule for the `rates` argument. This data frame must contain `age_start`
#' and `age_end` columns to define the half-open age intervals (inclusive at the
#' start and exclusive at the end), a `sex` column ("f" for female and "m" for
#' male) and a `rate` column for the age-specific fertility rates. The example
#' below shows the basic format:
#'
#' ```{r, echo = FALSE}
#' usa_2000_mortality
#' ```
#'
#' Rates are always specified as deaths over person-years lived. An open-ended
#' age group specified with `age_end` of `Inf` must be included. No gap in the
#' age band is allowed.
#'
#' Columns named by declared group dimensions (see [set_groups()]), such as
#' `race` or `region`, supply separate schedules for each combination of their
#' categories; dimensions without a column share the schedule. Every combination
#' of categories of the dimensions must be included in the rates if the
#' dimension is used.
#'
#' Optional `time_start`/`time_end` columns select half-open calendar intervals
#' in absolute years, including fractional years. Both columns are required
#' together; starts are finite and ends exceed starts (possibly Inf). Active
#' slices must satisfy the strata and coverage rules and cover the requested
#' simulation run.
#'
#' The [stratify_rates()] function provides a quick means to add stratification
#' to an existing rate schedule.
#'
#' @param rates A nonempty data frame defining the mortality rate schedules.
#'   Must contain at least numeric `age_start`, `age_end`, and `rate`, as well
#'   as character `sex` columns. See details for format and description of
#'   stratification variables that can be applied.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [weight_model()] specification. Weights choose which member of a cell
#'   experiences an event without changing the cell's input rate.
#' @param model A `hammel_model` object, used by `set_mortality()`.
#' @return `mortality()` returns an inspectable process inheriting from
#'   `hammel_mortality` and `hammel_process`, with `name`, `rates`,
#'   and `weights`. `set_mortality()` returns an updated `hammel_model`.
#' @seealso [set_process()], [fertility()], [union_formation()],
#'   [union_dissolution()], [group_transition()]
#' @export
#' @examples
#' mort_process <- mortality(usa_2000_mortality)
#' model <- hammel_model() |> set_process(mort_process)
#' # Equivalent convenience setter:
#' model <- hammel_model() |> set_mortality(usa_2000_mortality)
#' @order 1
mortality <- function(rates, weights = NULL) {
  validate_rates(rates, mortality = TRUE)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(
    weights, weight_strata("mortality", rates), "mortality"
  )
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


# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a mortality process
#'
#' @description
#' 
#' Summary method for class `"hammel_mortality"`
#'
#' @param object A `hammel_mortality` object returned by [mortality()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_morality` containing:
#' \describe{
#'   \item{rate_strata}{A list describing the input rate schedule, containing
#'     `stratified_by` (stratum column names), `strata` (a tibble of observed
#'     stratum combinations), `bands` (a tibble of exact band combinations),
#'     and `bands_by_stratum` (whether band combinations differ by stratum).
#'     Shared bands are stored once; otherwise `bands` includes the stratum
#'     columns. Band endpoints define half-open intervals `[start, end)`.}
#' }
#'
#' @seealso [mortality()], [calculate_life_table()]
#' @examples
#' process <- mortality(usa_2000_mortality)
#' summary(process)
#' @export
summary.hammel_mortality <- function(object, ...) {

  structure(
    list(
      rate_strata = summarize_rate_strata(object$rates)
    ),
    class = "summary.hammel_mortality"
  )
}

#' @rdname summary.hammel_fertility
#' @param x A `summary.hammel_fertility` object returned by
#'   `summary.hammel_fertility()`.
#' @export
print.summary.hammel_mortality <- function(x, ...) {

  cat(cli::style_bold("Mortality Process Summary", "\n\n"))

  print_rate_strata(x$rate_strata)
}
