# Create and validate model specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Create an empty Hammel model specification
#'
#' @description
#'
#' This function creates the empty `hammel_model` object which can then be built
#' out by adding populations and processes.
#'
#' @details
#'
#' The `hammel_model` is a container for a single `population` and a named list
#' of `processes`. These can all be added by the various setters through a
#' pipable interface or by feeding the model into [set_population()] or
#' [set_process()].
#'
#' @return A named list of class `hammel_model` containing `population`
#'   (initially `NULL`) and a named list of `processes`.
#' @export
#' @examples
#' model <- hammel_model() |>
#'   set_population(generate_uniform_population(1000)) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility)
hammel_model <- function() {
  structure(list(population = NULL, processes = list()), class = "hammel_model")
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a Hammel model
#'
#' @description
#'
#' Summary method for class `"hammel_model"`
#'
#' @details
#'
#' The summary combines the summaries of the starting population and of each
#' supported process in the model, and records which elements required by
#' [run_simulation()] are missing. Group transition processes are not yet
#' included in the summary.
#'
#' @param object A `hammel_model` object returned by [hammel_model()].
#' @param ... Additional arguments passed to [summary.hammel_population()].
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_model` containing:
#' \describe{
#'   \item{population}{The `summary.hammel_population` of the starting
#'     population, or `NULL` when no population is set.}
#'   \item{processes}{A named list of process summaries, in the order
#'     mortality, fertility, union formation, union dissolution, for the
#'     processes present.}
#'   \item{missing}{A character vector naming the required elements that are
#'     not set: any of `"population"`, `"mortality"`, `"fertility"`, and
#'     `"union_formation"`. Empty when the model can be simulated.}
#' }
#'
#' `print()` displays each component summary in turn, followed by whether the
#' model has the elements needed to run a simulation and, if not, which are
#' missing. It returns its input invisibly.
#'
#' @seealso [hammel_model()], [run_simulation()]
#' @examples
#' model <- hammel_model() |>
#'   set_population(generate_uniform_population(1000)) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility)
#' summary(model)
#'
#' model <- set_union_formation(model, "on_demand")
#' summary(model)
#' @export
summary.hammel_model <- function(object, ...) {
  population <- if (is.null(object$population)) {
    NULL
  } else {
    summary(object$population, ...)
  }
  processes <- intersect(summarized_processes, names(object$processes))

  structure(
    list(
      population = population,
      processes = lapply(object$processes[processes], summary),
      missing = missing_model_elements(object)
    ),
    class = "summary.hammel_model"
  )
}

#' @rdname summary.hammel_model
#' @param x A `summary.hammel_model` object returned by
#'   `summary.hammel_model()`.
#' @export
print.summary.hammel_model <- function(x, ...) {

  cli::cli_h1("Hammel model")

  if (!is.null(x$population)) {
    print(x$population)
  }
  for (process in x$processes) {
    print(process)
  }

  cli::cli_verbatim("")
  # Widen the gap after the status symbol, which some fonts render wide.
  status <- cli::cli_div(
    theme = list(
      ".alert-success" = list(
        before = function() {
          paste0(cli::col_green(cli::symbol$tick), "  ")
        }
      ),
      ".alert-danger" = list(
        before = function() {
          paste0(cli::col_red(cli::symbol$cross), "  ")
        }
      )
    )
  )
  if (!length(x$missing)) {
    cli::cli_alert_success(
      "The model has the elements needed to run a simulation."
    )
  } else {
    cli::cli_alert_danger(
      "The model needs the following before it can run a simulation:"
    )
    list <- open_summary_list()
    for (element in x$missing) {
      cli::cli_li(required_model_elements[[element]], .auto_close = FALSE)
    }
    cli::cli_end(list)
  }
  cli::cli_end(status)
  invisible(x)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Elements run_simulation() requires, with the setter that supplies each.
# Labels are cli format strings.
required_model_elements <- c(
  population = "A starting population, from {.fn set_population}",
  mortality = "A mortality process, from {.fn set_mortality}",
  fertility = "A fertility process, from {.fn set_fertility}",
  union_formation = "A union formation process, from {.fn set_union_formation}"
)

# Processes included in summary.hammel_model(), in display order.
summarized_processes <- c(
  "mortality", "fertility", "union_formation", "union_dissolution"
)

# Names of the required elements absent from a model.
missing_model_elements <- function(model) {
  present <- c(
    population = !is.null(model$population),
    mortality = !is.null(model$processes$mortality),
    fertility = !is.null(model$processes$fertility),
    union_formation = !is.null(model$processes$union_formation)
  )
  names(present)[!present]
}

check_model <- function(model) {
  if (!inherits(model, "hammel_model") || !is.list(model)) {
    stop("`model` must be a hammel_model object.", call. = FALSE)
  }
}
