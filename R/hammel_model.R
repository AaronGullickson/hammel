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
#' The summary combines the summaries of the starting population, of each
#' group dimension and its transition process, and of each other process in
#' the model, and records which elements required by [run_simulation()] are
#' missing.
#'
#' @param object A `hammel_model` object returned by [hammel_model()].
#' @param ... Additional arguments passed to [summary.hammel_population()].
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_model` containing:
#' \describe{
#'   \item{population}{The `summary.hammel_population` of the starting
#'     population, or `NULL` when no population is set.}
#'   \item{groups}{A named list with one element per group dimension, or
#'     `NULL` when the model has none. Each element contains `group` (the
#'     `summary.hammel_group` of the declared dimension, reporting the
#'     population's categories when a population is set, or `NULL` when the
#'     dimension is not declared) and `transition` (the
#'     `summary.hammel_group_transition` of the dimension's transition
#'     process, or `NULL` when it has none).}
#'   \item{processes}{A named list of process summaries, in the order
#'     mortality, fertility, union formation, union dissolution, for the
#'     processes present.}
#'   \item{missing}{A character vector naming the required elements that are
#'     not set: any of `"population"`, `"mortality"`, `"fertility"`, and
#'     `"union_formation"`, plus `"transition:<dimension>"` for each group
#'     transition process whose dimension is not declared
#'     (`"transition:.default"` for a process without `group` when the model
#'     does not declare exactly one dimension). Empty when the model has the
#'     elements needed to run a simulation.}
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
      groups = summarize_model_groups(object),
      processes = lapply(object$processes[processes], summary),
      missing = c(
        missing_model_elements(object),
        unresolved_transition_groups(object)
      )
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
  if (length(x$groups)) {
    print_model_groups(x$groups)
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
      print_missing_element(element)
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

# Per-dimension group information for summary.hammel_model(). Declared
# dimensions come first, in declaration order, followed by any dimension
# named only by a transition process.
summarize_model_groups <- function(model) {
  declaration <- model_groups(model)
  transitions <- transition_processes(model)
  dimensions <- names(declaration)
  keys <- names(transitions)
  # A process stored without a dimension applies to the only declared one.
  if (length(dimensions) == 1L) {
    keys[keys == ".default"] <- dimensions
  }
  names(transitions) <- keys
  dimensions <- union(dimensions, keys)
  if (!length(dimensions)) {
    return(NULL)
  }

  groups <- lapply(
    dimensions,
    function(dimension) {
      group <- declaration[[dimension]]
      categories <- model$population$groups[[dimension]]$categories
      if (is.null(categories)) {
        categories <- group$categories
      }
      transition <- transitions[[dimension]]
      list(
        group = if (is.null(group)) {
          NULL
        } else {
          summary(group, categories = categories)
        },
        transition = if (is.null(transition)) {
          NULL
        } else {
          summary(transition, categories = categories)
        }
      )
    }
  )
  names(groups) <- dimensions
  groups
}

# Print the Groups section of a model summary.
print_model_groups <- function(groups) {
  cli::cli_h2("Groups")
  for (dimension in names(groups)) {
    group <- groups[[dimension]]
    label <- if (dimension == ".default") "Unassigned dimension" else dimension
    cli::cli_h3(label)

    list <- open_summary_list()
    if (is.null(group$group)) {
      print_stat("Declared", "no")
    } else {
      print_group_details(group$group)
    }
    if (is.null(group$transition)) {
      print_stat("Transition process", "none")
    } else {
      transition <- open_nested_stat("Transition process")
      print_transition_details(group$transition)
      cli::cli_end(transition)
    }
    cli::cli_end(list)
  }
  invisible(groups)
}

# Transition process keys that run_simulation() cannot resolve to a declared
# group dimension, each prefixed "transition:". A key is the process's
# dimension name, or ".default" for a process stored without one.
unresolved_transition_groups <- function(model) {
  dimensions <- names(model_groups(model))
  keys <- names(transition_processes(model))
  unresolved <- keys[!keys %in% dimensions]
  if (length(dimensions) == 1L) {
    unresolved <- setdiff(unresolved, ".default")
  }
  if (!length(unresolved)) {
    return(character())
  }
  paste0("transition:", unresolved)
}

# Print one bullet describing an element returned by missing_model_elements()
# or unresolved_transition_groups().
print_missing_element <- function(element) {
  if (!startsWith(element, "transition:")) {
    cli::cli_li(required_model_elements[[element]], .auto_close = FALSE)
    return(invisible(element))
  }
  dimension <- sub("^transition:", "", element)
  if (dimension != ".default") {
    cli::cli_li(
      paste0(
        "A declaration of group dimension {.field {dimension}}, used by a ",
        "transition process, from {.fn set_groups}"
      ),
      .auto_close = FALSE
    )
  } else {
    cli::cli_li(
      paste0(
        "A single declared group dimension for the transition process ",
        "without {.arg group}, from {.fn set_groups}, or a {.arg group} ",
        "from {.fn set_group_transition}"
      ),
      .auto_close = FALSE
    )
  }
  invisible(element)
}

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
