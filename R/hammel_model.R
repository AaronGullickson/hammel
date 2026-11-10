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

# HELPER FUNCTIONS ---------------------------------------------------------

check_model <- function(model) {
  if (!inherits(model, "hammel_model") || !is.list(model)) {
    stop("`model` must be a hammel_model object.", call. = FALSE)
  }
}
