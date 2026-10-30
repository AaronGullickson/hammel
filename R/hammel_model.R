#' Create an empty Hammel model specification
#'
#' Set the starting population and demographic processes using pipeable setters.
#' Construction and configuration do not execute a simulation.
#'
#' @return A named list of class `hammel_model` containing `population`
#'   (initially `NULL`) and a named list of `processes`. Once set, `population`
#'   is a `hammel_population` object whose `end_time` is the model's starting time.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1:2, sex = c("m", "f"),
#'   birth_time = c(-25, -23), death_time = NA_real_,
#'   mother_id = NA_integer_, father_id = NA_integer_)
#' model <- hammel_model() |>
#'   set_population(hammel_population(persons, end_time = 0)) |>
#'   set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
#'     age_end = Inf, rate = c(0.01, 0.008)))
hammel_model <- function() {
  structure(
    list(population = NULL, processes = list()),
    class = "hammel_model"
  )
}

#' Set the starting population of a Hammel model
#'
#' Supply a validated population containing person and union histories and
#' their reference time. Tables are normalized to tibbles; column contents and
#' additional components are retained. Existing process configuration is retained.
#' The input model is not modified.
#' Grouped populations without ancestry proportions initialize each starting
#' person at 100% ancestry in their recorded group. Supplied ancestry is
#' retained and normalized to a tibble in person order; see [hammel_population()].
#'
#' @param model A `hammel_model` object.
#' @param population A `hammel_population` object. Construct it from tables with
#'   [hammel_population()]. Its `end_time` supplies the model's starting time.
#'
#' @return An updated `hammel_model`.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1:2, sex = c("m", "f"),
#'   birth_time = c(-25, -23), death_time = NA_real_,
#'   mother_id = NA_integer_, father_id = NA_integer_)
#' population <- hammel_population(persons, end_time = 0)
#' model <- hammel_model() |>
#'   set_population(population)
set_population <- function(model, population) {
  check_model(model)
  validate_population(population)
  population <- normalize_population_transitions(population)
  population <- normalize_population_ancestry(population)
  population$persons <- tibble::as_tibble(population$persons)
  population$unions <- tibble::as_tibble(population$unions)
  model$population <- population
  model
}

check_model <- function(model) {
  if (!inherits(model, "hammel_model") || !is.list(model)) {
    stop("`model` must be a hammel_model object.", call. = FALSE)
  }
}
