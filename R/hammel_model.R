#' Create an empty Hammel model specification
#'
#' Set the starting population and demographic processes using pipeable setters.
#' Construction and configuration do not execute a simulation.
#'
#' @return A named list of class `hammel_model` containing `persons`, `unions`,
#'   `start_time`, and a named list of `processes`. Population and starting
#'   time are initially `NULL`.
#' @export
#' @examples
#' model <- hammel_model() |>
#'   set_population(data.frame(person_id = 1:2, sex = c(0L, 1L),
#'     birth_time = c(-25, -23), death_time = NA_real_,
#'     mother_id = NA_integer_, father_id = NA_integer_)) |>
#'   set_mortality(data.frame(age = 0, rate = 0.01))
hammel_model <- function() {
  structure(
    list(persons = NULL, unions = NULL, start_time = NULL, processes = list()),
    class = "hammel_model"
  )
}

#' Set the starting population of a Hammel model
#'
#' Supply a person data frame with optional unions, or a `hammel_results` object.
#' Omitted or `NULL` unions represent no initial unions and are stored as an
#' empty data frame. Results supply
#' both histories and their endpoint time. Existing process configuration is
#' retained. The input model is not modified.
#'
#' @param model A `hammel_model` object.
#' @param persons A person data frame, or a `hammel_results` object.
#' @param unions A union data frame, or `NULL` (the default) for an unmarried
#'   initial population when `persons` is a data frame. Omit this argument when
#'   supplying results.
#' @param start_time Optional finite starting time in years for
#'   table inputs. Omit this argument when supplying results: their `end_time`
#'   is used. `NULL` means unknown; no starting time is inferred from dates.
#'
#' @return An updated `hammel_model`.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1:2, sex = c(0L, 1L),
#'   birth_time = c(-25, -23), death_time = NA_real_,
#'   mother_id = NA_integer_, father_id = NA_integer_)
#' unions <- NULL
#' model <- hammel_model() |>
#'   set_population(persons, unions, start_time = 0)
#' unmarried <- hammel_model() |>
#'   set_population(persons, start_time = 0)
#' results <- hammel_results(persons, unions, end_time = 120)
#' restarted <- hammel_model() |>
#'   set_population(results)
set_population <- function(model, persons, unions = NULL, start_time = NULL) {
  check_model(model)
  if (inherits(persons, "hammel_results")) {
    if (!missing(unions) || !missing(start_time)) {
      stop("When supplying results, omit `unions` and `start_time`.",
           call. = FALSE)
    }
    start_time <- persons$end_time
    unions <- persons$unions
    persons <- persons$persons
  } else if (is.null(unions)) {
    unions <- empty_unions()
  }
  check_time(start_time, "start_time")
  unions <- validate_population(persons, unions, start_time)
  model$persons <- persons
  model$unions <- unions
  model["start_time"] <- list(start_time)
  model
}

check_model <- function(model) {
  if (!inherits(model, "hammel_model") || !is.list(model)) {
    stop("`model` must be a hammel_model object.", call. = FALSE)
  }
}
