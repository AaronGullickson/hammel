#' Create an empty Hammel model specification
#'
#' Set the starting population and demographic processes using pipeable setters.
#' Construction and configuration do not execute a simulation.
#'
#' @return A named list of class `hammel_model` containing `persons`, `unions`,
#'   `start_month`, and a named list of `processes`. Population and starting
#'   month are initially `NULL`.
#' @export
#' @examples
#' model <- hammel_model() |>
#'   set_population(data.frame(person_id = 1:2), data.frame()) |>
#'   set_mortality(data.frame(age = 0, rate = 0.01))
hammel_model <- function() {
  structure(
    list(persons = NULL, unions = NULL, start_month = NULL, processes = list()),
    class = "hammel_model"
  )
}

#' Set the starting population of a Hammel model
#'
#' Supply a person data frame with optional unions, or a `hammel_results` object.
#' Omitted or `NULL` unions represent no initial unions and are stored as an
#' empty data frame. Results supply
#' both histories and their endpoint month. Existing process configuration is
#' retained. The input model is not modified.
#'
#' @param model A `hammel_model` object.
#' @param persons A person data frame, or a `hammel_results` object.
#' @param unions A union data frame, or `NULL` (the default) for an unmarried
#'   initial population when `persons` is a data frame. Omit this argument when
#'   supplying results.
#' @param start_month Optional finite, whole-number starting month index for
#'   table inputs. Omit this argument when supplying results: their `end_month`
#'   is used. `NULL` means unknown; no starting month is inferred from dates.
#'
#' @return An updated `hammel_model`.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1:2)
#' unions <- data.frame()
#' model <- hammel_model() |>
#'   set_population(persons, unions, start_month = 0)
#' unmarried <- hammel_model() |>
#'   set_population(persons, start_month = 0)
#' results <- hammel_results(persons, unions, end_month = 120)
#' restarted <- hammel_model() |>
#'   set_population(results)
set_population <- function(model, persons, unions = NULL, start_month = NULL) {
  check_model(model)
  if (inherits(persons, "hammel_results")) {
    if (!missing(unions) || !missing(start_month)) {
      stop("When supplying results, omit `unions` and `start_month`.",
           call. = FALSE)
    }
    start_month <- persons$end_month
    unions <- persons$unions
    persons <- persons$persons
  } else if (is.null(unions)) {
    unions <- data.frame()
  }
  if (!is.data.frame(persons) || !is.data.frame(unions)) {
    stop("`persons` and `unions` must be data frames.", call. = FALSE)
  }
  check_month(start_month, "start_month")
  model$persons <- persons
  model$unions <- unions
  model["start_month"] <- list(start_month)
  model
}

check_model <- function(model) {
  if (!inherits(model, "hammel_model") || !is.list(model)) {
    stop("`model` must be a hammel_model object.", call. = FALSE)
  }
}

check_month <- function(month, argument) {
  if (!is.null(month) &&
      (!is.numeric(month) || length(month) != 1L ||
       !is.finite(month) || month != trunc(month))) {
    stop(sprintf("`%s` must be NULL or a finite, whole-number month index.",
                 argument), call. = FALSE)
  }
}
