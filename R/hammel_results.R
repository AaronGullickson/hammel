#' Create a Hammel simulation result
#'
#' Collect person and union histories in an inspectable simulation result.
#' This constructor does not run a simulation. It validates
#' consistent tables against the population schema. Custom columns are retained.
#'
#' @param persons A data frame with one row per person, including deceased
#'   persons. Data frame subclasses, such as tibbles, are preserved.
#' @param unions A data frame with one row per union episode, or `NULL`
#'   for no unions. An empty data frame is normalized to the union schema.
#' @param end_time Optional finite simulation endpoint in years. Its calendar origin is supplied by the user. `NULL` means
#'   that the endpoint is unknown.
#'
#' @return A named list of class `hammel_results` with components `persons`,
#'   `unions`, and `end_time`.
#' @export
#'
#' @examples
#' sim <- hammel_results(
#'   persons = data.frame(person_id = 1:2, sex = c(0L, 1L),
#'     birth_time = c(-25, -23), death_time = NA_real_,
#'     mother_id = NA_integer_, father_id = NA_integer_),
#'   end_time = 0
#' )
#' sim$persons
#' sim$unions
hammel_results <- function(persons, unions = NULL, end_time = NULL) {
  check_time(end_time, "end_time")
  unions <- validate_population(persons, unions, end_time)

  structure(
    list(persons = persons, unions = unions, end_time = end_time),
    class = "hammel_results"
  )
}
