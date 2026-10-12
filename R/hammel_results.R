#' Create a Hammel simulation result
#'
#' Collect person and union histories in an inspectable simulation result.
#' This constructor does not run a simulation. It checks that both inputs are
#' data frames but does not yet validate their columns or demographic content.
#'
#' @param persons A data frame with one row per person, including deceased
#'   persons. Data frame subclasses, such as tibbles, are preserved.
#' @param unions A data frame with one row per union episode.
#' @param end_month Optional simulation endpoint as a finite, whole-number
#'   month index. Its calendar origin is supplied by the user. `NULL` means
#'   that the endpoint is unknown.
#'
#' @return A named list of class `hammel_results` with components `persons`,
#'   `unions`, and `end_month`.
#' @export
#'
#' @examples
#' sim <- hammel_results(
#'   persons = data.frame(person_id = c(1L, 2L)),
#'   unions = data.frame(
#'     union_id = 1L,
#'     person_id_1 = 1L,
#'     person_id_2 = 2L
#'   )
#' )
#' sim$persons
#' sim$unions
hammel_results <- function(persons, unions, end_month = NULL) {
  if (!is.data.frame(persons)) {
    stop("`persons` must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(unions)) {
    stop("`unions` must be a data frame.", call. = FALSE)
  }
  check_month(end_month, "end_month")

  structure(
    list(persons = persons, unions = unions, end_month = end_month),
    class = "hammel_results"
  )
}
