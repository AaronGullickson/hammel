#' Create a Hammel population
#'
#' Collect person and union histories in an inspectable population object.
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
#' @return A named list of class `hammel_population` with components `persons`,
#'   `unions`, and `end_time`.
#' @export
#'
#' @examples
#' sim <- hammel_population(
#'   persons = data.frame(person_id = 1:2, sex = c(0L, 1L),
#'     birth_time = c(-25, -23), death_time = NA_real_,
#'     mother_id = NA_integer_, father_id = NA_integer_),
#'   end_time = 0
#' )
#' sim$persons
#' sim$unions
hammel_population <- function(persons, unions = NULL, end_time = NULL) {
  results <- structure(
    list(persons = persons, unions = normalize_unions(unions), end_time = end_time),
    class = "hammel_population"
  )
  validate_population(results)
  results
}

#' Validate a Hammel population
#'
#' Check the population structure, population schemas, IDs, relationship references,
#' and date consistency. This function does not normalize or change the object.
#' Invalid objects raise an error; valid objects are returned invisibly.
#'
#' @param results A `hammel_population` object.
#' @return The supplied object, invisibly and unchanged.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1L, sex = 1L, birth_time = -20,
#'   death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
#' results <- hammel_population(persons, end_time = 0)
#' validate_population(results)
validate_population <- function(results) {
  if (!inherits(results, "hammel_population") || !is.list(results) ||
      anyDuplicated(names(results)) ||
      !all(c("persons", "unions", "end_time") %in% names(results))) {
    stop("`results` must be a hammel_population object with persons, unions, and end_time components.",
         call. = FALSE)
  }
  persons <- results$persons
  unions <- results$unions
  time <- results$end_time
  check_time(time, "end_time")
  check_columns(persons, c("person_id", "sex", "birth_time", "death_time",
                          "mother_id", "father_id"), "persons")
  check_columns(unions, names(empty_unions()), "unions")
  check_numeric_column(persons$person_id, "persons$person_id", id = TRUE)
  if (anyDuplicated(persons$person_id)) {
    stop("`persons$person_id` must be unique.", call. = FALSE)
  }
  check_numeric_column(persons$sex, "persons$sex")
  if (any(!persons$sex %in% c(0, 1))) {
    stop("`persons$sex` must be 0 (male) or 1 (female).", call. = FALSE)
  }
  for (field in c("birth_time", "death_time")) {
    check_numeric_column(persons[[field]], paste0("persons$", field),
                         missing = field == "death_time")
  }
  if (any(persons$death_time < persons$birth_time, na.rm = TRUE)) {
    stop("Death times must not precede birth times.", call. = FALSE)
  }
  for (field in c("mother_id", "father_id")) {
    x <- persons[[field]]
    check_numeric_column(x, paste0("persons$", field), missing = TRUE, id = TRUE)
    if (any(!is.na(x) & (!x %in% persons$person_id | x == persons$person_id))) {
      stop(sprintf("`persons$%s` must refer to another person or be NA.", field),
           call. = FALSE)
    }
  }
  for (field in c("union_id", "partner1_id", "partner2_id")) {
    check_numeric_column(unions[[field]], paste0("unions$", field), id = TRUE)
  }
  if (anyDuplicated(unions$union_id)) {
    stop("`unions$union_id` must be unique.", call. = FALSE)
  }
  if (any(!unions$partner1_id %in% persons$person_id |
          !unions$partner2_id %in% persons$person_id |
          unions$partner1_id == unions$partner2_id)) {
    stop("Union partners must refer to two distinct persons.", call. = FALSE)
  }
  check_numeric_column(unions$start_time, "unions$start_time")
  check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
  if (any(unions$end_time < unions$start_time, na.rm = TRUE)) {
    stop("Union end times must not precede start times.", call. = FALSE)
  }
  reason <- unions$end_reason
  if (!is.character(reason) || !is.null(dim(reason)) ||
      any(!is.na(reason) & !reason %in% c("death", "dissolution")) ||
      any(is.na(reason) != is.na(unions$end_time))) {
    stop("`unions$end_reason` must be death or dissolution for ended unions, and NA for ongoing unions.",
         call. = FALSE)
  }
  for (partner in c("partner1_id", "partner2_id")) {
    born <- persons$birth_time[match(unions[[partner]], persons$person_id)]
    if (any(unions$start_time < born)) {
      stop("Unions must not start before either partner's birth.", call. = FALSE)
    }
  }
  if (!is.null(time) &&
      any(c(persons$birth_time, persons$death_time, unions$start_time,
            unions$end_time) > time, na.rm = TRUE)) {
    stop("Recorded events must not occur after the population's reference time.",
         call. = FALSE)
  }
  invisible(results)
}

normalize_unions <- function(unions) {
  if (is.null(unions) ||
      (is.data.frame(unions) && nrow(unions) == 0L && ncol(unions) == 0L)) {
    return(empty_unions())
  }
  unions
}

empty_unions <- function() {
  data.frame(
    union_id = integer(), partner1_id = integer(), partner2_id = integer(),
    start_time = numeric(), end_time = numeric(), end_reason = character()
  )
}

check_time <- function(time, argument) {
  if (!is.null(time) &&
      (!is.numeric(time) || length(time) != 1L || !is.finite(time))) {
    stop(sprintf("`%s` must be NULL or a finite time in years.", argument),
         call. = FALSE)
  }
}

check_columns <- function(table, columns, name) {
  if (!is.data.frame(table)) {
    stop(sprintf("`%s` must be a data frame.", name), call. = FALSE)
  }
  if (anyDuplicated(names(table))) {
    stop(sprintf("`%s` must have unique column names.", name), call. = FALSE)
  }
  absent <- setdiff(columns, names(table))
  if (length(absent)) {
    stop(sprintf("`%s` is missing columns: %s.", name,
                 paste(absent, collapse = ", ")), call. = FALSE)
  }
}

check_numeric_column <- function(x, name, missing = FALSE, id = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) ||
      any(is.nan(x)) || (!missing && anyNA(x)) ||
      any(!is.finite(x[!is.na(x)])) ||
      (id && any(x[!is.na(x)] <= 0 | x[!is.na(x)] != trunc(x[!is.na(x)])))) {
    stop(sprintf("`%s` must contain %s%s.", name,
                 if (id) "positive whole-number IDs" else "finite numeric values",
                 if (missing) " or NA" else " without missing values"),
         call. = FALSE)
  }
}
