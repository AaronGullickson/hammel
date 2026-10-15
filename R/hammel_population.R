#' Create a Hammel population
#'
#' Collect person and union histories in an inspectable population object.
#' This constructor does not run a simulation. It validates
#' consistent tables against the population schema. Custom columns are retained.
#'
#' @param persons A data frame with one row per person, including deceased
#'   persons. Columns are preserved and the table is returned as a tibble.
#' @param unions A data frame with one row per union episode, or `NULL`
#'   for no unions. `wife_id` references the female partner and `husband_id`
#'   the male partner; these roles do not imply legal marriage. An empty data
#'   frame is normalized to the union schema.
#' @param end_time Optional finite simulation endpoint in years. Its calendar origin is supplied by the user. `NULL` means
#'   that the endpoint is unknown.
#'
#' @return A named list of class `hammel_population` with components `persons`,
#'   `unions`, and `end_time`.
#' @export
#'
#' @examples
#' sim <- hammel_population(
#'   persons = data.frame(person_id = 1:2, sex = c("m", "f"),
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
  results$persons <- tibble::as_tibble(results$persons)
  results$unions <- tibble::as_tibble(results$unions)
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
#' persons <- data.frame(person_id = 1L, sex = "f", birth_time = -20,
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
  if (!is.character(persons$sex) || !is.null(dim(persons$sex)) ||
      anyNA(persons$sex) || any(!persons$sex %in% c("m", "f"))) {
    stop('`persons$sex` must contain "m" (male) or "f" (female).', call. = FALSE)
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
  for (field in c("union_id", "wife_id", "husband_id")) {
    check_numeric_column(unions[[field]], paste0("unions$", field), id = TRUE)
  }
  if (anyDuplicated(unions$union_id)) {
    stop("`unions$union_id` must be unique.", call. = FALSE)
  }
  if (any(!unions$wife_id %in% persons$person_id |
          !unions$husband_id %in% persons$person_id |
          unions$wife_id == unions$husband_id)) {
    stop("Union partners must refer to two distinct persons.", call. = FALSE)
  }
  wife_sex <- persons$sex[match(unions$wife_id, persons$person_id)]
  husband_sex <- persons$sex[match(unions$husband_id, persons$person_id)]
  if (any(wife_sex != "f" | husband_sex != "m")) {
    stop('`wife_id` must reference a female person and `husband_id` a male person.',
         call. = FALSE)
  }
  check_numeric_column(unions$start_time, "unions$start_time")
  check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
  if (any(unions$end_time < unions$start_time, na.rm = TRUE)) {
    stop("Union end times must not precede start times.", call. = FALSE)
  }
  reason <- unions$end_reason
  if (!is.character(reason) || !is.null(dim(reason)) ||
      any(!is.na(reason) & !reason %in% c("death_husband", "death_wife", "dissolution")) ||
      any(is.na(reason) != is.na(unions$end_time))) {
    stop("`unions$end_reason` must be death_husband, death_wife, or dissolution for ended unions, and NA for ongoing unions.",
         call. = FALSE)
  }
  for (partner in c("wife_id", "husband_id")) {
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
  tibble::tibble(
    union_id = integer(), wife_id = integer(), husband_id = integer(),
    start_time = numeric(), end_time = numeric(), end_reason = character()
  )
}

#' Generate a simple baseline population
#'
#' Create a starter population with independent uniformly distributed ages,
#' randomly assigned sex, unknown parents, and no unions. Everyone is alive
#' at the reference time. This is not a stable or stationary population.
#'
#' Ages are sampled in continuous years without rounding to simulation steps.
#' Generation uses and advances R's current random-number state; use
#' [set.seed()] for reproducibility.
#'
#' @param size Number of persons, a nonnegative whole number.
#' @param max_age Maximum initial age in years, finite and nonnegative.
#' @param end_time Reference time in years, finite. Birth times are this value
#'   minus the sampled ages. Defaults to a relative time origin of zero.
#' @param female_probability Probability of assigning sex "f" (female).
#'   The other sex value is "m" (male). This does not enforce exact proportions.
#'
#' @return A validated `hammel_population` object.
#' @export
#' @examples
#' set.seed(42)
#' population <- create_baseline_population(size = 50000, max_age = 70)
#' model <- hammel_model() |>
#'   set_population(population)
create_baseline_population <- function(size, max_age = 70, end_time = 0,
                                female_probability = 0.5) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
      size < 0 || size != trunc(size) || size > .Machine$integer.max) {
    stop("`size` must be a nonnegative whole number no larger than R's integer limit.",
         call. = FALSE)
  }
  if (!is.numeric(max_age) || length(max_age) != 1L ||
      !is.finite(max_age) || max_age < 0) {
    stop("`max_age` must be finite and nonnegative, in years.", call. = FALSE)
  }
  check_time(end_time, "end_time")
  if (is.null(end_time)) {
    stop("`end_time` must be a finite time in years.", call. = FALSE)
  }
  if (!is.numeric(female_probability) || length(female_probability) != 1L ||
      !is.finite(female_probability) || female_probability < 0 ||
      female_probability > 1) {
    stop("`female_probability` must be a number between 0 and 1.", call. = FALSE)
  }

  persons <- data.frame(
    person_id = seq_len(size),
    sex = as.character(ifelse(stats::rbinom(size, size = 1, prob = female_probability) == 1,
                 "f", "m")),
    birth_time = end_time - stats::runif(size, min = 0, max = max_age),
    death_time = rep(NA_real_, size),
    mother_id = rep(NA_integer_, size),
    father_id = rep(NA_integer_, size)
  )
  hammel_population(persons, end_time = end_time)
}
