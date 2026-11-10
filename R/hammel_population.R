# Construct, validate, and attach populations to model specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Create a Hammel population
#'
#' @description
#'
#' Create a `hammel_population` object to serve as the starting population for a
#' simulation.
#'
#' @details
#'
#' In most cases, it will be more common to create a `hammel_population` object
#' by using [generate_uniform_population()] or to recieve a `hammel_population`
#' object as the output of [run_simulation()].
#'
#' The `hammel_population` object must always contain a `persons` data frame and
#' an `end_time` value. The `end_time` value defines the starting point for the
#' simulation and no dates in the `hammel_population` may exceed the `end_time`
#' value.
#'
#' @param persons A data frame with one row per person, including deceased
#'   persons. Must include a unique `person_id`, nonmissing `birth_time`,
#'   `death_time`, `mother_id`, and `father_id`.
#' @param unions A data frame with one row per union episode, or `NULL` for no
#'   unions. `wife_id` references the female partner and `husband_id` the male
#'   partner. Optional `union_type` is cohab or marriage.
#' @param end_time Required finite numeric reference time in years. Its calendar
#'   origin is supplied by the user. This supplies the next simulation's
#'   starting time. Recorded events must not occur after this time.
#' @param groups Optional group declaration: a character vector of dimension
#'   names (each with [group()] defaults) or a named list of [group()] objects.
#'   Persons columns that are not declared are ordinary custom columns. Models
#'   usually declare groups with [set_groups()].
#' @param transitions Optional data frame of individual category changes with
#'   `person_id`, `time`, `group` (the dimension; may be omitted when one is
#'   declared), `origin`, and `destination`. Each person's history in each
#'   dimension must chain to their current category and lie within life dates
#'   and the known endpoint. Grouped populations always include this tibble,
#'   even when empty, with character `group`, `origin` and `destination`.
#' @param ancestry Optional ancestry shares for dimensions declared with
#'   `group(ancestry = TRUE)`: a named list with one data frame per such group
#'   dimension (a bare data frame when only one dimension tracks ancestry), each
#'   with exactly one row per person, keyed by `person_id`, one column of shares
#'   per category (named by the category) and an `unknown` share. Shares must be
#'   finite, between zero and one, and sum to one per row (tolerance `1e-8`).
#'
#' @return A named list of class `hammel_population` with components `persons`,
#'   `unions`, and `end_time`, plus an `ancestry` tibble when supplied or when
#'   persons has groups, and a `transitions` tibble when groups are present.
#' @export
#' @examples
#' population <- hammel_population(
#'   persons = data.frame(
#'     person_id = 1:2, sex = c("m", "f"), birth_time = c(-25, -23),
#'     death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
#'   ),
#'   end_time = 0
#' )
hammel_population <- function(
  persons,
  unions = NULL,
  end_time,
  ancestry = NULL,
  transitions = NULL,
  groups = NULL
) {
  if (missing(end_time)) {
    stop("`end_time` is required, in years.", call. = FALSE)
  }
  check_time(end_time, "end_time")
  results <- structure(
    list(
      persons = persons,
      unions = normalize_unions(unions),
      end_time = end_time
    ),
    class = "hammel_population"
  )
  results$groups <- normalize_group_declaration(groups)
  if (!is.null(transitions)) {
    results$transitions <- transitions
  }
  validate_group_declaration(results$groups, results)
  if (!is.null(transitions)) {
    results$transitions <- transitions
  }
  results <- normalize_population_transitions(results)
  if (!is.null(ancestry)) {
    results$ancestry <- ancestry
  }
  validate_population(results)
  results <- normalize_population_ancestry(results)
  results$persons <- tibble::as_tibble(results$persons)
  results$unions <- tibble::as_tibble(results$unions)
  results
}

#' Set the starting population of a Hammel model
#'
#' @description
#'
#' Set the starting population for a `hammel_model` with a `hammel_population`
#' object.
#'
#' @details
#'
#' While a `hammel_population` object can be built from scratch with
#' [hammel_population()], this object will more commonly be created by
#' [generate_uniform_population()] or returned from [run_simulation()].
#'
#'
#' @param model A `hammel_model` object.
#' @param population A `hammel_population` object. Its `end_time` supplies the
#'   model's starting time.
#' @return An updated `hammel_model`.
#' @export
#' @examples
#' persons <- data.frame(
#'   person_id = 1:2, sex = c("m", "f"), birth_time = c(-25, -23),
#'   death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
#' )
#' population <- hammel_population(persons, end_time = 0)
#' model <- hammel_model() |>
#'   set_population(population)
set_population <- function(model, population) {
  check_model(model)
  validate_population(population)
  # A declaration made with set_groups() replaces the population's own.
  if (!is.null(model$groups)) {
    population <- apply_group_declaration(population, model$groups)
  }
  population <- normalize_population_transitions(population)
  population <- normalize_population_ancestry(population)
  population$persons <- tibble::as_tibble(population$persons)
  population$unions <- tibble::as_tibble(population$unions)
  model$population <- population
  model
}

#' Validate a Hammel population
#'
#' @description
#'
#' Check the population structure, population schemas, IDs, relationship
#' references, and date consistency. This function does not normalize or change
#' the object. Invalid objects raise an error; valid objects are returned
#' invisibly.
#'
#' @param results A `hammel_population` object.
#' @return The supplied object, invisibly and unchanged.
#' @export
validate_population <- function(results) {
  if (!inherits(results, "hammel_population") || !is.list(results) ||
        anyDuplicated(names(results)) ||
        !all(c("persons", "unions", "end_time") %in% names(results))) {
    stop(
      paste0(
        "`results` must be a hammel_population object with ",
        "persons, unions, and end_time components."
      ),
      call. = FALSE
    )
  }
  persons <- results$persons
  unions <- results$unions
  time <- results$end_time
  check_time(time, "end_time")
  check_columns(
    persons, c(
      "person_id", "sex", "birth_time", "death_time",
      "mother_id", "father_id"
    ), "persons"
  )
  validate_group_declaration(results$groups, results)
  validate_transition_history(results)
  check_columns(unions, names(empty_unions()), "unions")
  check_numeric_column(persons$person_id, "persons$person_id", id = TRUE)
  if (anyDuplicated(persons$person_id)) {
    stop("`persons$person_id` must be unique.", call. = FALSE)
  }
  if (!is.character(persons$sex) || !is.null(dim(persons$sex)) ||
        anyNA(persons$sex) ||
        any(!persons$sex %in% c("m", "f"))) {
    stop(
      '`persons$sex` must contain "m" (male) or "f" (female).',
      call. = FALSE
    )
  }
  for (field in c("birth_time", "death_time")) {
    check_numeric_column(
      persons[[field]], paste0("persons$", field),
      missing = field == "death_time"
    )
  }
  if (any(persons$death_time < persons$birth_time, na.rm = TRUE)) {
    stop("Death times must not precede birth times.", call. = FALSE)
  }
  for (field in c("mother_id", "father_id")) {
    x <- persons[[field]]
    check_numeric_column(
      x, paste0("persons$", field),
      missing = TRUE,
      id = TRUE
    )
    if (any(!is.na(x) & (!x %in% persons$person_id | x == persons$person_id))) {
      stop(
        sprintf("`persons$%s` must refer to another person or be NA.", field),
        call. = FALSE
      )
    }
  }
  for (field in c("union_id", "wife_id", "husband_id")) {
    check_numeric_column(unions[[field]], paste0("unions$", field), id = TRUE)
  }
  if (anyDuplicated(unions$union_id)) {
    stop("`unions$union_id` must be unique.", call. = FALSE)
  }
  if (any(
        !unions$wife_id %in% persons$person_id |
          !unions$husband_id %in% persons$person_id |
          unions$wife_id == unions$husband_id
      )) {
    stop("Union partners must refer to two distinct persons.", call. = FALSE)
  }
  wife_sex <- persons$sex[match(unions$wife_id, persons$person_id)]
  husband_sex <- persons$sex[match(unions$husband_id, persons$person_id)]
  if (any(wife_sex != "f" | husband_sex != "m")) {
    stop(
      paste0(
        "`wife_id` must reference a female person and ",
        "`husband_id` a male person."
      ),
      call. = FALSE
    )
  }
  check_numeric_column(unions$start_time, "unions$start_time")
  check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
  if (any(unions$end_time < unions$start_time, na.rm = TRUE)) {
    stop("Union end times must not precede start times.", call. = FALSE)
  }
  reason <- unions$end_reason
  if (!is.character(reason) || !is.null(dim(reason)) ||
        any(
          !is.na(reason) &
            !reason %in%
               c("death_husband", "death_wife", "dissolution", "conversion")
        ) ||
        any(is.na(reason) != is.na(unions$end_time))) {
    stop(
      paste0(
        "`unions$end_reason` must be death_husband, ",
        "death_wife, dissolution, or conversion for ended ",
        "unions, and NA for ongoing unions."
      ),
      call. = FALSE
    )
  }
  validate_union_types(unions)
  if (any(reason == "conversion", na.rm = TRUE) &&
        !"union_type" %in% names(unions)) {
    stop(
      "Conversion records require union_type and conversion links.",
      call. = FALSE
    )
  }
  for (partner in c("wife_id", "husband_id")) {
    born <- persons$birth_time[match(unions[[partner]], persons$person_id)]
    if (any(unions$start_time < born)) {
      stop(
        "Unions must not start before either partner's birth.",
        call. = FALSE
      )
    }
  }
  if (any(
        c(
          persons$birth_time,
          persons$death_time,
          unions$start_time,
          unions$end_time
        ) > time,
        na.rm = TRUE
      )) {
    stop(
      paste0(
        "Recorded events must not occur after the population's ",
        "reference time."
      ),
      call. = FALSE
    )
  }
  validate_ancestry(results[["ancestry"]], persons, results$groups)
  invisible(results)
}

# HELPER FUNCTIONS ---------------------------------------------------------

normalize_unions <- function(unions) {
  if (is.null(unions) ||
        (is.data.frame(unions) && nrow(unions) == 0L && ncol(unions) == 0L)) {
    return(empty_unions())
  }
  unions
}

empty_unions <- function() {
  tibble::tibble(
    union_id = integer(),
    wife_id = integer(),
    husband_id = integer(),
    start_time = numeric(),
    end_time = numeric(),
    end_reason = character()
  )
}

validate_union_types <- function(unions) {
  if (!"union_type" %in% names(unions)) {
    return(invisible(unions))
  }
  x <- unions$union_type
  if (!is.character(x) || !is.null(dim(x)) || anyNA(x) ||
        any(!x %in% c("cohab", "marriage"))) {
    stop(
      "`unions$union_type` must contain 'cohab' or 'marriage'.",
      call. = FALSE
    )
  }
  if (!"previous_cohab_id" %in% names(unions)) {
    if (any(unions$end_reason == "conversion", na.rm = TRUE)) {
      stop(
        "Conversion histories require `previous_cohab_id` links.",
        call. = FALSE
      )
    }
    return(invisible(unions))
  }
  links <- unions$previous_cohab_id
  check_numeric_column(links, "previous_cohab_id", missing = TRUE, id = TRUE)
  rows <- which(!is.na(links))
  parent <- match(links[rows], unions$union_id)
  if (anyNA(parent) || anyDuplicated(parent) ||
        any(unions$union_type[rows] != "marriage") ||
        any(unions$union_type[parent] != "cohab") ||
        any(
          is.na(unions$end_reason[parent]) |
            unions$end_reason[parent] != "conversion"
        ) ||
        any(
          is.na(unions$end_time[parent]) |
            unions$end_time[parent] != unions$start_time[rows]
        ) ||
        any(
          unions$wife_id[parent] != unions$wife_id[rows] |
            unions$husband_id[parent] != unions$husband_id[rows]
        ) ||
        any(!which(unions$end_reason == "conversion") %in% parent)) {
    stop(
      paste0(
        "Conversions require unique same-partner ",
        "cohab-to-marriage links at the same time."
      ),
      call. = FALSE
    )
  }
  invisible(unions)
}
