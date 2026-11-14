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
#' by using [generate_uniform_population()] or to receive a `hammel_population`
#' object as the output of [run_simulation()].
#'
#' The `hammel_population` object must always contain a `persons` data frame and
#' an `end_time` value. The `end_time` value defines the starting point for the
#' simulation and no dates in the `hammel_population` may exceed the `end_time`
#' value.
#'
#' @param persons A data frame with one row per person, including deceased
#'   persons. Must include a unique `person_id`, nonmissing `sex`, nonmissing
#'   `birth_time`, `death_time`, `mother_id`, and `father_id`.
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

# SUMMARY AND PRINT FUNCTIONS -------------------------------------------------

#' Summarize a Hammel population
#'
#' @description
#'
#' Summary method for class `"hammel_population"`
#'
#' @param object an object of class `"hammel_population"`, usually the result
#'   of a call to [run_simulation].
#' @param interval_pyramid length of age intervals used for drawing
#'   population pyramid
#' @param ... further arguments passed on from other methods.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_population` containing:
#' \describe{
#'   \item{first_time}{The earliest recorded birth time, in years.}
#'   \item{end_time}{The population's reference time, in years.}
#'   \item{persons_alive}{The number of persons alive at the reference time.}
#'   \item{n_births}{The number of recorded persons with a birth time,
#'     including persons in the starting population.}
#'   \item{n_births_u}{The number of retained unrealized birth records, or
#'     `NA` when no unrealized birth table is present.}
#'   \item{n_deaths}{The number of recorded persons with a death time.}
#'   \item{n_unions}{The number of recorded union episodes, including
#'     ended unions.}
#'   \item{n_dissolutions}{The number of union episodes with end reason
#'     `"dissolution"`.}
#'   \item{age_counts}{A tibble of living male and female counts by five-year
#'     age group, with columns `age`, `male`, and `female`.}
#'   \item{simulation}{A list describing the latest simulation run, with
#'     `start` and `end` times in years and living population counts
#'     `pop_size_start` and `pop_size_end`, or `NULL` when no simulation
#'     metadata is present.}
#'   \item{group_counts}{A named list of tibbles, one per declared group
#'     dimension, or `NULL` when no groups are declared. Each tibble contains
#'     `category`, `count`, and `percent` for the living population at the
#'     reference time, including declared categories with zero members.
#'     Percentages are on a 0--100 scale and are `NaN` when the dimension's
#'     total living count is zero.}
#' }
#'
#' `print()` displays the summary, group tables when present, and a text
#' population pyramid, returning its input invisibly. Pyramid bars share a
#' common scale across both sexes.
#'
#' @examples
#'
#' # get summary of starting population
#' population <- generate_uniform_population(1000)
#' summary(population)
#'
#' # get summary of simulatied population
#' model <- hammel_model() |>
#'   set_population(population) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_union_formation("on_demand")
#' result <- run_simulation(model, years = 50, seed = 42)
#' summary(result)
#' @export
summary.hammel_population <- function(object, interval_pyramid = 5, ...) {

  # get stats for drawing population pyramid
  alive <- is.na(object$persons$death_time)
  ages <- object$end_time - object$persons$birth_time[alive]

  breaks <- seq(0, max(c(5, ages)) + 5, by = 5)
  age_group <- cut(ages, breaks, right = FALSE)

  age_counts <- table(
    age_group,
    factor(object$persons$sex[alive], levels = c("m", "f"))
  )

  age_counts = tibble::tibble(
    age = paste0(
      head(breaks, -1), "\u002d", tail(breaks, -1) - 1
    ),
    male = as.integer(age_counts[, "m"]),
    female = as.integer(age_counts[, "f"])
  )

  if (!is.null(object$simulation)) {
    start <- object$simulation$start_time
    pop_init <- sum(
      object$persons$birth_time <= start &
        (is.na(object$persons$death_time) | object$persons$death_time > start)
    )
    simulation <- list(
      start = start,
      end = object$simulation$end_time,
      pop_size_start = pop_init,
      pop_size_end = sum(alive)
    )
  } else {
    simulation <- NULL
  }

  if (!is.null(object$groups)) {
    group_counts <- lapply(
      names(object$groups),
      function(group_name) {
        categories <- object$groups[[group_name]]$categories
        counts <- table(
          factor(
            object$persons[[group_name]][alive],
            levels = categories
          )
        ) |>
          as.numeric()
        tibble(
          category = categories,
          count = counts,
          percent = 100 * counts / sum(counts)
        )
      }
    )
    names(group_counts) <- names(object$groups)
  } else {
    group_counts <- NULL
  }

  structure(
    list(
      first_time = min(object$persons$birth_time),
      end_time = object$end_time,
      persons_alive = sum(alive),
      n_births = sum(!is.na(object$persons$birth_time)),
      n_births_u = ifelse(
        is.null(object$unrealized_births), NA, nrow(object$unrealized_births)
      ),
      n_deaths = sum(!is.na(object$persons$death_time)),
      n_unions = nrow(object$unions),
      n_dissolutions = (sum(!is.na(object$unions$end_reason) &
                              object$unions$end_reason == "dissolution")),
      age_counts = age_counts,
      simulation = simulation,
      group_counts = group_counts
    ),
    class = "summary.hammel_population"
  )
}

#' @rdname summary.hammel_population
#' @param x A `summary.hammel_population` object returned by
#'   `summary.hammel_population()`.
#' @export
print.summary.hammel_population <- function(x, ...) {

  cat(cli::style_bold("Hammel Population Summary"), "\n\n")

  print_stat("First date:", round(x$first_time, 3))
  print_stat("End time:", x$end_time)
  print_stat("Persons alive at end time:", x$persons_alive)
  print_stat("Total births:", x$n_births)
  if (!is.na(x$n_births_u)) {
    print_stat("Total unrealized births:", x$n_births_u)
  }
  print_stat("Total deaths", x$n_deaths)
  print_stat("Total unions", x$n_unions)
  print_stat("Total dissolutions", x$n_dissolutions)

  if(!is.null(x$simulation)) {
    cat("\n", cli::style_bold("Latest simulation run"), "\n", sep = "")

    print_stat(
      "Time period:", paste0(x$simulation$start, "\u002d", x$simulation$end)
    )
    print_stat("Starting population size:", x$simulation$pop_size_start)
    print_stat("Ending population size:", x$simulation$pop_size_end)
  }

  if(!is.null(x$group_counts)) {
    cat(
      "\n", cli::style_bold("Living population by group"), "\n", sep = ""
    )

    for (dimension in names(x$group_counts)) {
      counts <- x$group_counts[[dimension]]
      label_width <- max(nchar(c("Category", counts$category)))

      cat("\n", cli::style_bold(dimension), "\n", sep = "")
      cat(sprintf(
        "  %-*s %10s %9s\n",
        label_width, "Category", "Count", "Percent"
      ))

      for (i in seq_len(nrow(counts))) {
        cat(sprintf(
          "  %-*s %10.0f %9.1f\n",
          label_width,
          counts$category[i],
          counts$count[i],
          counts$percent[i]
        ))
      }
    }
  }

  counts <- x$age_counts
  width <- 25
  peak <- max(c(counts$male, counts$female, 1))
  male <- strrep("#", round(width * counts$male / peak))
  female <- strrep("#", round(width * counts$female / peak))
  age_labels <- format(counts$age, width = 7, justify = "centre")

  cat(
    "\n", cli::style_bold("Living population by age and sex"), "\n\n", sep = ""
  )
  cat(sprintf(
    "%25s   %s   %-25s\n",
    "Male",
    format("Age", width = 7, justify = "centre"),
    "Female"
  ))
  for (i in rev(seq_len(nrow(counts)))) {
    cat(sprintf(
      "%25s | %s | %-25s\n",
      male[i], age_labels[i], female[i]
    ))
  }

  invisible(x)
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

