# Generate populations with uniform ages and sampled group categories.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Generate a population with uniformly distributed ages
#'
#' @description
#'
#' Create a starter population with independent uniformly distributed ages,
#' randomly assigned sex, unknown parents, and no unions.
#'
#' @details
#'
#' Ages are sample from the uniform distribution between `min_age` and `max_age`
#' and only the living population is generated. This population will have no
#' father or mother ids. Sex and group categories follow their supplied
#' probabilities.
#'
#' The generated population will not have a realistic age distribution.
#' Typically, users will run a "burn in" simulation to produce a more realistic
#' population distribution and kinship history for the desired simulation.
#'
#' The `end_time` for the population will define the starting point for the
#' simulation. By default, this value is set to zero so that time is measured
#' since the beginning of the simulation. This means that birth times for the
#' existing population will be negative. Users who are drawing from input data
#' with time increments from real calendar time, should reset this starting
#' value as appropriate.
#'
#' If the `groups` argument is used, the population declares each dimension with
#' [group()] defaults, except that its `categories` are those given, in order:
#'
#' - the weights' names
#' - each dimension's values in order of first appearance in the table's rows.
#'
#' The first is the reference category, and zero-weight categories are declared,
#' allowing for groups that may emerge later. [set_groups()] can replace that
#' declaration; categories it omits are kept.
#'
#' @param size Number of living persons, a nonnegative whole number.
#' @param max_age Maximum initial age in years, finite and nonnegative.
#' @param end_time Reference time in years, finite. Birth times are this value
#'   minus the sampled ages. Defaults to a relative time origin of zero.
#' @param min_age Minimum initial age in years, finite and nonnegative, and no
#'   greater than `max_age`. Defaults to zero.
#' @param female_probability Probability of assigning sex "f" (female). The
#'   other sex value is "m" (male).
#' @param groups Optional group dimensions to sample, independently of age and
#'   sex. Either a named list with one named numeric vector of weights per
#'   dimension (categories as names; dimensions are drawn independently), or a
#'   data frame with one column per dimension and a numeric `weight` column,
#'   giving the joint distribution of category combinations. Weights are finite
#'   and nonnegative with at least one positive value, and are normalized.
#'   `NULL` (the default) adds no groups.
#'
#' @return A validated `hammel_population` object.
#' @export
#' @examples
#' set.seed(42)
#'
#' # generate 50000 people with a maximum age of 50
#' population <- generate_uniform_population(size = 50000, max_age = 50)
#'
#' # generate a population divided by race and region, with unequal division
#' # of race by region
#' region_race_dist <- tibble::tribble(
#'   ~race, ~region, ~weight,
#'  "W",   "NE",    0.145,
#'  "W",   "MW",    0.190,
#'  "W",   "S",     0.257,
#'  "W",   "W",     0.152,
#'  "B",   "NE",    0.022,
#'  "B",   "MW",    0.023,
#'  "B",   "S",     0.067,
#'  "B",   "W",     0.019,
#'  "O",   "NE",    0.021,
#'  "O",   "MW",    0.014,
#'  "O",   "S",     0.030,
#'  "O",   "W",     0.060,
#' )
#' population <- generate_uniform_population(50000, groups = region_race_dist)
#' prop.table(table(population$persons$race, population$persons$region))
generate_uniform_population <- function(
  size,
  max_age = 70,
  end_time = 0,
  min_age = 0,
  female_probability = 0.5,
  groups = NULL
) {
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size < 0 ||
        size != trunc(size) ||
        size > .Machine$integer.max) {
    stop(
      paste0(
        "`size` must be a nonnegative whole number no ",
        "larger than R's integer limit."
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(max_age) || is.complex(max_age) || !is.null(dim(max_age)) ||
        length(max_age) != 1L ||
        !is.finite(max_age) ||
        max_age < 0) {
    stop("`max_age` must be finite and nonnegative, in years.", call. = FALSE)
  }
  if (!is.numeric(min_age) || is.complex(min_age) || !is.null(dim(min_age)) ||
        length(min_age) != 1L ||
        !is.finite(min_age) ||
        min_age < 0) {
    stop("`min_age` must be finite and nonnegative, in years.", call. = FALSE)
  }
  if (min_age > max_age) {
    stop("`min_age` must not exceed `max_age`.", call. = FALSE)
  }
  check_time(end_time, "end_time")
  if (!is.numeric(female_probability) || length(female_probability) != 1L ||
        !is.finite(female_probability) ||
        female_probability < 0 ||
        female_probability > 1) {
    stop(
      "`female_probability` must be a number between 0 and 1.",
      call. = FALSE
    )
  }

  draw <- baseline_group_sampler(groups)
  persons <- data.frame(
    person_id = seq_len(size),
    sex = as.character(
      ifelse(
        stats::rbinom(size, size = 1, prob = female_probability) == 1,
        "f", "m"
      )
    ),
    birth_time = end_time - stats::runif(size, min = min_age, max = max_age),
    death_time = rep(NA_real_, size),
    mother_id = rep(NA_integer_, size),
    father_id = rep(NA_integer_, size)
  )
  if (!is.null(draw)) {
    for (column in names(drawn <- draw(size))) {
      persons[[column]] <- drawn[[column]]
    }
  }
  hammel_population(
    persons,
    end_time = end_time,
    groups = if (is.null(draw)) {
      NULL
    } else {
      lapply(attr(draw, "categories"), function(x) group(categories = x))
    }
  )
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Validate baseline group weights and return a function drawing `size` people.
baseline_group_sampler <- function(groups) {
  if (is.null(groups)) {
    return(NULL)
  }
  check_weights <- function(w, what) {
    if (!is.numeric(w) || !is.null(dim(w)) || !length(w) || anyNA(w) ||
          any(!is.finite(w)) ||
          any(w < 0) ||
          !any(w > 0)) {
      stop(
        sprintf(
          "%s must be finite and nonnegative, with a positive value.",
          what
        ),
        call. = FALSE
      )
    }
    # Scale before summation to avoid overflow for large finite weights.
    p <- w / max(w)
    p / sum(p)
  }
  if (is.data.frame(groups)) {
    dimensions <- setdiff(names(groups), "weight")
    if (!"weight" %in% names(groups) || !length(dimensions) || !nrow(groups)) {
      stop(
        paste0(
          "A `groups` data frame needs one column per ",
          "dimension and a `weight` column."
        ),
        call. = FALSE
      )
    }
    check_group_names(dimensions)
    for (d in dimensions) {
      check_groups(data.frame(group = groups[[d]]), "groups")
    }
    p <- check_weights(groups$weight, "`groups$weight`")
    draw <- function(size) {
      rows <- sample.int(nrow(groups), size, replace = TRUE, prob = p)
      lapply(
        stats::setNames(dimensions, dimensions), function(d) groups[[d]][rows]
      )
    }
    categories <- lapply(
      stats::setNames(dimensions, dimensions),
      function(d) unique(groups[[d]])
    )
  } else {
    if (!is.list(groups) || !length(groups)) {
      stop(
        "`groups` must be a named list of weight vectors or a data frame.",
        call. = FALSE
      )
    }
    dimensions <- names(groups)
    check_group_names(dimensions)
    probabilities <- lapply(
      dimensions, function(d) {
        labels <- names(groups[[d]])
        if (is.null(labels) || anyNA(labels) || any(!nzchar(trimws(labels))) ||
              anyDuplicated(labels)) {
          stop(
            sprintf(
              "`groups$%s` must be named by unique, nonempty categories.", d
            ),
            call. = FALSE
          )
        }
        check_weights(groups[[d]], sprintf("`groups$%s`", d))
      }
    )
    draw <- function(size) {
      stats::setNames(
        lapply(
          seq_along(dimensions),
          function(k) {
            names(groups[[k]])[sample.int(
              length(groups[[k]]), size,
              replace = TRUE,
              prob = probabilities[[k]]
            )]
          }
        ), dimensions
      )
    }
    categories <- lapply(
      stats::setNames(dimensions, dimensions),
      function(d) names(groups[[d]])
    )
  }
  attr(draw, "dimensions") <- dimensions
  # The declared categories, in the order given (the first is the reference).
  attr(draw, "categories") <- categories
  draw
}
