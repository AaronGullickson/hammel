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
#' @param groups Optional group dimensions to sample, independently of age and
#'   sex. Either a named list with one named numeric vector of weights per
#'   dimension (categories as names; dimensions are drawn independently), or a
#'   data frame with one column per dimension and a numeric `weight` column,
#'   giving the joint distribution of category combinations. Weights are finite
#'   and nonnegative with at least one positive value, and are normalized.
#'   Realized counts are random. The population declares each dimension with
#'   [group()] defaults, except that its `categories` are those given, in
#'   order: the weights' names, or each dimension's values in order of first
#'   appearance in the table's rows. The first is the reference category, and
#'   zero-weight categories are declared (they start empty). [set_groups()] can
#'   replace that declaration; categories it omits are kept. `NULL`
#'   (the default) adds no groups and consumes no extra random draws.
#'
#' @return A validated `hammel_population` object.
#' @export
#' @examples
#' set.seed(42)
#' population <- create_baseline_population(size = 50000, max_age = 70)
#' model <- hammel_model() |>
#'   set_population(population)
create_baseline_population <- function(size, max_age = 70, end_time = 0,
                                female_probability = 0.5, groups = NULL) {
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

  draw <- baseline_group_sampler(groups)
  persons <- data.frame(
    person_id = seq_len(size),
    sex = as.character(ifelse(stats::rbinom(size, size = 1, prob = female_probability) == 1,
                 "f", "m")),
    birth_time = end_time - stats::runif(size, min = 0, max = max_age),
    death_time = rep(NA_real_, size),
    mother_id = rep(NA_integer_, size),
    father_id = rep(NA_integer_, size)
  )
  if (!is.null(draw)) for (column in names(drawn <- draw(size))) persons[[column]] <- drawn[[column]]
  hammel_population(persons, end_time = end_time,
                    groups = if (is.null(draw)) NULL else
                      lapply(attr(draw, "categories"), function(x) group(categories = x)))
}

# Validate baseline group weights and return a function drawing `size` people.
baseline_group_sampler <- function(groups) {
  if (is.null(groups)) return(NULL)
  check_weights <- function(w, what) {
    if (!is.numeric(w) || !is.null(dim(w)) || !length(w) || anyNA(w) ||
        any(!is.finite(w)) || any(w < 0) || !any(w > 0))
      stop(sprintf("%s must be finite and nonnegative, with a positive value.", what), call. = FALSE)
    # Scale before summation to avoid overflow for large finite weights.
    p <- w / max(w)
    p / sum(p)
  }
  if (is.data.frame(groups)) {
    dimensions <- setdiff(names(groups), "weight")
    if (!"weight" %in% names(groups) || !length(dimensions) || !nrow(groups))
      stop("A `groups` data frame needs one column per dimension and a `weight` column.", call. = FALSE)
    check_group_names(dimensions)
    for (d in dimensions) check_groups(data.frame(group = groups[[d]]), "groups")
    p <- check_weights(groups$weight, "`groups$weight`")
    draw <- function(size) {
      rows <- sample.int(nrow(groups), size, replace = TRUE, prob = p)
      lapply(stats::setNames(dimensions, dimensions), function(d) groups[[d]][rows])
    }
    categories <- lapply(stats::setNames(dimensions, dimensions), function(d) unique(groups[[d]]))
  } else {
    if (!is.list(groups) || !length(groups))
      stop("`groups` must be a named list of weight vectors or a data frame.", call. = FALSE)
    dimensions <- names(groups)
    check_group_names(dimensions)
    probabilities <- lapply(dimensions, function(d) {
      labels <- names(groups[[d]])
      if (is.null(labels) || anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels))
        stop(sprintf("`groups$%s` must be named by unique, nonempty categories.", d), call. = FALSE)
      check_weights(groups[[d]], sprintf("`groups$%s`", d))
    })
    draw <- function(size) stats::setNames(lapply(seq_along(dimensions), function(k)
      names(groups[[k]])[sample.int(length(groups[[k]]), size, replace = TRUE,
                                    prob = probabilities[[k]])]), dimensions)
    categories <- lapply(stats::setNames(dimensions, dimensions), function(d) names(groups[[d]]))
  }
  attr(draw, "dimensions") <- dimensions
  # The declared categories, in the order given (the first is the reference).
  attr(draw, "categories") <- categories
  draw
}
