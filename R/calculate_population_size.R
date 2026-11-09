# Count living persons at requested times and within requested strata.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Calculate the living population at time points
#'
#' @description
#'
#'
#' @details
#'
#'
#' Count persons with `birth_time <= time` and missing `death_time` or
#' `death_time > time`. Counts include births at the requested time and exclude
#' deaths at that time, representing the population after recorded events there.
#' States in `by` are also those after all events at each time.
#'
#' @param population A `hammel_population` object either returned from
#'   [run_simulations()] or built as a baseline starting population.
#' @param times A numeric vector of finite time points in years. A population
#'   object cannot be queried beyond its known `end_time`.
#' @param by `NULL` (default) or a character vector of variables to count
#'   separately: `"sex"`, declared group dimensions (such as `"race"`),
#'   `"union_type"` and `"previous_union"`. See [calculate_mortality_rates()]
#'   for their definitions.
#' @param use_detailed_previous_union Logical, default `TRUE`, which counts
#'   `"previous_union"` as `"none"`, `"separated"` and `"widowed"`. With `FALSE`
#'   and `"previous_union"` in `by`, combine separated and widowed as
#'   `"dissolved"`.
#' @return A tibble with `time`, the `by` columns and `population`. Rows follow
#'   the requested times (repeats kept), with time varying fastest and earlier
#'   `by` variables varying faster than later ones.
#' @export
#' @examples
#' population <- generate_uniform_population(100, end_time = 0)
#' calculate_population_size(population, times = c(-10, -5, 0))
#' calculate_population_size(population, times = 0, by = "sex")
calculate_population_size <- function(
  population,
  times,
  by = NULL,
  use_detailed_previous_union = TRUE
) {
  if (!is.numeric(times) || !is.null(dim(times)) || any(!is.finite(times))) {
    stop(
      "`times` must be a numeric vector of finite times in years.",
      call. =FALSE
    )
  }
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by,
    c("sex", "group", "union_type", "previous_union"),
    "calculate_population_size",
    alias$dims
  )
  check_detailed_previous_union(use_detailed_previous_union, by)
  obs <- person_observation(population, by, if (length(times)) {
    max(times)
  } else {
    -Inf
  })
  previous <- "previous_union" %in% by
  utype <- "union_type" %in% by
  map <- if (use_detailed_previous_union) 0:2 else c(0L, 1L, 1L)
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  options <- list(
    typed = obs$typed, exclusive = obs$typed || previous,
    counts_cohab = obs$counts_cohab,
    sex = "sex" %in% by,
    group = "group" %in% by,
    union_type = utype,
    previous_union = previous,
    previous_map = map
  )
  counts <- person_snapshot_cpp(
    obs$persons, obs$unions, obs$group_indices, ng, as.numeric(times), options
  )
  sizes <- c(
    t = length(times), a = 1L, sex = if ("sex" %in% by) 2L else 1L,
    parity = 1L,
    group = if ("group" %in% by) ng else 1L,
    union_type = if (utype) length(union_type_labels(obs$typed)) else 1L,
    previous_union = if (previous) length(unique(map)) else 1L
  )
  strides <- cumprod(c(1, head(sizes, -1)))
  names(strides) <- names(sizes)
  dims <- c(
    list(dim_of(data.frame(time = times), stride = 1)),
    person_by_dims(
      by, obs, list(strides = strides), use_detailed_previous_union
    )
  )
  rows <- build_rows(dims)
  out <- finish_rows(rows[!grepl("^\\.key", names(rows))])
  out$population <- as.integer(counts[rows$.index])
  finish_strata(out, alias)
}
