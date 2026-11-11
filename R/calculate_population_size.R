# Count living persons at requested times and within requested strata.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Calculate the living population at time points
#'
#' @description
#'
#' Calculate population counts from a returned simulation result or a starting
#' population, at chosen time points. Counts can be disaggregated by stratifying
#' variables.
#'
#' @details
#'
#' The `by` argument can take the standard variables of `sex`, `union_type`, or
#' `previous_union` as well as any group variables defined in the `population`
#' input.
#'
#' The user can input any numeric value for `times` (including decimal values),
#' but no value may come after the `end_time` of the simulation.
#'
#' @param population A `hammel_population` object either returned from
#'   [run_simulation()] or built as a baseline starting population.
#' @param times A numeric vector of finite time points in years at which to
#'   calculate population size.
#' @param by Character vector of variables to disaggregate population counts by.
#'   Defaults to `NULL` for no disaggregation. Populations counts can be
#'   disaggregated by `"sex"`, `"union_type"`, `"previous_union"` and any
#'   declared group dimensions.
#' @param use_detailed_previous_union Logical value that indicates whether
#'   disaggregation by `"previous_union"` will use `"separated"` and `"widowed"`
#'   or a combined `"dissolved"` category.
#' @return A tibble with `time`, the `by` columns and `population`.
#' @export
#' @examples
#' # calculate population size on a generatd population
#' pop_start <- generate_uniform_population(1000, end_time = 0)
#' calculate_population_size(pop_start, times = c(-10, -5, 0))
#'
#' # run a simulation
#' model <- hammel_model() |>
#'   set_population(pop_start) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_union_formation("on_demand")
#' pop_end <- run_simulation(model, 100)
#'
#' # calculate at 20 year intervals
#' calculate_population_size(pop_end, times = seq(0, 100, 20))
#'
#' # caclulate at 20 year intervals by sex
#' calculate_population_size(pop_end, times = seq(0, 100, 20), by = "sex")
calculate_population_size <- function(
  population,
  times,
  by = NULL,
  use_detailed_previous_union = TRUE
) {
  if (!is.numeric(times) || !is.null(dim(times)) || any(!is.finite(times))) {
    stop(
      "`times` must be a numeric vector of finite times in years.",
      call. = FALSE
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
  last_time <- if (length(times)) max(times) else -Inf
  obs <- person_observation(population, by, last_time)
  previous <- "previous_union" %in% by
  utype <- "union_type" %in% by
  map <- if (use_detailed_previous_union) 0:2 else c(0L, 1L, 1L)
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  options <- list(
    typed = obs$typed,
    exclusive = obs$typed || previous,
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
    t = length(times),
    a = 1L,
    sex = if ("sex" %in% by) 2L else 1L,
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
