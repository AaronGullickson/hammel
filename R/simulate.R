#' Run a fertility and mortality simulation
#'
#' Advance a configured model and return a population that can initialize a
#' later run. The C++ engine samples candidate waiting times across age bands,
#' executes at continuous times, and prioritizes death over birth on exact ties.
#' Events execute chronologically. All living females face the
#' fertility schedule. A union-formation process is required. With on-demand
#' union formation, use the woman's
#' current partner, or sample men and choose one using conditional-logit
#' probabilities from the matching formula to form a union at the birth.
#' Pair eligibility checks are not yet applied.
#' If no men are available, preserve the birth with unknown father and no union.
#' Rate-based union initiation is not yet implemented.
#'
#' Minimum birth spacing is configured in [fertility()] or [set_fertility()]
#' and uses recorded maternal births, including on restart.
#' Midpoint eligible-risk adjustment uses the input fertility integral over
#' the preceding birth interval at each original age band's midpoint. Each band
#' has one adjusted hazard and shares the mortality waiting-time sampler.
#' Initial eligibility is not calibrated.
#' This is not a guarantee of exact ASFR matching. Fertility outside supplied
#' bands is zero. Mortality must cover all ages for both sexes, ending at `Inf`.
#'
#' Events are recorded at their sampled times in `[start, start + years)`. Rates and all
#' times remain in person-year and year units. Newborn custom columns are missing;
#' existing custom columns and histories are preserved. Restarting samples new
#' pending events and does not reproduce an uninterrupted run exactly.
#'
#' @param model A `hammel_model` containing a population, mortality, fertility,
#'   and union-formation processes. The population must have a finite `end_time`.
#' @param years Nonnegative simulation duration in years.
#' @param seed Optional integer seed. With a seed, the caller's RNG state is
#'   restored afterward. With `NULL`, the simulation advances the current RNG.
#' @return A `hammel_population` with updated histories and endpoint, plus a
#'   `simulation` component recording run settings and birth/death counts.
#' @export
#' @examples
#' population <- create_baseline_population(100, max_age = 50)
#' model <- hammel_model() |>
#'   set_population(population) |>
#'   set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
#'     age_end = Inf, rate = 0.01)) |>
#'   set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0.1)) |>
#'   set_union_formation("on_demand",
#'     create_matching_model(~ age_difference, c(age_difference = 0)))
#' result <- simulate(model, years = 2, seed = 42)
#' result$simulation
simulate <- function(model, years, seed = NULL) {
  check_model(model)
  if (is.null(model$population) || is.null(model$processes$mortality) ||
      is.null(model$processes$fertility) || is.null(model$processes$union_formation)) {
    stop("The model needs a population, mortality, fertility, and union-formation process.",
         call. = FALSE)
  }
  population <- model$population
  validate_population(population)
  start <- population$end_time
  if (is.null(start)) stop("The population needs a finite `end_time`.", call. = FALSE)
  if (!is.numeric(years) || length(years) != 1L || !is.finite(years) || years < 0) {
    stop("`years` must be finite and nonnegative.", call. = FALSE)
  }
  if (!is.finite(start + years) || (years > 0 && start + years == start)) {
    stop("Simulation times cannot be represented at this scale.", call. = FALSE)
  }
  # Recheck edited specifications, without introducing a model diagnostic API.
  set_process(model, model$processes$mortality)
  set_process(model, model$processes$fertility)
  death_rates <- model$processes$mortality$rates
  birth_rates <- model$processes$fertility$rates
  birth_interval <- model$processes$fertility$birth_interval
  for (sex in c("m", "f")) {
    schedule <- death_rates[death_rates$sex == sex, , drop = FALSE]
    schedule <- schedule[order(schedule$age_start), , drop = FALSE]
    if (schedule$age_start[1] != 0 || tail(schedule$age_end, 1) != Inf ||
        (nrow(schedule) > 1L &&
         any(head(schedule$age_end, -1) != tail(schedule$age_start, -1)))) {
      stop("Mortality schedules must cover [0, Inf) without gaps for each sex.",
           call. = FALSE)
    }
  }
  if (any(!names(model$processes) %in% c("mortality", "fertility", "union_formation"))) {
    stop("This engine supports mortality, fertility, and on-demand union formation only.",
         call. = FALSE)
  }
  union_process <- model$processes$union_formation
  set_process(model, union_process)
  if (union_process$initiation != "on_demand") {
    stop("Rate-based union formation is not yet implemented.", call. = FALSE)
  }
  sample_size <- as.integer(union_process$sample_size)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != trunc(seed) || abs(seed) > .Machine$integer.max) {
      stop("`seed` must be NULL or an integer within R's integer range.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  result <- population
  births <- deaths <- 0
  if (years > 0) {
    native <- simulate_cpp(population$persons, death_rates, birth_rates,
                           start, years, birth_interval, population$unions,
                           sample_size, prepare_matching_model(union_process$matching))
    persons <- population$persons
    added <- nrow(native$persons) - nrow(persons)
    if (added > 0L) {
      persons <- rbind(persons, persons[rep(NA_integer_, added), , drop = FALSE])
    }
    for (column in names(native$persons)) persons[[column]] <- native$persons[[column]]
    result$persons <- tibble::as_tibble(persons)
    unions <- result$unions
    added_unions <- nrow(native$unions) - nrow(unions)
    if (added_unions > 0L) {
      unions <- rbind(unions, unions[rep(NA_integer_, added_unions), , drop = FALSE])
    }
    for (column in names(native$unions)) unions[[column]] <- native$unions[[column]]
    result$unions <- tibble::as_tibble(unions)
    births <- native$births
    deaths <- native$deaths
  }
  result$end_time <- start + years
  result$simulation <- list(start_time = start, end_time = result$end_time,
    birth_interval = birth_interval,
    births = births, deaths = deaths)
  result$persons <- tibble::as_tibble(result$persons)
  result$unions <- tibble::as_tibble(result$unions)
  validate_population(result)
  result
}

validate_birth_interval <- function(duration) {
  if (!is.numeric(duration) || !is.null(dim(duration)) || length(duration) != 1L ||
      !is.finite(duration) || duration < 0) {
    stop("`birth_interval` must be a finite nonnegative number of years.",
         call. = FALSE)
  }
  invisible(duration)
}
