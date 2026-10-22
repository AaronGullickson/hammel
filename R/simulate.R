#' Run a fertility and mortality simulation
#'
#' Advance a configured model and return a population that can initialize a
#' later run. The C++ engine samples candidate waiting times across age bands,
#' executes at continuous times, and prioritizes death over birth on exact ties.
#' Events execute chronologically. In on-demand mode all living females face the
#' fertility schedule. Rate mode uses either overall rates for all women (with
#' a warning), or separate `union_status` schedules. A union-formation process is required. With on-demand
#' union formation, use the woman's
#' current partner, or sample men and choose one using conditional-logit
#' probabilities from the matching formula to form a union at the birth.
#' Sampled men must meet the configured minimum age (default 15 years).
#' Known kin are excluded within the configured generation bound on each side
#' of a common ancestor (default two). Hard group exclusions also apply.
#' If no sampled men are eligible, cancel the birth by default. Explicitly set
#' `cancel_birth_if_no_match = FALSE` to preserve it with unknown father.
#' Rate initiation adds age-specific search events for unpartnered women only.
#' Failed attempts draw another waiting time. Widowhood cancels pending births
#' and switches to unpartnered fertility; a new union switches to partnered
#' fertility, preserving birth spacing. Unpartnered births do not trigger searches
#' in rate mode. Instead, select a father using the same available male pool,
#' screening and `father_matching` model (falling back to `matching`), without
#' forming a union or reserving the man.
#' A failed match cancels the birth unless explicitly allowed by the union process.
#' Cancellation starts no new birth-spacing interval. On-demand mode rejects
#' status schedules.
#' Optional `matching_method = "two_sided"` on [union_formation()] uses
#' probabilistic sequential temporary markets instead of default ego-only
#' female choice. Only ego's outcome is retained; event initiation is unchanged.
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
#' times remain in person-year and year units. Newborn group identity uses
#' fertility's inheritance model (maternal by default), if groups are present;
#' other newborn custom columns are missing. Grouped rate
#' tables use each person's group (the woman's group for union searches).
#' Ungrouped tables apply to everyone. Missing population groups in a grouped
#' rate table are errors. Without population groups, each grouped table uses its
#' first alphabetical group with a warning; this does not create group identities.
#' Matching can use group predictors and pair effects. Existing custom columns and histories are preserved. Restarting samples new
#' pending events and does not reproduce an uninterrupted run exactly.
#' When ancestry is tracked, each newborn receives half of each parent's shares,
#' with a missing parent contributing its half to unknown ancestry. Initial
#' shares are retained, including on restart; grouped populations without shares
#' initialize from recorded identities. No additional random draws are made.
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
  population <- normalize_population_ancestry(population)
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
  model <- set_process(model, model$processes$fertility)
  death_rates <- model$processes$mortality$rates
  birth_rates <- model$processes$fertility$rates
  birth_interval <- model$processes$fertility$birth_interval
  ffab <- model$processes$fertility$ffab
  for (rows in split_rate_groups(death_rates)) {
    for (sex in c("m", "f")) {
      schedule <- death_rates[rows, , drop = FALSE]
      schedule <- schedule[schedule$sex == sex, , drop = FALSE]
      schedule <- schedule[order(schedule$age_start), , drop = FALSE]
      if (schedule$age_start[1] != 0 || tail(schedule$age_end, 1) != Inf ||
          (nrow(schedule) > 1L &&
           any(head(schedule$age_end, -1) != tail(schedule$age_start, -1)))) {
        stop("Mortality schedules must cover [0, Inf) without gaps for each sex and group.",
             call. = FALSE)
      }
    }
  }
  if (any(!names(model$processes) %in% c("mortality", "fertility", "union_formation"))) {
    stop("This engine supports mortality, fertility, and union formation only.",
         call. = FALSE)
  }
  union_process <- model$processes$union_formation
  set_process(model, union_process)
  status_fertility <- "union_status" %in% names(birth_rates)
  if (union_process$initiation == "on_demand" && status_fertility) {
    stop('On-demand union formation requires an overall fertility schedule without `union_status`.',
         call. = FALSE)
  }
  if (union_process$initiation == "rate" && !status_fertility) {
    warning(paste("Rate-based union formation is using an overall fertility schedule for both",
                  "partnered and unpartnered women. To make fertility depend on union status,",
                  "supply separate schedules using `union_status`."), call. = FALSE)
  }
  inheritance <- prepare_group_inheritance(model$processes$fertility$inheritance, population)
  group_values <- inheritance$values
  death_rates <- prepare_group_rates(death_rates, population$persons, "mortality", group_values)
  birth_rates <- prepare_group_rates(birth_rates, population$persons, "fertility", group_values)
  union_rates <- if (is.null(union_process$rates))
    data.frame(age_start = numeric(), age_end = numeric(), rate = numeric()) else
    prepare_group_rates(union_process$rates, population$persons, "union formation", group_values)
  death_rates <- apply_group_prop_hazard(death_rates, population$persons,
    model$processes$mortality$group_prop_hazard, "mortality", group_values)
  birth_rates <- apply_group_prop_hazard(birth_rates, population$persons,
    model$processes$fertility$group_prop_hazard, "fertility", group_values)
  if (union_process$initiation == "rate") {
    union_rates <- apply_group_prop_hazard(union_rates, population$persons,
      union_process$group_prop_hazard, "union formation", group_values)
  }
  validate_matching_groups(union_process$matching, group_values)
  if (!is.null(union_process$father_matching))
    validate_matching_groups(union_process$father_matching, group_values)
  engine_persons <- population$persons
  engine_persons$.hammel_group <- if (is.null(group_values)) rep(0L, nrow(engine_persons)) else
    match(as.character(engine_persons$group), as.character(group_values)) - 1L
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
    two_sided <- identical(union_process$matching_method, "two_sided")
    matching <- prepare_matching_model(union_process$matching, group_values, two_sided)
    father_matching <- if (is.null(union_process$father_matching)) matching else
      prepare_matching_model(union_process$father_matching, group_values, two_sided)
    native <- simulate_cpp(engine_persons, death_rates, birth_rates,
                           start, years, birth_interval, population$unions,
                           sample_size, matching,
                           union_process$min_male_age,
                           as.integer(union_process$kinship_generations),
                           union_process$initiation == "rate",
                           union_rates, father_matching,
                           union_process$cancel_birth_if_no_match, union_process$diagnostics,
                           if (is.null(population[["ancestry"]]))
                             matrix(numeric(), nrow = nrow(engine_persons), ncol = 0L) else
                             as.matrix(population[["ancestry"]][setdiff(names(population[["ancestry"]]), "person_id")]),
                           inheritance, ffab)
    persons <- population$persons
    added <- nrow(native$persons) - nrow(persons)
    if (added > 0L) {
      persons <- rbind(persons, persons[rep(NA_integer_, added), , drop = FALSE])
    }
    for (column in names(native$persons)) persons[[column]] <- native$persons[[column]]
    if (!is.null(group_values)) persons$group <- group_values[native$groups + 1L]
    result$persons <- tibble::as_tibble(persons)
    if (!is.null(population[["ancestry"]])) {
      ancestry <- tibble::as_tibble(native$ancestry, .name_repair = "minimal")
      names(ancestry) <- setdiff(names(population[["ancestry"]]), "person_id")
      result$ancestry <- tibble::tibble(person_id = persons$person_id)
      for (column in names(ancestry)) result$ancestry[[column]] <- ancestry[[column]]
    }
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
  if (union_process$diagnostics) {
    for (column in c(".n_eligible", ".utility_rank")) {
      if (!column %in% names(result$unions))
        result$unions[[column]] <- rep(NA_integer_, nrow(result$unions))
    }
  }
  result$end_time <- start + years
  result$simulation <- list(start_time = start, end_time = result$end_time,
    birth_interval = birth_interval, ffab = ffab,
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
