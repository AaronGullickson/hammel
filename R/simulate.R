#' Run a demographic simulation
#'
#' Advance a configured model and return a population that can initialize a
#' later run. The C++ engine samples candidate waiting times across age bands
#' (fertility cell clocks for births), executes at continuous times, and
#' prioritizes death over birth on exact ties.
#' Events execute chronologically. In on-demand mode all living females face the
#' fertility schedule. Rate mode uses either overall rates for all women (with
#' a warning), or separate `union_type` schedules. A union-formation process is required. With on-demand
#' union formation, use the woman's
#' current partner, or draw one from all eligible men using conditional-logit
#' probabilities from the matching formula to form a union at the birth.
#' Men must meet the configured minimum age (default 15 years).
#' Known kin are excluded within the configured generation bound on each side
#' of a common ancestor (default two). Hard group exclusions also apply.
#' If no eligible man exists, cancel the birth by default. Explicitly set
#' `cancel_birth_if_no_match = FALSE` to preserve it with unknown father.
#' Rate initiation adds age- or duration-specific search events for unpartnered
#' women only. Duration schedules for dissolved/widowed categories measure time
#' since the most recent union ended; `none` remains age-specific.
#' Failed attempts draw another waiting time. Widowhood cancels pending births
#' and switches to unpartnered fertility; a new union switches to partnered
#' fertility, preserving birth spacing. Unpartnered births do not trigger searches
#' in rate mode. Instead, select a father using the same available male pool,
#' screening and `father_matching` model (falling back to `matching`), without
#' forming a union or reserving the man.
#' A failed match cancels the birth unless explicitly allowed by the union process.
#' Cancellation starts no new birth-spacing interval. On-demand mode rejects
#' status schedules.
#' Optional [union_dissolution()] adds one duration-dependent annual union hazard
#' carried by the woman. It preserves the pending dissolution through births,
#' ends the episode with `end_reason = "dissolution"`, and restores unpartnered
#' availability and fertility/search eligibility. Exact same-person ties select
#' death, birth, dissolution, then search. Dissolution `group_pair` schedules and
#' multipliers use symmetric pairings of both partners' groups. Rate-based
#' formation optionally selects age/duration schedules using `previous_union`:
#' none, dissolved, or widowed. Duration origins and categories are recovered
#' from recorded histories.
#' Matching uses male-group weights recalibrated at each search; simulation
#' warns once if any search needed approximate weights. A separate father model
#' has its own preferences and delta target adjustments. Event initiation is
#' unchanged. Optional cohabitation uses rate mode only: type is selected from
#' `cohab_prob` before marriage/cohabitation-specific matching. Conversion retains
#' partners and creates linked episodes, with same-person ties ordered death,
#' birth, dissolution, conversion, search. Fertility input strata use `union_type`.
#'
#' Optional [group_transition()] adds directed competing age- or duration-based
#' moves, with optional sex and union-type strata. A move updates current group
#' and affected event/pool/cell state, retaining ancestry. Grouped output always
#' includes a transitions tibble; restarts preserve the full transition history.
#'
#' Optional fertility `hetfert = het_fert(cv, alpha)` weights eligible mother
#' selection within fixed-intensity cells. Supplied positive female weights are
#' retained; absent weights are drawn under the simulation seed, and daughters
#' inherit on the log scale. Output persons includes `hetfert` whenever enabled,
#' and omits it whenever disabled, even if the input supplied weights.
#'
#' Minimum birth spacing is configured in [fertility()] or [set_fertility()]
#' and uses recorded maternal births, including on restart. Births are generated
#' by cell-accounted fertility (D101; see [fertility()]): each one-year
#' age/parity/union-type/group cell produces births at its input rate times all
#' its members, assigned to an eligible member uniformly or proportionally to
#' `hetfert` weights when enabled. Births that
#' cannot occur because a cell has no eligible member are reported in
#' `simulation$unrealized_births`. Fertility outside supplied bands is zero. Mortality must cover all ages for both sexes, ending at `Inf`.
#'
#' Events are recorded at their sampled times in `[start, start + years)`. Rates and all
#' times remain in person-year and year units. Newborn group identity uses
#' fertility's inheritance model (random parent by default, with maternal fallback), if groups are present;
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
#' Calendar-banded rates change within the same run at `time_start`/`time_end`
#' boundaries, before same-time demographic events. Histories, birth spacing and
#' duration origins persist. Pending candidates are resampled at boundaries;
#' adding redundant rate boundaries need not preserve a fixed-seed trajectory.
#' Calendar-banded `cohab_prob` is evaluated at the search instant; probability-
#' only boundaries refresh role calibration while retaining pending event dates.
#' Matching and inheritance arguments may be direct models or calendar schedules
#' with `time_start`, `time_end` and list-column `model`. The current matching
#' model is selected at partner/father choice and inheritance at birth. Model-
#' only boundaries preserve pending events and earlier histories.
#'
#' @param model A `hammel_model` containing a population, mortality, fertility,
#'   and union-formation processes, optionally union dissolution and group transitions. The population
#'   must have a finite `end_time`.
#' @param years Nonnegative simulation duration in years; fractional years are
#'   supported, as are fractional population endpoints and calendar-rate bands.
#' @param seed Optional integer seed. With a seed, the caller's RNG state is
#'   restored afterward. With `NULL`, the simulation advances the current RNG.
#' @return A `hammel_population` with updated histories and endpoint, plus a
#'   `simulation` component recording run settings and birth/death counts, plus
#'   `dissolutions` when a dissolution process is configured.
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
  population <- normalize_population_transitions(population)
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
  attach_process(model, model$processes$mortality)
  model <- attach_process(model, model$processes$fertility)
  death_rates <- model$processes$mortality$rates
  birth_rates <- model$processes$fertility$rates
  birth_interval <- model$processes$fertility$birth_interval
  ffab <- model$processes$fertility$ffab
  hetfert <- model$processes$fertility$hetfert
  if (!is.null(hetfert)) validate_fertility_weights(population$persons)
  for (death_slice in calendar_rate_slices(death_rates)) {
    if (!nrow(death_slice)) next
    for (rows in split_rate_groups(death_slice)) {
      for (sex in c("m", "f")) {
        schedule <- death_slice[rows, , drop = FALSE]
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
  }
  if (any(!names(model$processes) %in% c("mortality", "fertility", "union_formation", "union_dissolution", "group_transition"))) {
    stop("Unsupported demographic process for this engine.",
         call. = FALSE)
  }
  transition_process <- model$processes$group_transition
  if (!is.null(transition_process)) {
    attach_process(model, transition_process)
    if (!"group" %in% names(population$persons))
      stop("Group transitions require population group identities.", call. = FALSE)
  }
  dissolution_process <- model$processes$union_dissolution
  if (!is.null(dissolution_process)) attach_process(model, dissolution_process)
  union_process <- model$processes$union_formation
  attach_process(model, union_process)
  calendar_inputs <- list(mortality = death_rates, fertility = birth_rates,
    search = union_process$rates, dissolution = dissolution_process$rates,
    conversion = union_process$union_convert_rate, transition = transition_process$rates,
    cohabitation_probability = union_process$cohab_prob)
  for (name in names(calendar_inputs))
    check_calendar_coverage(calendar_inputs[[name]], start, start + years, name)
  cohab_enabled <- !is.null(union_process$cohab_prob)
  counts_cohab <- isTRUE(union_process$cohabitation_counts_as_previous_union)
  if (!cohab_enabled && "union_type" %in% names(population$unions) &&
      any(population$unions$union_type == "cohab"))
    stop("A cohabitation history requires cohab_prob to continue with union types enabled.", call. = FALSE)
  if (!cohab_enabled && (("union_type" %in% names(birth_rates) && any(birth_rates$union_type %in% c("cohab", "marriage"))) ||
      (!is.null(dissolution_process) && "union_type" %in% names(dissolution_process$rates))))
    stop("Type-specific fertility/dissolution requires cohab_prob.", call. = FALSE)
  status_fertility <- "union_type" %in% names(birth_rates)
  if (union_process$initiation == "on_demand" && status_fertility) {
    stop('On-demand union formation requires an overall fertility schedule without `union_type`.',
         call. = FALSE)
  }
  if (union_process$initiation == "rate" && !status_fertility) {
    warning(paste("Rate-based union formation is using an overall fertility schedule for both",
                  "partnered and unpartnered women. To make fertility depend on union status,",
                  "supply separate schedules using `union_type`."), call. = FALSE)
  }
  registry_population <- population
  registry_population$.transition_groups <- unique(c(
    as.character(population$transitions$group_origin), as.character(population$transitions$group_dest),
    as.character(transition_process$rates$group_origin), as.character(transition_process$rates$group_dest)))
  model_arguments <- list(matching = union_process$matching,
    father_matching = union_process$father_matching, cohab_matching = union_process$cohab_matching,
    inheritance = model$processes$fertility$inheritance)
  for (name in names(model_arguments))
    check_model_coverage(model_arguments[[name]], start, start + years, name)
  scheduled_models <- any(vapply(model_arguments, is_model_schedule, logical(1)))
  model_context <- if (scheduled_models)
    prepare_model_schedules(model_arguments, registry_population, start, start + years) else NULL
  inheritance <- if (scheduled_models) model_context$inheritance else
    prepare_group_inheritance(model$processes$fertility$inheritance, registry_population)
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
  dissolution_rates <- prepare_dissolution_rates(dissolution_process,
    population$persons, group_values)
  cohab <- prepare_cohabitation(union_process, population, group_values)
  if (!is.null(union_process$cohab_matching))
    validate_matching_groups(union_process$cohab_matching, group_values)
  if (cohab_enabled && !"union_type" %in% names(population$unions))
    population$unions$union_type <- rep("marriage", nrow(population$unions))
  if (cohab_enabled && !"previous_cohab_id" %in% names(population$unions))
    population$unions$previous_cohab_id <- rep(NA_real_, nrow(population$unions))
  engine_persons <- population$persons
  engine_persons$.hammel_group <- if (is.null(group_values)) rep(0L, nrow(engine_persons)) else
    match(as.character(engine_persons$group), as.character(group_values)) - 1L
  transition_rates <- if (is.null(transition_process)) NULL else transition_process$rates
  if (!is.null(transition_rates)) {
    transition_rates$.origin <- match(as.character(transition_rates$group_origin), as.character(group_values)) - 1L
    transition_rates$.dest <- match(as.character(transition_rates$group_dest), as.character(group_values)) - 1L
    engine_persons$.group_entry <- engine_persons$birth_time
    tx <- population$transitions
    if (nrow(tx)) {
      latest <- tapply(tx$time, tx$person_id, max)
      i <- match(names(latest), as.character(engine_persons$person_id))
      engine_persons$.group_entry[i] <- as.numeric(latest)
    }
  }
  if ("previous_union" %in% names(union_rates) || cohab_enabled) {
    previous <- previous_union_info_cpp(population$persons, population$unions, counts_cohab)
    engine_persons$.hammel_previous_union <- previous$previous_union
    engine_persons$.hammel_previous_union_end <- previous$end_time
    engine_persons$.hammel_legal_previous <- previous$legal_previous
  }
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
  population$persons <- initialize_fertility_weights(population$persons, hetfert)
  engine_persons[["hetfert"]] <- population$persons[["hetfert"]]
  result <- population
  births <- deaths <- 0
  if (years > 0) {
    matching <- if (scheduled_models) model_context$matching else
      prepare_matching_model(union_process$matching, group_values)
    father_matching <- if (scheduled_models) model_context$father_matching else
      if (is.null(union_process$father_matching)) matching else
        prepare_matching_model(union_process$father_matching, group_values)
    cohab_matching <- if (scheduled_models) model_context$cohab_matching else
      if (is.null(union_process$cohab_matching)) matching else
        prepare_matching_model(union_process$cohab_matching, group_values)
    native <- simulate_cpp(engine_persons, death_rates, birth_rates,
                           start, years, birth_interval, population$unions,
                           matching,
                           union_process$min_male_age,
                           as.integer(union_process$kinship_generations),
                           union_process$initiation == "rate",
                           union_rates, father_matching,
                           union_process$cancel_birth_if_no_match, union_process$diagnostics,
                           if (is.null(population[["ancestry"]]))
                             matrix(numeric(), nrow = nrow(engine_persons), ncol = 0L) else
                             as.matrix(population[["ancestry"]][setdiff(names(population[["ancestry"]]), "person_id")]),
                           inheritance, ffab, dissolution_rates, cohab$prob, cohab$conversion,
                           cohab_matching, counts_cohab,
                           isTRUE(getOption("hammel.check_cells", FALSE)),
                           if (is.null(hetfert)) NULL else unclass(hetfert), transition_rates,
                           if (scheduled_models) model_context$config else NULL)
    if (!isTRUE(native$cells_consistent))
      stop("Internal error: fertility cell counts are inconsistent.", call. = FALSE)
    persons <- population$persons
    added <- nrow(native$persons) - nrow(persons)
    if (added > 0L) {
      persons <- rbind(persons, persons[rep(NA_integer_, added), , drop = FALSE])
    }
    for (column in names(native$persons)) persons[[column]] <- native$persons[[column]]
    if (!is.null(hetfert)) persons$hetfert <- native$hetfert
    if (!is.null(group_values)) persons$group <- group_values[native$groups + 1L]
    result$persons <- tibble::as_tibble(persons)
    if (!is.null(group_values)) {
      tx <- tibble::as_tibble(native$transitions)
      tx$group_origin <- group_values[tx$group_origin + 1L]
      tx$group_dest <- group_values[tx$group_dest + 1L]
      result$transitions <- tibble::as_tibble(rbind(population$transitions, tx))
    }
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
    if (native$weighted_infeasible > 0) {
      warning(sprintf(paste("weighted_choice could not meet male-group targets in %d search(es),",
                            "for example because no eligible woman may partner a group's men;",
                            "capped approximate weights were used."), as.integer(native$weighted_infeasible)),
              call. = FALSE)
    }
    if (native$draw_approximate > 0) {
      warning(sprintf(paste("A custom matching formula exceeded its within-bin bound in %d draw(s);",
                            "partner draws were approximate. The formula varies faster than the",
                            "1/32-year bound grid resolves."), as.integer(native$draw_approximate)),
              call. = FALSE)
    }
  }
  if (union_process$diagnostics) {
    # Weighted choice records the log group weights applied at each union's
    # search (0 = unit weight), one column per group; earlier unions are NA.
    if (!is.null(group_values)) {
      added <- if (years > 0) native$diagnostic_log_weights else matrix(numeric(), 0, length(group_values))
      new_rows <- seq_len(nrow(added)) + nrow(population$unions)
      for (g in seq_along(group_values)) {
        column <- paste0(".log_weight_", as.character(group_values)[g])
        if (!column %in% names(result$unions))
          result$unions[[column]] <- rep(NA_real_, nrow(result$unions))
        result$unions[[column]][new_rows] <- added[, g]
      }
    }
  }
  result$end_time <- start + years
  result$simulation <- list(start_time = start, end_time = result$end_time,
    birth_interval = birth_interval, ffab = ffab,
    births = births, deaths = deaths)
  if (!is.null(hetfert)) result$simulation$hetfert <- hetfert
  result$simulation$unrealized_births <- format_unrealized_births(
    if (years > 0) native$unrealized else NULL, birth_rates, group_values)
  if (cohab_enabled) {
    result$simulation$cohabitation_counts_as_previous_union <- counts_cohab
    result$simulation$conversions <- if (years > 0) native$conversions else 0
  }
  if (!is.null(transition_process))
    result$simulation$transitions <- if (years > 0) nrow(native$transitions) else 0L
  if (!is.null(dissolution_process))
    result$simulation$dissolutions <- if (years > 0) native$dissolutions else 0
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

# Unrealized expected births from cells with members but no eligible woman
# (D101), by one-year period from the run start and cell. Columns mirror the
# fertility table's strata: group, union_type and parity bounds appear only
# when the schedules use them. Always a tibble, possibly with zero rows.
format_unrealized_births <- function(native, rates, group_values) {
  if (is.null(native) || !length(native$period_start)) {
    native <- data.frame(period_start = numeric(), period_end = numeric(),
      group = integer(), union_type = character(), parity_start = numeric(),
      parity_end = numeric(), age_start = numeric(), age_end = numeric(),
      unrealized_births = numeric(), years_without_eligible = numeric())
  }
  out <- tibble::tibble(period_start = native$period_start, period_end = native$period_end)
  if (!is.null(group_values)) out$group <- group_values[native$group + 1L]
  if ("union_type" %in% names(rates)) out$union_type <- as.character(native$union_type)
  if ("parity_start" %in% names(rates)) {
    out$parity_start <- native$parity_start
    out$parity_end <- native$parity_end
  }
  out$age_start <- native$age_start
  out$age_end <- native$age_end
  out$unrealized_births <- native$unrealized_births
  out$years_without_eligible <- native$years_without_eligible
  out
}
