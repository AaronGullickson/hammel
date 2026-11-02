#' Specify a conditional-logit matching model
#'
#' Define relative partner utilities with an R formula and named coefficients.
#' Predictors include `age_difference`: male age minus female age,
#' in continuous years at the search time. Formula transformations such as
#' `I(age_difference^2)` and `abs(age_difference)` are supported through
#' [stats::model.matrix()]. An intercept is ignored because it cancels from
#' partner-only conditional-logit probabilities. No wait alternative is included.
#' `female_group` and `male_group` are factors with explicit `group_levels` and
#' fixed treatment contrasts (first level is reference). `same_group` is numeric
#' 1 for equal identities, 0 otherwise. Pair effects are added after the formula,
#' including custom formulas. Group effects require population identities.
#' These modify relative candidate weights, not search hazards or eligibility.
#'
#' Coefficient names must match the non-intercept model-matrix columns exactly.
#' Supply all coefficients explicitly; no demographic coefficients are defaulted.
#' Formula functions must be deterministic and must not mutate simulation state.
#' Data-dependent transformations can change with the candidate pool; use fixed
#' transformations when coefficients should have a fixed interpretation.
#'
#' @param formula A one-sided R formula using `age_difference`, `female_group`,
#'   `male_group`, and/or `same_group`. Use `~ 1` and `numeric()` coefficients
#'   for pair-effects-only matching.
#' @param coefficients A finite named numeric vector of utility coefficients.
#' @param group_exogamy Optional finite named numeric vector of symmetric
#'   group-pair effects: `"A:B"` applies the same log-weight effect to female
#'   A/male B and female B/male A. Specify each unordered pair once; reverse
#'   duplicates and `::` syntax are errors. Omitted pairs have zero effect;
#'   within-group entries are allowed. Alternatively provide a data frame with
#'   `female_group`, `male_group`, and `log_odds`; each row specifies an unordered
#'   pair and expands symmetrically, including for colon-containing labels.
#'   Effects add to formula scores and do not prohibit matches.
#' @param group_levels Optional ordered vector of group labels (character or
#'   numeric codes). Required with at least two levels when the formula uses
#'   `female_group` or `male_group`; must cover population groups. Sets fixed
#'   treatment reference/columns, independent of the available pool. Not required
#'   for `same_group` or pair effects alone.
#' @param prevent_exogamy `FALSE` (default) for no group restrictions; `TRUE`
#'   excludes all different-group candidates. A character vector such as
#'   `c("A::B", "A:C")` excludes selected pairs with female-first `:` and
#'   symmetric `::` notation. Alternatively supply a directional data frame
#'   with `female_group` and `male_group`. Pairs must join different groups.
#'   Restricted men are ineligible, alongside age and kinship exclusions. A
#'   search fails only when no eligible man exists in the whole pool, following
#'   existing failed-search and birth cancellation rules. Restrictions override
#'   finite preferences, require population identities, and do not invalidate
#'   existing partners.
#' @param delta Optional finite numeric vector with unique group names.
#'   Multiply mirrored male targets by `exp(delta)`
#'   and normalize. Omitted groups have zero adjustment. Positive values raise
#'   the male target relative to baseline and other groups; negative values
#'   lower it. Only differences matter. Does not set pair-specific ratios or
#'   change female search rates. Names must belong to the population's identity
#'   registry (including declared future groups).
#' @return An inspectable `hammel_matching_model` list containing `formula`
#'   and `coefficients`. This specifies a model; it does not estimate one.
#' @export
#' @examples
#' # Illustrative coefficients, not demographic defaults.
#' matching <- create_matching_model(
#'   ~ age_difference + I(age_difference^2),
#'   c(age_difference = 0.1, "I(age_difference^2)" = -0.02)
#' )
create_matching_model <- function(formula, coefficients, group_exogamy = NULL, group_levels = NULL, prevent_exogamy = FALSE,
                                  delta = NULL) {
  model <- structure(list(formula = formula, coefficients = coefficients,
                          group_exogamy = group_exogamy, group_levels = group_levels,
                          prevent_exogamy = prevent_exogamy,
                          delta = delta),
                     class = "hammel_matching_model")
  validate_matching_model(model)
  warn_absorbed_group_formula(model)
  model
}

# Male-group formula effects can be absorbed by calibrated target weights.
warn_absorbed_group_formula <- function(model) {
  if (!"male_group" %in% all.vars(model$formula)) return(invisible(FALSE))
  warning(paste("Formula terms using male_group are absorbed by calibrated group weights:",
    "main effects have no effect and interactions are only partly reflected.",
    "Use delta to adjust male-group targets."), call. = FALSE)
  invisible(TRUE)
}

#' Create a quadratic matching model from interpretable age-gap preferences
#'
#' Convert a preferred age gap and a symmetric half-weight distance into a
#' formula and coefficients for [create_matching_model()]. Age difference is
#' male age minus female age. Candidates at the preferred gap plus or minus
#' `half_weight_distance` receive half the selection weight of a candidate at
#' the preferred gap. Availability and eligibility also affect realized gaps.
#'
#' The quadratic coefficient is `-log(2) / half_weight_distance^2` and the linear
#' coefficient is `-2 * preferred_age_gap * quadratic`. This specifies relative
#' selection weights, not a halving of utility or a realized standard deviation.
#' @param preferred_age_gap Finite numeric scalar in years; negative values allowed.
#' @param group_exogamy Optional additive group-pair effects passed to
#'   [create_matching_model()].
#' @param prevent_exogamy Optional hard group restrictions passed to
#'   [create_matching_model()]. Defaults to `FALSE`.
#' @param half_weight_distance Finite positive numeric scalar in years from the
#'   preferred gap at which selection weight halves.
#' @inheritParams create_matching_model
#' @return A `hammel_matching_model` with linear and quadratic age-difference terms.
#' @export
#' @examples
#' # Illustrative preferences, not demographic defaults.
#' matching <- create_age_gap_matching_model(3, 5)
create_age_gap_matching_model <- function(preferred_age_gap, half_weight_distance, group_exogamy = NULL, prevent_exogamy = FALSE,
                                          delta = NULL) {
  for (name in c("preferred_age_gap", "half_weight_distance")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.null(dim(value)) ||
        !is.finite(value)) {
      stop("`", name, "` must be a finite numeric scalar.", call. = FALSE)
    }
  }
  if (half_weight_distance <= 0)
    stop("`half_weight_distance` must be positive.", call. = FALSE)
  quadratic <- -(sqrt(log(2)) / half_weight_distance)^2
  linear <- -2 * (preferred_age_gap * quadratic)
  if (!is.finite(quadratic) || quadratic >= 0 || !is.finite(linear))
    stop("Inputs must yield finite coefficients with a negative quadratic term.",
         call. = FALSE)
  create_matching_model(~ age_difference + I(age_difference^2),
    c(age_difference = unname(linear), "I(age_difference^2)" = unname(quadratic)),
    group_exogamy = group_exogamy, prevent_exogamy = prevent_exogamy,
    delta = delta)
}

validate_delta <- function(x) {
  if (is.null(x)) return(invisible(x))
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
      (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
       any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x)))))
    stop("`delta` must be a finite numeric vector with unique nonempty group names.", call. = FALSE)
  invisible(x)
}

validate_matching_model <- function(model) {
  if (!inherits(model, "hammel_matching_model") || !is.list(model)) {
    stop("`matching` must be created with create_matching_model().", call. = FALSE)
  }
  validate_delta(model$delta)
  formula <- model$formula
  if (!inherits(formula, "formula") || length(formula) != 2L ||
      any(!all.vars(formula) %in% c("age_difference", "female_group", "male_group", "same_group"))) {
    stop("Use a one-sided matching formula with `age_difference`, `female_group`, `male_group`, or `same_group`.", call. = FALSE)
  }
  if (length(attr(stats::terms(formula), "offset"))) {
    stop("Matching formula offsets are not supported.", call. = FALSE)
  }
  beta <- model$coefficients
  if (!is.numeric(beta) || is.complex(beta) || !is.null(dim(beta)) || any(!is.finite(beta)) ||
      (length(beta) && is.null(names(beta))) || anyNA(names(beta)) || any(!nzchar(names(beta))) ||
      anyDuplicated(names(beta))) {
    stop("`coefficients` must be a finite numeric vector with unique nonempty names.",
         call. = FALSE)
  }
  expand_group_exogamy(model$group_exogamy)
  expand_prevent_exogamy(model$prevent_exogamy)
  levels <- model$group_levels
  if (!is.null(levels)) {
    labels <- matching_group_labels(levels, "group_levels")
    if (!length(labels) || anyDuplicated(labels))
      stop("`group_levels` must be unique nonempty group labels.", call. = FALSE)
  }
  categorical <- any(all.vars(formula) %in% c("female_group", "male_group"))
  if (categorical && length(levels) < 2L)
    stop("Categorical matching formulas require at least two explicit `group_levels`.", call. = FALSE)
  # A fixed validation grid covers both sexes' group levels, equality, and gaps.
  data <- data.frame(age_difference = seq(-50, 50, length.out = 101))
  if (matching_has_group_formula(model)) {
    labels <- if (categorical) as.character(levels) else c("A", "B")
    size <- max(101L, 2L * length(labels))
    female <- rep(labels, length.out = size)
    male <- rep(c(labels[-1L], labels[1L]), length.out = size)
    male[seq_len(size) %% 2L == 0L] <- female[seq_len(size) %% 2L == 0L]
    data <- matching_predictor_data(model, rep(data$age_difference, length.out = size), female, male)
  }
  design <- matching_design_matrix(formula, data)
  if (!setequal(colnames(design), names(beta))) {
    stop("Coefficient names must match all non-intercept model-matrix columns.",
         call. = FALSE)
  }
  invisible(model)
}

matching_design_matrix <- function(formula, data) {
  design <- withCallingHandlers({
    frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass,
                                drop.unused.levels = FALSE)
    factors <- frame[vapply(frame, is.factor, logical(1))]
    contrasts <- if (length(factors)) lapply(factors, function(x)
      stats::contr.treatment(levels(x))) else NULL
    stats::model.matrix(formula, data = frame, contrasts.arg = contrasts)
  }, warning = function(w) stop("Matching formula evaluation: ", conditionMessage(w),
                                call. = FALSE))
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (nrow(design) != nrow(data) || any(!is.finite(design))) {
    stop("Matching formula must produce finite values for every candidate.", call. = FALSE)
  }
  design
}

# Evaluate the same formula for R diagnostics and native fallback scoring.
calculate_matching_utilities <- function(matching, age_difference, female_group = NULL,
                                         male_group = NULL, add_group_exogamy = TRUE) {
  check_numeric_column(age_difference, "age_difference")
  if (!length(age_difference)) return(numeric())
  data <- matching_predictor_data(matching, age_difference, female_group, male_group)
  design <- matching_design_matrix(matching$formula, data)
  if (!setequal(colnames(design), names(matching$coefficients))) {
    stop("Matching formula columns changed for the supplied candidates.", call. = FALSE)
  }
  utility <- drop(design %*% matching$coefficients[colnames(design)])
  if (add_group_exogamy && nrow(expand_group_exogamy(matching$group_exogamy))) {
    if (is.null(female_group) || is.null(male_group))
      stop("Pair effects require female and male group data.", call. = FALSE)
    utility <- utility + matching_pair_offsets(matching, female_group, male_group)
  }
  if (any(!is.finite(utility))) {
    stop("Matching utilities must be finite.", call. = FALSE)
  }
  unname(utility)
}

calculate_matching_probabilities <- function(matching, age_difference, female_group = NULL, male_group = NULL) {
  utility <- calculate_matching_utilities(matching, age_difference, female_group, male_group)
  if (!length(utility)) return(numeric())
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  allowed <- rep(TRUE, length(utility))
  if (restrictions$all || nrow(restrictions$pairs)) {
    if (is.null(female_group) || is.null(male_group))
      stop("Group restrictions require female and male group data.", call. = FALSE)
    allowed <- matching_group_allowed(matching, female_group, male_group)
  }
  weights <- numeric(length(utility))
  if (any(allowed)) {
    weights[allowed] <- exp(utility[allowed] - max(utility[allowed]))
    weights <- weights / sum(weights)
  }
  weights
}

# Use an exact native path for the common linear/raw quadratic specification.
# All other formulas retain model.matrix semantics through an R evaluator.
prepare_matching_model <- function(matching, population_groups = NULL) {
  validate_matching_groups(matching, population_groups)
  beta <- matching$coefficients
  columns <- gsub("[[:space:]]", "", names(beta))
  native <- all(columns %in% c("age_difference", "I(age_difference^2)"))
  if ("I(age_difference^2)" %in% columns) {
    native <- native && identical(
      get("I", envir = environment(matching$formula), inherits = TRUE), base::I)
  }
  quadratic <- numeric()
  if (native) {
    quadratic <- c(sum(beta[columns == "age_difference"]),
                   sum(beta[columns == "I(age_difference^2)"]))
  }
  pairs <- expand_group_exogamy(matching$group_exogamy)
  offsets <- matrix(numeric(), 0, 0)
  labels <- as.character(population_groups)
  if (nrow(pairs)) {
    offsets <- matrix(0, length(labels), length(labels))
    for (i in seq_len(nrow(pairs))) {
      female <- match(pairs$female_group[i], labels)
      male <- match(pairs$male_group[i], labels)
      if (!is.na(female) && !is.na(male)) offsets[female, male] <- pairs$log_odds[i]
    }
  }
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  blocked <- matrix(FALSE, 0, 0)
  if (nrow(restrictions$pairs)) {
    blocked <- matrix(FALSE, length(labels), length(labels))
    for (i in seq_len(nrow(restrictions$pairs))) {
      female <- match(restrictions$pairs$female_group[i], labels)
      male <- match(restrictions$pairs$male_group[i], labels)
      if (!is.na(female) && !is.na(male)) blocked[female, male] <- TRUE
    }
  }
  group_formula <- matching_has_group_formula(matching)
  list(delta = unname(vapply(labels, function(g) {
         if (g %in% names(matching$delta)) matching$delta[[g]] else 0
       }, numeric(1))),
       prevent_exogamy_all = restrictions$all, blocked_pairs = blocked,
quadratic = quadratic, has_squared = "I(age_difference^2)" %in% columns,
       group_formula = group_formula, needs_groups = group_formula || nrow(pairs) > 0,
       pair_offsets = offsets,
       evaluate = function(age_difference, female_group = NULL, male_group = NULL) {
         calculate_matching_utilities(matching, age_difference,
           if (is.null(female_group)) NULL else labels[female_group + 1L],
           if (is.null(male_group)) NULL else labels[male_group + 1L], add_group_exogamy = FALSE)
       })
}

#' Specify a union-formation process
#'
#' Construct or attach a union-formation specification without executing it.
#' In on-demand mode, at a birth use
#' the woman's current partner, or choose a partner from all living, unpartnered
#' men using conditional-logit probabilities from the matching formula to form
#' a union. Men below `min_male_age` are excluded. If no eligible man exists,
#' cancel the birth by default; `cancel_birth_if_no_match = FALSE`
#' preserves it with unknown father. Known kin within
#' `kinship_generations` on both sides of a common ancestor are excluded.
#' Each person may have only one ongoing union. Partnered men are not
#' available; surviving husbands return to the pool at widowhood. Death ends affected
#' unions. In rate mode unpartnered women face competing search events and
#' use their unpartnered fertility schedule and may give birth without initiating
#' a search. Failed searches draw another waiting time; formation and widowhood
#' switch fertility schedules while preserving spacing. An overall fertility
#' schedule applies to both statuses, with a simulation warning. On-demand mode
#' requires an overall fertility schedule. Unpartnered rate-mode births use the
#' same eligibility and exact draw to assign a father, using
#' `father_matching` when supplied or otherwise `matching`, without
#' forming a union or removing him from availability. A search with no
#' eligible man cancels the birth by default. No previous-father preference
#' is applied.
#'
#' Partners are drawn exactly from the full conditional-logit distribution over
#' all eligible men, without a candidate sample. Living, unpartnered men are
#' grouped by group and quarter-year birth bin; a bin is chosen in proportion
#' to its size times the largest weight any member could have, a member is
#' accepted with the ratio of his weight to that maximum, and only he is
#' screened for age and kinship. For the built-in linear/quadratic age model
#' this is exact. Custom R formulas are bounded on a 1/32-year grid; if a drawn
#' man's utility exceeds the bound, simulation warns that draws were
#' approximate. Existing unions are retained even if a new-search restriction
#' would exclude the partner.
#'
#' Matching always uses calibrated male-group weights: a woman chooses among
#' screened candidates with each man's selection weight multiplied by a
#' weight for his group. Weights are recalculated at every search, in continuous time,
#' so that expected male union flows by group match a mirrored target: each
#' unpartnered man has the per-capita intensity of an unpartnered woman of his
#' group, adjusted by the matching model's `delta`. Women's
#' intensity is the search hazard for rate-mode searches and the unpartnered
#' fertility hazard (zero during birth spacing) for on-demand searches and
#' nonmarital father selection, which has its own weights and uses the father
#' model's preferences and delta when one is supplied. Weights depend only
#' on the man's group, so they do not change the association between partners'
#' characteristics in the static group-only illustration. With one group,
#' all group weights are one. When no
#' exact weights exist, approximate weights are used with a warning.
#'
#' @param initiation Either `"on_demand"` for birth-triggered searches or
#'   `"rate"` for scheduled searches. No default is selected.
#' @param matching A specification created by [create_matching_model()], or a
#'   tibble with numeric `time_start`, `time_end`, and list-column `model` holding
#'   complete matching objects. Model intervals are half-open absolute years,
#'   cannot overlap, and must cover the requested run. The active model is used
#'   at partner selection, without changing pending search dates or existing
#'   unions. Base data frames with list columns are also accepted. Defaults
#'   to `create_age_gap_matching_model(0, 5)`: preferred age gap zero years,
#'   with selection weight halving five years on either side.
#' @param cancel_birth_if_no_match Logical, default `TRUE`. Cancel a birth when
#'   its on-demand partner search or rate-mode nonmarital father selection fails.
#'   Set `FALSE` to allow births with unknown fathers. Historical missing father
#'   IDs are retained. A cancellation starts no new spacing interval and draws
#'   fresh waiting times from the failed attempt. Partner scarcity can reduce
#'   realized fertility below inputs.
#' @param father_matching Optional matching specification for fathers of
#'   unpartnered births in rate mode, created by [create_matching_model()], or a
#'   `time_start`/`time_end`/`model` schedule as for `matching`.
#'   Defaults to `NULL`, using the general `matching` model active at that birth. Union
#'   searches always use `matching`. Not used by on-demand initiation. Process-level
#'   age and kinship screening settings are shared; NULL inherits the general
#'   model's settings.
#' @param min_male_age Minimum male age in continuous years at a new partner
#'   search, default `15`. Finite nonnegative scalar; zero disables this age
#'   exclusion. Existing unions are retained. No minimum is added for actual
#'   female event eligibility. Weighted choice counts unpartnered men and women
#'   at or above this same cutoff when calibrating its weights.
#' @param kinship_generations Generation bound on both sides of a known common
#'   ancestor, default `2`. Excludes parents, grandparents, siblings, aunts/uncles,
#'   nieces/nephews and first cousins, including half relationships, when known.
#'   Zero imposes no cross-person kin exclusion. Missing ancestry is not proof of
#'   unrelatedness; only demonstrated relationships are screened. Existing unions
#'   are retained. See [extract_kin()] for the shared distance definitions.
#' @section Calendar rates:
#' Search `rates` and `union_convert_rate` optionally accept `time_start` and
#' `time_end` for half-open absolute calendar-year bands, including fractional
#' years. Both columns are required together, with finite starts and greater
#' ends (possibly Inf). Every active slice must satisfy its existing rate rules
#' and cover the requested run. Clock choices stay consistent across periods;
#' changes precede same-time events and retain duration origins. `cohab_prob`
#' accepts the same calendar columns, with probabilities selected at the search
#' instant. Probability-only changes refresh matching calibration without
#' resampling pending event clocks. All three matching arguments also accept
#' whole-model schedule tibbles with `time_start`, `time_end`, and list-column
#' `model`. Model changes update role-specific preferences and calibration before
#' same-time choices without resampling event clocks. All supplied models,
#' including future intervals, are validated; gaps in a run are errors. NULL
#' father/cohabitation matching dynamically follows general matching. Other
#' non-rate arguments remain static.
#'
#' @param rates For `"rate"` initiation, an age-specific female search-hazard
#'   table with `age_start`, `age_end`, and `rate`, in years and per person-year.
#'   These specify search attempts, not guaranteed realized union formation.
#'   An optional `group` column selects the woman's group. Without it, all groups
#'   share the schedule. Group does not restrict partner matching. With no
#'   population groups, simulation warns and selects the first alphabetical group.
#'   See [mortality()] for group-label and rate-coverage rules.
#'   Nonoverlapping bands may have gaps or an open end; gaps and ages outside
#'   support have zero hazard. Optional `previous_union` stratifies by `"none"`,
#'   `"dissolved"`, and `"widowed"`. Each group requires `"none"` and
#'   `"dissolved"`; omitted `"widowed"` uses the complete dissolved schedule
#'   in that group and active calendar period. Explicit widowed bands override
#'   the fallback as a whole, including zero hazard in their gaps.
#'   Without it, the schedule applies to all previous union types. Categories
#'   describe the most recent recorded episode; complete histories are required.
#'   Optionally include both `dur_start` and `dur_end` for elapsed years since
#'   the most recent union ended. Each row supplies either age or duration bands
#'   with both unused boundaries NA. `none` must use age; each other category
#'   uses one clock consistently across all its bands and groups. Age and duration
#'   cannot apply together. All-age inputs may omit duration columns.
#'   Must be NULL for `"on_demand"`.
#' @param diagnostics Logical, default `FALSE`. When the general matching model
#'   uses weighted choice in a grouped population, record `.log_weight_<group>`
#'   columns on newly formed unions: the log group
#'   weights applied at the search that formed each union (0 for the reference
#'   group and when calibration cannot affect choices). Use the column for each alternative's
#'   group as a conditional-logit offset in wife-side choice sets.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   search-rate multipliers by the woman's group. Missing groups use 1.0;
#'   zero suppresses searches. Applies to shared or group-specific input rates.
#'   Without population groups, simulation warns and ignores multipliers.
#'   Requires `initiation = "rate"`; on-demand formation has no search hazard.
#'   Matching utilities and input rate tables are unchanged.
#' @param cohab_prob Optional age-specific probability table with `age_start`,
#'   `age_end`, `prob`, optionally `group` and `previous_union`. Enables typed
#'   unions only in rate mode. Age gaps mean zero. Optional `time_start` and
#'   `time_end` supply half-open calendar intervals; each active slice must satisfy
#'   the same age/group/previous-category rules and cover the run. Type uses the
#'   active probability at the search instant, before matching. Existing unions
#'   keep their types. NULL retains untyped partnerships and their RNG path.
#' @param cohab_matching Optional cohabitation matching model or whole-model
#'   schedule; NULL shares the general `matching` model active at that search. Marriage, cohabitation and unpartnered father searches calibrate
#'   weights using their own intensities and share the available male pool.
#' @param union_convert_rate Optional age- or duration-specific conversion table,
#'   optionally by woman's group. Duration uses the current cohabitation start.
#'   Uses `rate` and either age or duration boundaries, not both clocks, with no
#'   previous-union axis. NULL gives no conversions and warns at simulation.
#'   Conversion keeps partners, closes cohab with reason `conversion`, and starts
#'   a linked marriage with a new dissolution-duration origin.
#' @param cohabitation_counts_as_previous_union Logical, default FALSE. FALSE
#'   restores previous legal-marriage category after a cohabitation ends; TRUE
#'   classifies cohab breakup/death as dissolved/widowed. Search duration still
#'   starts at the latest partnership ending. Only applies to typed unions.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model. Rate tables are tibbles.
#' @export
union_formation <- function(initiation, matching = create_age_gap_matching_model(0, 5),
                            rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE,
                            diagnostics = FALSE, group_prop_hazard = NULL,
                            cohab_prob = NULL, cohab_matching = NULL, union_convert_rate = NULL,
                            cohabitation_counts_as_previous_union = FALSE) {
  process <- structure(list(
    name = "union_formation", initiation = initiation, matching = matching,
    min_male_age = min_male_age,
    kinship_generations = kinship_generations, father_matching = father_matching,
    cancel_birth_if_no_match = cancel_birth_if_no_match, diagnostics = diagnostics,
    rates = rates, rate_unit = "person_year", group_prop_hazard = group_prop_hazard,
    cohab_prob = cohab_prob, cohab_matching = cohab_matching, union_convert_rate = union_convert_rate,
    cohabitation_counts_as_previous_union = cohabitation_counts_as_previous_union
  ), class = c("hammel_union_formation", "hammel_process"))
  validate_union_formation(process)
  for (field in c("matching", "father_matching", "cohab_matching"))
    process[field] <- list(normalize_model_argument(process[[field]]))
  for (field in c("cohab_prob", "union_convert_rate"))
    if (!is.null(process[[field]])) process[[field]] <- tibble::as_tibble(process[[field]])
  if (!is.null(rates)) {
    process$rates <- tibble::as_tibble(rates)
    if (!is.null(attr(rates, "source"))) process$source <- attr(rates, "source")
  }
  process
}

#' @rdname union_formation
#' @export
set_union_formation <- function(model, initiation, matching = create_age_gap_matching_model(0, 5),
                            rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE,
                            diagnostics = FALSE, group_prop_hazard = NULL,
                            cohab_prob = NULL, cohab_matching = NULL, union_convert_rate = NULL,
                            cohabitation_counts_as_previous_union = FALSE) {
  set_process(model, union_formation(initiation, matching, rates,
                                     min_male_age, kinship_generations, father_matching,
                                     cancel_birth_if_no_match, diagnostics, group_prop_hazard,
                                     cohab_prob, cohab_matching, union_convert_rate,
                                     cohabitation_counts_as_previous_union))
}

validate_union_formation <- function(process) {
  validate_group_prop_hazard(process$group_prop_hazard)
  validate_cohabitation(process)
  if (!identical(process$name, "union_formation") ||
      !identical(process$rate_unit, "person_year")) {
    stop("Invalid union-formation specification or units.", call. = FALSE)
  }
  if (!is.character(process$initiation) || length(process$initiation) != 1L ||
      is.na(process$initiation) || !process$initiation %in% c("rate", "on_demand")) {
    stop('`initiation` must be "rate" or "on_demand".', call. = FALSE)
  }
  age <- process$min_male_age
  if (!is.numeric(age) || length(age) != 1L || !is.null(dim(age)) ||
      !is.finite(age) || age < 0) {
    stop("`min_male_age` must be a finite nonnegative number of years.", call. = FALSE)
  }
  diagnostic <- process$diagnostics
  if (!is.logical(diagnostic) || length(diagnostic) != 1L ||
      !is.null(dim(diagnostic)) || is.na(diagnostic)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  cancel <- process$cancel_birth_if_no_match
  if (!is.logical(cancel) || length(cancel) != 1L || !is.null(dim(cancel)) || is.na(cancel)) {
    stop("`cancel_birth_if_no_match` must be TRUE or FALSE.", call. = FALSE)
  }
  validate_kinship_generations(process$kinship_generations)
  validate_model_argument(process$matching, validate_matching_model, "matching")
  if (!is.null(process$father_matching))
    validate_model_argument(process$father_matching, validate_matching_model, "father_matching")
  if (process$initiation == "on_demand") {
    if (!is.null(process$group_prop_hazard)) {
      stop("`group_prop_hazard` requires rate-based union initiation; on-demand formation has no search hazard.",
           call. = FALSE)
    }
    if (!is.null(process$rates)) {
      stop("On-demand union formation does not take rates.", call. = FALSE)
    }
  } else {
    validate_union_search_rates(process$rates)
  }
  invisible(process)
}


expand_widowed_search_rates <- function(rates) {
  if (!"previous_union" %in% names(rates)) return(rates)
  changed <- FALSE
  expand_slice <- function(x) {
    copies <- lapply(split_rate_groups(x), function(rows) {
      if ("widowed" %in% x$previous_union[rows]) return(NULL)
      copy <- x[rows[x$previous_union[rows] == "dissolved"], , drop = FALSE]
      if (nrow(copy)) changed <<- TRUE
      copy$previous_union <- rep("widowed", nrow(copy))
      copy
    })
    do.call(rbind, c(list(x), copies))
  }
  if (!has_calendar_rates(rates)) return(expand_slice(rates))
  times <- sort(unique(c(rates$time_start, rates$time_end)))
  pieces <- lapply(seq_len(length(times) - 1L), function(i) {
    rows <- rates$time_start <= times[i] & times[i] < rates$time_end
    x <- rates[rows, , drop = FALSE]
    if (!nrow(x)) return(NULL)
    x$time_start <- times[i]
    x$time_end <- times[i + 1L]
    expand_slice(x)
  })
  if (changed) do.call(rbind, pieces) else rates
}

validate_union_search_rates <- function(rates, widow_fallback = TRUE) {
    if (has_calendar_rates(rates)) {
      validate_calendar_rates(rates, function(x) validate_union_search_rates(x, widow_fallback))
      effective <- if (widow_fallback) expand_widowed_search_rates(rates) else rates
      if (all(c("previous_union", "dur_start") %in% names(effective)) &&
          any(vapply(split(!is.na(effective$dur_start), effective$previous_union),
            function(x) length(unique(x)) != 1L, logical(1))))
        stop("Each previous_union category must use the same age or duration clock across calendar periods.", call. = FALSE)
      return(invisible(rates))
    }
    check_columns(rates, c("age_start", "age_end", "rate"), "rates")
    check_numeric_column(rates$rate, "rates$rate")
    if (!nrow(rates) || any(rates$rate < 0) ||
        any(c("sex", "time_start", "time_end") %in% names(rates)))
      stop("Use nonnegative female search rates in valid age or duration bands, without sex or calendar columns.",
           call. = FALSE)
    duration_columns <- c("dur_start", "dur_end") %in% names(rates)
    if (any(duration_columns) && !all(duration_columns))
      stop("Supply both `dur_start` and `dur_end`, or neither.", call. = FALSE)
    boundaries <- c("age_start", "age_end", if (all(duration_columns)) c("dur_start", "dur_end"))
    for (column in boundaries) {
      x <- rates[[column]]
      # Bare NA is logical in R; permit all-NA placeholders for an unused axis.
      if (!(is.numeric(x) || (is.logical(x) && all(is.na(x)))) ||
          !is.null(dim(x)) || any(is.nan(x)))
        stop(sprintf("`rates$%s` must contain numeric boundaries or NA for an unused clock.", column),
             call. = FALSE)
    }
    age <- !is.na(rates$age_start) & !is.na(rates$age_end)
    duration <- if (all(duration_columns))
      !is.na(rates$dur_start) & !is.na(rates$dur_end) else rep(FALSE, nrow(rates))
    age_empty <- is.na(rates$age_start) & is.na(rates$age_end)
    duration_empty <- if (all(duration_columns))
      is.na(rates$dur_start) & is.na(rates$dur_end) else rep(TRUE, nrow(rates))
    if (any(!((age & duration_empty) | (duration & age_empty))))
      stop("Each row requires exactly one complete age or duration band; unused boundaries must both be NA.",
           call. = FALSE)
    starts <- rates$age_start
    ends <- rates$age_end
    if (any(duration)) {
      if (!"previous_union" %in% names(rates))
        stop("Duration search bands require `previous_union` ('dissolved' or 'widowed').", call. = FALSE)
      starts[duration] <- rates$dur_start[duration]
      ends[duration] <- rates$dur_end[duration]
    }
    if (any(!is.finite(starts)) || any(starts < 0) || any(ends <= starts))
      stop("Use valid nonnegative age or duration bands with finite starts and ends greater than starts.",
           call. = FALSE)
    strata <- split_rate_groups(rates)
    if ("previous_union" %in% names(rates)) {
      previous <- rates$previous_union
      categories <- c("none", "dissolved", "widowed")
      if (!is.character(previous) || !is.null(dim(previous)) || anyNA(previous) ||
          any(!previous %in% categories))
        stop("`previous_union` must contain 'none', 'dissolved', or 'widowed'.", call. = FALSE)
      if (any(vapply(strata, function(rows)
          !(if (widow_fallback) all(c("none", "dissolved") %in% previous[rows]) else
            setequal(previous[rows], categories)), logical(1))))
        stop("Union-search schedules require 'none' and 'dissolved' (or all three `previous_union` categories) within each group; cohab_prob requires all three.",
             call. = FALSE)
      if (any(duration & previous == "none"))
        stop("`previous_union = 'none'` requires age-specific bands.", call. = FALSE)
      effective <- if (widow_fallback) expand_widowed_search_rates(rates) else rates
      effective_duration <- if (all(duration_columns)) !is.na(effective$dur_start) else
        rep(FALSE, nrow(effective))
      if (any(vapply(split(effective_duration, effective$previous_union), function(x) length(unique(x)) != 1L, logical(1))))
        stop("Each `previous_union` category must use the same age or duration clock for all bands and groups.",
             call. = FALSE)
      strata <- unlist(lapply(strata, function(rows) split(rows, previous[rows])), recursive = FALSE)
    }
    for (rows in strata) {
      rows <- rows[order(starts[rows])]
      if (length(rows) > 1L &&
          any(tail(starts[rows], -1) < head(cummax(ends[rows]), -1)))
        stop("Union-search age or duration bands must not overlap.", call. = FALSE)
    }
  invisible(rates)
}
