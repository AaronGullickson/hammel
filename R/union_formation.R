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
#' @param group_exogamy Optional finite named numeric vector: `"A:B"` means
#'   female A/male B; `"A::B"` supplies the same effect in both directions.
#'   Omitted pairs have odds ratio 1.0 (zero log effect); within-group entries are
#'   allowed. Overlapping
#'   directed pairs after expansion are errors. Alternatively provide a data
#'   frame with `female_group`, `male_group`, and `log_odds`, one directed pair
#'   per row, including for group labels containing colons. Effects add to final
#'   formula scores; they do not create a matching prohibition.
#' @param group_levels Optional ordered vector of group labels (character or
#'   numeric codes). Required with at least two levels when the formula uses
#'   `female_group` or `male_group`; must cover population groups. Sets fixed
#'   treatment reference/columns, independent of candidate samples. Not required
#'   for `same_group` or pair effects alone.
#' @param prevent_exogamy `FALSE` (default) for no group restrictions; `TRUE`
#'   excludes all different-group candidates. A character vector such as
#'   `c("A::B", "A:C")` excludes selected pairs with female-first `:` and
#'   symmetric `::` notation. Alternatively supply a directional data frame
#'   with `female_group` and `male_group`. Pairs must join different groups.
#'   Restrictions filter the original sample alongside age/kin eligibility,
#'   without refilling. Empty samples follow existing failed-search and birth
#'   cancellation rules. Restrictions override finite preferences, require
#'   population identities, and do not invalidate existing partners.
#' @param sample_size Nonnegative whole number of men to sample before
#'   pair-specific screening. Defaults to `250`; zero requests the whole
#'   available pool. Rejected candidates are not replaced. A separate
#'   `father_matching` model uses its own sample size.
#' @param imbalance_absorption Used only by `matching_method = "weighted_choice"`
#'   in [union_formation()]. Share, in `[0, 1]`, of a group's sex imbalance among
#'   unpartnered people that is absorbed through cross-group unions. Default `1`
#'   gives each unpartnered man the per-capita union intensity of an unpartnered
#'   woman of his group. `0` makes each male group's share of unions equal its
#'   female group's share, so the excess sex waits longer. Intermediate values
#'   interpolate the male group targets. Symmetric preferences with balanced sex
#'   composition give equal directional pairing counts at every value.
#' @return An inspectable `hammel_matching_model` list containing `formula`
#'   and `coefficients`. This specifies a model; it does not estimate one.
#' @export
#' @examples
#' # Illustrative coefficients, not demographic defaults.
#' matching <- create_matching_model(
#'   ~ age_difference + I(age_difference^2),
#'   c(age_difference = 0.1, "I(age_difference^2)" = -0.02)
#' )
create_matching_model <- function(formula, coefficients, group_exogamy = NULL, group_levels = NULL, prevent_exogamy = FALSE, sample_size = 250,
                                  imbalance_absorption = 1) {
  model <- structure(list(formula = formula, coefficients = coefficients,
                          group_exogamy = group_exogamy, group_levels = group_levels,
                          prevent_exogamy = prevent_exogamy, sample_size = sample_size,
                          imbalance_absorption = imbalance_absorption),
                     class = "hammel_matching_model")
  validate_matching_model(model)
  model
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
#' @param prevent_exogamy Optional post-sampling group restrictions passed to
#'   [create_matching_model()]. Defaults to `FALSE`.
#' @param half_weight_distance Finite positive numeric scalar in years from the
#'   preferred gap at which selection weight halves.
#' @inheritParams create_matching_model
#' @return A `hammel_matching_model` with linear and quadratic age-difference terms.
#' @export
#' @examples
#' # Illustrative preferences, not demographic defaults.
#' matching <- create_age_gap_matching_model(3, 5)
create_age_gap_matching_model <- function(preferred_age_gap, half_weight_distance, group_exogamy = NULL, prevent_exogamy = FALSE, sample_size = 250,
                                          imbalance_absorption = 1) {
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
    sample_size = sample_size, imbalance_absorption = imbalance_absorption)
}

# Older matching objects lack the field and use the default.
matching_imbalance_absorption <- function(model) {
  if ("imbalance_absorption" %in% names(model)) model$imbalance_absorption else 1
}

validate_imbalance_absorption <- function(x) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.null(dim(x)) ||
      !is.finite(x) || x < 0 || x > 1) {
    stop("`imbalance_absorption` must be a finite number between 0 and 1.", call. = FALSE)
  }
  invisible(x)
}

validate_sample_size <- function(n) {
  if (!is.numeric(n) || is.complex(n) || length(n) != 1L || !is.null(dim(n)) ||
      !is.finite(n) || n < 0 || n != floor(n) || n > .Machine$integer.max) {
    stop("`sample_size` must be a nonnegative whole number within integer range.",
         call. = FALSE)
  }
  invisible(NULL)
}

validate_matching_model <- function(model) {
  if (!inherits(model, "hammel_matching_model") || !is.list(model)) {
    stop("`matching` must be created with create_matching_model().", call. = FALSE)
  }
  validate_sample_size(model$sample_size)
  validate_imbalance_absorption(matching_imbalance_absorption(model))
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
prepare_matching_model <- function(matching, population_groups = NULL, two_sided = FALSE,
                                   weighted = FALSE) {
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
  list(sample_size = as.integer(matching$sample_size),
       two_sided = two_sided, weighted = weighted,
       imbalance_absorption = as.numeric(matching_imbalance_absorption(matching)),
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
#' the woman's current partner, or sample living, unpartnered men and choose a partner using
#' conditional-logit probabilities from the matching formula to form a union.
#' Sampled men below `min_male_age` are excluded. If no eligible sampled men
#' remain, cancel the birth by default; `cancel_birth_if_no_match = FALSE`
#' preserves it with unknown father. Known kin within
#' `kinship_generations` on both sides of a common ancestor are excluded.
#' Each person may have only one ongoing union. Partnered men are excluded
#' before sampling; surviving husbands return to the pool at widowhood. Death ends affected
#' unions. In rate mode unpartnered women face competing search events and
#' use their unpartnered fertility schedule and may give birth without initiating
#' a search. Failed searches draw another waiting time; formation and widowhood
#' switch fertility schedules while preserving spacing. An overall fertility
#' schedule applies to both statuses, with a simulation warning. On-demand mode
#' requires an overall fertility schedule. Unpartnered rate-mode births use the
#' same candidate sampling and eligibility to assign a father, using
#' `father_matching` when supplied or otherwise `matching`, without
#' forming a union or removing him from availability. Empty screened pools
#' cancel births by default. No previous-father preference is applied.
#'
#' Sampling is from the available male pool BEFORE pair-specific screening.
#' The matching model defaults to `sample_size = 250`; zero means the whole pool.
#' Otherwise sample without replacement
#' up to that many men. Age screening occurs afterward, without refilling the
#' sample. Kinship screening follows the same ordering. Existing unions are
#' retained even if a new-search restriction would exclude the partner.
#' Candidate sets will not be retained in simulation outputs.
#'
#' Explicit `matching_method = "weighted_choice"` keeps female choice among the
#' screened candidates but multiplies each man's selection weight by a weight
#' for his group. Weights are recalculated at every search, in continuous time,
#' so that expected male union flows by group match a mirrored target: each
#' unpartnered man has the per-capita intensity of an unpartnered woman of his
#' group, adjusted by the matching model's `imbalance_absorption`. Women's
#' intensity is the search hazard for rate-mode searches and the unpartnered
#' fertility hazard (zero during birth spacing) for on-demand searches and
#' nonmarital father selection, which has its own weights. Weights depend only
#' on the man's group, so they do not change the association between partners'
#' characteristics. With one group the method equals female choice. When no
#' exact weights exist, approximate weights are used with a warning.
#'
#' Default `matching_method = "two_sided"` adds alternate women to a
#' temporary market after ego's male eligibility screening. Each fixed female
#' row is scored once, including for R formulas; row and column probabilities
#' are combined for one partner draw.
#' Hard group exclusions apply to hypothetical pairs; kinship applies only to
#' ego. If fewer alternate women exist, randomly reduce the eligible men to
#' keep the initial market square. Only ego's final outcome changes real state.
#'
#' @param initiation Either `"on_demand"` for birth-triggered searches or
#'   `"rate"` for scheduled searches. No default is selected.
#' @param matching A specification created by [create_matching_model()]. Defaults
#'   to `create_age_gap_matching_model(0, 5)`: preferred age gap zero years,
#'   with selection weight halving five years on either side.
#' @param cancel_birth_if_no_match Logical, default `TRUE`. Cancel a birth when
#'   its on-demand partner search or rate-mode nonmarital father selection fails.
#'   Set `FALSE` to allow births with unknown fathers. Historical missing father
#'   IDs are retained. A cancellation starts no new spacing interval and draws
#'   fresh waiting times from the failed attempt. Partner scarcity can reduce
#'   realized fertility below inputs.
#' @param father_matching Optional matching specification for fathers of
#'   unpartnered births in rate mode, created by [create_matching_model()].
#'   Defaults to `NULL`, using the current general `matching` model. Union
#'   searches always use `matching`. Not used by on-demand initiation. Process-level
#'   age and kinship screening settings are shared; each matching model supplies its own
#'   candidate sample size. NULL inherits the general model's sample size.
#' @param min_male_age Minimum male age in continuous years at a new partner
#'   search, default `15`. Finite nonnegative scalar; zero disables this age
#'   exclusion. Existing unions are retained. No minimum is added for actual
#'   female event eligibility. Two-sided alternate women use this same cutoff.
#' @param kinship_generations Generation bound on both sides of a known common
#'   ancestor, default `2`. Excludes parents, grandparents, siblings, aunts/uncles,
#'   nieces/nephews and first cousins, including half relationships, when known.
#'   Zero imposes no cross-person kin exclusion. Missing ancestry is not proof of
#'   unrelatedness; only demonstrated relationships are screened. Existing unions
#'   are retained. See [extract_kin()] for the shared distance definitions.
#' @param rates For `"rate"` initiation, an age-specific female search-hazard
#'   table with `age_start`, `age_end`, and `rate`, in years and per person-year.
#'   These specify search attempts, not guaranteed realized union formation.
#'   An optional `group` column selects the woman's group. Without it, all groups
#'   share the schedule. Group does not restrict partner matching. With no
#'   population groups, simulation warns and selects the first alphabetical group.
#'   See [mortality()] for group-label and rate-coverage rules.
#'   Nonoverlapping bands may have gaps or an open end; gaps and ages outside
#'   support have zero hazard. Must be NULL for `"on_demand"`.
#' @param diagnostics Logical, default `FALSE`. Record `.n_eligible` (sampled men
#'   remaining after screening) and `.utility_rank` (1 = highest utility, ties
#'   share the minimum rank) on newly formed unions. Historical unknowns are NA.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   search-rate multipliers by the woman's group. Missing groups use 1.0;
#'   zero suppresses searches. Applies to shared or group-specific input rates.
#'   Without population groups, simulation warns and ignores multipliers.
#'   Requires `initiation = "rate"`; on-demand formation has no search hazard.
#'   Matching utilities and input rate tables are unchanged.
#' @param matching_method `"two_sided"` (default) for fixed-pool two-sided
#'   probabilities, `"female_choice"` for ego-only choice, or experimental
#'   `"weighted_choice"` for female choice with calibrated male-group weights.
#'   Female choice accounts for available men but not reciprocal competition
#'   among women. With unequal group sizes, symmetric group preferences can
#'   produce unequal female-A/male-B and female-B/male-A pairing counts.
#'   These imbalances can systematically affect group inheritance and subsequent
#'   group composition, especially when maternal and paternal inheritance differ.
#'   Two-sided matching samples alternate unpartnered women with the male age cutoff,
#'   combines ego-choice and male-choice probabilities without removing pairs.
#'   It applies to union formation and nonmarital father selection. Temporary
#'   pairs never change population state. Aggregate parameter recovery is under
#'   investigation; the method can cost quadratic work in candidate sample size.
#'   Older process objects without this field retain female choice. New calls
#'   omitting the argument use two-sided matching, changing seeded histories
#'   relative to the previous default.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model. Rate tables are tibbles.
#' @export
union_formation <- function(initiation, matching = create_age_gap_matching_model(0, 5),
                            rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE,
                            diagnostics = FALSE, group_prop_hazard = NULL,
                            matching_method = "two_sided") {
  process <- structure(list(
    name = "union_formation", initiation = initiation, matching = matching,
    min_male_age = min_male_age,
    kinship_generations = kinship_generations, father_matching = father_matching,
    cancel_birth_if_no_match = cancel_birth_if_no_match, diagnostics = diagnostics,
    rates = rates, rate_unit = "person_year", group_prop_hazard = group_prop_hazard,
    sampling_stage = "before_pair_screening", matching_method = matching_method
  ), class = c("hammel_union_formation", "hammel_process"))
  validate_union_formation(process)
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
                            matching_method = "two_sided") {
  set_process(model, union_formation(initiation, matching, rates,
                                     min_male_age, kinship_generations, father_matching,
                                     cancel_birth_if_no_match, diagnostics, group_prop_hazard, matching_method))
}

validate_union_formation <- function(process) {
  method <- if ("matching_method" %in% names(process)) process$matching_method else "female_choice"
  if (!is.character(method) || length(method) != 1L || !is.null(dim(method)) ||
      is.na(method) || !method %in% c("female_choice", "two_sided", "weighted_choice")) {
    stop('`matching_method` must be "female_choice", "two_sided", or "weighted_choice".', call. = FALSE)
  }
  validate_group_prop_hazard(process$group_prop_hazard)
  if (!identical(process$name, "union_formation") ||
      !identical(process$rate_unit, "person_year") ||
      !identical(process$sampling_stage, "before_pair_screening")) {
    stop("Invalid union-formation specification or units.", call. = FALSE)
  }
  if (!is.character(process$initiation) || length(process$initiation) != 1L ||
      is.na(process$initiation) || !process$initiation %in% c("rate", "on_demand")) {
    stop('`initiation` must be "rate" or "on_demand".', call. = FALSE)
  }
  if ("sample_size" %in% names(process))
    stop("Set `sample_size` on the matching model, not union formation; rebuild this process.",
         call. = FALSE)
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
  validate_matching_model(process$matching)
  if (!is.null(process$father_matching)) validate_matching_model(process$father_matching)
  if (process$initiation == "on_demand") {
    if (!is.null(process$group_prop_hazard)) {
      stop("`group_prop_hazard` requires rate-based union initiation; on-demand formation has no search hazard.",
           call. = FALSE)
    }
    if (!is.null(process$rates)) {
      stop("On-demand union formation does not take rates.", call. = FALSE)
    }
  } else {
    rates <- process$rates
    check_columns(rates, c("age_start", "age_end", "rate"), "rates")
    check_numeric_column(rates$age_start, "rates$age_start")
    check_numeric_column(rates$rate, "rates$rate")
    ends <- rates$age_end
    if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) ||
        anyNA(ends) || any(ends <= rates$age_start) ||
        any(rates$age_start < 0) || any(rates$rate < 0) ||
        any(c("sex", "time_start", "time_end") %in% names(rates))) {
      stop("Use nonnegative female search rates in valid age bands, without sex or calendar columns.",
           call. = FALSE)
    }
    for (rows in split_rate_groups(rates)) {
      ordered <- rates[rows, , drop = FALSE]
      ordered <- ordered[order(ordered$age_start), , drop = FALSE]
      if (nrow(ordered) > 1L &&
          any(tail(ordered$age_start, -1) < head(cummax(ordered$age_end), -1))) {
        stop("Union-search age bands must not overlap.", call. = FALSE)
      }
    }
  }
  invisible(process)
}
