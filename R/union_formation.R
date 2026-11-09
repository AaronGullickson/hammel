# Construct and validate union-formation processes and search rate schedules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

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
#' a search. Failed searches leave the woman searching; formation and widowhood
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
#' `union_formation()` creates the process specification;
#' `set_union_formation()` creates it and attaches it to a model through
#' [set_process()]. Neither function runs a simulation.
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
#'   IDs are retained. A cancellation starts no new spacing interval, and the
#'   attempt is not replaced. Partner scarcity can reduce
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
#' changing event clocks. All three matching arguments also accept
#' whole-model schedule tibbles with `time_start`, `time_end`, and list-column
#' `model`. Model changes update role-specific preferences and calibration before
#' same-time choices without changing event clocks. All supplied models,
#' including future intervals, are validated; gaps in a run are errors. NULL
#' father/cohabitation matching dynamically follows general matching. Other
#' non-rate arguments remain static.
#'
#' @param rates For `"rate"` initiation, an age-specific female search-hazard
#'   table with `age_start`, `age_end`, and `rate`, in years and per person-year.
#'   These specify search attempts, not guaranteed realized union formation.
#'   Columns named by declared group dimensions select the woman's categories;
#'   without them, all groups share the schedule. Groups do not restrict
#'   partner matching.
#'   See [mortality()] for group-label and rate-coverage rules.
#'   Nonoverlapping bands may have gaps or an open end; gaps and ages outside
#'   support have zero hazard. Optional `previous_union` supplies either
#'   `"none"`, `"separated"`, `"widowed"`, or pooled `"none"`, `"dissolved"`.
#'   Dissolved pools separation and widowhood. Each group and active calendar
#'   slice must supply one complete set; mixing the sets is not allowed.
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
#'   optionally by the partners' categories (`<group>_w`, `<group>_h` for a
#'   declared dimension). Duration uses the current cohabitation start.
#'   Uses `rate` and either age or duration boundaries, not both clocks, with no
#'   previous-union axis. NULL gives no conversions and warns at simulation.
#'   Conversion keeps partners, closes cohab with reason `conversion`, and starts
#'   a linked marriage with a new dissolution-duration origin.
#' @param weights NULL (default) for unweighted selection of searching women,
#'   or a [create_weight_model()] using person variables, which chooses which
#'   unpartnered woman in a cell searches without changing the cell's input
#'   search rate. Rate mode only. With weighted matching calibration the
#'   weights are kept current for every searching woman, so calibration uses
#'   the weighted search intensities.
#' @param weights_conversion NULL (default) for unweighted selection, or a
#'   [create_weight_model()] using union variables, which chooses which
#'   cohabitation in a cell converts without changing the cell's input rate.
#' @param cohabitation_counts_as_previous_union Logical, default FALSE. FALSE
#'   restores previous legal-marriage category after a cohabitation ends; TRUE
#'   classifies cohab breakup/death as separated/widowed. Search duration still
#'   starts at the latest partnership ending. Only applies to typed unions.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model. Rate tables are tibbles.
#' @export
union_formation <- function(initiation, matching = create_age_gap_matching_model(0, 5),
                            rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE,
                            diagnostics = FALSE,
                            cohab_prob = NULL, cohab_matching = NULL, union_convert_rate = NULL,
                            cohabitation_counts_as_previous_union = FALSE,
                            weights_conversion = NULL, weights = NULL) {
  process <- structure(list(
    name = "union_formation", initiation = initiation, matching = matching,
    min_male_age = min_male_age,
    kinship_generations = kinship_generations, father_matching = father_matching,
    cancel_birth_if_no_match = cancel_birth_if_no_match, diagnostics = diagnostics,
    rates = rates, rate_unit = "person_year",
    cohab_prob = cohab_prob, cohab_matching = cohab_matching, union_convert_rate = union_convert_rate,
    cohabitation_counts_as_previous_union = cohabitation_counts_as_previous_union,
    weights_conversion = normalize_model_argument(weights_conversion),
    weights = normalize_model_argument(weights)
  ), class = c("hammel_union_formation", "hammel_process"))
  validate_union_formation(process)
  validate_process_weights(weights_conversion, "union")
  validate_process_weights(weights, "person")
  if (!is.null(weights) && !identical(initiation, "rate"))
    stop("Search `weights` require rate-based union formation.", call. = FALSE)
  if (!is.null(weights) && !is.null(rates))
    warn_cancelled_weight_terms(weights, weight_strata("search", rates), "search")
  if (!is.null(weights_conversion) && !is.null(union_convert_rate))
    warn_cancelled_weight_terms(weights_conversion, weight_strata("conversion", union_convert_rate),
                                "conversion")
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
                            diagnostics = FALSE,
                            cohab_prob = NULL, cohab_matching = NULL, union_convert_rate = NULL,
                            cohabitation_counts_as_previous_union = FALSE,
                            weights_conversion = NULL, weights = NULL) {
  set_process(model, union_formation(initiation, matching, rates,
                                     min_male_age, kinship_generations, father_matching,
                                     cancel_birth_if_no_match, diagnostics,
                                     cohab_prob, cohab_matching, union_convert_rate,
                                     cohabitation_counts_as_previous_union, weights_conversion,
                                     weights))
}

# HELPER FUNCTIONS ---------------------------------------------------------

validate_union_formation <- function(process) {
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
    if (!is.null(process$rates)) {
      stop("On-demand union formation does not take rates.", call. = FALSE)
    }
  } else {
    validate_union_search_rates(process$rates)
  }
  invisible(process)
}

expand_previous_union_rates <- function(rates) {
  if (!"previous_union" %in% names(rates)) return(rates)
  changed <- FALSE
  expand_slice <- function(x) {
    copies <- lapply(split_rate_groups(x), function(rows) {
      copy <- x[rows[x$previous_union[rows] == "dissolved"], , drop = FALSE]
      if (nrow(copy)) changed <<- TRUE
      do.call(rbind, lapply(c("separated", "widowed"), function(label) {
        out <- copy
        out$previous_union <- rep(label, nrow(copy))
        out
      }))
    })
    do.call(rbind, c(list(x[x$previous_union != "dissolved", , drop = FALSE]), copies))
  }
  if (!has_calendar_rates(rates)) return(expand_slice(rates))
  pieces <- lapply(calendar_rate_segments(rates), function(segment) {
    x <- segment$rates
    if (!nrow(x)) return(NULL)
    x$time_start <- segment$start
    x$time_end <- segment$end
    x <- x[, names(rates), drop = FALSE]
    expand_slice(x)
  })
  if (changed) do.call(rbind, pieces) else rates
}

validate_union_search_rates <- function(rates) {
    check_stratified_columns(rates, "previous_union", "Union-search")
    if (has_calendar_rates(rates)) {
      validate_calendar_rates(rates, validate_union_search_rates)
      effective <- expand_previous_union_rates(rates)
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
        stop("Duration search bands require `previous_union` ('dissolved', 'separated', or 'widowed').", call. = FALSE)
      starts[duration] <- rates$dur_start[duration]
      ends[duration] <- rates$dur_end[duration]
    }
    if (any(!is.finite(starts)) || any(starts < 0) || any(ends <= starts))
      stop("Use valid nonnegative age or duration bands with finite starts and ends greater than starts.",
           call. = FALSE)
    strata <- split_rate_groups(rates, candidate_strata(rates))
    if ("previous_union" %in% names(rates)) {
      previous <- rates$previous_union
      categories <- c("none", "separated", "widowed", "dissolved")
      if (!is.character(previous) || !is.null(dim(previous)) || anyNA(previous) ||
          any(!previous %in% categories))
        stop("`previous_union` must contain 'none', 'separated', 'widowed', or pooled 'dissolved'.", call. = FALSE)
      if (any(vapply(strata, function(rows)
          !(setequal(previous[rows], c("none", "dissolved")) ||
            setequal(previous[rows], c("none", "separated", "widowed"))), logical(1))))
        stop("Supply either 'none' and 'dissolved', or all three detailed `previous_union` categories: 'none', 'separated', 'widowed', within each group and calendar period. Do not mix pooled and detailed categories.",
             call. = FALSE)
      if (any(duration & previous == "none"))
        stop("`previous_union = 'none'` requires age-specific bands.", call. = FALSE)
      effective <- expand_previous_union_rates(rates)
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
