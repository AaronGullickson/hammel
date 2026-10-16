#' Specify a conditional-logit matching model
#'
#' Define relative partner utilities with an R formula and named coefficients.
#' Currently the only predictor is `age_difference`: male age minus female age,
#' in continuous years at the search time. Formula transformations such as
#' `I(age_difference^2)` and `abs(age_difference)` are supported through
#' [stats::model.matrix()]. An intercept is ignored because it cancels from
#' partner-only conditional-logit probabilities. No wait alternative is included.
#'
#' Coefficient names must match the non-intercept model-matrix columns exactly.
#' Supply all coefficients explicitly; no demographic coefficients are defaulted.
#' Formula functions must be deterministic and must not mutate simulation state.
#' Data-dependent transformations can change with the candidate pool; use fixed
#' transformations when coefficients should have a fixed interpretation.
#'
#' @param formula A one-sided R formula using only `age_difference`.
#' @param coefficients A finite named numeric vector of utility coefficients.
#' @return An inspectable `hammel_matching_model` list containing `formula`
#'   and `coefficients`. This specifies a model; it does not estimate one.
#' @export
#' @examples
#' # Illustrative coefficients, not demographic defaults.
#' matching <- create_matching_model(
#'   ~ age_difference + I(age_difference^2),
#'   c(age_difference = 0.1, "I(age_difference^2)" = -0.02)
#' )
create_matching_model <- function(formula, coefficients) {
  model <- structure(list(formula = formula, coefficients = coefficients),
                     class = "hammel_matching_model")
  validate_matching_model(model)
  model
}

validate_matching_model <- function(model) {
  if (!inherits(model, "hammel_matching_model") || !is.list(model)) {
    stop("`matching` must be created with create_matching_model().", call. = FALSE)
  }
  formula <- model$formula
  if (!inherits(formula, "formula") || length(formula) != 2L ||
      any(!all.vars(formula) %in% "age_difference")) {
    stop("Use a one-sided matching formula with only `age_difference`.", call. = FALSE)
  }
  if (length(attr(stats::terms(formula), "offset"))) {
    stop("Matching formula offsets are not supported.", call. = FALSE)
  }
  beta <- model$coefficients
  if (!is.numeric(beta) || !is.null(dim(beta)) || any(!is.finite(beta)) ||
      is.null(names(beta)) || anyNA(names(beta)) || any(!nzchar(names(beta))) ||
      anyDuplicated(names(beta))) {
    stop("`coefficients` must be a finite numeric vector with unique nonempty names.",
         call. = FALSE)
  }
  # Validate expansion and coefficient names on representative age gaps.
  data <- data.frame(age_difference = seq(-50, 50, length.out = 101))
  design <- matching_design_matrix(formula, data)
  if (!ncol(design) || !setequal(colnames(design), names(beta))) {
    stop("Coefficient names must match all non-intercept model-matrix columns.",
         call. = FALSE)
  }
  invisible(model)
}

matching_design_matrix <- function(formula, data) {
  design <- withCallingHandlers(
    stats::model.matrix(formula, data = data, na.action = stats::na.pass),
    warning = function(w) stop("Matching formula evaluation: ", conditionMessage(w),
                               call. = FALSE)
  )
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (nrow(design) != nrow(data) || any(!is.finite(design))) {
    stop("Matching formula must produce finite values for every candidate.", call. = FALSE)
  }
  design
}

# Evaluate the same formula for R diagnostics and native fallback scoring.
calculate_matching_utilities <- function(matching, age_difference) {
  check_numeric_column(age_difference, "age_difference")
  if (!length(age_difference)) return(numeric())
  design <- matching_design_matrix(matching$formula,
                                  data.frame(age_difference = age_difference))
  if (!setequal(colnames(design), names(matching$coefficients))) {
    stop("Matching formula columns changed for the supplied candidates.", call. = FALSE)
  }
  utility <- drop(design %*% matching$coefficients[colnames(design)])
  if (any(!is.finite(utility))) {
    stop("Matching utilities must be finite.", call. = FALSE)
  }
  unname(utility)
}

calculate_matching_probabilities <- function(matching, age_difference) {
  utility <- calculate_matching_utilities(matching, age_difference)
  if (!length(utility)) return(numeric())
  weights <- exp(utility - max(utility))
  weights / sum(weights)
}

# Use an exact native path for the common linear/raw quadratic specification.
# All other formulas retain model.matrix semantics through an R evaluator.
prepare_matching_model <- function(matching) {
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
  list(quadratic = quadratic, has_squared = "I(age_difference^2)" %in% columns,
       evaluate = function(age_difference) {
         calculate_matching_utilities(matching, age_difference)
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
#' `sample_size = 0` means the whole pool. Otherwise sample without replacement
#' up to that many men. Age screening occurs afterward, without refilling the
#' sample. Kinship screening follows the same ordering. Existing unions are
#' retained even if a new-search restriction would exclude the partner.
#' Candidate sets will not be retained in simulation outputs.
#'
#' @param initiation Either `"on_demand"` for birth-triggered searches or
#'   `"rate"` for scheduled searches. No default is selected.
#' @param matching A specification created by [create_matching_model()].
#' @param cancel_birth_if_no_match Logical, default `TRUE`. Cancel a birth when
#'   its on-demand partner search or rate-mode nonmarital father selection fails.
#'   Set `FALSE` to allow births with unknown fathers. Historical missing father
#'   IDs are retained. A cancellation starts no new spacing interval and draws
#'   fresh waiting times from the failed attempt. Partner scarcity can reduce
#'   realized fertility below inputs.
#' @param father_matching Optional matching specification for fathers of
#'   unpartnered births in rate mode, created by [create_matching_model()].
#'   Defaults to `NULL`, using the current general `matching` model. Union
#'   searches always use `matching`. Not used by on-demand initiation. Sampling
#'   and eligibility settings are shared by both models.
#' @param sample_size Nonnegative whole number of men to sample before screening;
#'   zero requests the whole pool. Defaults to `1000`.
#' @param min_male_age Minimum male age in continuous years at a new partner
#'   search, default `15`. Finite nonnegative scalar; zero disables this age
#'   exclusion. Existing unions are retained. No female minimum age is added.
#' @param kinship_generations Generation bound on both sides of a known common
#'   ancestor, default `2`. Excludes parents, grandparents, siblings, aunts/uncles,
#'   nieces/nephews and first cousins, including half relationships, when known.
#'   Zero imposes no cross-person kin exclusion. Missing ancestry is not proof of
#'   unrelatedness; only demonstrated relationships are screened. Existing unions
#'   are retained. See [extract_kin()] for the shared distance definitions.
#' @param rates For `"rate"` initiation, an age-specific female search-hazard
#'   table with `age_start`, `age_end`, and `rate`, in years and per person-year.
#'   These specify search attempts, not guaranteed realized union formation.
#'   Nonoverlapping bands may have gaps or an open end; gaps and ages outside
#'   support have zero hazard. Must be NULL for `"on_demand"`.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model. Rate tables are tibbles.
#' @export
union_formation <- function(initiation, matching, sample_size = 1000, rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE) {
  process <- structure(list(
    name = "union_formation", initiation = initiation, matching = matching,
    sample_size = sample_size, min_male_age = min_male_age,
    kinship_generations = kinship_generations, father_matching = father_matching,
    cancel_birth_if_no_match = cancel_birth_if_no_match,
    rates = rates, rate_unit = "person_year",
    sampling_stage = "before_pair_screening"
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
set_union_formation <- function(model, initiation, matching, sample_size = 1000, rates = NULL,
                            min_male_age = 15, kinship_generations = 2,
                            father_matching = NULL, cancel_birth_if_no_match = TRUE) {
  set_process(model, union_formation(initiation, matching, sample_size, rates,
                                     min_male_age, kinship_generations, father_matching,
                                     cancel_birth_if_no_match))
}

validate_union_formation <- function(process) {
  if (!identical(process$name, "union_formation") ||
      !identical(process$rate_unit, "person_year") ||
      !identical(process$sampling_stage, "before_pair_screening")) {
    stop("Invalid union-formation specification or units.", call. = FALSE)
  }
  if (!is.character(process$initiation) || length(process$initiation) != 1L ||
      is.na(process$initiation) || !process$initiation %in% c("rate", "on_demand")) {
    stop('`initiation` must be "rate" or "on_demand".', call. = FALSE)
  }
  n <- process$sample_size
  if (!is.numeric(n) || length(n) != 1L || !is.null(dim(n)) ||
      !is.finite(n) || n < 0 || n != floor(n) || n > .Machine$integer.max) {
    stop("`sample_size` must be a nonnegative whole number within integer range.",
         call. = FALSE)
  }
  age <- process$min_male_age
  if (!is.numeric(age) || length(age) != 1L || !is.null(dim(age)) ||
      !is.finite(age) || age < 0) {
    stop("`min_male_age` must be a finite nonnegative number of years.", call. = FALSE)
  }
  cancel <- process$cancel_birth_if_no_match
  if (!is.logical(cancel) || length(cancel) != 1L || !is.null(dim(cancel)) || is.na(cancel)) {
    stop("`cancel_birth_if_no_match` must be TRUE or FALSE.", call. = FALSE)
  }
  validate_kinship_generations(process$kinship_generations)
  validate_matching_model(process$matching)
  if (!is.null(process$father_matching)) validate_matching_model(process$father_matching)
  if (process$initiation == "on_demand") {
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
    ordered <- rates[order(rates$age_start), ]
    if (nrow(ordered) > 1L &&
        any(tail(ordered$age_start, -1) < head(cummax(ordered$age_end), -1))) {
      stop("Union-search age bands must not overlap.", call. = FALSE)
    }
  }
  invisible(process)
}
