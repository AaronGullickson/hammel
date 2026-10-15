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
#' [simulate()] supports a temporary on-demand implementation: at a birth use
#' the woman's current partner, or sample living men and choose a partner using
#' conditional-logit probabilities from the matching formula to form a union. If no men exist, preserve the birth with unknown father.
#' Pair eligibility is not yet applied, so even
#' partnered men, relatives, and children can be selected. Death ends affected
#' unions. Rate-based initiation and conditional fertility are not implemented.
#'
#' Sampling is from the available male pool BEFORE pair-specific screening.
#' `sample_size = 0` means the whole pool. Otherwise sample without replacement
#' up to that many men. Pair screening is deferred in the temporary engine path.
#' Candidate sets will not be retained in simulation outputs.
#'
#' @param initiation Either `"on_demand"` for birth-triggered searches or
#'   `"rate"` for scheduled searches. No default is selected.
#' @param matching A specification created by [create_matching_model()].
#' @param sample_size Nonnegative whole number of men to sample before screening;
#'   zero requests the whole pool. Defaults to `1000`.
#' @param rates For `"rate"` initiation, an age-specific female search-hazard
#'   table with `age_start`, `age_end`, and `rate`, in years and per person-year.
#'   These specify search attempts, not guaranteed realized union formation.
#'   Nonoverlapping bands may have gaps or an open end; execution semantics for
#'   eligibility and gaps are not yet implemented. Must be NULL for `"on_demand"`.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model. Rate tables are tibbles.
#' @export
union_formation <- function(initiation, matching, sample_size = 1000, rates = NULL) {
  process <- structure(list(
    name = "union_formation", initiation = initiation, matching = matching,
    sample_size = sample_size, rates = rates, rate_unit = "person_year",
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
set_union_formation <- function(model, initiation, matching, sample_size = 1000, rates = NULL) {
  set_process(model, union_formation(initiation, matching, sample_size, rates))
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
  validate_matching_model(process$matching)
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
