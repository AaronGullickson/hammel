#' Set a demographic process specification
#'
#' Store a process specification under its name, replacing an existing process
#' with the same name. Other processes are retained. Built-in rate schedules
#' are validated again when attached. No simulation is executed.
#'
#' `set_mortality()` and `set_fertility()` construct the corresponding process
#' and delegate to `set_process()`.
#'
#' @param model A `hammel_model` object.
#' @param process A `hammel_process` specification.
#' @param rates An age-specific rate table passed to [mortality()] or [fertility()].
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative.
#' @param hetfert NULL (default) for uniform mother selection, or a
#'   [het_fert()] specification for persistent lognormal selection weights.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [create_weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate. In
#'   fertility they combine with `hetfert` on the log scale.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
#'
#' @return An updated `hammel_model`. The input model is not modified.
#' @export
#' @examples
#' mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008))
#' fertility_rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
#' model <- hammel_model() |>
#'   set_mortality(mortality_rates) |>
#'   set_fertility(fertility_rates)
#' model$processes
#' @order 0
set_process <- function(model, process) {
  model <- attach_process(model, process)
  model
}

# Validate and store a process without user-facing advisory warnings; used by
# set_process() and by simulate()'s re-validation of edited specifications.
attach_process <- function(model, process) {
  check_model(model)
  if (!inherits(process, "hammel_process") || !is.list(process)) {
    stop("`process` must be a hammel_process object.", call. = FALSE)
  }
  name <- process$name
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name))) {
    stop("`name` must be a single nonempty string.", call. = FALSE)
  }
  if (name %in% c("mortality", "fertility")) {
    if (!identical(process$rate_unit, "person_year")) {
      stop("Built-in process rates must use person-year units.", call. = FALSE)
    }
    validate_rates(process$rates, mortality = name == "mortality")
    validate_process_weights(process$weights, "person")
    if (name == "fertility") {
      validate_birth_interval(process$birth_interval)
      if (!"ffab" %in% names(process)) process$ffab <- 0.4886
      validate_ffab(process$ffab)
      validate_het_fert(process$hetfert)
      # Specifications saved while midpoint compensation existed (before D101
      # became the only method) may carry a spacing_method field.
      if (identical(process$spacing_method, "midpoint"))
        warning(paste('`spacing_method = "midpoint"` is no longer available;',
                      "cell-accounted fertility is used."), call. = FALSE)
      process$spacing_method <- NULL
    }
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_formation") {
    validate_union_formation(process)
    validate_process_weights(process$weights_conversion, "union")
    validate_process_weights(process$weights, "person")
    if (!is.null(process$rates)) process$rates <- tibble::as_tibble(process$rates)
    for (field in c("matching", "father_matching", "cohab_matching"))
      process[field] <- list(normalize_model_argument(process[[field]]))
  }
  if (name == "group_transition") {
    validate_group_transition(process)
    validate_process_weights(process$weights, "transition")
    process$rates <- tibble::as_tibble(process$rates)
    return(store_transition_process(model, process))
  }
  if (name == "union_dissolution") {
    validate_union_dissolution(process)
    validate_process_weights(process$weights, "union")
    process$rates <- tibble::as_tibble(process$rates)
  }
  model$processes[name] <- list(process)
  model
}

new_rate_process <- function(name, rates) {
  fields <- list(name = name, rates = tibble::as_tibble(rates), rate_unit = "person_year")
  if (!is.null(attr(rates, "source"))) fields$source <- attr(rates, "source")
  structure(fields,
            class = c(paste0("hammel_", name), "hammel_process"))
}
