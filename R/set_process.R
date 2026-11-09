# Validate and attach demographic process specifications to models.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Set a demographic process specification
#'
#' Store a process specification under its name, replacing an existing process
#' with the same name. Other processes are retained. Built-in rate schedules are
#' validated again when attached. No simulation is executed.
#'
#' Use this function to attach an already-created process. The convenience
#' setters [set_mortality()], [set_fertility()], [set_union_formation()],
#' [set_union_dissolution()], and [set_group_transition()] construct and attach
#' the corresponding process in one call.
#'
#' @param model A `hammel_model` object.
#' @param process A `hammel_process` specification.
#' @return An updated `hammel_model`. The input model is not modified.
#' @seealso [mortality()], [fertility()], [union_formation()],
#'   [union_dissolution()], [group_transition()]
#' @export
#' @examples
#' mortality_rates <- data.frame(
#'   sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008)
#' )
#' process <- mortality(mortality_rates)
#' model <- hammel_model() |> set_process(process)
#' model$processes
#' @order 0
set_process <- function(model, process) {
  model <- attach_process(model, process)
  model
}

# HELPER FUNCTIONS ---------------------------------------------------------

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
        warning(paste(
          '`spacing_method = "midpoint"` is no longer available;',
          "cell-accounted fertility is used."
        ), call. = FALSE)
      process$spacing_method <- NULL
    }
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_formation") {
    validate_union_formation(process)
    validate_process_weights(process$weights_conversion, "union")
    validate_process_weights(process$weights, "person")
    if (!is.null(process$rates)) process$rates <- tibble::as_tibble(
      process$rates
    )
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
  fields <- list(
    name = name, rates = tibble::as_tibble(rates),
    rate_unit = "person_year"
  )
  if (!is.null(attr(rates, "source"))) fields$source <- attr(rates, "source")
  structure(fields,
    class = c(paste0("hammel_", name), "hammel_process")
  )
}
