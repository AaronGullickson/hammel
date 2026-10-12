#' Set a demographic process specification
#'
#' Store a process specification under a name, replacing an existing process
#' with the same name. Other processes are retained. Specifications are stored
#' unchanged: rate schemas and the executable process contract are not yet
#' defined. No rates are converted and no simulation is executed.
#'
#' `set_mortality()` and `set_fertility()` delegate to `set_process()` using
#' the names `"mortality"` and `"fertility"`, respectively.
#'
#' @param model A `hammel_model` object.
#' @param name A single nonempty process name.
#' @param process A non-NULL process specification. For the demographic
#'   shortcuts, this can be a supplied rate table in person-year units.
#'
#' @return An updated `hammel_model`. The input model is not modified.
#' @export
#' @examples
#' mortality_rates <- data.frame(age = 0, rate = 0.01)
#' fertility_rates <- data.frame(age = 20, rate = 0.1)
#' model <- hammel_model() |>
#'   set_mortality(mortality_rates) |>
#'   set_fertility(fertility_rates)
#' model$processes
set_process <- function(model, name, process) {
  check_model(model)
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name))) {
    stop("`name` must be a single nonempty string.", call. = FALSE)
  }
  if (is.null(process)) {
    stop("`process` must not be NULL.", call. = FALSE)
  }
  model$processes[name] <- list(process)
  model
}

#' @rdname set_process
#' @export
set_mortality <- function(model, process) {
  set_process(model, "mortality", process)
}

#' @rdname set_process
#' @export
set_fertility <- function(model, process) {
  set_process(model, "fertility", process)
}
