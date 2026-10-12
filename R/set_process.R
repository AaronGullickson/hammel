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
set_process <- function(model, process) {
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
    process$rates <- tibble::as_tibble(process$rates)
  }
  model$processes[name] <- list(process)
  model
}

#' @rdname set_process
#' @export
set_mortality <- function(model, rates) {
  set_process(model, mortality(rates))
}

#' @rdname set_process
#' @export
set_fertility <- function(model, rates) {
  set_process(model, fertility(rates))
}

#' Specify age-specific demographic rates
#'
#' Construct a process specification without running a simulation. Rates are
#' piecewise constant over age intervals `[age_start, age_end)`, in years,
#' and expressed per person-year. Overlapping intervals are rejected; adjacent
#' intervals are allowed. Input order and custom columns are preserved; rate
#' schedules are returned as tibbles.
#'
#' Mortality requires separate schedules for both sexes. Gaps and uncovered
#' ages are not filled or interpreted as zero rates. Coverage and fertility
#' eligibility must be specified before simulation. Calendar-time-varying
#' schedules are not yet supported; time interval columns are rejected.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`. Starts must be finite and nonnegative; ends must exceed starts
#'   and may be `Inf`. Rates must be finite and nonnegative. Mortality additionally
#'   requires character `sex` values `"m"` and `"f"`, with both represented.
#' @return A named list containing `name`, `rates`, and `rate_unit`, inheriting
#'   from `hammel_mortality` or `hammel_fertility` and `hammel_process`.
#' @export
#' @examples
#' deaths <- mortality(data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008)))
#' births <- fertility(data.frame(age_start = 15, age_end = 50, rate = 0.1))
#' model <- hammel_model() |>
#'   set_process(deaths) |>
#'   set_process(births)
mortality <- function(rates) {
  validate_rates(rates, mortality = TRUE)
  new_rate_process("mortality", rates)
}

#' @rdname mortality
#' @export
fertility <- function(rates) {
  validate_rates(rates, mortality = FALSE)
  new_rate_process("fertility", rates)
}

new_rate_process <- function(name, rates) {
  structure(list(name = name, rates = tibble::as_tibble(rates), rate_unit = "person_year"),
            class = c(paste0("hammel_", name), "hammel_process"))
}

validate_rates <- function(rates, mortality) {
  columns <- c("age_start", "age_end", "rate")
  if (mortality) columns <- c("sex", columns)
  check_columns(rates, columns, "rates")
  if (nrow(rates) == 0L) {
    stop("`rates` must contain at least one age interval.", call. = FALSE)
  }
  if (any(c("time_start", "time_end") %in% names(rates))) {
    stop("Calendar-time-varying rate schedules are not yet supported.", call. = FALSE)
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  check_numeric_column(rates$rate, "rates$rate")
  ends <- rates$age_end
  if (!is.numeric(ends) || !is.null(dim(ends)) || anyNA(ends) ||
      any(ends <= rates$age_start) || any(rates$age_start < 0)) {
    stop("Age intervals must have nonnegative starts and greater ends (possibly Inf).",
         call. = FALSE)
  }
  if (any(rates$rate < 0)) {
    stop("`rates$rate` must be nonnegative.", call. = FALSE)
  }
  if (mortality) {
    if (!is.character(rates$sex) || !is.null(dim(rates$sex)) ||
        anyNA(rates$sex) || any(!rates$sex %in% c("m", "f")) ||
        !all(c("m", "f") %in% rates$sex)) {
      stop('Mortality rates must include both "m" and "f" sex schedules.',
           call. = FALSE)
    }
    schedules <- split(seq_len(nrow(rates)), rates$sex)
  } else {
    if ("sex" %in% names(rates)) {
      stop("Fertility schedules do not take a sex column in this initial interface.",
           call. = FALSE)
    }
    schedules <- list(seq_len(nrow(rates)))
  }
  for (rows in schedules) {
    rows <- rows[order(rates$age_start[rows])]
    if (length(rows) > 1L &&
        any(rates$age_start[rows[-1L]] <
            cummax(rates$age_end[rows])[-length(rows)])) {
      stop("Age intervals must not overlap within a rate schedule.", call. = FALSE)
    }
  }
  invisible(rates)
}
