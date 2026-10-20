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
    if (name == "fertility") validate_birth_interval(process$birth_interval)
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_formation") {
    validate_union_formation(process)
    if (!is.null(process$rates)) process$rates <- tibble::as_tibble(process$rates)
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
set_fertility <- function(model, rates, birth_interval = 0.75) {
  set_process(model, fertility(rates, birth_interval = birth_interval))
}

#' Specify age-specific demographic rates
#'
#' Construct a process specification without running a simulation. Rates are
#' piecewise constant over age intervals `[age_start, age_end)`, in years,
#' and expressed per person-year. Overlapping intervals are rejected; adjacent
#' intervals are allowed. Input order and custom columns are preserved; rate
#' schedules are returned as tibbles. A `source` attribute on the input table
#' is copied to the process's `source` field.
#'
#' Mortality requires separate schedules for both sexes. Gaps and uncovered
#' ages are not filled or interpreted as zero rates. Coverage and fertility
#' eligibility must be specified before simulation. Calendar-time-varying
#' schedules are not yet supported; time interval columns are rejected.
#'
#' Fertility optionally accepts `union_status` containing both `"partnered"`
#' and `"unpartnered"`. Rates then use woman-years in the indicated status.
#' Without that column, one overall schedule applies to all women. A partial
#' status table or any other label is rejected; use explicit zero rates to
#' suppress births in a status. Overlap is checked separately within each status.
#' On-demand union formation requires overall rates; rate-based formation accepts
#' either form and warns when an overall schedule is used.
#'
#' For fertility, one adjusted hazard is calculated at each original band's
#' midpoint: the input rate divided by one minus the input fertility integral
#' over the preceding `birth_interval`. The integral accounts for overlaps with
#' other input bands within the same status schedule. Execution rejects
#' nonpositive midpoint eligibility. Birth spacing persists across status changes;
#' compensation for the prior status history is not implemented.
#' This approximates the age-varying adjustment; initialization is not calibrated.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`. Starts must be finite and nonnegative; ends must exceed starts
#'   and must be finite for fertility; mortality ends may be `Inf`. Rates must
#'   be finite and nonnegative. Mortality additionally
#'   requires character `sex` values `"m"` and `"f"`, with both represented.
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative. Applies to fertility
#'   only; eligible hazards use midpoint birth-spacing compensation.
#' @return A named list containing `name`, `rates`, and `rate_unit`, plus
#'   `birth_interval` for fertility, inheriting
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
fertility <- function(rates, birth_interval = 0.75) {
  validate_rates(rates, mortality = FALSE)
  validate_birth_interval(birth_interval)
  process <- new_rate_process("fertility", rates)
  process$birth_interval <- birth_interval
  process
}

new_rate_process <- function(name, rates) {
  fields <- list(name = name, rates = tibble::as_tibble(rates), rate_unit = "person_year")
  if (!is.null(attr(rates, "source"))) fields$source <- attr(rates, "source")
  structure(fields,
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
  if (!mortality && any(!is.finite(ends))) {
    stop("Fertility age bands must have finite ends.", call. = FALSE)
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
    if ("union_status" %in% names(rates)) {
      status <- rates$union_status
      if (!is.character(status) || !is.null(dim(status)) || anyNA(status) ||
          !setequal(status, c("partnered", "unpartnered"))) {
        stop('Fertility `union_status` must contain both "partnered" and "unpartnered", with no other or missing values.',
             call. = FALSE)
      }
      schedules <- split(seq_len(nrow(rates)), status)
    } else {
      schedules <- list(seq_len(nrow(rates)))
    }
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
