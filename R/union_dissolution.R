# Construct and validate union-dissolution process specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify duration-dependent union dissolution
#'
#' @description
#'
#' Constructs a union dissolution process specification for a simulation model.
#' The convenience function `set_union_dissolution` is equivalent to
#' `set_process(model, union_dissolution(...))` and allows for direct piping
#' from a model.
#'
#' @details
#'
#' The primary and only required input to the function is a union dissolution
#' rate schedule for the `rates` argument. This data frame must contain
#' `dur_start` and `dur_end` columns to define the half-open age intervals
#' (inclusive at the start and exclusive at the end) and a `rate` column for the
#' duration-specific dissolution rates.
#'
#' Rates are always specified as union dissolutions over union-years lived. Gaps
#' in the duration band are treated as cases with zero rates.
#'
#' Columns named by declared group dimensions (see [set_groups()]), such as
#' `race` or `region`, supply separate schedules for each combination of their
#' categories, but they must be suffixed by `_h` or `_w` to indicate the group
#' of the husband or wife, respectively (e.g. `race_h` and `race_w`).
#'
#' Optional `time_start`/`time_end` columns select half-open calendar intervals
#' in absolute years, including fractional years. Both columns are required
#' together; starts are finite and ends exceed starts (possibly Inf). Active
#' slices must satisfy the strata and coverage rules and cover the requested
#' simulation run.
#'
#' If cohabitations are in use in the simulation, the rates can be stratified by
#' `union_type` to supply both `cohab` and `marriage` schedules; without it the
#' same duration risk applies to either type. The [stratify_rates()] function
#' provides a quick means to add stratification to an existing rate schedule.
#'
#' @param rates A nonempty data frame with at least numeric `dur_start`,
#'   `dur_end`, and `rate`. See details for format and description of
#'   stratification variables that can be applied.
#' @param weights NULL (default) for unweighted selection, or a [weight_model()]
#'   using union variables, which chooses which union in a cell dissolves
#'   without changing the cell's input rate.
#' @param model A `hammel_model` object.
#' @return `union_dissolution()` returns an inspectable
#'   `hammel_union_dissolution` process with a tibble of rates;
#'   `set_union_dissolution()` returns an updated model.
#' @export
#' @examples
#' # Illustrative rates
#' diss_rates <- data.frame(
#'   dur_start = c(0, 5), dur_end = c(5, Inf),
#'   rate = c(0.04, 0.02)
#' )
#' model <- hammel_model() |> set_union_dissolution(diss_rates)
union_dissolution <- function(rates, weights = NULL) {
  new_union_dissolution(rates, weights)
}

#' @rdname union_dissolution
#' @export
set_union_dissolution <- function(model, rates, weights = NULL) {
  attach_process(model, new_union_dissolution(rates, weights))
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a union dissolution process
#'
#' @description
#' 
#' Summary method for class `"hammel_union_dissolution"`
#'
#' @param object A `hammel_union_dissolution` object returned by [union_dissolution()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_union_dissolution` containing:
#' \describe{
#'   \item{rate_strata}{A list describing the input rate schedule, containing
#'     `stratified_by` (stratum column names), `strata` (a tibble of observed
#'     stratum combinations), `bands` (a tibble of exact band combinations),
#'     and `bands_by_stratum` (whether band combinations differ by stratum).
#'     Shared bands are stored once; otherwise `bands` includes the stratum
#'     columns. Band endpoints define half-open intervals `[start, end)`.}
#'   \item{weights}{A summary of the weight model used for selection within
#'     rate cells: a `summary.hammel_weight_model`, a list of `periods` and
#'     `models` for a calendar schedule, or `NULL` when selection is
#'     unweighted.}
#' }
#'
#' @seealso [union_dissolution()]
#' @examples
#' diss_rates <- tibble::tribble(
#'  ~union_type,  ~dur_start, ~dur_end, ~rate,
#'   "cohab",      0,          5,        0.08,
#'   "cohab",      5,          10,       0.04,
#'   "cohab",      10,         Inf,      0.02,
#'   "marriage",   0,          5,        0.06,
#'   "marriage",   5,          10,       0.02,
#'   "marriage",   10,         15,       0.01,
#'   "marriage",   15,         Inf,      0.015
#' )
#' process <- union_dissolution(diss_rates)
#' summary(process)
#' @export
summary.hammel_union_dissolution <- function(object, ...) {

  structure(
    list(
      rate_strata = summarize_rate_strata(object$rates),
      weights = summarize_model_argument(object$weights)
    ),
    class = "summary.hammel_union_dissolution"
  )
}

#' @rdname summary.hammel_union_dissolution
#' @param x A `summary.hammel_union_dissolution` object returned by
#'   `[summary.hammel_union_dissolution()]`.
#' @export
print.summary.hammel_union_dissolution <- function(x, ...) {

  cli::cli_h2("Union dissolution process")

  list <- open_summary_list()
  print_rate_strata(x$rate_strata)
  print_model_argument("Weight model", x$weights)
  cli::cli_end(list)
  invisible(x)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Build a validated dissolution process for union_dissolution() and
# set_union_dissolution(); errors name `call`.
new_union_dissolution <- function(rates, weights, call = rlang::caller_env()) {
  process <- new_rate_process("union_dissolution", rates)
  validate_union_dissolution(process, call = call)
  validate_process_weights(weights, "union", call = call)
  warn_cancelled_weight_terms(
    weights, weight_strata("dissolution", rates), "dissolution"
  )
  process$weights <- normalize_model_argument(weights)
  process
}

# Validate in two passes, as validate_rates(): row checks on the whole table,
# then schedule checks within each calendar period and group stratum.
validate_union_dissolution <- function(process, call = rlang::caller_env()) {
  if (!identical(process$name, "union_dissolution")) {
    abort_hammel(
      "{.arg process} must be a union dissolution process.",
      "specification",
      call = call
    )
  }
  rates <- process$rates
  check_stratified_columns(
    rates, "union_type", "Union dissolution",
    call = call
  )
  check_dissolution_rows(rates, call = call)
  validate_dissolution_schedules(rates, call = call)
  invisible(process)
}

check_dissolution_rows <- function(rates, call = rlang::caller_env()) {
  check_columns(rates, c("dur_start", "dur_end", "rate"), "rates", call = call)
  if (!nrow(rates)) {
    abort_hammel(
      c(
        "{.arg rates} must contain at least one duration band.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  calendar_rate_segments(rates, call = call)
  bad <- intersect(
    c("sex", "age_start", "age_end", "union_status"),
    names(rates)
  )
  if (length(bad)) {
    abort_hammel(
      c(
        "Union dissolution rates can't have {.field {bad}} columns.",
        "i" = paste(
          "Dissolution rates use union duration bands. Stratify by group",
          "with {.field <group>_w} and {.field <group>_h} columns."
        )
      ),
      "rates",
      argument = "rates",
      columns = bad,
      call = call
    )
  }
  for (column in candidate_strata(rates)) {
    check_group_labels(rates[[column]], paste0("rates$", column), call = call)
  }
  check_bands(rates, "dur", call = call)
  check_numeric_column(rates$rate, "rates$rate", call = call)
  check_nonnegative(rates$rate, "rates$rate", call = call)
  if ("union_type" %in% names(rates)) {
    check_labels(
      rates$union_type, "rates$union_type", c("cohab", "marriage"),
      call = call
    )
  }
  invisible(rates)
}

# Within each calendar period and group stratum: both union types when rates
# are typed, and duration bands that do not overlap.
validate_dissolution_schedules <- function(rates, call = rlang::caller_env()) {
  if (validate_calendar_rates(
        rates,
        function(x) {
          validate_dissolution_schedules(x, call = call)
        },
        call = call
      )) {
    return(invisible(rates))
  }
  # Possible group columns (`<group>_w`, `<group>_h`) define strata.
  strata <- candidate_strata(rates)
  pair_labels <- if (length(strata)) {
    do.call(paste, c(lapply(rates[strata], as.character), sep = ":"))
  } else {
    rep("", nrow(rates))
  }
  types <- if ("union_type" %in% names(rates)) {
    rates$union_type
  } else {
    rep("", nrow(rates))
  }
  stratum_text <- function(pair) {
    if (length(strata)) {
      paste0(" for group stratum ", paste(strata, collapse = ":"), " = ", pair)
    } else {
      ""
    }
  }
  if ("union_type" %in% names(rates)) {
    for (pair in unique(pair_labels)) {
      absent <- setdiff(c("cohab", "marriage"), types[pair_labels == pair])
      if (!length(absent)) {
        next
      }
      stratum <- stratum_text(pair)
      abort_hammel(
        c(
          paste(
            "Union dissolution rates with a {.field union_type} column need",
            "both {.val cohab} and {.val marriage} schedules."
          ),
          "x" = "No {.val {absent}} rates{stratum}."
        ),
        "rates",
        argument = "rates",
        call = call
      )
    }
  }
  keys <- paste(pair_labels, types, sep = "\r")
  for (rows in split(seq_len(nrow(rates)), keys)) {
    schedule <- stratum_text(pair_labels[rows[1]])
    if (nzchar(types[rows[1]])) {
      schedule <- paste0(
        if (nzchar(schedule)) paste0(schedule, ",") else " for",
        " union type \"", types[rows[1]], "\""
      )
    }
    check_band_overlap(
      rates, rows, "dur",
      schedule = schedule,
      call = call
    )
  }
  invisible(rates)
}
