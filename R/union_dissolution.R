# Construct and validate union-dissolution process specifications.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify duration-dependent union dissolution
#'
#' Construct or attach an optional dissolution process. One continuous hazard
#' applies per union, using years since the union began, independently of either
#' partner's age. Rates are per partnered woman-year (equivalently union-year),
#' not annual probabilities. Omitting this process retains death-only union
#' endings.
#'
#' Births do not affect a union's dissolution risk. Exact ties select death,
#' then birth, then dissolution; a tied birth occurs while partnered.
#' Dissolution ends the episode, restores both partners' availability, and moves
#' the woman to her unpartnered fertility and search risks, retaining birth
#' spacing. Optional age-specific search schedules by previous union are
#' specified through [union_formation()]. Restarts resume at attained union
#' duration rather than reset duration.
#'
#' @section Calendar rates:
#' Optional `time_start`/`time_end` select half-open absolute calendar-year
#' bands, including fractional years. Active slices retain the existing strata
#' rules and must cover the requested run. Changes precede same-time events
#' without resetting the union-duration origin. Other arguments remain static.
#'
#' `union_dissolution()` creates the process specification;
#' `set_union_dissolution()` creates it and attaches it to a model through
#' [set_process()]. Neither function runs a simulation.
#'
#' @param weights NULL (default) for unweighted selection, or a [weight_model()]
#'   using union variables, which chooses which union in a cell dissolves
#'   without changing the cell's input rate.
#' @param rates A nonempty data frame with numeric `dur_start`, `dur_end`, and
#'   `rate`. Bands are half-open `[dur_start, dur_end)` in years. Starts and
#'   rates must be finite and nonnegative; ends must exceed starts and may be
#'   `Inf`. Bands cannot overlap within a stratum. Gaps and durations outside
#'   support have zero hazard. Sex, age and union-status columns are not
#'   supported. Input order, custom columns and provenance are retained.
#'   Optional `<group>_w` and `<group>_h` columns (for a declared group
#'   dimension such as `race`: `race_w`, `race_h`) select the wife's and
#'   husband's categories; strata are directional, so A-wife/B-husband and
#'   B-wife/A-husband are separate. Stratified tables must cover every
#'   combination of the dimension's categories. Optional `union_type` supplies
#'   both `cohab` and `marriage` schedules within each stratum; without it the
#'   same duration risk applies to either type. Use [stratify_rates()] to add
#'   proportional differences, for example by listing both directions of a pair.
#' @param model A `hammel_model` object.
#' @return `union_dissolution()` returns an inspectable
#'   `hammel_union_dissolution` process with a tibble of rates;
#'   `set_union_dissolution()` returns an updated model.
#' @export
#' @examples
#' # Illustrative rates, not demographic defaults.
#' rates <- data.frame(
#'   dur_start = c(0, 5), dur_end = c(5, Inf),
#'   rate = c(0.04, 0.02)
#' )
#' model <- hammel_model() |> set_union_dissolution(rates)
union_dissolution <- function(rates, weights = NULL) {
  process <- new_rate_process("union_dissolution", rates)
  validate_union_dissolution(process)
  validate_process_weights(weights, "union")
  warn_cancelled_weight_terms(
    weights, weight_strata("dissolution", rates),
    "dissolution"
  )
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname union_dissolution
#' @export
set_union_dissolution <- function(model, rates, weights = NULL) {
  set_process(model, union_dissolution(rates, weights))
}

# HELPER FUNCTIONS ---------------------------------------------------------

validate_union_dissolution <- function(process, calendar_interval = NULL) {
  if (!identical(process$name, "union_dissolution") ||
    !identical(process$rate_unit, "person_year"))
    stop("Invalid union-dissolution specification or units.", call. = FALSE)
  rates <- process$rates
  check_stratified_columns(rates, "union_type", "Union-dissolution")
  if (validate_calendar_rates(rates, function(slice, interval) {
    copy <- process
    copy$rates <- slice
    validate_union_dissolution(copy, interval)
  }, with_interval = TRUE)) return(invisible(process))
  check_columns(rates, c("dur_start", "dur_end", "rate"), "rates")
  check_numeric_column(rates$dur_start, "rates$dur_start")
  check_numeric_column(rates$rate, "rates$rate")
  ends <- rates$dur_end
  if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) ||
    anyNA(ends) || any(ends <= rates$dur_start) ||
    any(rates$dur_start < 0) || any(rates$rate < 0) ||
    any(c(
      "sex", "age_start", "age_end", "union_status", "time_start",
      "time_end"
    ) %in% names(rates)))
    stop(
      paste0(
        "Use nonnegative dissolution rates in valid ",
        "duration bands; use `<group>_w`/`<group>_h` for ",
        "group strata, without sex, age, status or ",
        "calendar columns."
      ),
      call. = FALSE
    )
  # Possible group columns (`<group>_w`, `<group>_h`) define strata.
  strata <- candidate_strata(rates)
  pair_labels <- if (length(strata)) do.call(paste, c(lapply(
    rates[strata],
    as.character
  ), sep = ":")) else
    rep("shared", nrow(rates))
  types <- if ("union_type" %in% names(rates)) rates$union_type else rep(
    "shared", nrow(rates)
  )
  if ("union_type" %in% names(rates)) {
    if (!is.character(types) || anyNA(types) || any(!types %in% c(
      "cohab",
      "marriage"
    )))
      stop(paste0(
        "Dissolution `union_type` must use cohab and ",
        "marriage labels without missing values."
      ), call. = FALSE)
    for (pair in unique(pair_labels)) {
      missing <- setdiff(c("cohab", "marriage"), types[pair_labels == pair])
      if (!length(missing)) next
      context <- if (is.null(calendar_interval)) "" else sprintf(
        " during calendar interval [%s, %s)",
        format(calendar_interval[1], trim = TRUE, scientific = FALSE),
        format(calendar_interval[2], trim = TRUE, scientific = FALSE)
      )
      if (length(strata)) context <- paste0(context, sprintf(
        " for group stratum '%s'", pair
      ))
      stop(sprintf(
        paste0(
          "Dissolution rates are missing a '%s' schedule%s. ",
          paste0(
            "Supply both cohab and marriage schedules wherever ",
            "typed dissolution rates are supplied."
          )
        ),
        missing, context
      ), call. = FALSE)
    }
  }
  for (rows in split(seq_len(nrow(rates)), paste(pair_labels, types,
    sep =
      "\r"
  ))) {
    ordered <- rates[rows, , drop = FALSE]
    ordered <- ordered[order(ordered$dur_start), , drop = FALSE]
    if (nrow(ordered) > 1L &&
      any(tail(ordered$dur_start, -1) < head(cummax(ordered$dur_end), -1)))
      stop("Union-dissolution duration bands must not overlap.", call. = FALSE)
  }
  invisible(process)
}
