#' Specify duration-dependent union dissolution
#'
#' Construct or attach an optional dissolution process. One continuous hazard
#' applies per union, using years since the union
#' began, independently of either partner's age. Rates are per partnered
#' woman-year (equivalently union-year), not annual probabilities. Omitting
#' this process retains death-only union endings.
#'
#' Births do not affect a union's dissolution risk. Exact ties select death,
#' then birth, then dissolution; a tied birth occurs while partnered.
#' Dissolution ends the episode, restores both partners' availability, and moves
#' the woman to her unpartnered fertility and search risks, retaining birth
#' spacing. Optional age-specific search schedules by
#' previous union are specified through [union_formation()].
#' Restarts resume at attained union duration rather than reset duration.
#'
#' @section Calendar rates:
#' Optional `time_start`/`time_end` select half-open absolute calendar-year bands,
#' including fractional years. Active slices retain the existing strata rules
#' and must cover the requested run. Changes precede same-time events without
#' resetting the union-duration origin. Other arguments remain static.
#'
#' @param weights NULL (default) for unweighted selection, or a
#'   [create_weight_model()] using union variables, which chooses which union
#'   in a cell dissolves without changing the cell's input rate.
#' @param rates A nonempty data frame with numeric `dur_start`, `dur_end`, and
#'   `rate`. Bands are half-open `[dur_start, dur_end)` in years. Starts and
#'   rates must be finite and nonnegative; ends must exceed starts and may be
#'   `Inf`. Bands cannot overlap within an unordered group pair. Gaps and
#'   durations outside support have zero hazard. Sex, age and union-status columns
#'   are not supported. Input order, custom columns and provenance are retained.
#'   An optional `group_pair` column selects an unordered pair such as `"A:A"`,
#'   `"A:B"`, or `"B:B"`; `"A:B"` and `"B:A"` select the same schedule.
#'   Optional `union_type` supplies both `cohab` and `marriage` schedules within
#'   each pair stratum; without it the same duration risk applies to either type.
#'   Pair names require nonempty colon-free group labels. Pair-stratified tables
#'   must cover every unordered pair of possible identities, including future
#'   ones. Without population groups, warn and use the first canonical pair
#'   alphabetically. Omit `group_pair` for shared rates. `group`, `group_w` and
#'   `group_h` schedule columns are not supported.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   rate multipliers by unordered group pair, e.g. `c("A:B" = 1.2)` applies
#'   in both directions. Omitted pairs use 1; zero suppresses risk. Reversed
#'   duplicate keys are errors. With no population groups, simulation warns
#'   and ignores multipliers.
#' @param model A `hammel_model` object.
#' @return `union_dissolution()` returns an inspectable
#'   `hammel_union_dissolution` process with a tibble of rates;
#'   `set_union_dissolution()` returns an updated model.
#' @export
#' @examples
#' # Illustrative rates, not demographic defaults.
#' rates <- data.frame(dur_start = c(0, 5), dur_end = c(5, Inf),
#'                     rate = c(0.04, 0.02))
#' model <- hammel_model() |> set_union_dissolution(rates)
union_dissolution <- function(rates, group_prop_hazard = NULL, weights = NULL) {
  process <- new_rate_process("union_dissolution", rates)
  process$group_prop_hazard <- group_prop_hazard
  validate_union_dissolution(process)
  validate_process_weights(weights, "union")
  warn_cancelled_weight_terms(weights, weight_strata("dissolution", rates), "dissolution")
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname union_dissolution
#' @export
set_union_dissolution <- function(model, rates, group_prop_hazard = NULL, weights = NULL) {
  set_process(model, union_dissolution(rates, group_prop_hazard, weights))
}

validate_union_dissolution <- function(process, calendar_interval = NULL) {
  if (!identical(process$name, "union_dissolution") ||
      !identical(process$rate_unit, "person_year"))
    stop("Invalid union-dissolution specification or units.", call. = FALSE)
  validate_dissolution_multipliers(process$group_prop_hazard)
  rates <- process$rates
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
      any(c("sex", "age_start", "age_end", "union_status", "time_start", "time_end",
            "group", "group_w", "group_h") %in% names(rates)))
    stop("Use nonnegative dissolution rates in valid duration bands; use group_pair for symmetric strata, without sex, age, status or calendar columns.",
         call. = FALSE)
  pair_labels <- if ("group_pair" %in% names(rates))
    canonical_union_pairs(rates$group_pair) else rep("shared", nrow(rates))
  types <- if ("union_type" %in% names(rates)) rates$union_type else rep("shared", nrow(rates))
  if ("union_type" %in% names(rates)) {
    if (!is.character(types) || anyNA(types) || any(!types %in% c("cohab", "marriage")))
      stop("Dissolution `union_type` must use cohab and marriage labels without missing values.", call. = FALSE)
    for (pair in unique(pair_labels)) {
      missing <- setdiff(c("cohab", "marriage"), types[pair_labels == pair])
      if (!length(missing)) next
      context <- if (is.null(calendar_interval)) "" else sprintf(
        " during calendar interval [%s, %s)",
        format(calendar_interval[1], trim = TRUE, scientific = FALSE),
        format(calendar_interval[2], trim = TRUE, scientific = FALSE))
      if ("group_pair" %in% names(rates))
        context <- paste0(context, sprintf(" for group pair '%s'", pair))
      stop(sprintf(paste0("Dissolution rates are missing a '%s' schedule%s. ",
        "Supply both cohab and marriage schedules wherever typed dissolution rates are supplied."),
        missing, context), call. = FALSE)
    }
  }
  for (rows in split(seq_len(nrow(rates)), paste(pair_labels, types, sep = "\r"))) {
    ordered <- rates[rows, , drop = FALSE]
    ordered <- ordered[order(ordered$dur_start), , drop = FALSE]
    if (nrow(ordered) > 1L &&
        any(tail(ordered$dur_start, -1) < head(cummax(ordered$dur_end), -1)))
      stop("Union-dissolution duration bands must not overlap.", call. = FALSE)
  }
  invisible(process)
}

# Use the shared Band reader with private ordered-pair indices. The public
# strata remain unordered labels, and stored rates/provenance are never changed.
prepare_dissolution_rates <- function(process, persons, group_values, start = -Inf, end = Inf) {
  if (is.null(process))
    return(data.frame(age_start = numeric(), age_end = numeric(), rate = numeric()))
  rates <- process$rates
  rates$.hammel_group <- NULL
  grouped <- "group_pair" %in% names(rates)
  labels <- if (grouped) canonical_union_pairs(rates$group_pair) else NULL
  factors <- process$group_prop_hazard
  if (!is.null(factors)) names(factors) <- canonical_union_pairs(names(factors), name = "group_prop_hazard")
  if (!"group" %in% names(persons)) {
    if (grouped) {
      first <- sort(unique(labels), method = "radix")[1L]
    if (any(vapply(calendar_rate_slices(rates, start, end), function(x)
          nrow(x) > 0L && !first %in% canonical_union_pairs(x$group_pair), logical(1))))
        stop("Dissolution rates are missing schedules for the fallback pair in a calendar interval.", call. = FALSE)
      warning(sprintf("Pair rates were specified for union dissolution but no groups found in the population; using pair '%s' (first alphabetically).",
                      first), call. = FALSE)
      rates <- rates[labels == first, , drop = FALSE]
    }
    if (!is.null(factors))
      warning("`group_prop_hazard` was specified for union dissolution but no groups found in the population; ignoring multipliers.", call. = FALSE)
  } else if (grouped || !is.null(factors)) {
    pairs <- possible_union_pairs(group_values)
    if (grouped) {
      missing <- unique(unlist(lapply(calendar_rate_slices(rates, start, end), function(x)
        if (nrow(x)) setdiff(pairs, canonical_union_pairs(x$group_pair)) else character())))
      if (length(missing))
        stop(sprintf("Union dissolution rates are missing schedules for group pair(s): %s.",
                     paste(missing, collapse = ", ")), call. = FALSE)
    }
    groups <- as.character(group_values)
    copies <- list()
    for (male in seq_along(groups)) for (female in seq_along(groups)) {
      key <- canonical_union_pairs(groups[female], groups[male])
      selected <- if (grouped) rates[labels == key, , drop = FALSE] else rates
      factor <- if (key %in% names(factors)) unname(factors[[key]]) else 1
      selected$rate <- selected$rate * factor
      selected$.hammel_group <- (female - 1L) + length(groups) * (male - 1L)
      copies[[length(copies) + 1L]] <- selected
    }
    rates <- if (length(copies)) tibble::as_tibble(do.call(rbind, copies)) else {
      rates$.hammel_group <- rep(-1L, nrow(rates))
      rates
    }
  }
  if (any(!is.finite(rates$rate)))
    stop("Union dissolution rates multiplied by `group_prop_hazard` must remain finite.", call. = FALSE)
  rates$age_start <- rates$dur_start
  rates$age_end <- rates$dur_end
  rates
}
