# Construct fertility processes and validate fertility settings and feasibility.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify an age-specific fertility process
#'
#' `fertility()` constructs a process specification without running a simulation.
#' `set_fertility()` constructs the same specification and attaches it to a
#' model, equivalent to `set_process(model, fertility(...))` with the same
#' fertility arguments. The input model is not modified and neither function
#' runs a simulation.
#'
#' Rates are piecewise constant over half-open age intervals, in years, and
#' expressed per woman-year. Overlapping intervals are rejected; adjacent
#' intervals are allowed. Input order and custom columns are preserved; rate
#' schedules become tibbles. A `source` attribute on the input table is copied
#' into the process.
#'
#' Columns named by declared group dimensions (see [set_groups()]), such as
#' `race` or `region`, supply separate schedules for each combination of their
#' categories; dimensions without a column share the schedule. Each conditional
#' fertility stratum needs both statuses (or all three typed statuses).
#' Overlap is checked within strata. Simulation requires schedules for every
#' combination of categories of the dimensions a table uses; a table using
#' an undeclared persons column is an error.
#'
#' Optional `time_start`/`time_end` columns select half-open calendar intervals
#' in absolute years, including fractional years. Both columns are required
#' together; starts are finite and ends exceed starts (possibly Inf). Active
#' slices must satisfy the strata and coverage rules and cover the requested
#' run. Rates are constant within intervals, with no interpolation. Calendar
#' changes precede same-time events, preserving histories, birth spacing,
#' parity and duration origins. Inheritance also accepts whole-model calendar
#' schedules; other fertility arguments remain static.
#'
#' Fertility optionally accepts `union_type` containing both `"partnered"`
#' and `"unpartnered"`, or all three `"unpartnered"`, `"cohab"`, `"marriage"`.
#' Shared partnered rates cover either partnership type. Rates use woman-years
#' in the indicated type; the old fertility input `union_status` is rejected.
#' Without that column, one overall schedule applies to all women. A partial
#' status table or any other label is rejected; use explicit zero rates to
#' suppress births in a status. Overlap is checked separately within each status.
#' On-demand union formation requires overall rates; rate-based formation accepts
#' either form and warns when an overall schedule is used.
#'
#' Optional `parity_start`/`parity_end` use half-open current maternal parity
#' bands, covering zero through Inf without gaps or overlap in every group/type
#' stratum. Parity is counted from recorded maternal births (including deceased
#' children) and incremented only on successful births; no custom parity column
#' overrides history. A pooled parity band applies its rate to each contained
#' actual parity.
#'
#' Birth spacing uses cell-accounted fertility (D101). Women are grouped into
#' cells: one-year (or shorter) age sub-bands of each group/union-type schedule
#' after group multipliers, split at whole ages and band edges, and one actual
#' parity when rates are parity-specific. Each cell produces births at its input
#' rate times all its living members, including those inside the birth
#' interval, and each birth goes to an eligible member, uniformly by default
#' or proportional to persistent weights when `hetfert` is enabled. Observed
#' rates then equal inputs in expectation for any population composition,
#' including after restarts and union-type changes, while a cell has an eligible
#' member. A cell with members but none eligible keeps its intensity; when it
#' would produce a birth, no child is born and the event is counted in
#' `simulation$unrealized_births`. An individual woman's hazard
#' depends on recent births by others in her cell; this is a calibration
#' device. Upper-open age bands stay one cell when `birth_interval * rate` is
#' below 0.001 and are otherwise split yearly. Zero-rate bands form no cells.
#'
#' For parity-specific schedules with a positive `birth_interval`, `fertility()`
#' warns when positive rates fall in one-year age cells where no woman can be
#' eligible: births are singletons at least `birth_interval` apart, starting
#' from the schedule's first positive rate, so parity k is first eligible at the
#' age found by following each parity's first positive rate. Union types are
#' pooled and groups checked separately. Partly reachable cells and higher
#' actual parities within an open parity band (always unreachable at the
#' youngest ages) are not reported.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`. Starts must be finite and nonnegative; ends must exceed starts
#'   and may be `Inf` for upper-open bands. Rates must be finite and nonnegative,
#'   per woman-year. Use [stratify_rates()] to add proportional group differences.
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative; see cell-accounted
#'   fertility above.
#' @param hetfert NULL (default) for uniform mother selection, or a
#'   [het_fert()] specification for persistent lognormal selection weights.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [create_weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate. They
#'   combine with `hetfert` on the log scale.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`.
#' @param model A `hammel_model` object, used by `set_fertility()`.
#' @return `fertility()` returns an inspectable process inheriting from
#'   `hammel_fertility` and `hammel_process`, with `name`, `rates`, `rate_unit`,
#'   `birth_interval`, `ffab`, `hetfert`, and `weights`.
#'   `set_fertility()` returns an updated `hammel_model`.
#' @seealso [set_process()], [mortality()], [union_formation()],
#'   [het_fert()], [create_group_inheritance()]
#' @export
#' @examples
#' rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
#' process <- fertility(rates)
#' model <- hammel_model() |> set_process(process)
#' # Equivalent convenience setter:
#' model <- hammel_model() |> set_fertility(rates)
#' @order 1
fertility <- function(rates, birth_interval = 0.75, ffab = 0.4886,
                      hetfert = NULL, weights = NULL) {
  validate_rates(rates, mortality = FALSE)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(weights, weight_strata("fertility", rates), "fertility")
  validate_birth_interval(birth_interval)
  validate_ffab(ffab)
  validate_het_fert(hetfert)
  process <- new_rate_process("fertility", rates)
  process$birth_interval <- birth_interval
  process$ffab <- ffab
  process$hetfert <- hetfert
  process$weights <- normalize_model_argument(weights)
  warn_unreachable_parity_cells(process$rates, birth_interval)
  process
}

#' @rdname fertility
#' @export
#' @order 2
set_fertility <- function(model, rates, birth_interval = 0.75, ffab = 0.4886,
                          hetfert = NULL, weights = NULL) {
  set_process(model, fertility(rates, birth_interval = birth_interval,
                              ffab = ffab, hetfert = hetfert, weights = weights))
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Positive parity-specific rates in one-year age cells (the fertility cells of
# D101) that no woman can occupy while eligible, given the
# birth interval: singleton births, each at least `interval` apart, starting at
# the schedule's first positive rate. Union types are pooled, because women move
# between them; groups are checked separately. Pooled parity bands are checked
# at their first parity only: higher actual parities in a pooled band reuse one
# schedule and always have unreachable young cells. Cells reachable for part of
# the year are not reported. Returns one row per unreachable age range.
unreachable_parity_cells <- function(rates, interval) {
  empty <- data.frame(group = character(), parity = numeric(), age_start = numeric(),
                      age_end = numeric(), earliest_eligible = numeric())
  if (!all(c("parity_start", "parity_end") %in% names(rates)) || !(interval > 0)) return(empty)
  out <- list()
  # Earliest next birth from eligibility age `e` given positive rows `r`.
  next_birth <- function(r, e) {
    ok <- r$age_end > e
    if (any(ok)) min(pmax(r$age_start[ok], e)) else Inf
  }
  for (rows in split_rate_groups(rates)) {
    r <- rates[rows, , drop = FALSE]
    r <- r[r$rate > 0, , drop = FALSE]
    if (!nrow(r)) next
    group <- if ("group" %in% names(r)) as.character(r$group[1]) else NA_character_
    starts <- sort(unique(rates$parity_start[rows]))
    ends <- vapply(starts, function(k) rates$parity_end[rows][match(k, rates$parity_start[rows])], 0)
    eligible <- 0 # earliest age at which a woman at the current parity is eligible
    for (b in seq_along(starts)) {
      k <- starts[b]
      at_k <- r[r$parity_start == k, , drop = FALSE]
      if (k > 0) {
        # Report the band's first parity only (pooled bands apply one schedule).
        for (i in seq_len(nrow(at_k))) {
          lo <- at_k$age_start[i]
          if (!is.finite(eligible)) {
            out[[length(out) + 1L]] <- data.frame(group = group, parity = k,
              age_start = lo, age_end = at_k$age_end[i], earliest_eligible = Inf)
            next
          }
          # Cells split at whole ages; report the run of cells ending no later
          # than the earliest eligible age as one range.
          end <- if (at_k$age_end[i] <= eligible) at_k$age_end[i] else floor(eligible)
          if (end > lo) out[[length(out) + 1L]] <- data.frame(group = group,
            parity = k, age_start = lo, age_end = end, earliest_eligible = eligible)
        }
      }
      # Carry eligibility through every parity of this band to the next band.
      if (b < length(starts)) {
        steps <- ends[b] - k
        while (steps > 0 && is.finite(eligible)) {
          birth <- next_birth(at_k, eligible)
          if (!is.finite(birth)) { eligible <- Inf; break }
          inside <- at_k$age_start <= birth & birth < at_k$age_end
          # Within an upper-open positive row every later step adds the interval.
          if (any(inside & is.infinite(at_k$age_end))) {
            eligible <- birth + steps * interval
            break
          }
          eligible <- birth + interval
          steps <- steps - 1
        }
      }
    }
  }
  if (length(out)) do.call(rbind, out) else empty
}

warn_unreachable_parity_cells <- function(rates, interval) {
  # The static reachability heuristic assumes the same schedules throughout a
  # life. Applying it to calendar slices could falsely reject reachable histories.
  if (has_calendar_rates(rates)) return(invisible(NULL))
  cells <- unreachable_parity_cells(rates, interval)
  if (!nrow(cells)) return(invisible(cells))
  keys <- unique(cells[c("group", "parity")])
  parts <- vapply(seq_len(nrow(keys)), function(i) {
    d <- cells[(is.na(cells$group) & is.na(keys$group[i]) | cells$group %in% keys$group[i]) &
               cells$parity == keys$parity[i], ]
    sprintf("%sparity %g at ages [%g, %g)%s",
            if (is.na(keys$group[i])) "" else paste0("group ", keys$group[i], ", "),
            keys$parity[i], min(d$age_start), max(d$age_end),
            if (is.finite(d$earliest_eligible[1]))
              sprintf(" (earliest eligible age %g)", d$earliest_eligible[1])
            else " (unreachable at any age)")
  }, character(1))
  shown <- head(parts, 5)
  if (length(parts) > 5) shown <- c(shown, sprintf("and %d more", length(parts) - 5))
  warning(sprintf(paste0("With `birth_interval` = %g, some positive parity-specific ",
    "fertility rates fall in one-year age cells where no woman can be eligible, so ",
    "those births cannot occur: %s."), interval, paste(shown, collapse = "; ")),
    call. = FALSE)
  invisible(cells)
}

validate_ffab <- function(probability) {
  if (!is.numeric(probability) || !is.null(dim(probability)) ||
      length(probability) != 1L || !is.finite(probability) ||
      probability < 0 || probability > 1) {
    stop("`ffab` must be a finite numeric probability in [0, 1].", call. = FALSE)
  }
  invisible(probability)
}

validate_birth_interval <- function(duration) {
  if (!is.numeric(duration) || !is.null(dim(duration)) || length(duration) != 1L ||
      !is.finite(duration) || duration < 0) {
    stop("`birth_interval` must be a finite nonnegative number of years.",
         call. = FALSE)
  }
  invisible(duration)
}
