# Construct fertility processes and validate fertility settings and feasibility.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' @rdname mortality
#' @export
#' @order 2
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

#' @rdname set_process
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
