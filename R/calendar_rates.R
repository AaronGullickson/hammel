# Shared calendar validation for flat rate tables. Each active slice follows
# the existing demographic validator; time gaps are checked against the run.
has_calendar_rates <- function(rates) {
  any(c("time_start", "time_end") %in% names(rates))
}

calendar_rate_segments <- function(rates, start = -Inf, end = Inf) {
  if (!has_calendar_rates(rates))
    return(list(list(start = -Inf, end = Inf, rates = rates)))
  check_columns(rates, c("time_start", "time_end"), "rates")
  check_numeric_column(rates$time_start, "rates$time_start")
  ends <- rates$time_end
  if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) ||
      anyNA(ends) || any(ends <= rates$time_start))
    stop("Calendar intervals need finite starts and greater ends (possibly Inf).", call. = FALSE)
  times <- sort(unique(c(rates$time_start, ends)))
  indices <- seq_len(length(times) - 1L)
  keep <- if (start == end) times[indices] <= start & times[indices + 1L] > start else
    times[indices] < end & times[indices + 1L] > start
  lapply(indices[keep], function(i) {
    time <- times[i]
    rows <- rates$time_start <= time & time < rates$time_end
    list(start = time, end = times[i + 1L], rates =
      rates[rows, setdiff(names(rates), c("time_start", "time_end")), drop = FALSE])
  })
}

calendar_rate_slices <- function(rates, start = -Inf, end = Inf) {
  lapply(calendar_rate_segments(rates, start, end), `[[`, "rates")
}

validate_calendar_rates <- function(rates, validate, with_interval = FALSE) {
  if (!has_calendar_rates(rates)) return(FALSE)
  for (segment in calendar_rate_segments(rates)) if (nrow(segment$rates)) {
    if (with_interval) validate(segment$rates, c(segment$start, segment$end)) else
      validate(segment$rates)
  }
  TRUE
}

check_time_coverage <- function(starts, ends, start, end, description) {
  times <- sort(unique(c(start, starts, ends)))
  times <- times[times >= start & (times < end | times == start)]
  for (time in times) if (!any(starts <= time & time < ends))
    stop(sprintf("%s at time %g.", description, time), call. = FALSE)
  invisible(NULL)
}

check_calendar_coverage <- function(rates, start, end, name) {
  if (is.null(rates) || !has_calendar_rates(rates)) return(invisible(NULL))
  calendar_rate_slices(rates) # validate bounds, including edited inputs
  check_time_coverage(rates$time_start, rates$time_end, start, end,
    sprintf("%s rates have a calendar coverage gap", name))
}
