# Shared calendar validation for flat rate tables. Each active slice follows
# the existing demographic validator; time gaps are checked against the run.
has_calendar_rates <- function(rates) {
  any(c("time_start", "time_end") %in% names(rates))
}

calendar_rate_slices <- function(rates) {
  if (!has_calendar_rates(rates)) return(list(rates))
  check_columns(rates, c("time_start", "time_end"), "rates")
  check_numeric_column(rates$time_start, "rates$time_start")
  ends <- rates$time_end
  if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) ||
      anyNA(ends) || any(ends <= rates$time_start))
    stop("Calendar intervals need finite starts and greater ends (possibly Inf).", call. = FALSE)
  times <- sort(unique(c(rates$time_start, ends)))
  lapply(head(times, -1L), function(time) {
    rows <- rates$time_start <= time & time < rates$time_end
    rates[rows, setdiff(names(rates), c("time_start", "time_end")), drop = FALSE]
  })
}

validate_calendar_rates <- function(rates, validate) {
  if (!has_calendar_rates(rates)) return(FALSE)
  for (slice in calendar_rate_slices(rates)) if (nrow(slice)) validate(slice)
  TRUE
}

check_calendar_coverage <- function(rates, start, end, name) {
  if (is.null(rates) || !has_calendar_rates(rates)) return(invisible(NULL))
  calendar_rate_slices(rates) # validate bounds, including edited inputs
  times <- sort(unique(c(start, rates$time_start, rates$time_end)))
  times <- times[times >= start & (times < end | times == start)]
  for (time in times) {
    if (!any(rates$time_start <= time & time < rates$time_end))
      stop(sprintf("%s rates have a calendar coverage gap at time %g.", name, time), call. = FALSE)
  }
  invisible(NULL)
}
