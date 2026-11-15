# Identify strata and validate intervals for summaries of rate schedules.

# INTERNAL CONSTANTS -------------------------------------------------------

# Bands and measures of rate tables; summaries treat any other column (such
# as a group dimension, `race` or `region`) as identifying a stratum.
summary_measures <- c(
  "age_start", "age_end", "parity_start", "parity_end",
  "dur_start", "dur_end",
  "mx", "rate", "rate_mortality", "rate_fertility", "rate_fertility_adjusted",
  "person_years",
  "deaths", "births", "unrealized_births", "conf_low", "conf_high"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# A summary table's stratum columns, in table order: every column other than
# bands and measures, such as period, sex, union status and group dimensions.
summary_strata <- function(table) {
  setdiff(names(table), summary_measures)
}

summary_rows <- function(table) {
  columns <- summary_strata(table)
  if (!length(columns)) {
    return(list(seq_len(nrow(table))))
  }
  codes <- lapply(
    table[columns], function(column) match(column, unique(column))
  )
  keys <- do.call(paste, c(codes, sep = ":"))
  split(seq_len(nrow(table)), factor(keys, levels = unique(keys)))
}

# Validate that one stratum's age bands are valid, contiguous and
# nonoverlapping. `strata` names the stratum columns so that errors identify
# the failing stratum.
check_summary_intervals <- function(
  rates,
  strata = character(),
  call = rlang::caller_env()
) {
  where <- stratum_label(rates, strata)
  if (!nrow(rates)) {
    abort_no_bands(call = call)
  }
  check_numeric_column(rates$age_start, "rates$age_start", call = call)
  ends <- rates$age_end
  if (!is.numeric(ends) || anyNA(ends) || any(rates$age_start < 0) ||
        any(ends <= rates$age_start)) {
    abort_hammel(
      c(
        paste(
          "Age bands in {.arg rates} must have nonnegative starts and",
          "greater ends."
        ),
        "x" = "Some bands{where} are invalid.",
        "i" = "Use {.code Inf} for an open-ended band."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  ordered <- rates[order(rates$age_start), , drop = FALSE]
  if (nrow(ordered) > 1L &&
        any(head(ordered$age_end, -1) != tail(ordered$age_start, -1))) {
    bands <- format_band(ordered$age_start, ordered$age_end)
    shown <- cli::cli_vec(bands, list("vec-trunc" = 8))
    abort_hammel(
      c(
        "Age bands in {.arg rates} must meet without gaps or overlaps.",
        "x" = "The bands{where} are {shown}."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
}

# A summary rate table must have at least one band.
abort_no_bands <- function(call = rlang::caller_env()) {
  abort_hammel(
    c(
      "{.arg rates} must contain at least one age band.",
      "x" = "It has no rows."
    ),
    "rates",
    argument = "rates",
    call = call
  )
}

# Describe a stratum for error messages, as " in stratum sex = f, region = a
# (columns other than age bands, rates and counts define strata)", or "" when
# the table has no stratum columns.
stratum_label <- function(rates, strata) {
  strata <- intersect(strata, names(rates))
  if (!length(strata) || !nrow(rates)) {
    return("")
  }
  values <- vapply(
    strata,
    function(column) as.character(rates[[column]][1]),
    character(1)
  )
  paste0(
    " in stratum ", paste0(strata, " = ", values, collapse = ", "),
    " (columns other than age bands, rates and counts define strata)"
  )
}

# A summary function received a process of the wrong kind.
abort_process_kind <- function(kind, call = rlang::caller_env()) {
  abort_hammel(
    c(
      "{.arg rates} must be a {kind} rate table or {kind} process.",
      "i" = "Pass the rates themselves, or the process made by {.fn {kind}}."
    ),
    "argument",
    argument = "rates",
    call = call
  )
}

# A summary rate table needs one of the named rate columns.
abort_rate_column <- function(columns, call = rlang::caller_env()) {
  abort_hammel(
    "{.arg rates} needs a {.or {.field {columns}}} column.",
    "table",
    argument = "rates",
    columns = columns,
    call = call
  )
}

# TFR and PATFR exclude open-ended age bands.
warn_open_fertility_band <- function(measure) {
  warn_hammel(
    c(
      "{measure} omits an open-ended age band with nonzero fertility.",
      "i" = "{measure} is underestimated over the full age range."
    ),
    "estimate"
  )
}

abort_no_finite_bands <- function(measure, call = rlang::caller_env()) {
  abort_hammel(
    c(
      "{measure} needs at least one female age band with a finite end.",
      "i" = "Open-ended bands are excluded from {measure}."
    ),
    "rates",
    argument = "rates",
    call = call
  )
}
