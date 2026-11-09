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
  codes <- lapply(table[columns], function(column) match(column, unique(
    column
  )))
  keys <- do.call(paste, c(codes, sep = ":"))
  split(seq_len(nrow(table)), factor(keys, levels = unique(keys)))
}

# Validate that one stratum's age bands are valid, contiguous and
# nonoverlapping. `strata` names the stratum columns so that errors identify
# the failing stratum.
check_summary_intervals <- function(rates, strata = character()) {
  where <- stratum_label(rates, strata)
  if (!nrow(rates)) {
    stop("A rate table needs at least one age band.",
      call. =
        FALSE
    )
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  ends <- rates$age_end
  if (!is.numeric(ends) || anyNA(ends) || any(rates$age_start < 0) ||
    any(ends <= rates$age_start)) {
    stop(paste0("Invalid age intervals", where, "."), call. = FALSE)
  }
  ordered <- rates[order(rates$age_start), , drop = FALSE]
  if (nrow(ordered) > 1L && any(head(ordered$age_end, -1) != tail(
    ordered$age_start, -1
  ))) {
    stop(
      paste0(
        "Age intervals must be contiguous and nonoverlapping",
        where, "."
      ),
      call. = FALSE
    )
  }
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
