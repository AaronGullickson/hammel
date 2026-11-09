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

# A summary table's stratum columns, in table order: period, sex (when kept),
# union status (when kept) and every other non-measure column.
summary_strata <- function(table, by_sex = FALSE, by_union_status = FALSE) {
  status <- c("union_status", "union_type", "previous_union")
  keep <- c(
    "time_start", "time_end", if (by_sex) "sex", if (
      by_union_status) {
      status
    },
    setdiff(names(table), c(
      summary_measures, "time_start", "time_end",
      "sex", status
    ))
  )
  names(table)[names(table) %in% keep]
}

summary_rows <- function(table, by_sex = FALSE, by_union_status = FALSE) {
  columns <- summary_strata(table, by_sex, by_union_status)
  if (!length(columns)) {
    return(list(seq_len(nrow(table))))
  }
  codes <- lapply(table[columns], function(column) match(column, unique(
    column
  )))
  keys <- do.call(paste, c(codes, sep = ":"))
  split(seq_len(nrow(table)), factor(keys, levels = unique(keys)))
}

# Validate one stratum's age bands. `strata` names the stratum columns so that
# errors identify the failing stratum.
check_summary_intervals <- function(rates, open, strata = character()) {
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
  if (open && (ordered$age_start[1] != 0 || tail(ordered$age_end, 1) != Inf)) {
    stop(
      paste0(
        "Life tables require ages from zero through an ",
        "explicit open interval", where, "."
      ),
      call. = FALSE
    )
  }
  if (!open && any(!is.finite(ends))) {
    stop(paste0("TFR requires finite age intervals", where, "."),
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
