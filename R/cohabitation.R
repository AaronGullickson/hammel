# Validate and prepare cohabitation settings and conversion rate schedules.

# HELPER FUNCTIONS ---------------------------------------------------------

# Cohabitation input helpers: no demographic defaults beyond the accepted
# zero probability in gaps and shared matching/rates when explicitly
# unstratified.
validate_cohabitation <- function(process, call = rlang::caller_env()) {
  if (!is.null(process$cohab_as_previous_union)) {
    check_flag(
      process$cohab_as_previous_union, "cohab_as_previous_union",
      call = call
    )
  }
  if (!is.null(process$cohab_matching)) {
    validate_model_argument(
      process$cohab_matching, validate_matching_model, "cohab_matching",
      call = call
    )
  }
  enabled <- !is.null(process$cohab_prob)
  if (!enabled) {
    supplied <- c("cohab_matching", "union_convert_rate")[
      c(!is.null(process$cohab_matching), !is.null(process$union_convert_rate))
    ]
    if (length(supplied)) {
      abort_hammel(
        c(
          "{.arg {supplied}} need{?s/} cohabitation to be enabled.",
          "i" = "Enable cohabitation with {.arg cohab_prob}."
        ),
        "specification",
        argument = supplied,
        call = call
      )
    }
    return(invisible(process))
  }
  if (!identical(process$initiation, "rate")) {
    abort_hammel(
      c(
        "Cohabitation needs rate-based union formation.",
        "x" = "{.arg initiation} is {.val {process$initiation}}.",
        "i" = paste(
          "Use {.code initiation = \"rate\"} with search {.arg rates}, or",
          "remove {.arg cohab_prob}."
        )
      ),
      "specification",
      argument = "cohab_prob",
      call = call
    )
  }
  x <- process$cohab_prob
  check_columns(x, c("age_start", "age_end", "prob"), "cohab_prob", call = call)
  bad <- intersect(c("rate", "dur_start", "dur_end"), names(x))
  if (length(bad)) {
    abort_hammel(
      c(
        "{.arg cohab_prob} can't have {.field {bad}} columns.",
        "i" = paste(
          "Give age-specific probabilities in a {.field prob} column with",
          "{.field age_start} and {.field age_end}."
        )
      ),
      "rates",
      argument = "cohab_prob",
      columns = bad,
      call = call
    )
  }
  check_numeric_column(x$prob, "cohab_prob$prob", call = call)
  rows <- which(x$prob < 0 | x$prob > 1)
  if (length(rows)) {
    abort_hammel(
      c(
        "Column {.field prob} of {.arg cohab_prob} must be between 0 and 1.",
        "x" = "Out of range at {format_rows(rows)}."
      ),
      "table",
      argument = "cohab_prob$prob",
      rows = rows,
      call = call
    )
  }
  x$rate <- x$prob
  validate_union_search_rates(x, "cohab_prob", "prob", call = call)
  if (!is.null(process$union_convert_rate)) {
    validate_conversion_rates(process$union_convert_rate, call = call)
  }
  invisible(process)
}

# Conversion rates use one clock (age or union duration) throughout, with
# optional `<group>_w`/`<group>_h` strata. Row checks run on the whole table,
# then overlap checks within each calendar period.
validate_conversion_rates <- function(rates, call = rlang::caller_env()) {
  name <- "union_convert_rate"
  check_stratified_columns(rates, character(), "Conversion", call = call)
  check_columns(rates, "rate", name, call = call)
  if (!nrow(rates)) {
    abort_hammel(
      c(
        "{.arg {name}} must contain at least one band.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = name,
      call = call
    )
  }
  calendar_rate_segments(rates, name = name, call = call)
  for (column in candidate_strata(rates)) {
    check_group_labels(
      rates[[column]], paste0(name, "$", column),
      call = call
    )
  }
  check_numeric_column(rates$rate, paste0(name, "$rate"), call = call)
  check_nonnegative(rates$rate, paste0(name, "$rate"), call = call)
  duration <- check_clock_rows(rates, name, call = call)
  if (length(unique(duration)) != 1L) {
    abort_hammel(
      c(
        paste(
          "{.arg {name}} must use one clock throughout: age bands or",
          "duration bands."
        ),
        "x" = "Age bands at {format_rows(which(!duration))}.",
        "x" = "Duration bands at {format_rows(which(duration))}."
      ),
      "rates",
      argument = name,
      call = call
    )
  }
  validate_conversion_schedules(
    rates, if (duration[1]) "dur" else "age",
    call = call
  )
  invisible(rates)
}

validate_conversion_schedules <- function(
  rates,
  axis,
  call = rlang::caller_env()
) {
  if (validate_calendar_rates(
        rates,
        function(x) {
          validate_conversion_schedules(x, axis, call = call)
        },
        call = call
      )) {
    return(invisible(rates))
  }
  strata <- candidate_strata(rates)
  for (rows in split_rate_groups(rates, strata, call = call)) {
    schedule <- if (length(strata)) {
      values <- vapply(
        strata,
        function(column) {
          as.character(rates[[column]][rows[1]])
        },
        character(1)
      )
      paste0(" for ", paste(strata, "=", values, collapse = ", "))
    } else {
      ""
    }
    check_band_overlap(
      rates, rows, axis, "union_convert_rate", schedule,
      call = call
    )
  }
  invisible(rates)
}

# Cohabitation probabilities (by the woman) and conversion rates (by both
# partners) with their group strata; `strata` prepares one rate table.
prepare_cohabitation <- function(process, strata) {
  empty <- data.frame(
    age_start = numeric(),
    age_end = numeric(),
    rate = numeric()
  )
  none <- strata(data.frame(rate = numeric()), "cohabitation rates")$spec
  if (is.null(process$cohab_prob)) {
    return(
      list(
        enabled = FALSE,
        prob = empty,
        conversion = empty,
        prob_spec = none,
        conversion_spec = none
      )
    )
  }
  prob <- strata(
    expand_previous_union_rates(process$cohab_prob),
    "cohabitation probabilities"
  )
  prob$rates$rate <- prob$rates$prob # private Band convention
  conversion <- process$union_convert_rate
  if (is.null(conversion)) {
    warn_hammel(
      c(
        "Cohabitation is enabled without {.arg union_convert_rate}.",
        "i" = "No cohabitations will convert to marriage."
      ),
      "no_effect"
    )
    conversion <- list(rates = empty, spec = none)
  } else {
    conversion <- strata(
      conversion, "cohabitation conversion rates",
      union_level = TRUE
    )
    if (!"age_start" %in% names(conversion$rates)) {
      conversion$rates$age_start <- NA_real_
    }
    if (!"age_end" %in% names(conversion$rates)) {
      conversion$rates$age_end <- NA_real_
    }
  }
  list(
    enabled = TRUE,
    prob = prob$rates,
    conversion = conversion$rates,
    prob_spec = prob$spec,
    conversion_spec = conversion$spec
  )
}
