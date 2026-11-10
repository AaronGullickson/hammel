# Validate and prepare cohabitation settings and conversion rate schedules.

# HELPER FUNCTIONS ---------------------------------------------------------

# Cohabitation input helpers: no demographic defaults beyond the accepted
# zero probability in gaps and shared matching/rates when explicitly
# unstratified.
validate_cohabitation <- function(process) {
  flag <- process$cohabitation_counts_as_previous_union
  if (!is.null(flag) &&
        (!is.logical(flag) || length(flag) != 1L || !is.null(dim(flag)) ||
           is.na(flag))) {
    stop(
      "`cohabitation_counts_as_previous_union` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  if (!is.null(process$cohab_matching)) {
    validate_model_argument(
      process$cohab_matching, validate_matching_model,
      "cohab_matching"
    )
  }
  enabled <- !is.null(process$cohab_prob)
  if (!enabled &&
        (!is.null(process$cohab_matching) ||
           !is.null(process$union_convert_rate))) {
    stop(
      "`cohab_matching` and `union_convert_rate` require `cohab_prob`.",
      call. = FALSE
    )
  }
  if (enabled) {
    if (!identical(process$initiation, "rate")) {
      stop(
        "Cohabitation is available only with rate-based initiation.",
        call. = FALSE
      )
    }
    x <- process$cohab_prob
    check_columns(x, c("age_start", "age_end", "prob"), "cohab_prob")
    check_numeric_column(x$prob, "cohab_prob$prob")
    if (any(x$prob < 0 | x$prob > 1) ||
          any(c("rate", "dur_start", "dur_end") %in% names(x))) {
      stop(
        paste0(
          "`cohab_prob` needs age-specific probabilities in ",
          "[0, 1], with `prob`, not `rate`."
        ),
        call. = FALSE
      )
    }
    x$rate <- x$prob
    validate_union_search_rates(x)
  }
  if (!is.null(process$union_convert_rate)) {
    validate_conversion_rates(process$union_convert_rate)
  }
  invisible(process)
}

validate_conversion_rates <- function(rates) {
  check_stratified_columns(rates, character(), "Conversion")
  if (has_calendar_rates(rates)) {
    if (all(c("age_start", "dur_start") %in% names(rates)) &&
          length(unique(!is.na(rates$dur_start))) != 1L) {
      stop(
        paste0(
          "Conversion rates must use one age or duration ",
          "clock across calendar periods."
        ),
        call. = FALSE
      )
    }
    validate_calendar_rates(rates, validate_conversion_rates)
    return(invisible(rates))
  }
  check_columns(rates, "rate", "union_convert_rate")
  if (any(
    c("previous_union", "sex", "time_start", "time_end", "union_type") %in%
      names(rates)
  )) {
    stop(
      paste0(
        "Conversion rates allow age or duration and ",
        "optional `<group>_w`/`<group>_h` strata, without ",
        "previous_union, sex, type or calendar strata."
      ),
      call. = FALSE
    )
  }
  # Reuse search validation with a synthetic previous category for each row;
  # conversion uses one clock across the entire table, not a previous axis.
  x <- rates
  for (field in c("age_start", "age_end", "dur_start", "dur_end")) {
    if (!field %in% names(x)) {
      x[[field]] <- NA_real_
    }
  }
  duration <- !is.na(x$dur_start)
  if (length(unique(duration)) != 1L) {
    stop(
      paste0(
        "Conversion rates must use one age or duration clock ",
        "across all groups."
      ),
      call. = FALSE
    )
  }
  x$previous_union <- "separated"
  # Fill unused synthetic categories with copies to exercise existing coverage.
  x <- do.call(
    rbind, lapply(
      c("none", "separated", "widowed"), function(
        label
      ) {
        copy <- x
        copy$previous_union <- label
        if (label == "none" && all(duration)) {
          copy$age_start <- copy$dur_start
          copy$age_end <- copy$dur_end
          copy$dur_start <- copy$dur_end <- NA_real_
        }
        copy
      }
    )
  )
  validate_union_search_rates(x)
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
  none <- strata(data.frame(rate = numeric()), "cohabitation")$spec
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
    "cohabitation probability"
  )
  prob$rates$rate <- prob$rates$prob # private Band convention
  conversion <- process$union_convert_rate
  if (is.null(conversion)) {
    warning(
      paste0(
        "Cohabitation is enabled without ",
        "`union_convert_rate`; no conversions will occur."
      ),
      call. = FALSE
    )
    conversion <- list(rates = empty, spec = none)
  } else {
    conversion <- strata(
      conversion, "cohabitation conversion",
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
