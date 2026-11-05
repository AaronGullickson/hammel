# Cohabitation input helpers: no demographic defaults beyond the accepted
# zero probability in gaps and shared matching/rates when explicitly unstratified.
validate_cohabitation <- function(process) {
  flag <- process$cohabitation_counts_as_previous_union
  if (!is.null(flag) && (!is.logical(flag) || length(flag) != 1L ||
      !is.null(dim(flag)) || is.na(flag)))
    stop("`cohabitation_counts_as_previous_union` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(process$cohab_matching))
    validate_model_argument(process$cohab_matching, validate_matching_model, "cohab_matching")
  enabled <- !is.null(process$cohab_prob)
  if (!enabled && (!is.null(process$cohab_matching) || !is.null(process$union_convert_rate)))
    stop("`cohab_matching` and `union_convert_rate` require `cohab_prob`.", call. = FALSE)
  if (enabled) {
    if (!identical(process$initiation, "rate"))
      stop("Cohabitation is available only with rate-based initiation.", call. = FALSE)
    x <- process$cohab_prob
    check_columns(x, c("age_start", "age_end", "prob"), "cohab_prob")
    check_numeric_column(x$prob, "cohab_prob$prob")
    if (any(x$prob < 0 | x$prob > 1) ||
        any(c("rate", "dur_start", "dur_end") %in% names(x)))
      stop("`cohab_prob` needs age-specific probabilities in [0, 1], with `prob`, not `rate`.", call. = FALSE)
    x$rate <- x$prob
    validate_union_search_rates(x)
  }
  if (!is.null(process$union_convert_rate)) validate_conversion_rates(process$union_convert_rate)
  invisible(process)
}

validate_conversion_rates <- function(rates) {
  check_stratified_columns(rates, "group", "Conversion")
  if (has_calendar_rates(rates)) {
    if (all(c("age_start", "dur_start") %in% names(rates)) &&
        length(unique(!is.na(rates$dur_start))) != 1L)
      stop("Conversion rates must use one age or duration clock across calendar periods.", call. = FALSE)
    validate_calendar_rates(rates, validate_conversion_rates)
    return(invisible(rates))
  }
  check_columns(rates, "rate", "union_convert_rate")
  if (any(c("previous_union", "sex", "time_start", "time_end", "union_type") %in% names(rates)))
    stop("Conversion rates allow age or duration and optional group, without previous_union, sex, type or calendar strata.", call. = FALSE)
  # Reuse search validation with a synthetic previous category for each row;
  # conversion uses one clock across the entire table, not a previous axis.
  x <- rates
  for (field in c("age_start", "age_end", "dur_start", "dur_end"))
    if (!field %in% names(x)) x[[field]] <- NA_real_
  duration <- !is.na(x$dur_start)
  if (length(unique(duration)) != 1L)
    stop("Conversion rates must use one age or duration clock across all groups.", call. = FALSE)
  x$previous_union <- "separated"
  # Fill unused synthetic categories with copies to exercise existing coverage.
  x <- do.call(rbind, lapply(c("none", "separated", "widowed"), function(label) {
    copy <- x; copy$previous_union <- label
    if (label == "none" && all(duration)) {
      copy$age_start <- copy$dur_start; copy$age_end <- copy$dur_end
      copy$dur_start <- copy$dur_end <- NA_real_
    }
    copy
  }))
  validate_union_search_rates(x)
  invisible(rates)
}

prepare_cohabitation <- function(process, population, groups, start = -Inf, end = Inf) {
  empty <- data.frame(age_start = numeric(), age_end = numeric(), rate = numeric())
  if (is.null(process$cohab_prob)) return(list(enabled = FALSE, prob = empty, conversion = empty))
  prob <- prepare_group_rates(expand_previous_union_rates(process$cohab_prob),
    population$persons, "cohabitation probability", groups, start, end)
  prob$rate <- prob$prob # private Band convention
  conversion <- process$union_convert_rate
  if (is.null(conversion)) {
    warning("Cohabitation is enabled without `union_convert_rate`; no conversions will occur.", call. = FALSE)
    conversion <- empty
  } else {
    conversion <- prepare_group_rates(conversion, population$persons, "cohabitation conversion", groups, start, end)
    if (!"age_start" %in% names(conversion)) conversion$age_start <- NA_real_
    if (!"age_end" %in% names(conversion)) conversion$age_end <- NA_real_
  }
  list(enabled = TRUE, prob = prob, conversion = conversion)
}

validate_union_types <- function(unions) {
  if (!"union_type" %in% names(unions)) return(invisible(unions))
  x <- unions$union_type
  if (!is.character(x) || !is.null(dim(x)) || anyNA(x) || any(!x %in% c("cohab", "marriage")))
    stop("`unions$union_type` must contain 'cohab' or 'marriage'.", call. = FALSE)
  if (!"previous_cohab_id" %in% names(unions)) {
    if (any(unions$end_reason == "conversion", na.rm = TRUE))
      stop("Conversion histories require `previous_cohab_id` links.", call. = FALSE)
    return(invisible(unions))
  }
  links <- unions$previous_cohab_id
  check_numeric_column(links, "previous_cohab_id", missing = TRUE, id = TRUE)
  rows <- which(!is.na(links)); parent <- match(links[rows], unions$union_id)
  if (anyNA(parent) || anyDuplicated(parent) || any(unions$union_type[rows] != "marriage") ||
      any(unions$union_type[parent] != "cohab") ||
      any(is.na(unions$end_reason[parent]) | unions$end_reason[parent] != "conversion") ||
      any(is.na(unions$end_time[parent]) | unions$end_time[parent] != unions$start_time[rows]) ||
      any(unions$wife_id[parent] != unions$wife_id[rows] | unions$husband_id[parent] != unions$husband_id[rows]) ||
      any(!which(unions$end_reason == "conversion") %in% parent))
    stop("Conversions require unique same-partner cohab-to-marriage links at the same time.", call. = FALSE)
  invisible(unions)
}
