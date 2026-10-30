# Shared group validation and execution-only schedule indexing. Public labels
# and their types are preserved; the engine uses compact integer indices.
validate_group_prop_hazard <- function(x) {
  if (is.null(x)) return(invisible(x))
  labels <- names(x)
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) ||
      anyNA(x) || any(!is.finite(x)) || any(x < 0) ||
      is.null(labels) || anyNA(labels) || any(!nzchar(trimws(labels))) ||
      anyDuplicated(labels)) {
    stop("`group_prop_hazard` must be NULL or a named numeric vector of finite nonnegative multipliers with unique nonempty group names.",
         call. = FALSE)
  }
  invisible(x)
}

validate_fertility_group_prop_hazard <- function(x, rates) {
  if (!is.list(x)) return(validate_group_prop_hazard(x))
  labels <- names(x)
  if (!length(x) || is.null(labels) || anyNA(labels) ||
      anyDuplicated(labels) || any(!labels %in% c("partnered", "unpartnered", "cohab", "marriage"))) {
    stop("A fertility `group_prop_hazard` list must have unique union_type names: unpartnered, partnered, cohab or marriage.",
         call. = FALSE)
  }
  if (!"union_type" %in% names(rates)) {
    stop("A fertility `group_prop_hazard` status list requires rates with `union_type`; use a named vector for overall fertility.",
         call. = FALSE)
  }
  if (any(!labels %in% rates$union_type))
    stop("Fertility multiplier list names must occur in union_type schedules.", call. = FALSE)
  for (factors in x) {
    if (is.null(factors)) {
      stop("Each fertility `group_prop_hazard` list entry must be a named numeric vector; omit an unchanged status.",
           call. = FALSE)
    }
    validate_group_prop_hazard(factors)
  }
  invisible(x)
}

check_groups <- function(table, name) {
  if (!"group" %in% names(table)) return(invisible(table))
  x <- table$group
  if (!(is.character(x) || is.numeric(x)) ||
      !is.null(dim(x)) || anyNA(x) ||
      any(!nzchar(trimws(as.character(x)))) ||
      (is.numeric(x) && any(!is.finite(x)))) {
    stop(sprintf("`%s$group` must contain nonmissing, nonempty group labels.", name),
         call. = FALSE)
  }
  invisible(table)
}

split_rate_groups <- function(rates) {
  check_groups(rates, "rates")
  if ("group" %in% names(rates))
    split(seq_len(nrow(rates)), as.character(rates$group), drop = TRUE)
  else list(seq_len(nrow(rates)))
}

prepare_group_rates <- function(rates, persons, process,
                                group_values = if ("group" %in% names(persons)) unique(persons$group) else NULL) {
  rates$.hammel_group <- NULL
  if (!"group" %in% names(rates)) return(rates)
  labels <- as.character(rates$group)
  if (!"group" %in% names(persons)) {
    first <- sort(unique(labels), method = "radix")[1L]
    if (any(vapply(calendar_rate_slices(rates), function(x)
        nrow(x) > 0L && !first %in% as.character(x$group), logical(1))))
      stop(sprintf("%s rates are missing schedules for fallback group '%s' in a calendar interval.",
                   process, first), call. = FALSE)
    warning(sprintf("Group rates were specified for %s but no groups found in the population; using group '%s' (first alphabetically).",
                    process, first), call. = FALSE)
    return(rates[labels == first, , drop = FALSE])
  }
  groups <- as.character(group_values)
  missing <- unique(unlist(lapply(calendar_rate_slices(rates), function(x)
    if (nrow(x)) setdiff(groups, as.character(x$group)) else character())))
  if (length(missing)) {
    stop(sprintf("%s rates are missing schedules for population group(s): %s.",
                 process, paste(missing, collapse = ", ")), call. = FALSE)
  }
  rates$.hammel_group <- match(labels, groups) - 1L
  rates
}

# Rates have already passed group coverage/fallback checks. Expand shared
# schedules only on this execution copy; preserve the user's process inputs.
apply_group_prop_hazard <- function(rates, persons, multipliers, process,
                                    group_values = if ("group" %in% names(persons)) unique(persons$group) else NULL) {
  if (is.null(multipliers)) return(rates)
  if (!"group" %in% names(persons)) {
    warning(sprintf("`group_prop_hazard` was specified for %s but no groups found in the population; ignoring multipliers.",
                    process), call. = FALSE)
    return(rates)
  }
  groups <- as.character(group_values)
  if (!length(groups)) return(rates)
  if (!".hammel_group" %in% names(rates)) {
    rows <- nrow(rates)
    rates <- rates[rep(seq_len(rows), times = length(groups)), , drop = FALSE]
    rates$.hammel_group <- rep(seq_along(groups) - 1L, each = rows)
  }
  labels <- groups[rates$.hammel_group + 1L]
  if (is.list(multipliers)) {
    factors <- rep(1, nrow(rates))
    for (status in names(multipliers)) {
      rows <- rates$union_type == status
      vector <- multipliers[[status]]
      factors[rows] <- vector[match(labels[rows], names(vector))]
    }
  } else {
    factors <- multipliers[match(labels, names(multipliers))]
  }
  factors[is.na(factors)] <- 1
  rates$rate <- rates$rate * factors
  if (any(!is.finite(rates$rate))) {
    stop(sprintf("%s rates multiplied by `group_prop_hazard` must remain finite.", process),
         call. = FALSE)
  }
  rates
}
