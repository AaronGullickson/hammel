# Shared group validation and execution-only schedule indexing. Public labels
# and their types are preserved; the engine uses compact integer indices.
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

prepare_group_rates <- function(rates, persons, process) {
  rates$.hammel_group <- NULL
  if (!"group" %in% names(rates)) return(rates)
  labels <- as.character(rates$group)
  if (!"group" %in% names(persons)) {
    first <- sort(unique(labels), method = "radix")[1L]
    warning(sprintf("Group rates were specified for %s but no groups found in the population; using group '%s' (first alphabetically).",
                    process, first), call. = FALSE)
    return(rates[labels == first, , drop = FALSE])
  }
  groups <- unique(as.character(persons$group))
  missing <- setdiff(groups, labels)
  if (length(missing)) {
    stop(sprintf("%s rates are missing schedules for population group(s): %s.",
                 process, paste(missing, collapse = ", ")), call. = FALSE)
  }
  rates$.hammel_group <- match(labels, groups) - 1L
  rates
}
