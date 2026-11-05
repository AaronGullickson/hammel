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

prepare_group_rates <- function(rates, persons, process,
                                group_values = if ("group" %in% names(persons)) unique(persons$group) else NULL,
                                start = -Inf, end = Inf) {
  rates$.hammel_group <- NULL
  if (!"group" %in% names(rates)) return(rates)
  # Simulation has already rejected group rates without a declared group.
  labels <- as.character(rates$group)
  groups <- as.character(group_values)
  missing <- unique(unlist(lapply(calendar_rate_slices(rates, start, end), function(x)
    if (nrow(x)) setdiff(groups, as.character(x$group)) else character())))
  if (length(missing)) {
    stop(sprintf("%s rates are missing schedules for population group(s): %s.",
                 process, paste(missing, collapse = ", ")), call. = FALSE)
  }
  rates$.hammel_group <- match(labels, groups) - 1L
  rates
}
