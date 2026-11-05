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

# Row sets of each stratum of a rate table, by the given group columns.
split_rate_groups <- function(rates, columns = "group") {
  columns <- intersect(columns, names(rates))
  for (column in columns) check_groups(tibble::tibble(group = rates[[column]]), paste0("rates$", column))
  if (!length(columns)) return(list(seq_len(nrow(rates))))
  key <- do.call(paste, c(lapply(rates[columns], as.character), sep = "\r"))
  split(seq_len(nrow(rates)), key, drop = TRUE)
}

