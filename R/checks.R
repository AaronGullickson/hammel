# Provide general checks for table columns, numeric values, and times.

# HELPER FUNCTIONS ---------------------------------------------------------

# General validation helpers shared across package components.
# Domain-specific checks stay with the API they validate.

check_time <- function(time, argument) {
  if (!is.null(time) &&
    (!is.numeric(time) || length(time) != 1L || !is.finite(time))) {
    stop(sprintf("`%s` must be NULL or a finite time in years.", argument),
      call. = FALSE
    )
  }
}

check_columns <- function(table, columns, name) {
  if (!is.data.frame(table)) {
    stop(sprintf("`%s` must be a data frame.", name), call. = FALSE)
  }
  if (anyDuplicated(names(table))) {
    stop(sprintf("`%s` must have unique column names.", name), call. = FALSE)
  }
  absent <- setdiff(columns, names(table))
  if (length(absent)) {
    stop(sprintf(
      "`%s` is missing columns: %s.", name,
      paste(absent, collapse = ", ")
    ), call. = FALSE)
  }
}

check_numeric_column <- function(x, name, missing = FALSE, id = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) ||
    any(is.nan(x)) || (!missing && anyNA(x)) ||
    any(!is.finite(x[!is.na(x)])) ||
    (id && any(x[!is.na(x)] <= 0 | x[!is.na(x)] != trunc(x[!is.na(x)])))) {
    stop(
      sprintf(
        "`%s` must contain %s%s.", name,
        if (id) "positive whole-number IDs" else
          "finite numeric values",
        if (missing) " or NA" else " without missing values"
      ),
      call. = FALSE
    )
  }
}
