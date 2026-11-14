# Provide general checks for table columns, numeric values, and times.

# HELPER FUNCTIONS ---------------------------------------------------------

# General validation helpers shared across package components.
# Domain-specific checks stay with the API they validate. Each check raises
# through abort_hammel() and names the function in `call` (D183).

check_time <- function(time, argument, call = rlang::caller_env()) {
  if (!is.numeric(time) || is.complex(time) || !is.null(dim(time)) ||
        length(time) != 1L ||
        !is.finite(time)) {
    abort_hammel(
      c(
        "{.arg {argument}} must be a single finite number of years.",
        "x" = "It is {describe_value(time)}."
      ),
      "argument",
      argument = argument,
      call = call
    )
  }
}

check_columns <- function(table, columns, name, call = rlang::caller_env()) {
  if (!is.data.frame(table)) {
    abort_hammel(
      c(
        "{.arg {name}} must be a data frame.",
        "x" = "It is {describe_value(table)}."
      ),
      "table",
      argument = name,
      call = call
    )
  }
  duplicated_names <- unique(names(table)[duplicated(names(table))])
  if (length(duplicated_names)) {
    abort_hammel(
      c(
        "{.arg {name}} must have unique column names.",
        "x" = "Duplicated: {.field {duplicated_names}}."
      ),
      "table",
      argument = name,
      columns = duplicated_names,
      call = call
    )
  }
  absent <- setdiff(columns, names(table))
  if (length(absent)) {
    abort_hammel(
      c(
        "{.arg {name}} is missing required columns.",
        "x" = "Missing: {.field {absent}}."
      ),
      "table",
      argument = name,
      columns = absent,
      call = call
    )
  }
}

# `name` is an argument (`end_time`) or a column written as `table$column`.
check_numeric_column <- function(
  x,
  name,
  missing = FALSE,
  id = FALSE,
  call = rlang::caller_env()
) {
  observed <- x[!is.na(x)]
  if (!is.numeric(x) || !is.null(dim(x)) || any(is.nan(x)) ||
        (!missing && anyNA(x)) ||
        any(!is.finite(observed)) ||
        (id && any(observed <= 0 | observed != trunc(observed)))) {
    requirement <- if (id) "positive whole-number IDs" else "finite numbers"
    if (missing) {
      requirement <- paste(requirement, "or NA")
    }
    parts <- strsplit(name, "$", fixed = TRUE)[[1]]
    subject <- if (length(parts) == 2L) {
      "Column {.field {parts[2]}} of {.arg {parts[1]}}"
    } else {
      "{.arg {name}}"
    }
    rows <- integer()
    problem <- if (!is.numeric(x) || !is.null(dim(x))) {
      "It is {describe_value(x)}."
    } else {
      bad <- is.nan(x) | (!missing & is.na(x)) | (!is.na(x) & !is.finite(x))
      if (id) {
        bad <- bad | (!is.na(x) & is.finite(x) & (x <= 0 | x != trunc(x)))
      }
      rows <- which(bad)
      values <- x[rows]
      shown_rows <- cli::cli_vec(rows, list("vec-trunc" = 5))
      shown_values <- cli::cli_vec(values, list("vec-trunc" = 5))
      row_noun <- if (length(rows) == 1L) "row" else "rows"
      value_noun <- if (length(rows) == 1L) "value" else "values"
      if (all(is.na(values) & !is.nan(values))) {
        "Missing at {row_noun} {shown_rows}."
      } else {
        "Invalid {value_noun} {.val {shown_values}} at {row_noun} {shown_rows}."
      }
    }
    abort_hammel(
      c(paste(subject, "must contain {requirement}."), "x" = problem),
      "table",
      argument = name,
      rows = rows,
      call = call
    )
  }
}
