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
    subject <- describe_column(name)
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
      shown_values <- cli::cli_vec(values, list("vec-trunc" = 5))
      value_noun <- if (length(rows) == 1L) "value" else "values"
      if (all(is.na(values) & !is.nan(values))) {
        "Missing at {format_rows(rows)}."
      } else {
        "Invalid {value_noun} {.val {shown_values}} at {format_rows(rows)}."
      }
    }
    abort_hammel(
      c("{subject} must contain {requirement}.", "x" = problem),
      "table",
      argument = name,
      rows = rows,
      call = call
    )
  }
}

# Numeric values (already checked) must not be negative; NA is ignored.
check_nonnegative <- function(x, name, call = rlang::caller_env()) {
  rows <- which(!is.na(x) & x < 0)
  if (length(rows)) {
    subject <- describe_column(name)
    shown <- cli::cli_vec(x[rows], list("vec-trunc" = 5))
    abort_hammel(
      c(
        "{subject} must be nonnegative.",
        "x" = "Negative {.val {shown}} at {format_rows(rows)}."
      ),
      "table",
      argument = name,
      rows = rows,
      call = call
    )
  }
  invisible(x)
}

# Character labels drawn from `allowed`, without missing values.
check_labels <- function(x, name, allowed, call = rlang::caller_env()) {
  subject <- describe_column(name)
  if (!is.character(x) || !is.null(dim(x))) {
    abort_hammel(
      c(
        "{subject} must contain {.or {.val {allowed}}}.",
        "x" = "It is {describe_value(x)}."
      ),
      "table",
      argument = name,
      call = call
    )
  }
  missing_rows <- which(is.na(x))
  rows <- if (length(missing_rows)) missing_rows else which(!x %in% allowed)
  if (length(rows)) {
    problem <- if (length(missing_rows)) {
      "Missing at {format_rows(rows)}."
    } else {
      shown <- cli::cli_vec(unique(x[rows]), list("vec-trunc" = 5))
      "Invalid {.val {shown}} at {format_rows(rows)}."
    }
    abort_hammel(
      c("{subject} must contain {.or {.val {allowed}}}.", "x" = problem),
      "table",
      argument = name,
      rows = rows,
      call = call
    )
  }
  invisible(x)
}

# A single TRUE or FALSE setting.
check_flag <- function(x, name, call = rlang::caller_env()) {
  if (!is.logical(x) || length(x) != 1L || !is.null(dim(x)) || is.na(x)) {
    abort_hammel(
      c(
        "{.arg {name}} must be {.code TRUE} or {.code FALSE}.",
        "x" = "It is {describe_value(x)}."
      ),
      "argument",
      argument = name,
      call = call
    )
  }
  invisible(x)
}
