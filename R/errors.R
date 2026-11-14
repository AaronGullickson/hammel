# Raise classified errors and warnings with cli messages (D183).

# INTERNAL CONSTANTS -------------------------------------------------------

# Specific condition kinds; see design/error-messages.md. Each condition is
# classed "hammel_error_<kind>" or "hammel_warning_<kind>" plus the base
# "hammel_error" or "hammel_warning".
error_kinds <- c(
  "argument", "table", "rates", "coverage", "groups", "population",
  "specification", "database", "runtime", "internal"
)

warning_kinds <- c(
  "no_effect", "estimate", "simulation", "performance", "deprecated"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Raise a classified error. `message` is a cli character vector: a header
# followed by optional "x" (problem) and "i" (information) bullets, with
# inline markup interpolated in the caller's environment. `kind` is one of
# error_kinds. Named arguments in `...` become fields of the condition.
# `call` is the frame of the user-facing function to name in the error;
# helpers that raise on behalf of another function forward their own `call`.
abort_hammel <- function(
  message,
  kind,
  ...,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_abort(
    message,
    class = condition_classes(kind, "error"),
    ...,
    call = call,
    .envir = .envir
  )
}

# Raise a classified warning; arguments as for abort_hammel(). Warnings do
# not name a calling function.
warn_hammel <- function(message, kind, ..., .envir = parent.frame()) {
  cli::cli_warn(
    message,
    class = condition_classes(kind, "warning"),
    ...,
    .envir = .envir
  )
}

# Raise an internal consistency error: a failure that valid inputs cannot
# cause, so it indicates a bug in Hammel.
abort_internal <- function(
  message,
  ...,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  url <- utils::packageDescription("hammel")$BugReports
  abort_hammel(
    c(
      paste("Internal error:", message),
      "i" = paste0(
        "This is likely a bug in Hammel. Please report it at {.url ",
        url, "}."
      )
    ),
    "internal",
    ...,
    call = call,
    .envir = .envir
  )
}

# Evaluate `expr`. A Hammel error it raises gains an "i" bullet saying where
# the problem was found, such as a calendar period or a group stratum, and is
# re-raised with its class, fields and call unchanged.
with_error_context <- function(expr, context) {
  tryCatch(
    expr,
    hammel_error = function(error) {
      error$message <- c(error$message, i = context)
      stop(error)
    }
  )
}

# "row 3" or "rows 3, 7, 12, ..., 31, and 40" for a message; long lists show
# the first and last rows.
format_rows <- function(rows) {
  noun <- if (length(rows) == 1L) "row" else "rows"
  paste(noun, cli::ansi_collapse(rows, trunc = 5, style = "both-ends"))
}

# A half-open band, as in "[15, 20)".
format_band <- function(start, end) {
  paste0("[", start, ", ", end, ")")
}

# The subject of a message about an argument (`end_time`) or a column given
# as `table$column`: "`end_time`" or "Column age_start of `rates`", styled.
describe_column <- function(name) {
  parts <- strsplit(name, "$", fixed = TRUE)[[1]]
  if (length(parts) == 2L) {
    cli::format_inline("Column {.field {parts[2]}} of {.arg {parts[1]}}")
  } else {
    cli::format_inline("{.arg {name}}")
  }
}

# Describe a value's type and size for an "x" bullet, as in "It is a
# character vector of length 2." Single atomic values are shown as code.
describe_value <- function(x) {
  if (is.null(x)) {
    return("`NULL`")
  }
  if (is.function(x)) {
    return("a function")
  }
  if (is.data.frame(x)) {
    return("a data frame")
  }
  type <- if (is.factor(x)) {
    "factor"
  } else if (is.atomic(x)) {
    typeof(x)
  } else if (is.list(x)) {
    "list"
  } else {
    class(x)[1]
  }
  type <- sub("^(double|integer)$", "numeric", type)
  if (!is.null(dim(x))) {
    return(paste("a", type, "matrix"))
  }
  if (is.atomic(x) && !is.factor(x) && length(x) == 1L) {
    value <- if (is.na(x)) "NA" else deparse(unname(x))
    return(paste0("`", value, "`"))
  }
  if (!length(x)) {
    return(paste0("an empty ", type, if (type == "list") "" else " vector"))
  }
  article <- if (grepl("^[aeiou]", type)) "an" else "a"
  noun <- if (type == "list") "" else " vector"
  paste0(article, " ", type, noun, " of length ", length(x))
}

# The class vector of a condition of `type` ("error" or "warning"). An
# unknown kind is a programming error in Hammel.
condition_classes <- function(kind, type) {
  kinds <- if (type == "error") error_kinds else warning_kinds
  if (!is.character(kind) || length(kind) != 1L || !kind %in% kinds) {
    stop(
      sprintf("Internal error: unknown %s kind `%s`.", type, format(kind)),
      call. = FALSE
    )
  }
  base <- paste0("hammel_", type)
  c(paste0(base, "_", kind), base)
}

# Capitalize the first letter, as at the start of a message.
upper_first <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
