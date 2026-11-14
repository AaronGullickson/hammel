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
