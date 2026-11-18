# Package help page, options, native-library registration and imports.

# PACKAGE DOCUMENTATION ----------------------------------------------------

#' @section Options:
#' hammel reads the following option, which can be set with [options()],
#' for example in a session or an `.Rprofile`:
#'
#' * `hammel.progress`: the default for the `progress` argument of
#'   [run_simulation()]. When unset, a progress bar is shown in
#'   interactive sessions ([interactive()]) and hidden otherwise, such as
#'   in scripts, rendered documents and tests. Set it to `TRUE` or `FALSE`
#'   to override this. An explicit `progress` argument takes precedence.
"_PACKAGE"

#' @useDynLib hammel, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head tail
NULL
