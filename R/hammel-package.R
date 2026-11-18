# Package help page, options, native-library registration and imports.

# PACKAGE DOCUMENTATION ----------------------------------------------------

#' @section Options:
#' hammel reads the following option, which can be set with [options()],
#' or in an `.Rprofile`:
#'
#' * `hammel.progress`: the default for the `progress` argument of
#'   [run_simulation()]. When unset, a progress bar is shown in
#'   interactive sessions and hidden otherwise. An explicit `progress` 
#'   argument takes precedence.
"_PACKAGE"

#' @useDynLib hammel, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head tail
NULL
