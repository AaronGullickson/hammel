# Expect a Hammel error of `kind` whose call names the exported function
# `fn`, as in expect_hammel_error(set_mortality(m, bad), "table",
# "set_mortality"). Returns the condition for further checks.
expect_hammel_error <- function(object, kind, fn) {
  error <- expect_error(object, class = paste0("hammel_error_", kind))
  called <- if (is.null(error$call)) "" else rlang::call_name(error$call)
  expect_identical(called, fn)
  expect_true(fn %in% getNamespaceExports("hammel"))
  invisible(error)
}
