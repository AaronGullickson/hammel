# Keep test output independent of whether the session is interactive. Tests
# of the hammel.progress option set it explicitly.
withr::local_options(
  hammel.progress = FALSE,
  .local_envir = testthat::teardown_env()
)
