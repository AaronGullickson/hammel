test_that("abort_hammel classes errors and names the calling function", {
  user_function <- function(columns) {
    abort_hammel(
      c(
        "{.arg rates} is missing required columns.",
        "x" = "Missing: {.field {columns}}."
      ),
      "table",
      columns = columns
    )
  }
  error <- expect_error(
    user_function(c("age_end", "rate")),
    class = "hammel_error_table"
  )
  expect_s3_class(error, "hammel_error")
  expect_identical(error$columns, c("age_end", "rate"))
  expect_identical(
    deparse(error$call),
    "user_function(c(\"age_end\", \"rate\"))"
  )
  expect_match(conditionMessage(error), "missing required columns")
  expect_match(conditionMessage(error), "age_end")
})

test_that("helpers raising for another function forward its call", {
  check_value <- function(x, call = rlang::caller_env()) {
    if (x < 0) {
      abort_hammel("{.arg x} must be nonnegative.", "argument", call = call)
    }
    invisible(x)
  }
  set_value <- function(x) {
    check_value(x)
  }
  error <- expect_error(set_value(-1), class = "hammel_error_argument")
  expect_identical(deparse(error$call), "set_value(-1)")
})

test_that("interpolated values print literally", {
  category <- "a{b}"
  error <- expect_error(
    abort_hammel("Unknown category {.val {category}}.", "groups"),
    class = "hammel_error_groups"
  )
  expect_match(conditionMessage(error), "a{b}", fixed = TRUE)
})

test_that("warn_hammel classes warnings and stores fields", {
  term <- "sex"
  warning <- expect_warning(
    warn_hammel(
      "The weight term {.field {term}} has no effect.",
      "no_effect",
      term = term
    ),
    class = "hammel_warning_no_effect"
  )
  expect_s3_class(warning, "hammel_warning")
  expect_identical(warning$term, "sex")
})

test_that("every condition kind has its classes", {
  for (kind in error_kinds) {
    expect_error(
      abort_hammel("Problem.", kind),
      class = paste0("hammel_error_", kind)
    )
  }
  for (kind in warning_kinds) {
    expect_warning(
      warn_hammel("Problem.", kind),
      class = paste0("hammel_warning_", kind)
    )
  }
})

test_that("unknown condition kinds are internal errors", {
  expect_error(abort_hammel("Problem.", "unknown"), "unknown error kind")
  expect_error(warn_hammel("Problem.", "table"), "unknown warning kind")
  expect_error(abort_hammel("Problem.", c("table", "rates")), "unknown")
})

test_that("messages keep their bullet layout", {
  user_function <- function() {
    abort_hammel(
      c(
        "{.arg rates} is missing required columns.",
        "x" = "Missing: {.field age_end} and {.field rate}.",
        "i" = "Mortality rates need {.field age_start} and {.field rate}."
      ),
      "table"
    )
  }
  expect_snapshot(user_function(), error = TRUE)
})
