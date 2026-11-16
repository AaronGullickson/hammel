test_that("check_columns reports type, duplicates and missing columns", {
  set_rates <- function(rates) {
    check_columns(rates, c("age_start", "age_end", "rate"), "rates")
  }
  error <- expect_error(set_rates(list(rate = 1)), class = "hammel_error_table")
  expect_identical(error$argument, "rates")
  expect_identical(deparse(error$call), "set_rates(list(rate = 1))")
  expect_match(conditionMessage(error), "It is a list of length 1.")

  rates <- data.frame(age_start = 0, rate = 1, rate = 2, check.names = FALSE)
  error <- expect_error(set_rates(rates), class = "hammel_error_table")
  expect_identical(error$columns, "rate")
  expect_match(conditionMessage(error), "unique column names")

  rates <- tibble::tibble(age_start = 0)
  error <- expect_error(set_rates(rates), class = "hammel_error_table")
  expect_identical(error$columns, c("age_end", "rate"))
  expect_null(set_rates(tibble::tibble(age_start = 0, age_end = 1, rate = 1)))
})

test_that("check_numeric_column reports invalid rows", {
  x <- c(1.5, NA, 3, NA)
  error <- expect_error(
    check_numeric_column(x, "rates$age_start"),
    class = "hammel_error_table"
  )
  expect_identical(error$argument, "rates$age_start")
  expect_identical(error$rows, c(2L, 4L))
  expect_match(conditionMessage(error), "Column age_start of `rates`")
  expect_match(conditionMessage(error), "Missing at rows 2 and 4.")
  expect_null(check_numeric_column(x, "rates$age_start", missing = TRUE))

  error <- expect_error(
    check_numeric_column(c(1, NaN, Inf), "rates$rate", missing = TRUE),
    class = "hammel_error_table"
  )
  expect_identical(error$rows, c(2L, 3L))
  expect_match(conditionMessage(error), "Invalid values NaN and Inf")

  error <- expect_error(
    check_numeric_column(c(1, 2.5, -1, NA), "persons$person_id", id = TRUE),
    class = "hammel_error_table"
  )
  expect_identical(error$rows, 2:4)
  expect_match(conditionMessage(error), "positive whole-number IDs")

  error <- expect_error(
    check_numeric_column(letters, "age_difference"),
    class = "hammel_error_table"
  )
  expect_identical(error$rows, integer())
  expect_match(conditionMessage(error), "`age_difference`")
  expect_match(conditionMessage(error), "character vector of length 26")

  many <- c(NA, NA, 1, NA, NA, NA, NA, NA)
  error <- expect_error(check_numeric_column(many, "rates$rate"))
  expect_length(error$rows, 7L)
})

test_that("check_time requires one finite time", {
  for (time in list(NULL, NA_real_, Inf, c(0, 1), "0", 0 + 0i, matrix(0))) {
    error <- expect_error(
      check_time(time, "end_time"),
      class = "hammel_error_argument"
    )
    expect_identical(error$argument, "end_time")
  }
  expect_null(check_time(2000.5, "end_time"))
})

test_that("checks format their messages consistently", {
  set_rates <- function(rates) {
    check_columns(rates, c("age_start", "age_end", "rate"), "rates")
    check_numeric_column(rates$rate, "rates$rate")
  }
  expect_snapshot(
    set_rates(tibble::tibble(age_start = 0, rate = 1)),
    error = TRUE
  )
  expect_snapshot(
    set_rates(tibble::tibble(age_start = 0, age_end = 1, rate = c(1, NA))),
    error = TRUE
  )
  expect_snapshot(check_time("2000", "end_time"), error = TRUE)
})
