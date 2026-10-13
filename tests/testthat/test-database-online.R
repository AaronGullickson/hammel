test_that("live HMD and HFD downloads satisfy the rate and provenance contracts", {
  skip_if_not(identical(Sys.getenv("HAMMEL_RUN_ONLINE_TESTS"), "true"), "Set HAMMEL_RUN_ONLINE_TESTS=true to enable live database tests")
  if (!requireNamespace("HMDHFDplus", quietly = TRUE)) {
    stop("The online test requires HMDHFDplus.", call. = FALSE)
  }
  credentials <- c("HMD_USERNAME", "HMD_PASSWORD", "HFD_USERNAME", "HFD_PASSWORD")
  if (any(!nzchar(Sys.getenv(credentials)))) {
    stop("The online test requires HMD and HFD credential environment variables.", call. = FALSE)
  }
  model <- set_fertility(set_mortality(hammel_model(), get_mortality_hmd("USA", year = 2000)),
    get_fertility_hfd("USA", year = 2000))
  deaths <- model$processes$mortality
  births <- model$processes$fertility
  expect_true(inherits(deaths, "hammel_mortality"))
  expect_true(inherits(births, "hammel_fertility"))
  expect_identical(deaths$rate_unit, "person_year")
  expect_identical(births$rate_unit, "person_year")
  expect_true(setequal(deaths$rates$sex, c("m", "f")))
  expect_true(any(is.infinite(deaths$rates$age_end[deaths$rates$sex == "m"])))
  expect_true(any(is.infinite(deaths$rates$age_end[deaths$rates$sex == "f"])))
  expect_true(all(is.finite(deaths$rates$rate)))
  expect_true(all(deaths$rates$rate >= 0))
  expect_true(all(is.finite(births$rates$rate)))
  expect_true(all(births$rates$rate >= 0))
  expect_true(any(deaths$rates$rate > 0))
  expect_true(any(births$rates$rate > 0))
  expect_true(all(births$rates$age_end - births$rates$age_start == 1))
  expect_equal(deaths$source$country, "USA", tolerance = 0)
  expect_equal(births$source$country, "USA", tolerance = 0)
  expect_equal(deaths$source$year, 2000, tolerance = 0)
  expect_equal(births$source$year, 2000, tolerance = 0)
  expect_equal(deaths$source$item, "Mx_1x1", tolerance = 0)
  expect_equal(births$source$item, "asfrRR", tolerance = 0)
  expect_equal(births$source$open_ages, "omit", tolerance = 0)
  expect_true(!any(c("username", "password") %in% names(deaths$source)))
  expect_true(!any(c("username", "password") %in% names(births$source)))
})

