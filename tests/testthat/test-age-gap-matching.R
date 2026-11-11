test_that("age-gap preferences give symmetric half-peak selection weights", {
  for (preferred in c(-3, 0, 3)) {
    matching <- age_gap_matching_model(preferred, 5)
    expect_s3_class(matching, "hammel_matching_model")
    probabilities <- calculate_matching_probabilities(matching,
      preferred + c(-10, -5, 0, 5, 10))
    expect_equal(probabilities / probabilities[3], c(1/16, 1/2, 1, 1/2, 1/16))
    expect_equal(prepare_matching_model(matching)$quadratic,
      c(2 * preferred * log(2) / 25, -log(2) / 25))
  }
})

test_that("age-gap helper rejects invalid or unrepresentable preferences", {
  for (bad in list(NA_real_, Inf, NaN, numeric(), c(1, 2), "3", TRUE, matrix(3))) {
    expect_error(age_gap_matching_model(bad, 5), "preferred_age_gap")
    expect_error(age_gap_matching_model(3, bad), "half_weight_distance")
  }
  for (bad in c(0, -1))
    expect_error(age_gap_matching_model(3, bad), "positive")
  for (bad in c(1e-300, 1e300))
    expect_error(age_gap_matching_model(3, bad), "finite coefficients")
  expect_error(age_gap_matching_model(1e308, 0.1), "finite coefficients")
})
test_that("union constructors default to zero-gap five-year half-weight matching", {
  expected <- age_gap_matching_model(0, 5)
  rates <- data.frame(age_start = 18, age_end = Inf, rate = 0.1)
  for (mode in c("on_demand", "rate")) {
    schedule <- if (mode == "rate") rates else NULL
    process <- union_formation(mode, rates = schedule)
    expect_equal(process$matching$coefficients, expected$coefficients)
    expect_equal(calculate_matching_probabilities(process$matching, c(-5, 0, 5)),
      c(0.25, 0.5, 0.25))
    model <- set_union_formation(hammel_model(), mode, rates = schedule)
    expect_equal(model$processes$union_formation$matching$coefficients,
      expected$coefficients)
    custom <- age_gap_matching_model(3, 7)
    expect_identical(union_formation(mode, custom, rates = schedule)$matching, custom)
    expect_identical(set_union_formation(hammel_model(), mode, custom,
      rates = schedule)$processes$union_formation$matching, custom)
  }
})

