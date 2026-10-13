test_that("annual midpoints weight the current and preceding input bands", {
  rates <- data.frame(age_start = c(24, 25, 26), age_end = c(25, 26, 27),
                      rate = c(0.1, 0.2, 0.3))
  adjusted <- hammel:::adjusted_fertility_cpp(rates, 0.75)
  expect_identical(adjusted$age_start, rates$age_start)
  expect_identical(adjusted$age_end, rates$age_end)
  # At 25.5: 0.25 years at 0.1 plus 0.5 years at 0.2 -> eligibility 0.875.
  expect_equal(adjusted$rate, rates$rate / c(0.95, 0.875, 0.8), tolerance = 1e-14)
  expect_equal(hammel:::adjusted_fertility_cpp(rates[3:1, ], 0.75), adjusted)
})

test_that("five-year band midpoints can have constant-rate renewal compensation", {
  rates <- data.frame(age_start = c(20, 25), age_end = c(25, 30),
                      rate = c(0.1, 0.2))
  adjusted <- hammel:::adjusted_fertility_cpp(rates, 0.75)
  expect_equal(adjusted$rate, rates$rate / (1 - 0.75 * rates$rate), tolerance = 1e-14)
  expect_equal(1 / (0.75 + 1 / adjusted$rate), rates$rate, tolerance = 1e-14)
})

test_that("long spacing windows include multiple original bands and gaps", {
  rates <- data.frame(age_start = c(20, 21, 23), age_end = c(21, 22, 24),
                      rate = c(0.1, 0.2, 0.3))
  adjusted <- hammel:::adjusted_fertility_cpp(rates, 2.5)
  # At 23.5, [21,23.5) contains 1 year at 0.2, a gap, and 0.5 at 0.3.
  expect_equal(adjusted$rate, rates$rate / c(0.95, 0.8, 0.65), tolerance = 1e-14)
  # The calculation must never substitute already-adjusted earlier rates.
  expect_equal(adjusted$rate[3], 0.3 / (1 - 0.2 - 0.5 * 0.3), tolerance = 1e-14)
})

test_that("adjusted fertility uses the shared constant-hazard sampler", {
  rates <- data.frame(age_start = c(24, 25), age_end = c(25, 26),
                      rate = c(0.1, 0.2))
  adjusted <- hammel:::adjusted_fertility_cpp(rates, 0.75)
  threshold <- 0.5 * adjusted$rate[1] + 0.1
  expect_equal(hammel:::hazard_wait_cpp(adjusted, 24.5, threshold),
               0.5 + 0.1 / adjusted$rate[2], tolerance = 1e-14)
  expect_true(is.infinite(hammel:::hazard_wait_cpp(adjusted, 26, 0.1)))
})

test_that("zero spacing and zero rates retain original bands and risks", {
  rates <- data.frame(age_start = c(15.13, 17.21),
                      age_end = c(16.07, 18.01), rate = c(0.2, 0.4))
  expect_equal(hammel:::adjusted_fertility_cpp(rates, 0), rates)
  rates$rate <- 0
  adjusted <- hammel:::adjusted_fertility_cpp(rates, 0.37)
  expect_equal(adjusted, rates)
  expect_true(all(is.infinite(hammel:::hazard_wait_cpp(adjusted,
    c(0, 15.5, 16.5, 17.8), rep(0.2, 4)))))
})

test_that("nonpositive midpoint eligibility is rejected", {
  for (rate in c(1, 1.1)) {
    rates <- data.frame(age_start = 15, age_end = 50, rate = rate)
    expect_error(hammel:::adjusted_fertility_cpp(rates, 1),
                 "nonpositive midpoint eligibility")
    expect_error(simulate(make_model(sex = "f", birth_rate = rate,
                                    birth_interval = 1), 1),
                 "nonpositive midpoint eligibility")
  }
  rates <- data.frame(age_start = 15, age_end = 50, rate = 2)
  expect_equal(hammel:::adjusted_fertility_cpp(rates, 0), rates)
})

test_that("open-ended fertility bands are rejected at construction and revalidation", {
  model <- make_model()
  for (rate in c(0, 0.2)) {
    for (b in c(0, 0.75)) {
      rates <- data.frame(age_start = 15, age_end = Inf, rate = rate)
      expect_error(fertility(rates, birth_interval = b), "finite ends")
      expect_error(set_fertility(model, rates, birth_interval = b), "finite ends")
      expect_error(hammel:::adjusted_fertility_cpp(rates, b), "finite ends")
      edited <- model$processes$fertility
      edited$rates$age_end <- Inf
      expect_error(set_process(model, edited), "finite ends")
      invalid <- model
      invalid$processes$fertility <- edited
      expect_error(simulate(invalid, 0), "finite ends")
    }
  }
  expect_true(is.infinite(mortality(model$processes$mortality$rates)$rates$age_end[1]))
})

test_that("independent mortality preserves constant fertility after burn-in", {
  n <- 10000
  model <- make_model(n, sex = "f", death_rate = 0.1, birth_rate = 0.2,
                      birth_interval = 0.37)
  result <- simulate(model, 20, seed = 456)
  deaths <- result$persons$death_time[seq_len(n)]
  deaths[is.na(deaths)] <- 20
  exposure <- sum(pmax(0, pmin(deaths, 20) - 5))
  births <- sum(result$persons$mother_id <= n &
                result$persons$birth_time >= 5, na.rm = TRUE)
  expect_lt(abs(births / exposure - 0.2), 0.01)
})
