test_that("hazard inversion crosses rate bands and handles zero hazards", {
  rates <- data.frame(age_start = c(0, 1), age_end = c(1, Inf), rate = c(0.12, 0.24))
  wait <- hammel:::hazard_wait_cpp(rates, 10/12, 0.05)
  expect_equal(as.numeric(wait), 3.5/12)
  constant <- data.frame(age_start = 0, age_end = Inf, rate = 0.2)
  expect_equal(as.numeric(hammel:::hazard_wait_cpp(constant, c(0, 10), c(0.2, 0.4))), c(1,
    2))
  zero <- rates
  zero$rate <- 0
  expect_true(is.infinite(hammel:::hazard_wait_cpp(zero, 0, 0.5)))
  zero$rate[2] <- 1
  expect_equal(hammel:::hazard_wait_cpp(zero, 0, 0.5), 1.5, tolerance = 0)
})

test_that("sampled survival matches integrated piecewise hazards", {
  rates <- data.frame(age_start = c(0, 1), age_end = c(1, Inf), rate = c(0.12, 0.24))
  set.seed(42)
  n <- 1e+05
  waits <- hammel:::hazard_wait_cpp(rates, rep(0, n), rexp(n))
  expected <- exp(-0.12 - 0.24)
  expect_true(abs(mean(waits > 2) - expected) < 6 * sqrt(expected * (1 - expected)/n))
})

