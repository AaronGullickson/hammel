make_rate_model <- function(n = 2, search_rate = 1, birth_rate = 0, birth_interval = 0) {
  model <- make_model(n, birth_rate = birth_rate, birth_interval = birth_interval)
  if (n) model$population$persons$sex[1] <- "f"
  set_union_formation(model, "rate",
    create_matching_model(~ age_difference, c(age_difference = 0)),
    rates = data.frame(age_start = 15, age_end = Inf, rate = search_rate),
    sample_size = 0)
}

test_that("rate searches form unions independently of births and never originate from men", {
  model <- make_rate_model(search_rate = 20)
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$wife_id, 1)
  expect_equal(result$unions$husband_id, 2)
  expect_equal(result$simulation$births, 0)
  expect_identical(result, simulate(model, 1, seed = 52))
  model$population$persons$sex[] <- "m"
  expect_equal(nrow(simulate(model, 1, seed = 52)$unions), 0)
})

test_that("unpartnered women cannot give birth in rate mode", {
  for (n in c(1, 2)) {
    model <- make_rate_model(n, search_rate = 0, birth_rate = 20)
    result <- simulate(model, 1, seed = 52)
    expect_equal(result$simulation$births, 0)
    expect_equal(nrow(result$unions), 0)
  }
  model <- make_rate_model(search_rate = 20, birth_rate = 20, birth_interval = 0.02)
  result <- simulate(model, 1, seed = 52)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 1)
  expect_true(all(children$birth_time > result$unions$start_time[1]))
  expect_true(all(children$father_id == 2))
  expect_true(all(diff(children$birth_time) >= 0.02))
  continued <- simulate(set_population(model, result), 0.2, seed = 53)
  expect_identical(continued$unions, result$unions)
})

test_that("search waiting times agree with constant hazards and competing mortality", {
  n <- 3000
  model <- make_rate_model(2 * n, search_rate = 0.4)
  model$population$persons$sex[1:n] <- "f"
  model <- set_mortality(model, data.frame(sex = c("f", "m"),
    age_start = 0, age_end = Inf, rate = c(0.2, 0)))
  result <- simulate(model, 2, seed = 502)
  # First search precedes death and the endpoint with this probability.
  p <- 0.4 / 0.6 * (1 - exp(-0.6 * 2))
  expect_lt(abs(nrow(result$unions) - n * p), 6 * sqrt(n * p * (1 - p)))
  expect_equal(anyDuplicated(result$unions$wife_id), 0L)
  # Husbands can legitimately remarry after a wife's death, never concurrently.
  episodes <- split(result$unions, result$unions$husband_id)
  no_overlap <- vapply(episodes, function(x) {
    x <- x[order(x$start_time), ]
    nrow(x) == 1L || all(head(x$end_time, -1) <= tail(x$start_time, -1))
  }, logical(1))
  expect_true(all(no_overlap))
})

test_that("failed searches retry and age-band gaps carry zero search hazard", {
  model <- make_rate_model(search_rate = 100)
  # Initially the only man is too young. Repeated attempts can match after 0.5.
  model$population$persons$birth_time[2] <- -14.5
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 1)
  expect_gte(result$unions$start_time, 0.5)
  model$processes$union_formation$rates <- data.frame(
    age_start = c(0, 20.5), age_end = c(20, 20.75), rate = c(0, 1e6))
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 1)
  expect_gte(result$unions$start_time, 0.5)
  expect_lt(result$unions$start_time, 0.75)
  model$population$persons$birth_time[1] <- -21
  expect_equal(nrow(simulate(model, 1, seed = 52)$unions), 0)
})

test_that("husband death cancels pending births and widowhood permits another search", {
  model <- make_rate_model(3, search_rate = 20, birth_rate = 20, birth_interval = 0.02)
  model$population$persons$birth_time[2] <- -30
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  model <- set_mortality(model, data.frame(sex = c("f", "m", "m"),
    age_start = c(0, 0, 30.3), age_end = c(Inf, 30.3, Inf), rate = c(0, 0, 1e8)))
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 2)
  expect_equal(result$unions$husband_id, c(2, 3))
  death <- result$persons$death_time[2]
  expect_equal(result$unions$end_time[1], death)
  expect_gt(result$unions$start_time[2], death)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_true(all(children$birth_time < death |
                  children$birth_time > result$unions$start_time[2]))
  expect_true(all(children$father_id[children$birth_time < death] == 2))
  expect_true(all(children$father_id[children$birth_time > death] == 3))
  # No replacement partner: all post-widowhood births must be suppressed,
  # including any candidate already queued while the husband was alive.
  model <- set_population(model, hammel_population(model$population$persons[1:2, ],
    model$population$unions, end_time = 0))
  result <- simulate(model, 1, seed = 52)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_true(all(children$birth_time < result$persons$death_time[2]))
})
