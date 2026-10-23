test_that("minimum male age excludes boys while preserving births and unknown fathers", {
  model <- make_model(2, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time[2] <- -14
  result <- simulate(model, 0.5, seed = 52)
  expect_gt(result$simulation$births, 0)
  expect_equal(nrow(result$unions), 0)
  expect_true(all(is.na(result$persons$father_id)))
  expect_identical(result, simulate(model, 0.5, seed = 52))
  model$processes$union_formation$min_male_age <- NA_real_
  expect_error(simulate(model, 1), "min_male_age")
})

test_that("exact and fractional age thresholds govern searches without rounding", {
  # Set the first birth to a known time using the engine's exponential draws.
  set.seed(52)
  waits <- rexp(2)
  model <- make_model(3, birth_rate = waits[2] / 0.25, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time <- c(-20, -16.25, -16.249)
  model$processes$union_formation$min_male_age <- 16.5
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$unions$start_time, 0.25)
  expect_equal(result$unions$husband_id, 2)
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 1] == 2))
})

test_that("boys age into eligibility during a run", {
  model <- make_model(2, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time[2] <- -14.5
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 1)
  expect_gte(result$unions$start_time, 0.5)
  children <- result$persons[result$persons$mother_id %in% 1, ]
  expect_true(any(children$birth_time < 0.5))
  expect_true(all(is.na(children$father_id[children$birth_time < 0.5])))
  expect_true(all(children$father_id[children$birth_time >= result$unions$start_time] == 2))
})

test_that("sampling precedes age screening and does not refill rejected samples", {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time <- c(-20, -1, -25)
  model$processes$union_formation$matching$sample_size <- 1
  results <- lapply(1:50, function(seed) simulate(model, 0.1, seed = seed))
  first_fathers <- vapply(results, function(result) {
    children <- result$persons[result$persons$mother_id %in% 1, ]
    if (!nrow(children)) return(-1)
    children$father_id[1]
  }, numeric(1))
  # Missing fathers occur despite an eligible adult outside the drawn sample.
  expect_true(anyNA(first_fathers))
  expect_true(any(first_fathers == 3, na.rm = TRUE))
  expect_true(all(is.na(first_fathers) | first_fathers %in% c(-1, 3)))
})

test_that("existing partners are retained and women have no additional age restriction", {
  model <- make_model(2, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time <- c(-14, -13)
  model <- set_fertility(model, data.frame(age_start = 14, age_end = 15, rate = 20),
                          birth_interval = 0)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -0.1, end_time = NA_real_, end_reason = NA_character_)
  result <- simulate(model, 0.1, seed = 52)
  expect_gt(result$simulation$births, 0)
  expect_identical(result$unions, model$population$unions)
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 1] == 2))
  # Without a supplied union, the same young woman can search for an adult man.
  model$population$unions <- create_baseline_population(0)$unions
  model$population$persons$birth_time[2] <- -20
  result <- simulate(model, 0.1, seed = 52)
  expect_equal(result$unions$wife_id, 1)
  expect_equal(result$unions$husband_id, 2)
})
