test_that("diagnostics validate and pass through the setter", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0))
  expect_false(union_formation("on_demand", matching)$diagnostics)
  expect_true(set_union_formation(hammel_model(), "on_demand", matching,
    diagnostics = TRUE)$processes$union_formation$diagnostics)
  for (bad in list(NULL, NA, 1, c(TRUE, FALSE), "yes")) {
    expect_error(union_formation("on_demand", matching, diagnostics = bad), "diagnostics")
    model <- make_model(0)
    model$processes$union_formation$diagnostics <- bad
    expect_error(simulate(model, 0), "diagnostics")
  }
})

test_that("diagnostics count screened candidates and rank utilities without changing draws", {
  for (mode in c("on_demand", "rate")) {
    model <- if (mode == "rate") make_rate_model(5, search_rate = 100) else
      make_model(5, birth_rate = 100, birth_interval = 0)
    model$population$persons$sex[1] <- "f"
    model$population$persons$birth_time <- c(-25, -20, -25, -30, -5)
    model$processes$union_formation$sample_size <- 0
    for (tied in c(FALSE, TRUE)) {
      model$processes$union_formation$matching <- create_matching_model(
        ~ age_difference, c(age_difference = if (tied) 0 else 0.1))
      plain <- simulate(model, 0.1, seed = 92)
      model$processes$union_formation$diagnostics <- TRUE
      observed <- simulate(model, 0.1, seed = 92)
      expect_equal(nrow(observed$unions), 1)
      expect_identical(observed$persons, plain$persons)
      expect_identical(observed$unions[names(plain$unions)], plain$unions)
      expect_identical(observed$unions$.n_eligible, 3L)
      husband <- observed$unions$husband_id
      scores <- if (tied) c(0, 0, 0) else c(-5, 0, 5)
      expected <- 1L + sum(scores > scores[match(husband, 2:4)])
      expect_equal(observed$unions$.utility_rank, expected)
      restarted <- set_population(model, observed)
      expect_identical(simulate(restarted, 0.1, seed = 93)$unions, observed$unions)
      restarted$processes$union_formation$diagnostics <- FALSE
      expect_identical(simulate(restarted, 0.1, seed = 93)$unions, observed$unions)
      model$processes$union_formation$diagnostics <- FALSE
    }
  }
})

test_that("diagnostic columns exist for empty runs and historical unknowns", {
  model <- make_rate_model(2, search_rate = 100)
  old <- simulate(model, 0.1, seed = 92)
  model$processes$union_formation$diagnostics <- TRUE
  empty <- simulate(model, 0, seed = 92)
  expect_identical(empty$unions$.n_eligible, integer())
  expect_identical(empty$unions$.utility_rank, integer())
  single <- simulate(model, 0.1, seed = 92)
  expect_identical(single$unions$.n_eligible, 1L)
  expect_identical(single$unions$.utility_rank, 1L)
  model <- set_population(model, old)
  for (duration in c(0, 0.1)) {
    historical <- simulate(model, duration, seed = 93)
    expect_identical(historical$unions$.n_eligible, NA_integer_)
    expect_identical(historical$unions$.utility_rank, NA_integer_)
  }
})
