exogamy_test_model <- function(group_exogamy = NULL) {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m", "m")
  model$population$persons$group <- c("A", "A", "B")
  model$population$groups <- list(group = group())
  identity(set_union_formation(model, "on_demand", diagnostics = TRUE,
    matching = age_gap_matching_model(0, 5, group_exogamy = group_exogamy)))
}

test_that("NA prohibitions exclude candidates before scoring and diagnostics", {
  # A whole-dimension NA, or the pair, prohibits the B man; each overrides any
  # preference (here a strong one for the other pair, A:A).
  for (prohibit in list(NA, c("A:B" = NA), list(group = NA), list(group = c("A:B" = NA, "A:A" = 1000)))) {
    model <- exogamy_test_model(prohibit)
    whole <- identical(prohibit, NA) || identical(prohibit, list(group = NA))
    if (whole) {
      expect_no_warning(result <- simulate(model, .3, seed = 52))
    } else {
      result <- suppressWarnings(simulate(model, .3, seed = 52))
    }
    expect_equal(result$unions$husband_id, 2)
  }
  matching <- as_internal_matching(matching_model(~ age_difference, c(age_difference = 0),
    group_exogamy = NA), c("A", "B"))
  expect_equal(calculate_matching_probabilities(matching, c(0, 0), "A", c("A", "B")), c(1, 0))
  expect_equal(calculate_matching_probabilities(matching, 0, "A", "B"), 0)
})

test_that("prohibitions are symmetric and validated with pair effects", {
  # "B:A" prohibits A women with B men too.
  model <- exogamy_test_model(c("B:A" = NA))
  result <- suppressWarnings(simulate(model, .3, seed = 52))
  expect_equal(result$unions$husband_id, 2)
  internal <- as_internal_matching(age_gap_matching_model(0, 5, group_exogamy = c("B:A" = NA)),
                                   c("A", "B"))
  expect_false(matching_group_allowed(internal, "A", "B"))
  expect_false(matching_group_allowed(internal, "B", "A"))
  expect_true(matching_group_allowed(internal, "A", "A"))
  for (bad in list(c("A:A" = NA), c("A:B" = NA, "B:A" = -1), c("A:B" = NA, "A:B" = NA),
                   c("A::B" = NA), c("A:B" = -Inf), "A:B", c(TRUE, FALSE), list(NA), c(NA, NA)))
    expect_error(age_gap_matching_model(0, 5, group_exogamy = bad))
  expect_error(age_gap_matching_model(0, 5, group_exogamy = c("A:A" = NA)), "within a category")
  expect_error(age_gap_matching_model(0, 5, group_exogamy = c("A:B" = NA, "B:A" = -1)),
               "Conflicting")
  # The argument is gone.
  expect_error(age_gap_matching_model(0, 5, prevent_exogamy = TRUE), "unused argument")
  model$processes$union_formation$matching$group_exogamy <- "x"
  expect_error(simulate(model, 0), "group_exogamy")
})

test_that("empty filtered samples use existing birth cancellation and unknown-father rules", {
  model <- exogamy_test_model(NA)
  model$population$persons$group <- c("A", "B", "B")
  model$population$groups <- list(group = group())
  model$processes$union_formation$cancel_birth_if_no_match <- TRUE
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$simulation$births, 0)
  expect_equal(nrow(result$unions), 0)
  model$processes$union_formation$cancel_birth_if_no_match <- FALSE
  result <- simulate(model, 0.3, seed = 52)
  expect_gt(result$simulation$births, 0)
  expect_equal(nrow(result$unions), 0)
  children <- result$persons$birth_time >= 0
  expect_true(all(is.na(result$persons$father_id[children])))
  expect_true(all(result$ancestry$unknown[children] == 0.5))
  # A finite penalty alone still selects a partner from a cross-group-only pool.
  model$processes$union_formation$matching <- age_gap_matching_model(0, 5,
    group_exogamy = c("A:B" = -10))
  expect_equal(nrow(simulate(model, 0.3, seed = 52)$unions), 1)
})

test_that("a compatible man anywhere in the pool is always found under restrictions", {
  model <- exogamy_test_model(NA)
  model$processes$union_formation$cancel_birth_if_no_match <- TRUE
  husbands <- vapply(1:30, function(seed) {
    unions <- simulate(model, 0.1, seed = seed)$unions
    if (nrow(unions)) unions$husband_id[1] else NA_real_
  }, numeric(1))
  births <- vapply(1:30, function(seed) simulate(model, 0.1, seed = seed)$simulation$births, numeric(1))
  expect_true(all(husbands[births > 0] == 2))
})

test_that("rate searches retry and father overrides carry independent restrictions", {
  model <- exogamy_test_model(NA)
  model <- identity(set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0))
  matching <- model$processes$union_formation$matching
  model <- identity(set_union_formation(model, "rate", matching = matching,
    rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
    cancel_birth_if_no_match = FALSE))
  result <- simulate(model, 0.3, seed = 52)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 2))
  model$processes$union_formation$father_matching <- age_gap_matching_model(0, 5,
    group_exogamy = c("A:B" = 1000))
  expect_warning(result <- simulate(model, .3, seed = 52), "could not meet male-group targets")
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model$processes$union_formation$rates$rate <- 100
  # The union may form before any unpartnered birth needs a father, so the
  # father-calibration warning is incidental here.
  result <- suppressWarnings(simulate(model, .3, seed = 52))
  expect_equal(result$unions$husband_id, 2)
  # No compatible men: repeated searches neither form unions nor reserve men.
  model$population$persons$group <- c("A", "B", "B")
  model$population$groups <- list(group = group())
  model$processes$union_formation$father_matching <- NULL
  model$processes$union_formation$cancel_birth_if_no_match <- TRUE
  expect_equal(nrow(simulate(model, 0.3, seed = 52)$unions), 0)
})

test_that("existing partners remain and disabled restrictions preserve RNG and histories", {
  model <- exogamy_test_model(NA)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1, husband_id = 3,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model <- exogamy_test_model()
  set.seed(432)
  original <- simulate(model, 1)
  state <- .Random.seed
  model$processes$union_formation$matching$group_exogamy <- numeric()
  set.seed(432)
  expect_identical(simulate(model, 1), original)
  expect_identical(.Random.seed, state)
})

test_that("restriction tables and missing or empty group populations validate", {
  matching <- as_internal_matching(age_gap_matching_model(0, 5, group_exogamy = data.frame(
    group_x = c("A:B", "A:B"), group_y = c("C::D", "E"), log_odds = c(NA, -1))), c("A:B", "C::D", "E"))
  expect_false(matching_group_allowed(matching, "A:B", "C::D"))
  expect_false(matching_group_allowed(matching, "C::D", "A:B"))
  expect_true(matching_group_allowed(matching, "A:B", "E"))
  model <- make_model(0)
  model$processes$union_formation$matching <- age_gap_matching_model(0, 5, group_exogamy = NA)
  expect_error(simulate(model, 0), "declared group dimension")
  model$population$persons$group <- character()
  model$population$groups <- list(group = group())
  expect_no_error(simulate(model, 1))
})
