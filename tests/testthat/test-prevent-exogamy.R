exogamy_test_model <- function(prevent = FALSE, group_exogamy = NULL) {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m", "m")
  model$population$persons$group <- c("A", "A", "B")
  set_union_formation(model, "on_demand", sample_size = 0, diagnostics = TRUE,
    matching = create_age_gap_matching_model(0, 5, group_exogamy = group_exogamy,
                                             prevent_exogamy = prevent))
}

test_that("hard restrictions exclude candidates before scoring and diagnostics", {
  for (prevent in list(TRUE, "A:B", "A::B")) {
    model <- exogamy_test_model(prevent, c("A:B" = 1000))
    result <- simulate(model, 0.3, seed = 52)
    expect_equal(result$unions$husband_id, 2)
    expect_equal(result$unions$.n_eligible, 1L)
    expect_equal(result$unions$.utility_rank, 1L)
  }
  matching <- create_matching_model(~ age_difference, c(age_difference = 0),
    group_exogamy = c("A:B" = 1000), prevent_exogamy = TRUE)
  expect_equal(calculate_matching_probabilities(matching, c(0, 0), "A", c("A", "B")), c(1, 0))
  expect_equal(calculate_matching_probabilities(matching, 0, "A", "B"), 0)
})

test_that("directional and symmetric prohibitions use female-first identities", {
  model <- exogamy_test_model("B:A", c("A:B" = 1000))
  expect_equal(simulate(model, 0.3, seed = 52)$unions$husband_id, 3)
  model <- exogamy_test_model("B::A", c("A:B" = 1000))
  expect_equal(simulate(model, 0.3, seed = 52)$unions$husband_id, 2)
  for (bad in list(NA, c(TRUE, FALSE), 1, "A:::B", "A:A", c("A:B", "A::B"), c(NA_character_)))
    expect_error(create_age_gap_matching_model(0, 5, prevent_exogamy = bad))
  model$processes$union_formation$matching$prevent_exogamy <- NA
  expect_error(simulate(model, 0), "prevent_exogamy")
})

test_that("empty filtered samples use existing birth cancellation and unknown-father rules", {
  model <- exogamy_test_model(TRUE)
  model$population$persons$group <- c("A", "B", "B")
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
  model$processes$union_formation$matching <- create_age_gap_matching_model(0, 5,
    group_exogamy = c("A::B" = -10))
  expect_equal(nrow(simulate(model, 0.3, seed = 52)$unions), 1)
})

test_that("sampling remains broad with no refilling after group exclusions", {
  small <- exogamy_test_model(TRUE)
  small$processes$union_formation$sample_size <- 1L
  small$processes$union_formation$cancel_birth_if_no_match <- TRUE
  full <- small
  full$processes$union_formation$sample_size <- 0L
  # Some one-man encounter samples miss the compatible man although he exists.
  outcomes <- vapply(1:30, function(seed) {
    nrow(simulate(full, 0.1, seed = seed)$unions) > 0 &&
      nrow(simulate(small, 0.1, seed = seed)$unions) == 0
  }, logical(1))
  expect_true(any(outcomes))
})

test_that("rate searches retry and father overrides carry independent restrictions", {
  model <- exogamy_test_model(TRUE)
  model <- set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0)
  model <- set_union_formation(model, "rate", matching = model$processes$union_formation$matching,
    rates = data.frame(age_start = 15, age_end = Inf, rate = 0), sample_size = 0,
    cancel_birth_if_no_match = FALSE)
  result <- simulate(model, 0.3, seed = 52)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 2))
  model$processes$union_formation$father_matching <- create_age_gap_matching_model(0, 5,
    group_exogamy = c("A:B" = 1000))
  result <- simulate(model, 0.3, seed = 52)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model$processes$union_formation$rates$rate <- 100
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$unions$husband_id, 2)
  # No compatible men: repeated searches neither form unions nor reserve men.
  model$population$persons$group <- c("A", "B", "B")
  model$processes$union_formation$father_matching <- NULL
  model$processes$union_formation$cancel_birth_if_no_match <- TRUE
  expect_equal(nrow(simulate(model, 0.3, seed = 52)$unions), 0)
})

test_that("existing partners remain and disabled restrictions preserve RNG and histories", {
  model <- exogamy_test_model(TRUE)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1, husband_id = 3,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model <- exogamy_test_model(FALSE)
  set.seed(432)
  original <- simulate(model, 1)
  state <- .Random.seed
  model$processes$union_formation$matching$prevent_exogamy <- character()
  set.seed(432)
  expect_identical(simulate(model, 1), original)
  expect_identical(.Random.seed, state)
})

test_that("restriction tables and missing or empty group populations validate", {
  matching <- create_age_gap_matching_model(0, 5, prevent_exogamy = data.frame(
    female_group = "A:B", male_group = "C::D"))
  expect_false(matching_group_allowed(matching, "A:B", "C::D"))
  expect_true(matching_group_allowed(matching, "C::D", "A:B"))
  model <- make_model(0)
  model$processes$union_formation$matching <- create_age_gap_matching_model(0, 5, prevent_exogamy = TRUE)
  expect_error(simulate(model, 0), "group identities")
  model$population$persons$group <- character()
  expect_no_error(simulate(model, 1))
})
