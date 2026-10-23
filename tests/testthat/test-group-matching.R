test_that("pair vectors expand directional and symmetric effects unambiguously", {
  pairs <- expand_group_exogamy(c("A:B" = -1, "B:A" = -0.5, "A::C" = -0.8, "A::A" = 0.2))
  expect_s3_class(pairs, "tbl_df")
  expect_equal(pairs$female_group, c("A", "B", "A", "C", "A"))
  expect_equal(pairs$male_group, c("B", "A", "C", "A", "A"))
  expect_equal(pairs$log_odds, c(-1, -0.5, -0.8, -0.8, 0.2))
  for (bad in list(c(1), c("A:B" = NA_real_), c("A:B" = Inf), c("A:::B" = 0),
                   c("A:" = 0), c(":B" = 0), c("A/B" = 0), c("A: :B" = 0),
                   c("A:B" = 0, "A::B" = 1), c("A::B" = 0, "B::A" = 0),
                   c("A:B" = 0, "A:B" = 1)))
    expect_error(expand_group_exogamy(bad))
  expect_equal(nrow(expand_group_exogamy(NULL)), 0)
  expect_equal(nrow(expand_group_exogamy(numeric())), 0)
  table <- data.frame(female_group = "A:B", male_group = "C::D", log_odds = 1)
  expect_identical(expand_group_exogamy(table)$female_group, "A:B")
  expect_error(expand_group_exogamy(rbind(table, table)), "overlapping")
})

test_that("native pair offsets match analytical probabilities and R scores", {
  model <- create_matching_model(~ age_difference, c(age_difference = 0),
    group_exogamy = c("A:B" = log(2), "A::C" = log(3)), matching_method = "female_choice")
  labels <- c("A", "B", "C")
  prepared <- prepare_matching_model(model, labels)
  result <- inspect_group_matching_cpp(c(0, 0, 0), prepared, 0L, 0:2)
  expect_equal(result$utilities, log(c(1, 2, 3)))
  expect_equal(result$probabilities, c(1, 2, 3)/6)
  expect_equal(result$utilities, calculate_matching_utilities(model, c(0, 0, 0), "A", labels))
  expect_equal(inspect_group_matching_cpp(c(0, 0), prepared, 1L, 0:1)$probabilities, c(0.5, 0.5))
  expect_equal(inspect_group_matching_cpp(c(0, 0), prepared, 2L, 0:1)$probabilities, c(3, 1)/4)
  only <- create_matching_model(~ 1, numeric(), group_exogamy = c("A:B" = log(2)), matching_method = "female_choice")
  expect_equal(calculate_matching_probabilities(only, c(0, 0), "A", c("A", "B")), c(1, 2)/3)
  helper <- create_age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1), matching_method = "female_choice")
  expect_identical(helper$group_exogamy, c("A:B" = -1))
})

test_that("same_group formula and additive effects combine for custom formulas", {
  matching <- create_matching_model(~ abs(age_difference) + same_group,
    c("abs(age_difference)" = -0.1, same_group = 0.3), group_exogamy = c("A:B" = 0.7), matching_method = "female_choice")
  gaps <- c(-2, 1, 0)
  expected <- c(-0.2 + 0.3, -0.1 + 0.7, 0)
  result <- inspect_group_matching_cpp(gaps, prepare_matching_model(matching, c("A", "B", "C")),
                                       0L, 0:2)
  expect_equal(result$utilities, expected)
  expect_equal(result$utilities, calculate_matching_utilities(matching, gaps, "A", c("A", "B", "C")))
  expect_equal(result$probabilities, exp(expected) / sum(exp(expected)))
})

test_that("categorical group predictors have stable treatment contrasts and interactions", {
  formula <- ~ age_difference + male_group + female_group:age_difference
  # Derive explicit coefficient names from a reference design, using fixed levels.
  data <- data.frame(age_difference = c(0, 1), male_group = factor(c("A", "B"), levels = c("A", "B")),
    female_group = factor(c("A", "B"), levels = c("A", "B")))
  names <- colnames(stats::model.matrix(formula, data))[-1]
  beta <- setNames(seq_along(names)/10, names)
  matching <- create_matching_model(formula, beta, group_levels = c("A", "B"), matching_method = "female_choice")
  gaps <- c(1, 3)
  supplied <- data.frame(age_difference = gaps,
    female_group = factor(c("B", "B"), levels = c("A", "B")),
    male_group = factor(c("A", "A"), levels = c("A", "B")))
  expected <- drop(stats::model.matrix(formula, supplied)[, -1, drop = FALSE] %*% beta)
  result <- inspect_group_matching_cpp(gaps, prepare_matching_model(matching, c("A", "B")), 1L, c(0L, 0L))
  expect_equal(result$utilities, unname(expected))
  # Session contrast settings must not change these specified coefficients.
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old))
  expect_equal(calculate_matching_utilities(matching, gaps, "B", c("A", "A")), unname(expected))
  expect_error(create_matching_model(~ male_group, c(male_groupB = 1)), "group_levels")
  expect_error(create_matching_model(~ male_group, c(male_groupB = 1), group_levels = "A"), "two")
  expect_error(prepare_matching_model(matching, c("A", "C")), "covered")
  expect_error(calculate_matching_utilities(matching, 1, "A", "C"), "undeclared")
})

group_choice_model <- function(initiation = "on_demand", matching, father_matching = NULL) {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m", "m")
  model$population$persons$group <- c("A", "A", "B")
  if (initiation == "rate") {
    model <- quiet_inheritance(set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
      age_start = 15, age_end = 50, rate = c(0, 20)), 0))
  }
  matching$sample_size <- 0
  matching$matching_method <- "female_choice"
  if (!is.null(father_matching)) {
    father_matching$sample_size <- 0
    father_matching$matching_method <- "female_choice"
  }
  quiet_inheritance(set_union_formation(model, initiation, matching,
    rates = if (initiation == "rate") data.frame(age_start = 15, age_end = Inf, rate = 0) else NULL,
    father_matching = father_matching, diagnostics = TRUE))
}

test_that("pair effects drive both union modes, father choice, and final utility ranks", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0), group_exogamy = c("A:B" = 1000), matching_method = "female_choice")
  result <- simulate(group_choice_model(matching = matching), 0.3, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  expect_equal(result$unions$.utility_rank, 1L)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model <- group_choice_model("rate", matching)
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(nrow(result$unions), 0)
  expect_true(any(result$persons$birth_time >= 0))
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 3))
  model$processes$union_formation$rates$rate <- 1e8
  result <- simulate(model, 0.3, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  opposite <- create_matching_model(~ age_difference, c(age_difference = 0), group_exogamy = c("A:A" = 1000))
  result <- simulate(group_choice_model("rate", matching, father_matching = opposite), 0.3, seed = 52)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 2))
})

test_that("group matching validates edited models and missing identities even at zero duration", {
  matching <- create_matching_model(~ same_group, c(same_group = 1))
  model <- make_model(0)
  model$processes$union_formation$matching <- matching
  expect_error(simulate(model, 0), "group identities")
  model$population$persons$group <- character()
  expect_no_error(simulate(model, 1))
  matching <- create_matching_model(~ age_difference, c(age_difference = 0), group_exogamy = c("A:B" = 0))
  matching$group_exogamy <- c("A:B" = 0, "A::B" = 0)
  expect_error(set_union_formation(make_model(0), "on_demand", matching), "overlapping")
  model$processes$union_formation$matching <- matching
  expect_error(simulate(model, 0), "overlapping")
})

test_that("zero or absent pair effects preserve histories and RNG consumption", {
  matching <- create_age_gap_matching_model(0, 5)
  model <- group_choice_model(matching = matching)
  set.seed(654)
  original <- simulate(model, 1)
  state <- .Random.seed
  model$processes$union_formation$matching$group_exogamy <- c("A::B" = 0, "absent:C" = 5)
  set.seed(654)
  expect_identical(simulate(model, 1), original)
  expect_identical(.Random.seed, state)
})

test_that("custom group comparisons work with absent candidate levels", {
  matching <- create_matching_model(~ age_difference + I(male_group == "B"),
    c(age_difference = 0.1, "I(male_group == \"B\")TRUE" = log(2)),
    group_levels = c("A", "B"), matching_method = "female_choice")
  prepared <- prepare_matching_model(matching, c("A", "B"))
  expect_equal(inspect_group_matching_cpp(c(1, 3), prepared, 0L, c(0L, 0L))$utilities,
               c(0.1, 0.3))
  expect_equal(inspect_group_matching_cpp(c(1, 3), prepared, 0L, c(0L, 1L))$utilities,
               c(0.1, 0.3 + log(2)))
})

test_that("pair tables support numeric and colon-containing group labels", {
  matching <- create_matching_model(~ same_group, c(same_group = 0.2),
    group_exogamy = data.frame(female_group = 20, male_group = 10, log_odds = log(2)), matching_method = "female_choice")
  prepared <- prepare_matching_model(matching, c(20, 10))
  expect_equal(inspect_group_matching_cpp(c(0, 0), prepared, 0L, 0:1)$utilities, c(0.2, log(2)))
  expect_equal(calculate_matching_utilities(matching, c(0, 0), 20, c(20, 10)), c(0.2, log(2)))
  matching <- create_matching_model(~ 1, numeric(), group_exogamy = data.frame(
    female_group = "A:B", male_group = "C::D", log_odds = -1), matching_method = "female_choice")
  prepared <- prepare_matching_model(matching, c("A:B", "C::D"))
  expect_equal(inspect_group_matching_cpp(c(0, 0), prepared, 0L, 0:1)$utilities, c(0, -1))
})
