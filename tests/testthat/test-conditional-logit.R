test_that("requested quadratic model peaks at 2.57 years and matches analytical probabilities", {
  matching <- create_matching_model(~ age_difference + I(age_difference^2),
    c(age_difference = 0.072, "I(age_difference^2)" = -0.014))
  peak <- -0.072 / (2 * -0.014)
  expect_equal(peak, 2.57142857142857)
  gaps <- c(peak - 1, peak, peak + 1, 0, 3, 10)
  native <- inspect_matching_cpp(gaps, prepare_matching_model(matching))
  expected <- 0.072 * gaps - 0.014 * gaps^2
  expect_equal(native, expected)
  expect_equal(which.max(native), 2L)
  expect_equal(native[1], native[3])
  expect_equal(calculate_matching_probabilities(matching, gaps), exp(expected) / sum(exp(expected)))
  # Coefficient order and formula whitespace do not change scoring.
  reversed <- create_matching_model(~ 0 + age_difference + I(age_difference ^ 2),
    c("I(age_difference^2)" = -0.014, age_difference = 0.072))
  expect_equal(inspect_matching_cpp(gaps, prepare_matching_model(reversed)), native)
})

# A frozen pool: a woman born at -20 and men at given male-minus-female gaps.
exact_draws <- function(matching, gaps, draws, groups = 1L, male_group = integer(length(gaps)),
                        female_group = 0L, log_weights = numeric(groups), excluded = logical(length(gaps)),
                        labels = NULL, female_born = -20) {
  inspect_exact_draw_cpp(female_born - gaps, male_group, excluded, female_born, female_group, log_weights,
                         prepare_matching_model(as_internal_matching(matching, labels), labels),
                         groups, 0, 15, draws)
}

test_that("exact full-pool draws follow conditional-logit probabilities, including rare men", {
  matching <- create_matching_model(~ age_difference + I(age_difference^2),
    c(age_difference = 0.072, "I(age_difference^2)" = -0.014))
  set.seed(8)
  n <- 50000
  result <- exact_draws(matching, c(0, 3, 10), n)
  probability <- exp(c(0, 0.09, -0.68)) / sum(exp(c(0, 0.09, -0.68)))
  expect_equal(result$probability, probability)
  frequency <- result$counts / n
  expect_true(all(abs(frequency - probability) < 6 * sqrt(probability * (1 - probability) / n)))
  expect_true(all(frequency > 0))
  expect_equal(result$failures, 0L)
  expect_equal(result$approximate, 0)
  uniform <- create_matching_model(~ age_difference, c(age_difference = 0))
  expect_equal(exact_draws(uniform, c(-4, 2, 30), 0)$probability, rep(1/3, 3))
})

test_that("exact draws handle pools, age edges, kin, group weights and custom formulas", {
  set.seed(9)
  matching <- create_age_gap_matching_model(1, 5)
  # A large pool across many bins, including men near and under the age cutoff.
  gaps <- c(runif(300, -6, 40), -5.05, -4.95, -5.2)   # woman aged 20: men aged 14.8 to 60
  n <- 200000
  result <- exact_draws(matching, gaps, n)
  under <- 20 + gaps < 15
  expect_true(all(result$counts[under] == 0) && all(result$probability[under] == 0))
  cell <- cut(20 + gaps, c(14, 15, 16, 18, 20, 22, 25, 30, 40, 70))
  observed <- tapply(result$counts, cell, sum)
  expected <- tapply(result$probability, cell, sum) * n
  keep <- expected > 50
  expect_lt(sum((observed[keep] - expected[keep])^2 / expected[keep]),
            qchisq(0.999, sum(keep) - 1))
  # Excluded (kin) men are never drawn; with all but one excluded the fallback
  # enumeration returns him, and with all excluded the draw fails.
  excluded <- seq_along(gaps) != 5L
  only <- exact_draws(matching, gaps, 200, excluded = excluded)
  expect_true(all(only$counts[-5] == 0) && only$counts[5] == 200)
  expect_equal(exact_draws(matching, gaps, 10, excluded = rep(TRUE, length(gaps)))$failures, 10L)
  expect_equal(exact_draws(matching, numeric(), 3)$failures, 3L)
  # Group weights and hard restrictions enter the exact probabilities.
  pair <- create_age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1, "A:C" = NA))
  groups <- rep(0:2, 100)
  result <- exact_draws(pair, gaps[1:300], 150000, groups = 3L, male_group = groups,
                        log_weights = c(0, 0.7, 0.2), labels = c("A", "B", "C"))
  expect_true(all(result$counts[groups == 2] == 0))
  share <- tapply(result$counts, groups, sum) / 150000
  exact_share <- tapply(result$probability, groups, sum)
  allowed <- exact_share > 0
  expect_equal(unname(exact_share[3]), 0)
  expect_true(all(abs(share - exact_share)[allowed] <
                  6 * sqrt(exact_share * (1 - exact_share) / 150000)[allowed]))
  # A custom formula uses its fine-grid bound; smooth formulas stay exact.
  custom <- create_matching_model(~ abs(age_difference - 2), c("abs(age_difference - 2)" = -0.3))
  result <- exact_draws(custom, gaps, n)
  expect_equal(result$approximate, 0)
  observed <- tapply(result$counts, cell, sum)
  expected <- tapply(result$probability, cell, sum) * n
  keep <- expected > 50
  expect_lt(sum((observed[keep] - expected[keep])^2 / expected[keep]),
            qchisq(0.999, sum(keep) - 1))
})

test_that("formula fallback has the same utilities as R model matrices", {
  matching <- create_matching_model(~ abs(age_difference - 2.5),
    c("abs(age_difference - 2.5)" = -0.2))
  expect_length(prepare_matching_model(matching)$quadratic, 0)
  gaps <- c(-10, -0.5, 2.5, 9.7)
  expect_equal(inspect_matching_cpp(gaps, prepare_matching_model(matching)), -0.2 * abs(gaps - 2.5))
  # A caller-defined deterministic transformation retains its formula semantics.
  transform_gap <- function(x) (x - 3)^2
  custom <- create_matching_model(~ transform_gap(age_difference),
    c("transform_gap(age_difference)" = -0.01))
  expect_equal(inspect_matching_cpp(gaps, prepare_matching_model(custom)),
               calculate_matching_utilities(custom, gaps))
})

test_that("extreme utilities draw the best man and invalid scores are errors", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 1000))
  prepared <- prepare_matching_model(matching)
  result <- exact_draws(matching, c(0, 1, 2), 50)
  expect_equal(result$probability, c(0, 0, 1))
  expect_true(all(result$counts == c(0, 0, 50)))
  expect_equal(exact_draws(matching, 2.5, 10)$counts, 10L)
  expect_error(inspect_matching_cpp(Inf, prepared), "finite")
  expect_error(inspect_matching_cpp(1e308, prepared), "finite")
  bad <- prepared
  bad$quadratic <- numeric()
  bad$evaluate <- function(x) rep(NA_real_, length(x))
  expect_error(inspect_matching_cpp(c(1, 2), bad), "finite")
  bad$evaluate <- function(x) 1
  expect_error(inspect_matching_cpp(c(1, 2), bad), "wrong number")
})

test_that("birth-time selection uses formula scores and male-minus-female gap", {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time <- c(-20, -20, -30)
  positive <- create_matching_model(~ age_difference, c(age_difference = 1000))
  model <- set_union_formation(model, "on_demand", positive)
  result <- simulate(model, 0.5, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 1] == 3))
  expect_identical(result, simulate(model, 0.5, seed = 52))
  negative <- create_matching_model(~ age_difference, c(age_difference = -1000))
  changed <- simulate(set_union_formation(model, "on_demand", negative), 0.5, seed = 52)
  expect_equal(changed$unions$husband_id, 2)
  fallback <- create_matching_model(~ abs(age_difference - 10),
    c("abs(age_difference - 10)" = -1000))
  expect_equal(simulate(set_union_formation(model, "on_demand", fallback),
                        0.5, seed = 52)$unions$husband_id, 3)
})
