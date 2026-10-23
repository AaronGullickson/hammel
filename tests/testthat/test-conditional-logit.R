test_that("requested quadratic model peaks at 2.57 years and matches analytical probabilities", {
  matching <- create_matching_model(~ age_difference + I(age_difference^2),
    c(age_difference = 0.072, "I(age_difference^2)" = -0.014))
  peak <- -0.072 / (2 * -0.014)
  expect_equal(peak, 2.57142857142857)
  gaps <- c(peak - 1, peak, peak + 1, 0, 3, 10)
  native <- inspect_matching_cpp(gaps, prepare_matching_model(matching), 0)
  expected <- 0.072 * gaps - 0.014 * gaps^2
  expect_equal(native$utilities, expected)
  expect_equal(which.max(native$utilities), 2L)
  expect_equal(native$utilities[1], native$utilities[3])
  expect_equal(native$probabilities, exp(expected) / sum(exp(expected)))
  expect_equal(native$probabilities, calculate_matching_probabilities(matching, gaps))
  # Coefficient order and formula whitespace do not change scoring.
  reversed <- create_matching_model(~ 0 + age_difference + I(age_difference ^ 2),
    c("I(age_difference^2)" = -0.014, age_difference = 0.072))
  expect_equal(inspect_matching_cpp(gaps, prepare_matching_model(reversed), 0)$utilities,
               native$utilities)
})

test_that("partner draws follow softmax probabilities, including rare alternatives", {
  matching <- create_matching_model(~ age_difference + I(age_difference^2),
    c(age_difference = 0.072, "I(age_difference^2)" = -0.014))
  set.seed(8)
  n <- 50000
  result <- inspect_matching_cpp(c(0, 3, 10), prepare_matching_model(matching), n)
  probability <- exp(c(0, 0.09, -0.68)) / sum(exp(c(0, 0.09, -0.68)))
  frequency <- tabulate(result$choices, nbins = 3) / n
  expect_true(all(abs(frequency - probability) < 6 * sqrt(probability * (1 - probability) / n)))
  expect_true(all(frequency > 0))
  uniform <- create_matching_model(~ age_difference, c(age_difference = 0))
  expect_equal(inspect_matching_cpp(c(-10, 2, 80), prepare_matching_model(uniform), 0)$probabilities,
               rep(1/3, 3))
})

test_that("formula fallback has the same utilities and probabilities as R model matrices", {
  matching <- create_matching_model(~ abs(age_difference - 2.5),
    c("abs(age_difference - 2.5)" = -0.2))
  expect_length(prepare_matching_model(matching)$quadratic, 0)
  gaps <- c(-10, -0.5, 2.5, 9.7)
  native <- inspect_matching_cpp(gaps, prepare_matching_model(matching), 0)
  expect_equal(native$utilities, -0.2 * abs(gaps - 2.5))
  expect_equal(native$probabilities, calculate_matching_probabilities(matching, gaps))
  # A caller-defined deterministic transformation retains its formula semantics.
  transform_gap <- function(x) (x - 3)^2
  custom <- create_matching_model(~ transform_gap(age_difference),
    c("transform_gap(age_difference)" = -0.01))
  expect_equal(inspect_matching_cpp(gaps, prepare_matching_model(custom), 0)$utilities,
               calculate_matching_utilities(custom, gaps))
})

test_that("stable softmax handles extreme utilities, empty pools, and invalid scores", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 1000))
  prepared <- prepare_matching_model(matching)
  result <- inspect_matching_cpp(c(0, 1, 2), prepared, 20)
  expect_equal(result$probabilities, c(0, 0, 1))
  expect_true(all(result$choices == 3L))
  expect_true(all(inspect_matching_cpp(2.5, prepared, 10)$choices == 1L))
  empty <- inspect_matching_cpp(numeric(), prepared, 2)
  expect_length(empty$probabilities, 0)
  expect_true(all(is.na(empty$choices)))
  expect_error(inspect_matching_cpp(Inf, prepared, 0), "finite")
  expect_error(inspect_matching_cpp(1e308, prepared, 0), "finite")
  bad <- prepared
  bad$quadratic <- numeric()
  bad$evaluate <- function(x) rep(NA_real_, length(x))
  expect_error(inspect_matching_cpp(c(1, 2), bad, 0), "finite")
  bad$evaluate <- function(x) 1
  expect_error(inspect_matching_cpp(c(1, 2), bad, 0), "wrong number")
})

test_that("female-choice birth-time selection uses formula scores and male-minus-female gap", {
  model <- make_model(3, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex[1] <- "f"
  model$population$persons$birth_time <- c(-20, -20, -30)
  positive <- create_matching_model(~ age_difference, c(age_difference = 1000))
  model <- set_union_formation(model, "on_demand", positive, sample_size = 0,
    matching_method = "female_choice")
  result <- simulate(model, 0.5, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 1] == 3))
  expect_identical(result, simulate(model, 0.5, seed = 52))
  negative <- create_matching_model(~ age_difference, c(age_difference = -1000))
  changed <- simulate(set_union_formation(model, "on_demand", negative, 0,
    matching_method = "female_choice"), 0.5, seed = 52)
  expect_equal(changed$unions$husband_id, 2)
  fallback <- create_matching_model(~ abs(age_difference - 10),
    c("abs(age_difference - 10)" = -1000))
  expect_equal(simulate(set_union_formation(model, "on_demand", fallback, 0,
                        matching_method = "female_choice"),
                        0.5, seed = 52)$unions$husband_id, 3)
})
