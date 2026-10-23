test_that("formula matching gives expected conditional-logit probabilities", {
  matching <- create_matching_model(~ age_difference,
                                    c(age_difference = log(2)))
  expect_equal(calculate_matching_probabilities(matching, c(0, 1, 2)),
               c(1, 2, 4) / 7)
  expect_equal(calculate_matching_probabilities(matching, numeric()), numeric())
  expect_equal(calculate_matching_probabilities(matching, 3), 1)
  quadratic <- create_matching_model(~ age_difference + I(age_difference^2),
    c("I(age_difference^2)" = -1, age_difference = 2))
  expect_equal(calculate_matching_probabilities(quadratic, c(0, 1, 2)),
               exp(c(0, 1, 0)) / sum(exp(c(0, 1, 0))))
  no_intercept <- create_matching_model(~ 0 + age_difference,
                                        c(age_difference = log(2)))
  expect_equal(calculate_matching_probabilities(no_intercept, c(0, 1)), c(1, 2)/3)
  large <- create_matching_model(~ age_difference, c(age_difference = 1000))
  expect_equal(calculate_matching_probabilities(large, c(1, 2)), c(0, 1))
  expect_error(calculate_matching_probabilities(large, NA_real_))
})

test_that("matching specifications reject ambiguous or invalid formulas and coefficients", {
  expect_error(create_matching_model(age ~ age_difference, c(age_difference = 1)))
  expect_error(create_matching_model(~ age, c(age = 1)))
  expect_error(create_matching_model(~ age_difference, 1))
  expect_error(create_matching_model(~ age_difference, c(wrong_name = 1)))
  expect_error(create_matching_model(~ age_difference, c(age_difference = Inf)))
  expect_error(create_matching_model(~ age_difference,
                                    c(age_difference = 1, age_difference = 2)))
  expect_error(create_matching_model(~ offset(age_difference), c(age_difference = 1)))
  expect_error(create_matching_model(~ log(age_difference), c("log(age_difference)" = 1)))
})

test_that("union processes construct and replace without modifying other processes", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0))
  rates <- data.frame(age_start = 18, age_end = Inf, rate = 0.1)
  attr(rates, "source") <- list(database = "synthetic")
  process <- union_formation("rate", matching, rates = rates)
  expect_s3_class(process, "hammel_process")
  expect_s3_class(process, "hammel_union_formation")
  expect_s3_class(process$rates, "tbl_df")
  expect_identical(process$source, attr(rates, "source"))
  model <- set_fertility(hammel_model(), data.frame(age_start = 20, age_end = 30, rate = 0.1))
  updated <- set_process(model, process)
  expect_equal(union_formation("on_demand", matching)$matching$sample_size, 250)
  expect_false("sample_size" %in% names(process))
  expect_identical(set_union_formation(model, "rate", matching, rates = rates),
                   set_process(model, process))
  expect_null(model$processes$union_formation)
  expect_identical(updated$processes$fertility, model$processes$fertility)
  expect_identical(updated, set_union_formation(model, "rate", matching, rates = rates))
  full_pool <- matching
  full_pool$sample_size <- 0
  replaced <- set_union_formation(updated, "on_demand", full_pool)
  expect_null(replaced$processes$union_formation$rates)
  expect_identical(replaced$processes$fertility, model$processes$fertility)
  edited <- process
  edited$matching$coefficients <- c(other = 1)
  expect_error(set_process(model, edited))
  edited <- process
  edited$matching$sample_size <- -1
  expect_error(set_process(model, edited))
  edited <- process
  edited$rates$rate <- -1
  expect_error(set_process(model, edited))
})

test_that("union process validation and execution boundaries are explicit", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0))
  rates <- data.frame(age_start = 18, age_end = Inf, rate = 0.1)
  expect_error(union_formation("other", matching))
  expect_error(union_formation("rate", matching))
  expect_error(union_formation("on_demand", matching, rates = rates))
  expect_error(union_formation("rate", matching, rates = rbind(rates, rates)))
  model <- hammel_model() |>
    set_population(create_baseline_population(0)) |>
    set_mortality(data.frame(sex = c("f", "m"), age_start = 0,
                            age_end = Inf, rate = 0.01)) |>
    set_fertility(data.frame(age_start = 20, age_end = 30, rate = 0.1)) |>
    set_union_formation("rate", matching, rates = rates)
  expect_warning(result <- simulate(model, years = 1), "overall fertility schedule")
  expect_equal(nrow(result$unions), 0)
})

test_that("minimum male age defaults to 15 and is configurable and revalidated", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0))
  expect_equal(union_formation("on_demand", matching)$min_male_age, 15)
  model <- set_union_formation(hammel_model(), "on_demand", matching, min_male_age = 16.5)
  expect_equal(model$processes$union_formation$min_male_age, 16.5)
  expect_equal(union_formation("on_demand", matching, min_male_age = 0)$min_male_age, 0)
  for (age in list(NULL, NA_real_, Inf, -1, "15", c(15, 16), matrix(15))) {
    expect_error(union_formation("on_demand", matching, min_male_age = age), "min_male_age")
    edited <- model$processes$union_formation
    edited$min_male_age <- age
    expect_error(set_process(hammel_model(), edited), "min_male_age")
  }
})
