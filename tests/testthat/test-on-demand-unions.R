make_on_demand_model <- function(n = 2, sample_size = 1000) {
  model <- make_model(n, sex = "m", birth_rate = 20, birth_interval = 0)
  if (n) model$population$persons$sex[1] <- "f"
  set_union_formation(model, "on_demand",
    create_matching_model(~ age_difference, c(age_difference = 100)), sample_size)
}

test_that("male pool sampling has correct size, uniqueness and inclusion probabilities", {
  expect_length(sample_male_positions_cpp(0, 1000), 0)
  expect_identical(sample_male_positions_cpp(3, 1000), 0:2)
  expect_identical(sample_male_positions_cpp(3, 0), 0:2)
  expect_identical(sample_male_positions_cpp(1, 1), 0L)
  set.seed(47)
  samples <- replicate(4000, sample_male_positions_cpp(8, 3))
  expect_true(all(apply(samples, 2, function(x) length(unique(x)) == 3)))
  expect_true(all(samples >= 0 & samples < 8))
  inclusion <- tabulate(as.vector(samples) + 1L, nbins = 8)
  expect_true(all(abs(inclusion - 1500) < 6 * sqrt(4000 * 3/8 * 5/8)))
})

test_that("on-demand birth creates one union and keeps the partner for later births", {
  model <- make_on_demand_model()
  model$population$persons$person_id <- c(10, 40)
  model$population$unions$note <- character()
  baseline <- simulate(make_model(), 1, seed = 7)
  set.seed(888)
  rng <- .Random.seed
  result <- simulate(model, 1, seed = 52)
  expect_identical(.Random.seed, rng)
  expect_identical(result, simulate(model, 1, seed = 52))
  expect_identical(simulate(make_model(), 1, seed = 7), baseline)
  children <- result$persons[result$persons$mother_id %in% 10, ]
  expect_gt(nrow(children), 2)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$partner1_id, 10)
  expect_equal(result$unions$partner2_id, 40)
  expect_equal(result$unions$start_time, min(children$birth_time))
  expect_true(all(children$father_id == 40))
  expect_true(is.na(result$unions$end_time))
  expect_true(is.na(result$unions$note))
  expect_s3_class(result$unions, "tbl_df")
  expect_equal(nrow(model$population$unions), 0)
  later <- simulate(set_population(model, result), 1, seed = 53)
  expect_equal(nrow(later$unions), 1)
  expect_true(all(later$persons$father_id[later$persons$mother_id %in% 10] == 40))
})

test_that("existing unions preserve partner identity, orientation and custom history", {
  model <- make_on_demand_model(3)
  model$population$unions <- tibble::tibble(union_id = 9, partner1_id = 3,
    partner2_id = 1, start_time = -1, end_time = NA_real_,
    end_reason = NA_character_, note = "retained")
  result <- simulate(model, 1, seed = 5)
  expect_identical(result$unions, model$population$unions)
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 1] == 3))
  model$population$unions$end_time <- -0.5
  model$population$unions$end_reason <- "dissolution"
  ended <- simulate(model, 1, seed = 5)
  expect_equal(ended$unions$union_id, c(9, 10))
  expect_equal(ended$unions$end_reason[1], "dissolution")
  expect_equal(ended$unions$note[1], "retained")
  expect_true(is.na(ended$unions$note[2]))
})

test_that("no men preserves the first birth with unknown father and no union", {
  model <- make_on_demand_model(1)
  result <- simulate(model, 0.1, seed = 3)
  expect_gt(result$simulation$births, 0)
  first <- which(!is.na(result$persons$mother_id))[1]
  expect_true(is.na(result$persons$father_id[first]))
  expect_true(all(result$unions$start_time > result$persons$birth_time[first]))
  expect_equal(nrow(simulate(make_on_demand_model(0), 1, seed = 3)$unions), 0)
})

test_that("widowhood closes new unions and subsequent births search again", {
  model <- make_on_demand_model()
  model <- set_mortality(model, data.frame(
    sex = c("f", "m", "m", "m"), age_start = c(0, 0, 20, 20.5),
    age_end = c(Inf, 20, 20.5, Inf), rate = c(0, 0, 0, 1e8)))
  result <- simulate(model, 1, seed = 52)
  expect_gt(nrow(result$unions), 1)
  death <- result$persons$death_time[result$persons$person_id == 2]
  expect_equal(result$unions$end_time[1], death)
  expect_equal(result$unions$end_reason[1], "death")
  expect_gt(result$unions$start_time[2], death)
  children <- result$persons[result$persons$mother_id %in% 1, ]
  expect_true(all(children$father_id[children$birth_time < death] == 2))
  expect_true(all(children$father_id[children$birth_time > death] != 2))
  continued <- simulate(set_population(model, result), 0.1, seed = 53)
  expect_identical(continued$unions[1, ], result$unions[1, ])
})

test_that("screening and utility matching are deliberately not applied yet", {
  model <- make_on_demand_model(3, 1)
  model$population$persons$sex[2] <- "f"
  model$population$persons$birth_time[3] <- -1 # only man is a child
  result <- simulate(model, 1, seed = 19)
  expect_equal(nrow(result$unions), 2)
  expect_true(all(result$unions$partner2_id == 3))
  expect_equal(sort(result$unions$partner1_id), c(1, 2))
  model$processes$union_formation$matching$coefficients[] <- -100
  expect_identical(simulate(model, 1, seed = 19), result)
})

test_that("on-demand input and edited process errors are explicit", {
  model <- make_on_demand_model(3)
  model$population$persons$sex[3] <- "f"
  model$population$unions <- tibble::tibble(union_id = 1, partner1_id = 1,
    partner2_id = 3, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  expect_error(simulate(model, 1), "living female/male")
  model$population$persons$sex[3] <- "m"
  model$population$unions <- rbind(model$population$unions,
    transform(model$population$unions, union_id = 2, partner2_id = 2))
  expect_error(simulate(model, 1), "multiple ongoing unions")
  model <- make_on_demand_model()
  model$processes$union_formation$sample_size <- -1
  expect_error(simulate(model, 1), "sample_size")
})

test_that("simulations require an explicit union process even for zero duration", {
  model <- make_on_demand_model()
  model$processes$union_formation <- NULL
  expect_error(simulate(model, 1), "union-formation process")
  expect_error(simulate(model, 0), "union-formation process")
})
