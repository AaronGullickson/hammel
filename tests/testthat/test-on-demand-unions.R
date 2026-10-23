make_on_demand_model <- function(n = 2, sample_size = 250) {
  model <- make_model(n, sex = "m", birth_rate = 20, birth_interval = 0)
  if (n) model$population$persons$sex[1] <- "f"
  set_union_formation(model, "on_demand",
    create_matching_model(~ age_difference, c(age_difference = 100), sample_size = sample_size),
    cancel_birth_if_no_match = FALSE)
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
  expect_equal(result$unions$wife_id, 10)
  expect_equal(result$unions$husband_id, 40)
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

test_that("existing unions preserve sex-specific partner identity and custom history", {
  model <- make_on_demand_model(3)
  model$population$unions <- tibble::tibble(union_id = 9, wife_id = 1,
    husband_id = 3, start_time = -1, end_time = NA_real_,
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
  model <- make_on_demand_model(3)
  model$population$persons$birth_time[3] <- -19
  model <- set_mortality(model, data.frame(
    sex = c("f", "m", "m", "m"), age_start = c(0, 0, 20, 20.5),
    age_end = c(Inf, 20, 20.5, Inf), rate = c(0, 0, 0, 1e8)))
  result <- simulate(model, 1, seed = 52)
  expect_gt(nrow(result$unions), 1)
  death <- result$persons$death_time[result$persons$person_id == 2]
  expect_equal(result$unions$end_time[1], death)
  expect_equal(result$unions$end_reason[1], "death_husband")
  expect_gt(result$unions$start_time[2], death)
  children <- result$persons[result$persons$mother_id %in% 1, ]
  expect_true(all(children$father_id[children$birth_time < death] == 2))
  expect_true(all(children$father_id[children$birth_time > death] != 2))
  continued <- simulate(set_population(model, result), 0.1, seed = 53)
  expect_identical(continued$unions[1, ], result$unions[1, ])
})

test_that("zero minimum age allows a child and singleton choices do not depend on utilities", {
  model <- make_on_demand_model(3, 1)
  model$processes$union_formation$min_male_age <- 0
  model$population$persons$sex[2] <- "f"
  model$population$persons$birth_time[3] <- -1 # only man is a child
  result <- simulate(model, 1, seed = 19)
  expect_equal(nrow(result$unions), 2)
  expect_true(3 %in% result$unions$husband_id)
  expect_equal(anyDuplicated(result$unions$husband_id), 0L)
  expect_equal(sort(result$unions$wife_id), c(1, 2))
  model$processes$union_formation$matching$coefficients[] <- -100
  expect_identical(simulate(model, 1, seed = 19), result)
})

test_that("on-demand input and edited process errors are explicit", {
  model <- make_on_demand_model(3)
  model$population$persons$sex[3] <- "f"
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 3, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  expect_error(simulate(model, 1), "wife_id.*female")
  model$population$persons$sex[3] <- "m"
  model$population$unions <- rbind(model$population$unions,
    transform(model$population$unions, union_id = 2, husband_id = 2))
  expect_error(simulate(model, 1), "multiple ongoing unions")
  model <- make_on_demand_model()
  model$processes$union_formation$matching$sample_size <- -1
  expect_error(simulate(model, 1), "sample_size")
})

test_that("simulations require an explicit union process even for zero duration", {
  model <- make_on_demand_model()
  model$processes$union_formation <- NULL
  expect_error(simulate(model, 1), "union-formation process")
  expect_error(simulate(model, 0), "union-formation process")
})

test_that("union death reasons identify the wife or husband for input and new unions", {
  for (dead_sex in c("f", "m")) {
    for (initial_union in c(FALSE, TRUE)) {
      model <- make_on_demand_model()
      if (initial_union) {
        model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
          husband_id = 2, start_time = -1, end_time = NA_real_,
          end_reason = NA_character_)
      }
      # Allow a union to form before the selected partner dies after age 20.5.
      model <- set_mortality(model, data.frame(
        sex = c(dead_sex, dead_sex, if (dead_sex == "f") "m" else "f"),
        age_start = c(0, 20.5, 0), age_end = c(20.5, Inf, Inf),
        rate = c(0, 1e8, 0)))
      result <- simulate(model, 1, seed = 52)
      expect_gt(nrow(result$unions), 0)
      expect_equal(result$unions$end_reason[1],
                   if (dead_sex == "f") "death_wife" else "death_husband")
      deceased <- if (dead_sex == "f") 1 else 2
      expect_equal(result$unions$end_time[1], result$persons$death_time[deceased])
      later <- simulate(set_population(model, result), 0.1, seed = 53)
      expect_identical(later$unions[1, ], result$unions[1, ])
    }
  }
})

test_that("partnered men are absent before even a singleton sample", {
  model <- make_on_demand_model(4, 1)
  model$population$persons$sex <- c("f", "m", "f", "m")
  model$population$persons$birth_time <- c(-60, -25, -20, -23)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  for (seed in 1:10) {
    result <- simulate(model, 0.5, seed = seed)
    expect_equal(result$unions$husband_id, c(2, 4))
    expect_equal(result$unions$wife_id, c(1, 3))
  }
  model$population$unions <- rbind(model$population$unions,
    transform(model$population$unions, union_id = 2, wife_id = 3))
  expect_error(simulate(model, 1), "multiple ongoing unions for a man")
})

test_that("scarcity preserves births without reusing partnered men", {
  model <- make_on_demand_model(150)
  model$population$persons$sex[1:100] <- "f"
  result <- simulate(model, 1, seed = 19)
  expect_equal(nrow(result$unions), 50)
  expect_equal(anyDuplicated(result$unions$husband_id), 0L)
  expect_equal(anyDuplicated(result$unions$wife_id), 0L)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_true(any(is.na(children$father_id)))
  continued <- simulate(set_population(model, result), 0.5, seed = 20)
  expect_identical(continued$unions, result$unions)
})

test_that("a surviving husband returns immediately after his wife's death", {
  model <- make_on_demand_model(3, 1)
  model$population$persons$sex <- c("f", "m", "f")
  model$population$persons$birth_time <- c(-60, -25, -20)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  model <- set_mortality(model, data.frame(sex = c("f", "f", "m"),
    age_start = c(0, 60.25, 0), age_end = c(60.25, Inf, Inf),
    rate = c(0, 1e8, 0)))
  # Use the full available pool so infant boys cannot mask the returned widower.
  model$processes$union_formation$matching$sample_size <- 0
  result <- simulate(model, 1, seed = 52)
  expect_equal(result$unions$husband_id, c(2, 2))
  expect_equal(result$unions$wife_id, c(1, 3))
  death <- result$persons$death_time[1]
  expect_equal(result$unions$end_time[1], death)
  expect_equal(result$unions$end_reason[1], "death_wife")
  expect_gt(result$unions$start_time[2], death)
  children <- result$persons[result$persons$mother_id %in% 3, ]
  expect_true(any(children$birth_time < death))
  expect_true(all(is.na(children$father_id[children$birth_time < death])))
  expect_true(all(children$father_id[children$birth_time > death] == 2))
})
