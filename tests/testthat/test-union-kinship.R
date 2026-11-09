make_kinship_pair_model <- function(ego, relative, generations = 2) {
  population <- make_kinship_population()
  population$persons$death_time <- -0.1
  population$persons$death_time[population$persons$person_id %in% c(ego, relative)] <- NA_real_
  age <- -population$persons$birth_time[match(ego, population$persons$person_id)]
  hammel_model() |>
    set_population(population) |>
    set_mortality(data.frame(sex = c("f", "m"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = age, age_end = age + 1, rate = 20), birth_interval = 0) |>
    set_union_formation("on_demand",
      matching_model(~ age_difference, c(age_difference = 0)),
      min_male_age = 0, kinship_generations = generations,
      cancel_birth_if_no_match = FALSE)
}

test_that("R extraction and native eligibility agree on known relatives", {
  cases <- data.frame(ego = c(rep(80, 9), 110, 110, 110),
    relative = c(70, 85, 86, 60, 90, 95, 20, 150, 120, 120, 120, 20),
    generations = c(rep(2, 10), 3, 3))
  for (j in seq_len(nrow(cases))) {
    model <- make_kinship_pair_model(cases$ego[j], cases$relative[j], cases$generations[j])
    known <- nrow(extract_kin(model$population, cases$ego[j],
      generations = cases$generations[j], relative_ids = cases$relative[j])) > 0
    result <- simulate(model, 0.2, seed = 52)
    expect_gt(result$simulation$births, 0)
    if (known) {
      expect_equal(nrow(result$unions), 0)
      newborn <- result$persons$birth_time >= 0
      expect_true(all(is.na(result$persons$father_id[newborn])))
    } else {
      expect_equal(result$unions$husband_id, cases$relative[j])
    }
  }
})

test_that("newborn parent links enter screening and survive restart reconstruction", {
  model <- make_model(1, sex = "f", birth_rate = 20, birth_interval = 0)
  model$processes$union_formation$min_male_age <- 0
  result <- simulate(model, 0.5, seed = 52)
  expect_true(any(result$persons$sex[-1] == "m"))
  expect_gt(result$simulation$births, 0)
  expect_equal(nrow(result$unions), 0)
  expect_true(all(is.na(result$persons$father_id)))
  continued <- simulate(set_population(model, result), 0.5, seed = 53)
  expect_equal(nrow(continued$unions), 0)
  expect_true(all(is.na(continued$persons$father_id)))
  kin <- extract_kin(continued, 1)
  expect_true(all(kin$ego_distance == 0 & kin$relative_distance == 1))
  expect_equal(sort(unique(kin$relative_id)), sort(continued$persons$person_id[-1]))
  expect_identical(result, simulate(model, 0.5, seed = 52))
  # Zero generations disables cross-person exclusions, retaining the old case.
  model$processes$union_formation$kinship_generations <- 0
  unrestricted <- simulate(model, 0.5, seed = 52)
  expect_equal(nrow(unrestricted$unions), 1)
})

test_that("existing unions remain intact and only demonstrated relationships are screened", {
  model <- make_kinship_pair_model(80, 85)
  model$population$unions <- tibble::tibble(union_id = 5, wife_id = 80,
    husband_id = 85, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  result <- simulate(model, 0.2, seed = 52)
  expect_identical(result$unions, model$population$unions)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 85))
  founders <- make_model(2, birth_rate = 20, birth_interval = 0)
  founders$population$persons$sex[1] <- "f"
  expect_equal(nrow(simulate(founders, 0.2, seed = 52)$unions), 1)
})

test_that("kinship settings default to two and reject invalid or edited depths", {
  matching <- matching_model(~ age_difference, c(age_difference = 0))
  expect_equal(union_formation("on_demand", matching)$kinship_generations, 2)
  expect_equal(set_union_formation(hammel_model(), "on_demand", matching,
    kinship_generations = 3)$processes$union_formation$kinship_generations, 3)
  for (value in list(NULL, NA_real_, Inf, -1, 1.5, "2", c(1, 2), matrix(2))) {
    expect_error(union_formation("on_demand", matching, kinship_generations = value), "kinship_generations")
    model <- make_model(0)
    model$processes$union_formation$kinship_generations <- value
    expect_error(simulate(model, 1), "kinship_generations")
  }
})

test_that("switching between women does not retain another woman's kin markers", {
  model <- make_model(4, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m", "f", "m")
  model$population$persons$father_id <- c(2, NA, NA, 2)
  model$population$persons$birth_time[2] <- -60
  model$population$persons$death_time[2] <- -0.1
  result <- simulate(model, 0.5, seed = 52)
  expect_equal(result$unions$wife_id, 3)
  expect_equal(result$unions$husband_id, 4)
  expect_true(all(is.na(result$persons$father_id[result$persons$mother_id %in% 1])))
  expect_true(all(result$persons$father_id[result$persons$mother_id %in% 3] == 4))
})
