make_nonmarital_model <- function(n = 2, birth_rate = 20, search_rate = 0) {
  model <- make_rate_model(n, search_rate = search_rate)
  set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(birth_rate, birth_rate)), 0)
}

test_that("nonmarital fathers stay available without forming unions and enter kinship", {
  model <- make_nonmarital_model(3)
  model$population$persons$sex <- c("f", "m", "f")
  model$population$persons$person_id <- c(10, 40, 70)
  set.seed(123)
  rng <- .Random.seed
  result <- simulate(model, 0.5, seed = 52)
  expect_identical(.Random.seed, rng)
  expect_identical(result, simulate(model, 0.5, seed = 52))
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 2)
  expect_true(all(c(10, 70) %in% children$mother_id))
  expect_true(all(children$father_id == 40))
  expect_equal(nrow(result$unions), 0)
  kin <- extract_kin(result, children$person_id[1], relative_ids = 40)
  expect_equal(kin$common_ancestor_id, 40)
  expect_equal(kin$ego_distance, 1)
  expect_equal(kin$relative_distance, 0)
  continued <- simulate(set_population(model, result), 0.5, seed = 53)
  expect_equal(nrow(continued$unions), 0)
  expect_true(all(continued$persons$father_id[continued$persons$birth_time >= 0.5] == 40))
})

test_that("father searches exclude partnered men, boys and demonstrated kin", {
  for (reason in c("partnered", "young", "kin", "dead")) {
    model <- make_nonmarital_model(3)
    model$population$persons$sex[3] <- "f"
    model$population$persons$birth_time[3] <- -60 # no births from this woman
    if (reason == "partnered") {
      model$population$unions <- tibble::tibble(union_id = 1, wife_id = 3,
        husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
    }
    if (reason == "young") model$population$persons$birth_time[2] <- -14
    if (reason == "kin") {
      model$population$persons$birth_time[2] <- -50
      model$population$persons$father_id[1] <- 2
    }
    if (reason == "dead") model$population$persons$death_time[2] <- -0.1
    result <- simulate(model, 0.5, seed = 52)
    children <- result$persons[result$persons$birth_time >= 0, ]
    expect_gt(nrow(children), 0)
    expect_true(all(is.na(children$father_id)), info = reason)
    expect_equal(result$unions, model$population$unions)
  }
})

test_that("female-choice father selection uses the entered conditional-logit probabilities", {
  model <- make_nonmarital_model(3, birth_rate = 500)
  model$population$persons$birth_time[3] <- -30
  model$processes$union_formation$matching <- create_matching_model(
    ~ age_difference, c(age_difference = log(3) / 10), matching_method = "female_choice")
  result <- simulate(model, 1, seed = 72)
  children <- result$persons[result$persons$birth_time >= 0, ]
  n <- nrow(children)
  expect_gt(n, 400)
  expect_true(all(children$father_id %in% c(2, 3)))
  expect_lt(abs(sum(children$father_id == 3) - 0.75 * n),
            6 * sqrt(n * 0.75 * 0.25))
  expect_equal(nrow(result$unions), 0)
})

test_that("a nonmarital father may later be selected for a union", {
  model <- make_nonmarital_model()
  # No union hazard before age 20.5. The sole eligible man remains available.
  model$processes$union_formation$rates <- data.frame(
    age_start = 20.5, age_end = Inf, rate = 100)
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$husband_id, 2)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_true(any(children$birth_time < result$unions$start_time[1]))
  expect_true(any(children$birth_time > result$unions$start_time[1]))
  expect_true(all(children$father_id == 2))
})

test_that("a father's death removes him from later nonmarital choices", {
  model <- make_nonmarital_model()
  model$population$persons$birth_time[2] <- -30
  model <- set_mortality(model, data.frame(sex = c("f", "m", "m"),
    age_start = c(0, 0, 30.3), age_end = c(Inf, 30.3, Inf), rate = c(0, 0, 1e8)))
  result <- simulate(model, 1, seed = 52)
  children <- result$persons[result$persons$birth_time >= 0, ]
  death <- result$persons$death_time[2]
  expect_true(any(children$birth_time < death))
  expect_true(any(children$birth_time > death))
  expect_true(all(children$father_id[children$birth_time < death] == 2))
  expect_true(all(is.na(children$father_id[children$birth_time > death])))
})

test_that("father_matching defaults to the general model and validates edited objects", {
  general <- create_matching_model(~ age_difference, c(age_difference = 0.1))
  process <- union_formation("on_demand", general)
  expect_null(process$father_matching)
  expect_identical(union_formation("on_demand", general, father_matching = general)$father_matching,
                   general)
  expect_error(union_formation("on_demand", general, father_matching = "invalid"), "matching")
  model <- make_nonmarital_model(3)
  model$population$persons$birth_time[3] <- -30
  model$processes$union_formation$matching <- general
  baseline <- simulate(model, 0.5, seed = 52)
  explicit <- model
  explicit$processes$union_formation$father_matching <- general
  expect_identical(simulate(explicit, 0.5, seed = 52), baseline)
  # Resolve NULL against the current general model, including edited models.
  model$processes$union_formation$matching$coefficients[] <- -100
  explicit$processes$union_formation$father_matching <- model$processes$union_formation$matching
  expect_identical(simulate(model, 0.5, seed = 52), simulate(explicit, 0.5, seed = 52))
  explicit$processes$union_formation$father_matching$coefficients <- c(other = 1)
  expect_error(set_process(hammel_model(), explicit$processes$union_formation), "Coefficient names")
  expect_error(simulate(explicit, 0), "Coefficient names")
  updated <- set_union_formation(model, "rate", general,
    rates = model$processes$union_formation$rates, father_matching = general)
  expect_identical(updated$processes$union_formation$father_matching, general)
})

test_that("father_matching changes nonmarital fathers but not later union choices", {
  model <- make_nonmarital_model(3)
  model$population$persons$birth_time[3] <- -30
  model$processes$union_formation$matching <- create_matching_model(
    ~ age_difference, c(age_difference = 100), matching_method = "female_choice")
  model$processes$union_formation$father_matching <- create_matching_model(
    ~ age_difference, c(age_difference = -100), matching_method = "female_choice")
  model$processes$union_formation$rates <- data.frame(
    age_start = 20.5, age_end = Inf, rate = 100)
  result <- simulate(model, 1, seed = 52)
  expect_equal(result$unions$husband_id, 3)
  children <- result$persons[result$persons$birth_time >= 0, ]
  before <- children$birth_time < result$unions$start_time[1]
  expect_true(any(before))
  expect_true(any(!before))
  expect_true(all(children$father_id[before] == 2))
  expect_true(all(children$father_id[!before] == 3))
})

test_that("father_matching supports the R formula fallback and is unused on demand", {
  model <- make_nonmarital_model(3)
  model$population$persons$birth_time[3] <- -30
  model$processes$union_formation$father_matching <- create_matching_model(
    ~ age_difference, c(age_difference = -0.1))
  native <- simulate(model, 1, seed = 72)
  model$processes$union_formation$father_matching <- create_matching_model(
    ~ I(age_difference^3), c("I(age_difference^3)" = -0.001))
  custom <- simulate(model, 1, seed = 72)
  # With available gaps 0 and 10 both specifications give utilities 0 and -1.
  # Native and custom formulas use different within-bin bounds, so they agree
  # in distribution rather than draw by draw.
  expected <- exp(-1) / (1 + exp(-1))
  for (result in list(native, custom)) {
    children <- result$persons[result$persons$birth_time >= 0, ]
    expect_true(all(children$father_id %in% c(2, 3)))
    n <- nrow(children)
    if (n) expect_lt(abs(mean(children$father_id == 3) - expected), 6 * sqrt(expected * (1 - expected) / n) + 1e-9)
  }
  model$processes$union_formation$initiation <- "on_demand"
  model$processes$union_formation$rates <- NULL
  model <- set_fertility(model, data.frame(age_start = 15, age_end = 50, rate = 20), 0)
  alternative <- simulate(model, 0.5, seed = 52)
  model$processes$union_formation$father_matching <- NULL
  expect_identical(simulate(model, 0.5, seed = 52), alternative)
})
