test_that("failed-match cancellation defaults on and validates logical scalars", {
  matching <- create_matching_model(~ age_difference, c(age_difference = 0))
  expect_true(union_formation("on_demand", matching)$cancel_birth_if_no_match)
  expect_false(union_formation("on_demand", matching,
    cancel_birth_if_no_match = FALSE)$cancel_birth_if_no_match)
  for (bad in list(NULL, NA, 0, 1, "TRUE", c(TRUE, FALSE), matrix(TRUE))) {
    expect_error(union_formation("on_demand", matching, cancel_birth_if_no_match = bad),
                 "cancel_birth_if_no_match")
    model <- make_model(0)
    process <- model$processes$union_formation
    process$cancel_birth_if_no_match <- bad
    expect_error(set_process(model, process), "cancel_birth_if_no_match")
    model$processes$union_formation <- process
    expect_error(simulate(model, 0), "cancel_birth_if_no_match")
  }
})

test_that("both initiation modes cancel unmatched births unless explicitly allowed", {
  for (mode in c("on_demand", "rate")) {
    model <- make_model(1, sex = "f", birth_rate = 20, birth_interval = 0)
    if (mode == "rate") {
      model <- set_union_formation(model, mode, model$processes$union_formation$matching,
        rates = data.frame(age_start = 15, age_end = Inf, rate = 1))
      model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
        age_start = 15, age_end = 50, rate = 20), 0)
    }
    model$processes$union_formation$cancel_birth_if_no_match <- TRUE
    set.seed(381)
    rng <- .Random.seed
    result <- simulate(model, 1, seed = 52)
    expect_identical(.Random.seed, rng)
    expect_equal(result$persons, tibble::as_tibble(model$population$persons))
    expect_equal(result$simulation$births, 0)
    expect_equal(nrow(result$unions), 0)
    expect_identical(result, simulate(model, 1, seed = 52))
    model$processes$union_formation$cancel_birth_if_no_match <- FALSE
    allowed <- simulate(model, 1, seed = 52)
    expect_gt(allowed$simulation$births, 0)
    expect_true(all(is.na(allowed$persons$father_id)))
  }
})

test_that("cancelled attempts do not create spacing, IDs or kinship and can later succeed", {
  for (mode in c("on_demand", "rate")) {
    model <- make_model(2, birth_rate = 1, birth_interval = 0.9)
    model$population$persons$sex[1] <- "f"
    model$population$persons$birth_time[2] <- -14.5
    if (mode == "rate") {
      model <- set_union_formation(model, "rate", model$processes$union_formation$matching,
        rates = data.frame(age_start = 0, age_end = Inf, rate = 0))
      model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
        age_start = 15, age_end = 50, rate = 1), 0.9)
    }
    model$processes$union_formation$cancel_birth_if_no_match <- TRUE
    # Compensated birth hazard is 10/year. Many early failed attempts; the man
    # becomes eligible at .5. Applying .9 spacing after failure would delay success.
    result <- simulate(model, 0.9, seed = 52)
    children <- result$persons[result$persons$birth_time >= 0, ]
    expect_equal(nrow(children), 1)
    expect_equal(children$person_id, 3)
    expect_gte(children$birth_time, 0.5)
    expect_lt(children$birth_time, 0.9)
    expect_equal(children$father_id, 2)
    kin <- extract_kin(result, 2, relative_ids = 3)
    expect_equal(kin$relative_distance, 1)
    continued <- simulate(set_population(model, result), 0.2, seed = 53)
    expect_equal(continued$simulation$births, 0)
  }
})

test_that("default cancellation covers samples rejected by eligibility in both modes", {
  for (mode in c("on_demand", "rate")) {
    for (restriction in c("young", "kin", "partnered")) {
      model <- make_model(3, birth_rate = 20, birth_interval = 0)
      model$population$persons$sex <- c("f", "m", "f")
      model$population$persons$birth_time[3] <- -60
      if (restriction == "young") model$population$persons$birth_time[2] <- -14
      if (restriction == "kin") {
        model$population$persons$birth_time[2] <- -50
        model$population$persons$father_id[1] <- 2
      }
      if (restriction == "partnered") {
        model$population$unions <- tibble::tibble(union_id = 1, wife_id = 3,
          husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
      }
      matching <- model$processes$union_formation$matching
      if (mode == "rate") {
        model <- set_union_formation(model, "rate", matching,
          rates = data.frame(age_start = 0, age_end = Inf, rate = 0))
        model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
          age_start = 15, age_end = 50, rate = 20), 0)
      } else model <- set_union_formation(model, "on_demand", matching)
      expect_true(model$processes$union_formation$cancel_birth_if_no_match)
      result <- simulate(model, 0.5, seed = 52)
      expect_equal(result$simulation$births, 0)
      expect_equal(nrow(result$unions), nrow(model$population$unions))
      model$processes$union_formation$cancel_birth_if_no_match <- FALSE
      result <- simulate(model, 0.5, seed = 52)
      expect_gt(result$simulation$births, 0)
    }
  }
})
