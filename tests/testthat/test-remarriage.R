previous_rates <- function(none = 0, dissolved = 0, widowed = 0) {
  data.frame(age_start = 15, age_end = Inf, rate = c(none, dissolved, widowed),
    previous_union = c("none", "dissolved", "widowed"))
}

test_that("previous-union schedules require complete valid strata within groups", {
  rates <- previous_rates(.1, .2, .3)
  expect_s3_class(union_formation("rate", rates = rates)$rates, "tbl_df")
  expect_error(union_formation("rate", rates = rates[-1, ]), "all three")
  for (value in list(c("none", "dissolved", NA), c("none", "divorced", "widowed"),
                     factor(rates$previous_union), matrix(rates$previous_union))) {
    bad <- rates; bad$previous_union <- value
    expect_error(union_formation("rate", rates = bad), "previous_union")
  }
  grouped <- rbind(transform(rates, group = "A"), transform(rates, group = "B"))
  expect_error(union_formation("rate", rates = grouped[-1, ]), "all three")
  expect_error(union_formation("rate", rates = rbind(rates, rates[1, ])), "overlap")
  # Gaps and different age bands in different categories are valid.
  rates$age_start <- c(15, 20, 30)
  expect_silent(union_formation("rate", rates = rates))
  model <- make_rate_model()
  model$processes$union_formation$rates <- rates[-1, ]
  expect_error(simulate(model, 0), "all three")
})

test_that("latest recorded union selects the search schedule on restart", {
  model <- make_rate_model(4)
  model$population$unions <- tibble::tibble(union_id = c(2, 1), wife_id = 1,
    husband_id = c(3, 2), start_time = c(-4, -8), end_time = c(-2, -6),
    end_reason = c("dissolution", "death_husband"))
  model$population$persons$death_time[2] <- -6
  for (category in c("none", "dissolved", "widowed")) {
    rates <- previous_rates()
    rates$rate[rates$previous_union == category] <- 1e4
    m <- set_union_formation(model, "rate", rates = rates, kinship_generations = 0)
    result <- simulate(m, .1, 42)
    expect_equal(nrow(result$unions), if (category == "dissolved") 3 else 2)
    if (category == "dissolved") expect_lt(tail(result$unions$start_time, 1), .1)
    # Reversing the input rows leaves execution unchanged.
    m$population$unions <- m$population$unions[2:1, ]
    other <- simulate(m, .1, 42)
    expect_equal(result$persons, other$persons)
    expect_equal(result$unions[order(result$unions$union_id), ],
                 other$unions[order(other$unions$union_id), ])
  }
})

test_that("dissolution and widowhood enable their own remarriage schedules", {
  model <- make_rate_model(3)
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  model <- set_union_formation(model, "rate", rates = previous_rates(0, 100, 0))
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = 10))
  result <- simulate(model, 1, 12)
  expect_gt(nrow(result$unions), 2)
  expect_true(all(head(result$unions$end_reason, -1) == "dissolution"))
  expect_equal(previous_union_state_cpp(result$persons, result$unions)[1], 1L)
  restarted <- simulate(set_population(model, result), 1, 13)
  expect_gt(nrow(restarted$unions), nrow(result$unions))
  model$processes$union_dissolution <- NULL
  model <- set_mortality(model, data.frame(sex = c("f", "m", "m"),
    age_start = c(0, 0, 30), age_end = c(Inf, 30, Inf), rate = c(0, 0, 1e5)))
  model$population$persons$birth_time[2] <- -29.5
  model <- set_union_formation(model, "rate", rates = previous_rates(0, 0, 1e5))
  widowed <- simulate(model, .75, 42)
  expect_equal(nrow(widowed$unions), 2)
  expect_equal(widowed$unions$end_reason[1], "death_husband")
  expect_gte(widowed$unions$start_time[2], widowed$unions$end_time[1])
  expect_equal(widowed$unions$husband_id[2], 3)
})

test_that("each previous-category hazard agrees with analytical incidence", {
  n <- 1000L
  model <- make_rate_model(6 * n)
  model$population$persons$sex[seq_len(3 * n)] <- "f"
  previous <- rep(0:2, each = n)
  women <- seq.int(n + 1L, 3L * n)
  model$population$unions <- tibble::tibble(union_id = seq_len(2 * n),
    wife_id = women, husband_id = women + 3L * n,
    start_time = -5, end_time = -1,
    end_reason = rep(c("dissolution", "death_husband"), each = n))
  model$population$persons$death_time[seq.int(5L * n + 1L, 6L * n)] <- -1
  rates <- c(.1, .3, .7)
  model <- set_union_formation(model, "rate",
    rates = previous_rates(rates[1], rates[2], rates[3]), kinship_generations = 0)
  result <- simulate(model, 2, 502)
  starts <- subset(result$unions, start_time >= 0)
  for (category in 0:2) {
    p <- 1 - exp(-rates[category + 1L] * 2)
    actual <- sum(previous[starts$wife_id] == category)
    expect_lt(abs(actual - n * p), 6 * sqrt(n * p * (1 - p)))
  }
})

test_that("equal previous-union schedules preserve ungrouped execution exactly", {
  model <- make_rate_model(20, search_rate = 2)
  model$population$persons$sex[1:10] <- "f"
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .5))
  shared <- simulate(model, 3, 71)
  model <- set_union_formation(model, "rate", matching = model$processes$union_formation$matching,
    rates = previous_rates(2, 2, 2), cancel_birth_if_no_match = FALSE)
  expect_identical(simulate(model, 3, 71), shared)
})

test_that("previous categories contribute the explicit calibrated intensity", {
  women <- data.frame(born = -30.125, group = rep(0:1, each = 6), fertile = TRUE)
  previous <- c(rep(0L, 6), rep(1L, 3), rep(2L, 3))
  men <- data.frame(born = -30.125, group = rep(0:1, each = 8))
  rates <- previous_rates(.1, .2, .4)
  matching <- prepare_matching_model(create_matching_model(~ 1, numeric(),
    group_exogamy = c("A:B" = -2)), c("A", "B"))
  got <- inspect_weighted_calibration_cpp(women$born, women$group, women$fertile,
    men$born, men$group, rates, matching, 2L, 0, 15, .25, FALSE, previous)
  # A: six * .1, B: three * .2 + three * .4; both groups have six women.
  expect_equal(got$target, c(.25, .75), tolerance = 1e-12)
  expect_equal(got$flow, got$target, tolerance = 1e-9)
  effective <- data.frame(age_start = 15, age_end = Inf, rate = c(.1, .3), .hammel_group = 0:1)
  expected <- inspect_weighted_calibration_cpp(women$born, women$group, women$fertile,
    men$born, men$group, effective, matching, 2L, 0, 15, .25, FALSE)
  expect_equal(got, expected, tolerance = 1e-9)
  # Zero previous-specific risk in B removes that group's female intensity.
  rates$rate <- c(.1, 0, 0)
  zero <- inspect_weighted_calibration_cpp(women$born, women$group, women$fertile,
    men$born, men$group, rates, matching, 2L, 0, 15, .25, FALSE, previous)
  expect_equal(zero$target, c(1, 0))
})

test_that("grouped schedules and multipliers preserve the explicit hazards", {
  model <- make_rate_model(80, search_rate = 2)
  model$population$persons$sex[1:40] <- "f"
  model$population$persons$group <- rep(c("A", "B"), 40)
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .8))
  rates <- rbind(transform(previous_rates(.8, .4, .2), group = "A"),
    transform(previous_rates(.5, .3, .1), group = "B"))
  matching <- create_matching_model(~ 1, numeric(), group_exogamy = c("A:B" = -1))
  scaled <- set_union_formation(model, "rate", matching, rates = rates,
    group_prop_hazard = c(B = 2), diagnostics = TRUE, kinship_generations = 0)
  explicit <- rates
  explicit$rate[explicit$group == "B"] <- 2 * explicit$rate[explicit$group == "B"]
  reference <- set_union_formation(model, "rate", matching, rates = explicit,
    diagnostics = TRUE, kinship_generations = 0)
  a <- simulate(scaled, 5, 73)
  expect_identical(a, simulate(reference, 5, 73))
  expect_true(any(a$unions$end_reason == "dissolution", na.rm = TRUE))
  counts <- calculate_marriage_counts(a, c(0, 5), c(0, Inf), TRUE, by_previous_union = TRUE)
  expect_gt(sum(counts$unions[counts$previous_union == "dissolved"]), 0)
  # Shared hazards versus three equal category hazards give the same execution.
  plain <- set_union_formation(model, "rate", matching,
    rates = data.frame(age_start = 15, age_end = Inf, rate = 2),
    diagnostics = TRUE, kinship_generations = 0)
  stratified <- set_union_formation(model, "rate", matching,
    rates = previous_rates(2, 2, 2), diagnostics = TRUE, kinship_generations = 0)
  expect_equal(simulate(plain, 5, 74), simulate(stratified, 5, 74), tolerance = 1e-10)
  # History overlap is ambiguous only for the new stratification contract.
  bad <- stratified
  bad$population$unions <- data.frame(union_id = 1:2, wife_id = 1,
    husband_id = c(41, 43), start_time = c(-2, -1.5), end_time = c(-1, -.5),
    end_reason = "dissolution")
  expect_error(simulate(bad, 1, 42), "nonoverlapping")
})
