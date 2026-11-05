duration_search_rates <- function(rate = .2) {
  data.frame(previous_union = c("none", "separated", "widowed"),
    age_start = c(15, NA, NA), age_end = c(Inf, NA, NA),
    dur_start = c(NA, 0, 0), dur_end = c(NA, Inf, Inf), rate = c(0, rate, rate))
}

ended_union_model <- function(n = 1L, ended = -2, reason = "dissolution") {
  model <- make_rate_model(2L * n)
  model$population$persons$sex[seq_len(n)] <- "f"
  model$population$unions <- tibble::tibble(union_id = seq_len(n),
    wife_id = seq_len(n), husband_id = n + seq_len(n),
    start_time = ended - 5, end_time = ended, end_reason = reason)
  set_union_formation(model, "rate", rates = duration_search_rates(), kinship_generations = 0L)
}

test_that("single search tables permit one clock per previous category", {
  rates <- duration_search_rates()
  expect_s3_class(union_formation("rate", rates = rates)$rates, "tbl_df")
  mixed <- rates
  mixed$dur_start[3] <- mixed$dur_end[3] <- NA
  mixed$age_start[3] <- 30; mixed$age_end[3] <- Inf
  expect_silent(union_formation("rate", rates = mixed))
  expect_error(union_formation("rate", rates = rates[, setdiff(names(rates), "dur_end")]), "both")
  for (column in c("age_start", "age_end")) {
    bad <- rates; bad[[column]][2] <- 20
    expect_error(union_formation("rate", rates = bad), "exactly one")
  }
  bad <- rates; bad$age_start[2] <- 20; bad$age_end[2] <- 50
  expect_error(union_formation("rate", rates = bad), "exactly one")
  bad <- rates; bad$dur_start[2] <- bad$dur_end[2] <- NA
  expect_error(union_formation("rate", rates = bad), "exactly one")
  bad <- rates; bad$age_start[1] <- bad$age_end[1] <- NA
  bad$dur_start[1] <- 0; bad$dur_end[1] <- Inf
  expect_error(union_formation("rate", rates = bad), "none.*age")
  # A category cannot switch clocks across either bands or groups.
  extra <- mixed[2, ]; extra$dur_start <- extra$dur_end <- NA
  extra$age_start <- 50; extra$age_end <- Inf
  expect_error(union_formation("rate", rates = rbind(mixed, extra)), "same.*clock")
  group_a <- transform(rates, group = "A")
  group_b <- transform(mixed, group = "B")
  expect_error(union_formation("rate", rates = rbind(group_a, group_b)), "same.*clock")
  bad <- rates; bad$previous_union <- NULL
  expect_error(union_formation("rate", rates = bad), "require.*previous_union")
  for (value in c(-1, Inf, NaN)) {
    bad <- rates; bad$dur_start[2] <- value
    expect_error(union_formation("rate", rates = bad))
  }
  bad <- rates; bad$dur_end[2] <- 0
  expect_error(union_formation("rate", rates = bad), "greater")
  expect_error(union_formation("rate", rates = rbind(rates, rates[2, ])), "overlap")
  all_age <- data.frame(age_start = 15, age_end = Inf, rate = .1)
  expect_silent(union_formation("rate", rates = transform(all_age, dur_start = NA, dur_end = NA)))
  expect_identical(simulate(set_union_formation(make_rate_model(), "rate", rates = all_age), 2, 42),
    simulate(set_union_formation(make_rate_model(), "rate", rates =
      transform(all_age, dur_start = NA, dur_end = NA)), 2, 42))
})

test_that("duration search clocks resume attained duration rather than age", {
  model <- ended_union_model(1, ended = -2)
  rates <- duration_search_rates(1e5)
  rates$dur_start[2:3] <- 5
  model <- set_union_formation(model, "rate", rates = rates)
  early <- simulate(model, 2, 42)
  expect_equal(nrow(early$unions), 1)
  late <- simulate(set_population(model, early), 2, 43)
  expect_equal(nrow(late$unions), 2)
  expect_gte(late$unions$start_time[2], 3)
  expect_lt(late$unions$start_time[2], 4)
  # Changing age leaves a duration-specific waiting time unchanged.
  older <- model; older$population$persons$birth_time <- -80
  expect_identical(simulate(model, 4, 43)$unions, simulate(older, 4, 43)$unions)
  # Input history order does not change either the previous category or origin.
  more <- model
  more$population$unions <- rbind(more$population$unions,
    transform(more$population$unions, union_id = 2, start_time = -12,
      end_time = -10, end_reason = "dissolution"))
  a <- simulate(more, 4, 43)
  more$population$unions <- more$population$unions[2:1, ]
  b <- simulate(more, 4, 43)
  expect_equal(tail(a$unions, 1), tail(b$unions, 1))
})

test_that("duration hazards agree with analytical cumulative hazard across gaps", {
  n <- 2500L
  model <- ended_union_model(n, ended = -2)
  rates <- duration_search_rates()
  dissolved <- rates[2, ]
  dissolved <- dissolved[rep(1L, 3), ]
  dissolved$dur_start <- c(0, 2, 4); dissolved$dur_end <- c(1, 3, Inf)
  dissolved$rate <- c(.1, .2, .4)
  model <- set_union_formation(model, "rate", rates = rbind(rates[c(1, 3), ], dissolved),
    kinship_generations = 0L)
  result <- simulate(model, 3, 813)
  # Starting at duration 2: .2 over [2,3), zero over [3,4), .4 over [4,5).
  probability <- 1 - exp(-.6)
  starts <- subset(result$unions, start_time >= 0)
  expect_lt(abs(nrow(starts) - n * probability), 6 * sqrt(n * probability * (1 - probability)))
  expect_false(any(starts$start_time >= 1 & starts$start_time < 2))
})

test_that("each new ending resets duration and widowhood starts at partner death", {
  model <- ended_union_model(1, ended = -2)
  rates <- duration_search_rates(1e5)
  rates$dur_start[2:3] <- 1
  model <- set_union_formation(model, "rate", rates = rates)
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = 1e5))
  result <- simulate(model, 3.5, 71)
  # The input episode ended two years ago, so first remarriage has no delay.
  # Every subsequent ending imposes a fresh year before the next search.
  expect_gte(nrow(result$unions), 4)
  starts <- result$unions$start_time
  ends <- result$unions$end_time
  expect_true(all(starts[3:length(starts)] - ends[2:(length(ends)-1)] >= 1))
  widow <- make_rate_model(3)
  widow$population$persons$birth_time[2] <- -29.5
  widow$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -2, end_time = NA_real_, end_reason = NA_character_)
  widow <- set_union_formation(widow, "rate", rates = rates)
  widow <- set_mortality(widow, data.frame(sex = c("f", "m", "m"),
    age_start = c(0, 0, 30), age_end = c(Inf, 30, Inf), rate = c(0, 0, 1e5)))
  result <- simulate(widow, 2, 71)
  expect_equal(nrow(result$unions), 2)
  expect_gte(result$unions$start_time[2] - result$unions$end_time[1], 1)
  expect_lt(result$unions$start_time[2], 2)
  info <- previous_union_info_cpp(result$persons, result$unions)
  expect_equal(info$end_time[1], result$unions$end_time[1])
})

test_that("grouped duration calibration matches the equivalent age hazards", {
  model <- ended_union_model(30, ended = -2)
  model$population$persons$group <- rep(c("A", "B"), 30)
  duration <- duration_search_rates(.4)
  duration$dur_start[2:3] <- 3
  age <- duration
  age$age_start[2:3] <- 21
  age$age_end[2:3] <- Inf
  age$dur_start <- age$dur_end <- NULL
  matching <- create_matching_model(~ 1, numeric(), group_exogamy = c("A:B" = -2))
  configure <- function(rates) set_union_formation(model, "rate", matching,
    rates = stratify_rates(rates, list(group = c(A = 1, B = 1.3))), diagnostics = TRUE,
    kinship_generations = 0L)
  set.seed(87); rng <- .Random.seed
  a <- simulate(configure(duration), 4, 71)
  b <- simulate(configure(age), 4, 71)
  expect_identical(a, b)
  expect_identical(.Random.seed, rng)
  expect_true(all(a$unions$start_time[a$unions$start_time >= 0] >= 1))
  expect_true(all(is.finite(a$unions$.log_weight_B[-seq_len(30)])))
  expect_identical(simulate(configure(duration), 4, 71), a)
})

test_that("duration searches calibrate on their ending origins", {
  rates <- duration_search_rates()
  women <- c(-20, -50)
  prepared <- prepare_matching_model(create_matching_model(~ 1, numeric(),
    group_exogamy = c("A:B" = -1)), c("A", "B"))
  # Durations 1 and 10, both in a constant .2 band: target splits evenly.
  result <- inspect_weighted_calibration_cpp(women, 0:1, c(TRUE, TRUE),
    c(-30, -30), 0:1, rates, prepared, 2L, 0, 15, .25, FALSE, c(1L, 2L), c(-1, -10))
  expect_equal(result$target, c(.5, .5), tolerance = 1e-12)
  expect_error(inspect_weighted_calibration_cpp(women, 0:1, c(TRUE, TRUE),
    c(-30, -30), 0:1, rates, prepared, 2L, 0, 15, .25, FALSE, c(1L, 2L)), "ending time")
})

test_that("unpartnered births do not reset the remarriage duration clock", {
  model <- ended_union_model(1, ended = 0)
  rates <- duration_search_rates(1e5)
  rates$dur_start[2:3] <- 1
  model <- set_union_formation(model, "rate", rates = rates, kinship_generations = 0L)
  model <- set_fertility(model, data.frame(age_start = 15, age_end = 50,
    union_type = c("partnered", "unpartnered"), rate = c(0, 20)), birth_interval = 0)
  result <- simulate(model, 1.1, 317)
  children <- subset(result$persons, birth_time >= 0)
  expect_gt(nrow(children), 0)
  expect_equal(nrow(result$unions), 2)
  expect_gte(result$unions$start_time[2], 1)
  expect_lt(result$unions$start_time[2], 1.01)
  expect_true(all(children$birth_time < result$unions$start_time[2]))
})

test_that("failed searches retain duration until an eligible man becomes available", {
  model <- make_rate_model(3)
  model$population$persons$death_time[2] <- 0
  model$population$persons$birth_time[3] <- -14.5
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -5, end_time = 0, end_reason = "death_husband")
  rates <- duration_search_rates(100)
  rates$dur_start[2:3] <- .1
  rates$dur_end[2:3] <- .75
  model <- set_union_formation(model, "rate", rates = rates, kinship_generations = 0L)
  result <- simulate(model, 1, 617)
  expect_equal(nrow(result$unions), 2)
  expect_equal(result$unions$husband_id[2], 3)
  expect_gte(result$unions$start_time[2], .5)
  expect_lt(result$unions$start_time[2], .75)
})
