transition_model <- function(n = 1, sex = "m", age = 20) {
  model <- make_model(n, sex = sex, age = age)
  model$population$persons$group <- rep("A", n)
  model
}
transition_rates <- function(rate = 1, origin = "A", dest = "B") {
  data.frame(group_origin = origin, group_dest = dest,
    age_start = 0, age_end = Inf, rate = rate)
}

test_that("transition specifications validate common bands and optional strata", {
  r <- transition_rates()
  attr(r, "source") <- list(example = "synthetic")
  p <- group_transition(r)
  expect_s3_class(p, "hammel_group_transition")
  expect_s3_class(p$rates, "tbl_df")
  expect_identical(p$source, attr(r, "source"))
  expect_identical(set_group_transition(hammel_model(), r)$processes$group_transition, p)
  for (bad in list(transform(r, rate = -1), transform(r, rate = Inf),
                   transform(r, age_end = 0), transform(r, group_dest = "A"),
                   transform(r, group_dest = NA_character_), transform(r, sex = "x"),
                   transform(r, union_type = "single"),
                   transform(r, dur_start = 0, dur_end = Inf), rbind(r, r)))
    expect_error(group_transition(bad))
  expect_s3_class(group_transition(transform(r, sex = "f", union_type = "cohab")), "hammel_process")
  expect_error(simulate(set_group_transition(make_model(1), r), 1), "group identities")
})

test_that("competing outgoing hazards recover analytical destinations and survival", {
  model <- transition_model(12000)
  model <- set_group_transition(model, rbind(transition_rates(.2), transition_rates(.3, dest = "C")))
  set.seed(982)
  before <- .Random.seed
  out <- simulate(model, 2, seed = 124)
  expect_identical(.Random.seed, before)
  counts <- table(factor(out$persons$group, levels = c("A", "B", "C"))) / 12000
  expected <- c(exp(-1), .4 * (1 - exp(-1)), .6 * (1 - exp(-1)))
  expect_lt(max(abs(as.numeric(counts) - expected) / sqrt(expected * (1 - expected) / 12000)), 5)
  expect_equal(nrow(out$transitions), sum(out$persons$group != "A"))
  expect_identical(out$simulation$transitions, nrow(out$transitions))
  expect_equal(out$ancestry$group_A, rep(1, 12000))
  expect_identical(simulate(model, 2, seed = 124), out)
})

test_that("grouped outputs always have typed histories, including zero and empty runs", {
  for (n in c(0, 2)) {
    model <- transition_model(n)
    out <- simulate(model, 0, seed = 1)
    expect_s3_class(out$transitions, "tbl_df")
    expect_identical(names(out$transitions), c("person_id", "time", "group_origin", "group_dest"))
    expect_equal(nrow(out$transitions), 0)
    expect_type(out$transitions$group_origin, "character")
    expect_s3_class(simulate(model, 1, seed = 1)$transitions, "tbl_df")
  }
  model <- transition_model(2)
  model$population$persons$group <- c(10L, 10L)
  out <- simulate(set_group_transition(model, transition_rates(1e5, 10L, 20L)), .1, seed = 1)
  expect_identical(out$persons$group, c(20L, 20L))
  expect_type(out$transitions$group_origin, "integer")
})

test_that("age support, sex and union type select actual applicable schedules", {
  model <- transition_model(2)
  model$population$persons$sex <- c("m", "f")
  r <- transform(transition_rates(1e5), age_start = 21, sex = "f", union_type = "unpartnered")
  out <- simulate(set_group_transition(model, r), 1.5, seed = 1)
  expect_equal(out$transitions$person_id, 2)
  expect_gte(out$transitions$time, 1)
  expect_identical(out$persons$group, c("A", "B"))
  expect_equal(nrow(simulate(set_group_transition(model, transform(r, union_type = "partnered")), 2, seed = 1)$transitions), 0)
})

test_that("duration resets at moves, survives restarts and falls back to age", {
  r <- rbind(transition_rates(1e5), transition_rates(1e5, "B", "A"))
  names(r)[names(r) == "age_start"] <- "dur_start"
  names(r)[names(r) == "age_end"] <- "dur_end"
  r$dur_start <- 2
  model <- set_group_transition(transition_model(), r)
  out <- simulate(model, 1, seed = 2)
  expect_equal(nrow(out$transitions), 1)
  expect_identical(out$persons$group, "B")
  second <- simulate(set_population(model, out), 2, seed = 3)
  expect_equal(nrow(second$transitions), 2)
  expect_gte(diff(second$transitions$time), 2)
  expect_identical(second$transitions[1, ], out$transitions)
  expect_identical(second$persons$group, "A")
  expect_error(hammel_population(second$persons, end_time = 3,
    transitions = transform(second$transitions, group_dest = "C")), "chain")
})

test_that("group changes select new mortality and update only the individual", {
  model <- transition_model(2)
  model$population$persons$sex <- c("f", "m")
  model$population$unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  death <- rbind(transform(model$processes$mortality$rates, group = "A"),
    transform(model$processes$mortality$rates, group = "B", rate = 1e5))
  model <- set_mortality(model, death)
  out <- simulate(set_group_transition(model, transform(transition_rates(1e5), sex = "f")), .1, seed = 1)
  expect_equal(out$transitions$person_id, 1)
  expect_identical(out$persons$group, c("B", "A"))
  expect_gt(out$persons$death_time[1], out$transitions$time)
  expect_true(is.na(out$persons$death_time[2]))
  expect_identical(out$ancestry$group_A, c(1, 1))
})

test_that("a moving man invalidates pair dissolution and union changes refresh his risks", {
  model <- transition_model(2)
  model$population$persons$sex <- c("f", "m")
  model$population$unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  dissolution <- data.frame(group_pair = c("A:A", "A:B", "B:B"),
    dur_start = 0, dur_end = Inf, rate = c(0, 1e5, 0))
  model <- set_union_dissolution(model, dissolution)
  transitions <- rbind(transform(transition_rates(1e5), sex = "m", union_type = "partnered"),
    transform(transition_rates(1e5, "B", "C"), sex = "m", union_type = "unpartnered"))
  # C must also have all dissolution pair schedules.
  dissolution <- rbind(dissolution, data.frame(group_pair = c("A:C", "B:C", "C:C"),
    dur_start = 0, dur_end = Inf, rate = 0))
  model <- set_union_dissolution(model, dissolution)
  out <- simulate(set_group_transition(model, transitions), .1, seed = 15)
  expect_equal(out$transitions$group_dest, c("B", "C"))
  expect_equal(out$unions$end_reason, "dissolution")
  expect_gt(out$unions$end_time, out$transitions$time[1])
  expect_lt(out$unions$end_time, out$transitions$time[2])
  expect_identical(out$persons$group, c("A", "C"))
})

# Known histories isolate observation rules from stochastic simulation.
transition_observation_population <- function() {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"),
    birth_time = c(-20, -20, 1), death_time = c(NA, 4, NA),
    mother_id = c(NA, NA, 1), father_id = c(NA, NA, 2), group = c("B", "B", "A"))
  unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = 0, end_time = 3, end_reason = "dissolution")
  tx <- data.frame(person_id = 1:2, time = c(2, 1), group_origin = "A", group_dest = "B")
  hammel_population(persons, unions, end_time = 5, transitions = tx)
}

test_that("observations split exposure and attribute events using group history", {
  pop <- transition_observation_population()
  snapshot <- calculate_population_size(pop, c(0, 1, 2, 5), by_group = TRUE)
  expect_equal(snapshot$population[snapshot$group == "B"], c(0, 1, 2, 1))
  mortality <- calculate_mortality_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE)
  expect_equal(mortality$person_years[mortality$sex == "m" & mortality$group == "A"], 1)
  expect_equal(mortality$person_years[mortality$sex == "m" & mortality$group == "B"], 3)
  expect_equal(mortality$deaths[mortality$sex == "m" & mortality$group == "B"], 1)
  fertility <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE)
  expect_equal(fertility$births[fertility$group == "A"], 1)
  expect_equal(fertility$births[fertility$group == "B"], 0)
  expect_equal(fertility$person_years[fertility$group == "A"], 6)
  typed <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE, by_union_type = TRUE)
  parity <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE, parity_breaks = c(0, 1, Inf))
  expect_equal(sum(typed$births), 1)
  expect_equal(sum(typed$person_years), sum(fertility$person_years))
  expect_equal(sum(parity$person_years), sum(fertility$person_years))
  expect_equal(sum(parity$births[parity$group == "A"]), 1)
  unions <- calculate_union_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE)
  expect_equal(sum(unions$unions[unions$group == "A"]), 2)
  dissolution <- calculate_dissolution_rates(pop, c(0, 5), c(0, Inf), by_group = TRUE)
  expect_equal(dissolution$union_years[dissolution$group_pair == "A:A"], 1)
  expect_equal(dissolution$union_years[dissolution$group_pair == "A:B"], 1)
  expect_equal(dissolution$union_years[dissolution$group_pair == "B:B"], 1)
  expect_equal(dissolution$dissolutions[dissolution$group_pair == "B:B"], 1)
})

test_that("moves refresh available male pools and active fertility/calibration cells", {
  old <- options(hammel.check_cells = TRUE)
  on.exit(options(old))
  model <- make_rate_model(2, search_rate = 0)
  model$population$persons$group <- c("B", "A")
  model <- set_union_formation(model, "rate",
    matching = create_age_gap_matching_model(0, 5, prevent_exogamy = TRUE),
    rates = data.frame(age_start = 21, age_end = Inf, rate = 1e5))
  out <- simulate(set_group_transition(model, transform(transition_rates(1e5), sex = "m")), 1.1, seed = 33)
  expect_equal(nrow(out$unions), 1)
  expect_equal(out$transitions$person_id, 2)
  expect_gt(out$unions$start_time, out$transitions$time)

  model <- transition_model(2, sex = "f")
  rates <- rbind(data.frame(group = "A", age_start = 15, age_end = 50, rate = 0),
                 data.frame(group = "B", age_start = 15, age_end = 50, rate = 30))
  model <- set_fertility(model, rates, birth_interval = 0, inheritance = "mother")
  out <- simulate(set_group_transition(model, transform(transition_rates(1e5), age_end = 21)), .3, seed = 51)
  expect_gt(out$simulation$births, 0)
  expect_equal(nrow(out$transitions), 2)
  expect_true(all(out$persons$group == "B"))
  expect_true(all(out$persons$birth_time[out$persons$birth_time > 0] > min(out$transitions$time)))
  f <- calculate_fertility_rates(out, c(0, .3), c(0, Inf), by_group = TRUE)
  expect_equal(f$births[f$group == "A"], 0)
  expect_equal(sum(f$births), out$simulation$births)
})

test_that("cohabitation conversion refreshes both partners' transition strata", {
  conversion <- data.frame(group = c("A", "B", "C"), dur_start = 0, dur_end = Inf,
    rate = c(0, 1e5, 0))
  model <- make_rate_model(2, search_rate = 0)
  model$population$persons$group <- c("A", "A")
  model$population$unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_, union_type = "cohab")
  model <- set_union_formation(model, "rate", rates = data.frame(age_start = 0, age_end = Inf, rate = 0),
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1), union_convert_rate = conversion)
  rates <- rbind(transform(transition_rates(1e5), sex = "f", union_type = "cohab"),
    transform(transition_rates(1e5, "A", "C"), sex = "m", union_type = "marriage"))
  out <- simulate(set_group_transition(model, rates), .1, seed = 1)
  expect_identical(out$transitions$group_dest, c("B", "C"))
  expect_identical(out$unions$end_reason[1], "conversion")
  expect_gt(out$unions$end_time[1], out$transitions$time[1])
  expect_lt(out$unions$end_time[1], out$transitions$time[2])
  counts <- calculate_cohab_convert_rate(out, c(0, .1), dur_breaks = c(0, Inf), by_group = TRUE)
  expect_equal(counts$conversions[counts$group == "B"], 1)
  expect_equal(counts$union_years[counts$group == "A"], out$transitions$time[1])
  expect_equal(sum(counts$union_years), out$unions$end_time[1])
})

test_that("history-only restarts and absent transition processes preserve provenance", {
  model <- transition_model()
  out <- simulate(set_group_transition(model, transition_rates(1e5)), .1, seed = 9)
  resumed <- simulate(set_population(model, out), .1, seed = 10)
  expect_identical(resumed$transitions, out$transitions)
  expect_identical(resumed$ancestry, out$ancestry)
  expect_identical(resumed$persons$group, "B")
  expect_null(model$population$transitions)
  expect_identical(model$population$persons$group, "A")
  zero <- simulate(set_group_transition(model, transition_rates(0)), 1, seed = 2)
  expect_equal(nrow(zero$transitions), 0)
  bad <- set_group_transition(model, transition_rates())
  bad$processes$group_transition$rates$rate <- -1
  expect_error(simulate(bad, 0), "nonnegative")
})
