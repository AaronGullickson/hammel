transition_model <- function(n = 1, sex = "m", age = 20) {
  model <- make_model(n, sex = sex, age = age)
  model$population$persons$group <- rep("A", n)
  model$population$groups <- list(group = group(ancestry = TRUE))
  model
}
transition_rates <- function(rate = 1, origin = "A", dest = "B") {
  data.frame(origin = origin, destination = dest,
    age_start = 0, age_end = Inf, rate = rate)
}

test_that("transition specifications validate common bands and optional strata", {
  r <- transition_rates()
  attr(r, "source") <- list(example = "synthetic")
  p <- group_transition(r)
  expect_s3_class(p, "hammel_group_transition")
  expect_s3_class(p$rates, "tbl_df")
  expect_identical(p$source, attr(r, "source"))
  # Without a declaration the process waits to take the only dimension.
  expect_identical(set_group_transition(hammel_model(), r)$processes$group_transition, list(.default = p))
  for (bad in list(transform(r, rate = -1), transform(r, rate = Inf),
                   transform(r, age_end = 0), transform(r, destination = "A"),
                   transform(r, destination = NA_character_), transform(r, sex = "x"),
                   transform(r, union_type = "single"),
                   transform(r, dur_start = 0, dur_end = Inf), rbind(r, r)))
    expect_error(group_transition(bad))
  expect_s3_class(group_transition(transform(r, sex = "f", union_type = "cohab")), "hammel_process")
  expect_error(simulate(set_group_transition(make_model(1), r), 1), "declared group")
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
  expect_equal(out$ancestry$group$A, rep(1, 12000))
  expect_identical(simulate(model, 2, seed = 124), out)
})

test_that("grouped outputs always have typed histories, including zero and empty runs", {
  for (n in c(0, 2)) {
    model <- transition_model(n)
    out <- simulate(model, 0, seed = 1)
    expect_s3_class(out$transitions, "tbl_df")
    expect_identical(names(out$transitions), c("person_id", "time", "group", "origin", "destination"))
    expect_equal(nrow(out$transitions), 0)
    expect_type(out$transitions$origin, "character")
    expect_s3_class(simulate(model, 1, seed = 1)$transitions, "tbl_df")
  }
  model <- transition_model(2)
  model$population$persons$group <- c(10L, 10L)
  model$population$groups <- list(group = group())
  out <- simulate(set_group_transition(model, transition_rates(1e5, 10L, 20L)), .1, seed = 1)
  expect_identical(out$persons$group, c(20L, 20L))
  # History labels are text in every dimension; persons keep their type.
  expect_identical(out$transitions$origin, c("10", "10"))
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
    transitions = transform(second$transitions, destination = "C"), groups = if ("group" %in% names(second$persons)) "group"), "chain")
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
  expect_identical(out$ancestry$group$A, c(1, 1))
})

test_that("a moving man invalidates pair dissolution and union changes refresh his risks", {
  model <- transition_model(2)
  model$population$persons$sex <- c("f", "m")
  model$population$unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  dissolution <- directional_pairs(data.frame(group_pair = c("A:A", "A:B", "B:B"),
    dur_start = 0, dur_end = Inf, rate = c(0, 1e5, 0)))
  model <- set_union_dissolution(model, dissolution)
  transitions <- rbind(transform(transition_rates(1e5), sex = "m", union_type = "partnered"),
    transform(transition_rates(1e5, "B", "C"), sex = "m", union_type = "unpartnered"))
  # C must also have all dissolution pair schedules.
  dissolution <- rbind(dissolution, directional_pairs(data.frame(group_pair = c("A:C", "B:C", "C:C"),
    dur_start = 0, dur_end = Inf, rate = 0)))
  model <- set_union_dissolution(model, dissolution)
  out <- simulate(set_group_transition(model, transitions), .1, seed = 15)
  expect_equal(out$transitions$destination, c("B", "C"))
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
  tx <- data.frame(person_id = 1:2, time = c(2, 1), origin = "A", destination = "B")
  hammel_population(persons, unions, end_time = 5, transitions = tx, groups = if ("group" %in% names(persons)) "group")
}

test_that("observations split exposure and attribute events using group history", {
  pop <- transition_observation_population()
  snapshot <- calculate_population_size(pop, c(0, 1, 2, 5), by = "group")
  expect_equal(snapshot$population[snapshot$group == "B"], c(0, 1, 2, 1))
  mortality <- calculate_mortality_rates(pop, c(0, 5), c(0, Inf), by = "group")
  expect_equal(mortality$person_years[mortality$sex == "m" & mortality$group == "A"], 1)
  expect_equal(mortality$person_years[mortality$sex == "m" & mortality$group == "B"], 3)
  expect_equal(mortality$deaths[mortality$sex == "m" & mortality$group == "B"], 1)
  fertility <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by = "group")
  expect_equal(fertility$births[fertility$group == "A"], 1)
  expect_equal(fertility$births[fertility$group == "B"], 0)
  expect_equal(fertility$person_years[fertility$group == "A"], 6)
  typed <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by = c("group", "union_type"))
  parity <- calculate_fertility_rates(pop, c(0, 5), c(0, Inf), by = "group", parity_breaks = c(0, 1, Inf))
  expect_equal(sum(typed$births), 1)
  expect_equal(sum(typed$person_years), sum(fertility$person_years))
  expect_equal(sum(parity$person_years), sum(fertility$person_years))
  expect_equal(sum(parity$births[parity$group == "A"]), 1)
  unions <- calculate_union_rates(pop, c(0, 5), c(0, Inf), by = "group")
  expect_equal(sum(unions$unions[unions$group == "A"]), 2)
  dissolution <- pair_rates(calculate_dissolution_rates(pop, c(0, 5), c(0, Inf), by = c("group_w", "group_h")))
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
  model$population$groups <- list(group = group())
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
  model <- set_inheritance(set_fertility(model, rates, birth_interval = 0), "mother")
  out <- simulate(set_group_transition(model, transform(transition_rates(1e5), age_end = 21)), .3, seed = 51)
  expect_gt(out$simulation$births, 0)
  expect_equal(nrow(out$transitions), 2)
  expect_true(all(out$persons$group == "B"))
  expect_true(all(out$persons$birth_time[out$persons$birth_time > 0] > min(out$transitions$time)))
  f <- calculate_fertility_rates(out, c(0, .3), c(0, Inf), by = "group")
  expect_equal(f$births[f$group == "A"], 0)
  expect_equal(sum(f$births), out$simulation$births)
})

test_that("cohabitation conversion refreshes both partners' transition strata", {
  conversion <- data.frame(group_w = c("A", "B", "C"), dur_start = 0, dur_end = Inf,
    rate = c(0, 1e5, 0))
  model <- make_rate_model(2, search_rate = 0)
  model$population$persons$group <- c("A", "A")
  model$population$groups <- list(group = group())
  model$population$unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_, union_type = "cohab")
  model <- set_union_formation(model, "rate", rates = data.frame(age_start = 0, age_end = Inf, rate = 0),
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1), union_convert_rate = conversion)
  rates <- rbind(transform(transition_rates(1e5), sex = "f", union_type = "cohab"),
    transform(transition_rates(1e5, "A", "C"), sex = "m", union_type = "marriage"))
  out <- simulate(set_group_transition(model, rates), .1, seed = 1)
  expect_identical(out$transitions$destination, c("B", "C"))
  expect_identical(out$unions$end_reason[1], "conversion")
  expect_gt(out$unions$end_time[1], out$transitions$time[1])
  expect_lt(out$unions$end_time[1], out$transitions$time[2])
  counts <- calculate_cohab_convert_rate(out, c(0, .1), dur_breaks = c(0, Inf), by = "group_w")
  expect_equal(counts$conversions[counts$group_w == "B"], 1)
  expect_equal(counts$union_years[counts$group_w == "A"], out$transitions$time[1])
  expect_equal(sum(counts$union_years), out$unions$end_time[1])
})

test_that("history-only restarts and absent transition processes preserve provenance", {
  model <- transition_model()
  out <- simulate(set_group_transition(model, transition_rates(1e5)), .1, seed = 9)
  resumed <- simulate(set_population(model, out), .1, seed = 10)
  expect_identical(resumed$transitions, out$transitions)
  expect_identical(resumed$ancestry$group, out$ancestry$group)
  expect_identical(resumed$persons$group, "B")
  expect_null(model$population$transitions)
  expect_identical(model$population$persons$group, "A")
  zero <- simulate(set_group_transition(model, transition_rates(0)), 1, seed = 2)
  expect_equal(nrow(zero$transitions), 0)
  bad <- set_group_transition(model, transition_rates())
  bad$processes$group_transition$group$rates$rate <- -1
  expect_error(simulate(bad, 0), "nonnegative")
})

test_that("snapshot event accounting handles repeated moves, deaths and arbitrary time order", {
  persons <- data.frame(person_id = c(30, 10, 20), sex = c("f", "m", "f"),
    birth_time = c(0, -1, 3), death_time = c(5, NA, NA),
    mother_id = NA_real_, father_id = NA_real_, group = c("A", "B", "C"))
  tx <- data.frame(person_id = c(30, 10, 30, 10), time = c(4, 2, 1, 4),
    origin = c("B", "C", "A", "A"), destination = c("A", "A", "B", "B"))
  pop <- hammel_population(persons, end_time = 6, transitions = tx, groups = if ("group" %in% names(persons)) "group")
  times <- c(5, -2, 0, 4, 1, 3, 4, 2, 6)
  out <- calculate_population_size(pop, times, by = "group")
  # Independently reconstruct identity from each person's ordered history.
  expected <- lapply(times, function(time) {
    living <- which(persons$birth_time <= time & (is.na(persons$death_time) | persons$death_time > time))
    labels <- vapply(living, function(i) {
      h <- tx[tx$person_id == persons$person_id[i], ]
      if (!nrow(h)) return(persons$group[i])
      h <- h[order(h$time), ]
      before <- which(h$time <= time)
      if (length(before)) tail(h$destination[before], 1) else h$origin[1]
    }, character(1))
    tabulate(match(labels, c("A", "B", "C")), nbins = 3)
  })
  # Rows vary by time fastest, then group.
  expect_equal(out$population, as.integer(t(matrix(unlist(expected), nrow = 3))))
  expect_identical(out$time, rep(times, 3))
  expect_equal(as.vector(tapply(out$population, rep(seq_along(times), 3), sum)),
    calculate_population_size(pop, times)$population)
  expect_equal(nrow(calculate_population_size(pop, numeric(), by = "group")), 0)
})

test_that("edited transition specifications reject incompatible units and unsupported processes", {
  model <- set_group_transition(transition_model(), transition_rates())
  model$processes$group_transition$group$rate_unit <- "annual_probability"
  expect_error(simulate(model, 0), "specification or units")
  model <- make_model(1)
  model$processes$unimplemented <- structure(list(name = "unimplemented"), class = "hammel_process")
  expect_error(simulate(model, 0), "Unsupported demographic process")
})

test_that("transition histories reject invalid references, life dates and missing identities", {
  original <- transition_observation_population()
  for (tx in list(
    transform(original$transitions, person_id = c(999, 2)),
    transform(original$transitions, time = c(-21, 1)),
    transform(original$transitions, time = c(2, 4)),
    transform(original$transitions, time = c(6, 1)),
    transform(original$transitions, origin = destination)
  )) {
    pop <- original
    pop$transitions <- tx
    expect_error(validate_population(pop), "invalid people, dates or self-transitions")
  }
  original$persons$group <- NULL
  original$groups <- NULL
  expect_error(validate_population(original), "declared group")
})

test_that("historical identity lookup retains unknown references and nonmoving people", {
  pop <- transition_observation_population()
  original <- pop
  # At the first woman's move, a left-limit query sees A and a snapshot sees B.
  expect_identical(hammel:::groups_at(pop, c(1, NA_integer_, 3, 2), c(2, 2, 2, 0)),
    c("B", NA_character_, "A", "A"))
  expect_identical(hammel:::groups_at(pop, 1, 2, left = TRUE), "A")
  expect_identical(hammel:::groups_at(pop, integer(), numeric()), character())
  pop$transitions <- pop$transitions[0, ]
  expect_identical(hammel:::groups_at(pop, 1:3, rep(0, 3)), pop$persons$group)
  expect_equal(nrow(original$transitions), 2)
})

test_that("dissolution observations retain unions whose partners never move", {
  pop <- transition_observation_population()
  # Only the child moves; neither partner in the recorded union changes group.
  pop$transitions <- tibble::tibble(person_id = 3, time = 2,
    origin = "B", destination = "A")
  out <- pair_rates(calculate_dissolution_rates(pop, c(0, 5), c(0, Inf), by = c("group_w", "group_h")))
  expect_equal(out$union_years[out$group_pair == "B:B"], 3)
  expect_equal(out$dissolutions[out$group_pair == "B:B"], 1)
  expect_equal(sum(out$union_years), 3)
})

test_that("union observations combine simultaneous partner moves and typed numeric pairs", {
  persons <- data.frame(person_id = c(4, 3, 2, 1), sex = c("m", "f", "m", "f"),
    birth_time = -20, death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_,
    group = c(10L, 10L, 10L, 20L))
  unions <- data.frame(union_id = c(2, 1), wife_id = c(3, 1), husband_id = c(4, 2),
    start_time = 0, end_time = 4, end_reason = "dissolution", union_type = "marriage")
  tx <- data.frame(person_id = c(1, 2), time = 2, origin = c(10L, 20L), destination = c(20L, 10L))
  pop <- hammel_population(persons, unions, end_time = 5, transitions = tx, groups = if ("group" %in% names(persons)) "group")
  out <- pair_rates(calculate_dissolution_rates(pop, c(0, 5), c(0, Inf),
    by = c("group_w", "group_h", "union_type")))
  expect_equal(subset(out, group_pair == "10:20" & union_type == "marriage")$union_years, 4)
  expect_equal(subset(out, group_pair == "10:10" & union_type == "marriage")$union_years, 4)
  expect_equal(sum(out$dissolutions), 2)
  expect_equal(sum(out$union_years), 8)
  expect_equal(sum(out$union_years[out$union_type == "cohab"]), 0)

  # For conversion exposure only the woman's moves matter; his simultaneous
  # move must not introduce a second boundary or change her duration origin.
  pop$unions <- transform(pop$unions, union_type = "cohab")
  conversion <- calculate_cohab_convert_rate(pop, c(0, 5), dur_breaks = c(0, Inf), by = "group_w")
  expect_identical(conversion$group_w, c(10L, 20L))
  expect_equal(conversion$union_years, c(6, 2))
  expect_equal(sum(conversion$conversions), 0)
})
