# Household moves (D140). One region dimension, N and S. A family: W (woman,
# 34), her husband H (36), their child C (10), H's son S (12) by his first
# wife F (dead), and W's daughter D (19, unpartnered; father unknown).
family <- function(h_region = "N", union_end = NA_real_, w_death = NA_real_, d_age = 19) {
  persons <- data.frame(
    person_id = 1:6, sex = c("f", "m", "f", "m", "f", "f"),
    birth_time = c(-34, -36, -10, -12, -d_age, -40),
    death_time = c(w_death, NA, NA, NA, NA, -5),
    mother_id = c(NA, NA, 1, 6, 1, NA), father_id = c(NA, NA, 2, 2, NA, NA),
    region = c("N", h_region, "N", "N", "N", "N"))
  unions <- data.frame(union_id = 1:2, wife_id = c(1, 6), husband_id = c(2, 2),
    start_time = c(-11, -14), end_time = c(union_end, -13),
    end_reason = c(if (is.na(union_end)) NA else "dissolution", "dissolution"))
  if (!is.na(w_death)) { unions$end_time[1] <- w_death; unions$end_reason[1] <- "death_wife" }
  hammel_population(persons, unions = unions, end_time = 0, groups = "region")
}

household_model <- function(population, rates, rule = household()) {
  hammel_model() |>
    set_population(population) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0) |>
    set_group_transition(rates, group = "region", household = rule)
}

# Only people aged 30-40 in N move (to S), almost at once.
adult_moves <- data.frame(origin = "N", destination = "S", age_start = 30, age_end = 40, rate = 1e8)

moves <- function(out) {
  tx <- out$transitions
  tx[order(tx$move_id, tx$person_id), ]
}

test_that("the woman carries the household's clock and her household moves with her", {
  withr::local_options(hammel.check_cells = TRUE)
  out <- simulate(household_model(family(), adult_moves), .5, seed = 1)
  tx <- moves(out)
  # W (34) carries; H (36) has no clock; C and S are dependents; D (19) is
  # independent with no risk at her age and stays.
  expect_identical(tx$person_id, c(1, 2, 3, 4))
  expect_identical(unique(tx$move_id), 1)
  expect_identical(unique(tx$initiator_id), 1)
  expect_identical(out$persons$region, c("S", "S", "S", "S", "N", "N"))
  # Without a rule, each adult moves on their own clock and children stay.
  individual <- moves(simulate(household_model(family(), adult_moves, rule = NULL), .5, seed = 1))
  expect_setequal(individual$person_id, c(1, 2))
  expect_identical(individual$initiator_id, individual$person_id)
})

test_that("members already at the destination stay and are not recorded", {
  # Example 4: H already lives in S.
  tx <- moves(simulate(household_model(family(h_region = "S"), adult_moves), .5, seed = 2))
  expect_identical(tx$person_id, c(1, 3, 4))
  expect_identical(tx$origin, rep("N", 3))
})

test_that("separation gives the man his own clock, with his own dependents", {
  # Example 5: W and H separated; C stays attached to W, S to H.
  tx <- moves(simulate(household_model(family(union_end = -1), adult_moves), .5, seed = 3))
  by_move <- split(tx$person_id, tx$move_id)
  expect_length(by_move, 2)
  expect_setequal(lapply(by_move, sort), list(c(1, 3), c(2, 4)))
})

test_that("children follow the father when the mother has died, and orphans are independent", {
  # Example 6: W died; H (a widower, so single) carries his own clock and
  # takes C (his child, mother dead) and S.
  tx <- moves(simulate(household_model(family(w_death = -1), adult_moves), .5, seed = 4))
  expect_identical(tx$person_id, c(2, 3, 4))
  expect_identical(unique(tx$initiator_id), 2)
  # With the father dead too, the children are independent and move on their
  # own clocks (rates for children here).
  pop <- family(w_death = -1)
  pop$persons$death_time[2] <- -1
  pop$unions$end_time[1] <- -1; pop$unions$end_reason[1] <- "death_wife"
  pop <- hammel_population(pop$persons, unions = pop$unions, end_time = 0, groups = "region")
  child_moves <- data.frame(origin = "N", destination = "S", age_start = 0, age_end = 15, rate = 1e8)
  tx <- moves(simulate(household_model(pop, child_moves), .5, seed = 5))
  expect_setequal(tx$person_id, c(3, 4))
  expect_identical(tx$initiator_id, tx$person_id)
})

test_that("dependents get their own clock at the age of independence", {
  # D is 17.5: a dependent of W until 18, then her own clock at the rate for 18.
  rates <- data.frame(origin = "N", destination = "S", age_start = 15, age_end = 20, rate = 1e8)
  out <- simulate(household_model(family(d_age = 17.5), rates), 1, seed = 6)
  tx <- out$transitions
  expect_identical(tx$person_id, 5)
  expect_equal(tx$time, .5, tolerance = 1e-6)
  expect_identical(tx$initiator_id, 5)
  # A different age of independence moves the boundary.
  out <- simulate(household_model(family(d_age = 17.5), rates, household(age_independence = 17.8)), 1, seed = 6)
  expect_equal(out$transitions$time, .3, tolerance = 1e-6)
})

test_that("household rules are validated and histories continue", {
  expect_error(household(-1), "nonnegative")
  expect_error(household(c(18, 19)), "finite")
  expect_error(group_transition(adult_moves, household = TRUE), "created with household")
  expect_error(group_transition(data.frame(origin = "N", destination = "S", dur_start = 0, dur_end = Inf,
    rate = 1), household = household()), "age bands")
  # Restarts continue move identifiers; transition rates count every mover.
  model <- household_model(family(), adult_moves)
  out <- simulate(model, .5, seed = 7)
  again <- simulate(set_population(model, out), .5, seed = 8)
  expect_identical(again$transitions[seq_len(nrow(out$transitions)), ], out$transitions)
  counts <- calculate_transition_rates(out, c(0, .5), age_breaks = c(0, Inf))
  expect_equal(sum(counts$transitions), 4)
  expect_equal(sum(out$transitions$person_id == out$transitions$initiator_id), 1)
})

test_that("with households living together, everyone moves at the input rate", {
  set.seed(11)
  pop <- create_baseline_population(8000, groups = list(region = c(N = .5, S = .5)))
  rates <- data.frame(origin = c("N", "S"), destination = c("S", "N"), age_start = 0, age_end = Inf, rate = .05)
  model <- hammel_model() |> set_population(pop) |>
    set_groups(region = group(inheritance = "mother")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .01)) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = .1)) |>
    set_union_formation("on_demand", kinship_generations = 0,
      matching = create_age_gap_matching_model(2, 5, group_exogamy = list(region = NA))) |>
    set_group_transition(rates, group = "region", household = household())
  out <- simulate(model, 20, seed = 12)
  r <- calculate_transition_rates(out, c(5, 20), age_breaks = c(0, 18, 40, Inf), by = "sex")
  agg <- aggregate(cbind(transitions, person_years) ~ sex + age_start, data = r, sum)
  # Partnered men and children have no clock of their own, yet move at 0.05.
  expect_equal(agg$transitions / agg$person_years, rep(.05, 6), tolerance = .2)
  expect_lt(length(unique(out$transitions$move_id)), nrow(out$transitions))
})
