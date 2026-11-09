transition_rate_population <- function() {
  persons <- data.frame(person_id = c(11, 29), sex = c("f", "m"), birth_time = 0,
    death_time = c(NA, 4), mother_id = NA_real_, father_id = NA_real_, group = c("A", "C"))
  tx <- data.frame(person_id = c(11, 11), time = c(5, 2),
    origin = c("B", "A"), destination = c("A", "B"))
  unions <- data.frame(union_id = 1, wife_id = 11, husband_id = 29,
    start_time = 1, end_time = 3, end_reason = "dissolution")
  hammel_population(persons, unions, end_time = 6, transitions = tx, groups = if ("group" %in% names(persons)) "group")
}

test_that("directed age rates use origin exposure and retain zero-event directions", {
  pop <- transition_rate_population()
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 6)
  expect_false(any(out$origin == out$destination))
  expect_false(any(c("sex", "union_type") %in% names(out)))
  expect_equal(out$person_years[out$origin == "A"], c(3, 3))
  expect_equal(out$person_years[out$origin == "B"], c(3, 3))
  expect_equal(out$person_years[out$origin == "C"], c(4, 4))
  expect_equal(out$transitions[out$origin == "A" & out$destination == "B"], 1)
  expect_equal(out$transitions[out$origin == "B" & out$destination == "A"], 1)
  expect_equal(sum(out$transitions), 2)
  expect_equal(out$rate_transition[out$transitions > 0], c(1/3, 1/3))
  clipped <- calculate_transition_rates(pop, c(3, 6), age_breaks = c(0, Inf))
  expect_equal(clipped$person_years[clipped$origin == "A"], c(1, 1))
  expect_equal(clipped$person_years[clipped$origin == "B"], c(2, 2))
  expect_equal(clipped$person_years[clipped$origin == "C"], c(1, 1))
  expect_equal(sum(clipped$transitions), 1)
})

test_that("duration counts reset at group entry, retain initial age and ignore union changes", {
  out <- calculate_transition_rates(transition_rate_population(), c(0, 6), dur_breaks = c(0, 2))
  ab <- subset(out, origin == "A" & destination == "B")
  ba <- subset(out, origin == "B" & destination == "A")
  expect_equal(ab$dur_end, c(2, Inf))
  expect_equal(ab$person_years, c(3, 0))
  expect_equal(ab$transitions, c(0, 1))
  expect_equal(ab$rate_transition[1], 0)
  expect_true(is.na(ab$rate_transition[2]))
  expect_equal(ba$person_years, c(2, 1))
  expect_equal(ba$transitions, c(0, 1))
  c_group <- subset(out, origin == "C" & destination == "A")
  expect_equal(c_group$person_years, c(2, 2))
  late <- calculate_transition_rates(transition_rate_population(), c(3, 6), dur_breaks = c(0, 2, Inf))
  expect_equal(subset(late, origin == "B" & destination == "A")$person_years, c(1, 1))
})

test_that("sex and union-type components recover pooled counts and exposure", {
  pop <- transition_rate_population()
  original <- pop
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf),
    by = c("sex", "union_type"))
  ab <- subset(out, origin == "A" & destination == "B" & sex == "f")
  ba <- subset(out, origin == "B" & destination == "A" & sex == "f")
  expect_equal(ab$person_years, c(2, 1))
  expect_equal(ab$transitions, c(0, 1))
  expect_equal(ba$person_years, c(2, 1))
  expect_equal(ba$transitions, c(1, 0))
  expect_true(all(out$person_years[out$sex == "m" & out$origin == "A"] == 0))
  expect_true(all(is.na(out$rate_transition[out$person_years == 0])))
  pooled <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  for (j in seq_len(nrow(pooled))) {
    rows <- out$origin == pooled$origin[j] & out$destination == pooled$destination[j]
    expect_equal(sum(out$person_years[rows]), pooled$person_years[j])
    expect_equal(sum(out$transitions[rows]), pooled$transitions[j])
  }
  expect_identical(pop, original)
})

test_that("typed union boundaries classify moves after conversion and dissolution", {
  pop <- transition_rate_population()
  pop$unions <- tibble::tibble(union_id = c(1, 2), wife_id = 11, husband_id = 29,
    start_time = c(1, 2), end_time = c(2, 3), end_reason = c("conversion", "dissolution"),
    union_type = c("cohab", "marriage"), previous_cohab_id = c(NA, 1))
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = "union_type")
  ab <- subset(out, origin == "A" & destination == "B")
  expect_equal(ab$union_type, c("unpartnered", "cohab", "marriage"))
  expect_equal(ab$person_years, c(2, 1, 0))
  expect_equal(ab$transitions, c(0, 0, 1))
  pop$transitions$time[1] <- 3 # B -> A exactly at dissolution
  pop$transitions$time[2] <- 1 # A -> B exactly at union formation
  out <- calculate_transition_rates(pop, c(0, 6), dur_breaks = c(0, Inf), by = "union_type")
  expect_equal(subset(out, origin == "A" & destination == "B" & union_type == "cohab")$transitions, 1)
  expect_equal(subset(out, origin == "B" & destination == "A" & union_type == "unpartnered")$transitions, 1)
})

test_that("calendar and age boundaries use half-open event cells", {
  out <- calculate_transition_rates(transition_rate_population(), c(0, 2, 5), age_breaks = c(0, 2, Inf))
  expect_equal(sum(out$transitions), 1) # move at endpoint 5 excluded
  event <- out[out$transitions > 0, ]
  expect_equal(event$time_start, 2)
  expect_equal(event$age_start, 2)
  expect_equal(event$person_years, 0)
  expect_true(is.na(event$rate_transition))
})

test_that("empty histories, empty populations and group label types are retained", {
  pop <- transition_rate_population()
  pop$transitions <- pop$transitions[0, ]
  pop$persons$group <- c(10L, 20L)
  pop$groups <- list(group = group())
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  expect_identical(out$origin, c(10L, 20L))
  expect_equal(out$person_years, c(6, 4))
  expect_equal(out$transitions, c(0, 0))
  expect_equal(out$rate_transition, c(0, 0))
  pop$persons <- pop$persons[0, ]
  pop$unions <- pop$unions[0, ]
  pop$ancestry <- pop$ancestry[0, ]
  empty <- calculate_transition_rates(pop, c(0, 6), dur_breaks = c(0, Inf), by = c("sex", "union_type"))
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 0)
  expect_type(empty$origin, "integer")
  expect_type(empty$rate_transition, "double")
  one <- hammel_population(data.frame(person_id = 1, sex = "m", birth_time = 0,
    death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_, group = "A"), end_time = 6, groups = "group")
  expect_equal(nrow(calculate_transition_rates(one, c(0, 6), age_breaks = c(0, Inf))), 0)
})

test_that("transition counter rejects invalid axes, endpoints and flags", {
  pop <- transition_rate_population()
  expect_error(calculate_transition_rates(pop, c(0, 6)), "exactly one")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), dur_breaks = c(0, Inf)), "exactly one")
  expect_error(calculate_transition_rates(pop, c(0, 7), age_breaks = c(0, Inf)), "endpoint")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, 0)), "increasing")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(-1, Inf)), "nonnegative")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = NA), "`by` must")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = 1), "`by` must")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = "group"),
               "not by `group` in `by`")
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = "region"), "cannot split")
  pop$persons$group <- NULL
  pop$transitions <- NULL
  pop$groups <- list(group = group())
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf)), "group")
})

test_that("simulated competing hazards recover directed rates with shared origin exposure", {
  model <- make_model(10000)
  model$population$persons$group <- rep("A", 10000)
  model$population$groups <- list(group = group())
  model <- set_group_transition(model, data.frame(origin = "A", destination = c("B", "C"),
    age_start = 0, age_end = Inf, rate = c(.2, .3)))
  pop <- run_simulation(model, 2, seed = 317)
  out <- calculate_transition_rates(pop, c(0, 2), age_breaks = c(0, Inf))
  a <- subset(out, origin == "A")
  expect_equal(a$person_years[1], a$person_years[2])
  targets <- c(B = .2, C = .3)[a$destination]
  expect_lt(max(abs(a$rate_transition - targets) / sqrt(targets / a$person_years)), 5)
  expect_equal(sum(out$transitions), nrow(pop$transitions))
  # No reproductive or death events: origin exposures sum to 20,000 years,
  # after removing the intentional denominator duplication by destination.
  expect_equal(sum(out$person_years) / 2, 20000)
})

test_that("transition counters accept legacy grouped populations without history components", {
  pop <- transition_rate_population()
  pop$transitions <- NULL
  expected <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  expect_equal(expected$transitions, c(0, 0))
  expect_equal(expected$person_years, c(6, 4))
  expect_error(calculate_transition_rates(pop$persons, c(0, 6), age_breaks = c(0, Inf)),
    "grouped hammel_population")
})

test_that("transition-rate errors identify the requested observation", {
  pop <- transition_rate_population()
  pop$persons$group <- NULL
  pop$groups <- NULL
  pop$transitions <- NULL
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf)),
    "calculate_transition_rates.*requires a declared group")
  pop <- transition_rate_population()
  pop$unions <- rbind(pop$unions, transform(pop$unions, union_id = 2, start_time = 2))
  expect_error(calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by = "union_type"),
    "requires nonoverlapping")
})
