transition_rate_population <- function() {
  persons <- data.frame(person_id = c(11, 29), sex = c("f", "m"), birth_time = 0,
    death_time = c(NA, 4), mother_id = NA_real_, father_id = NA_real_, group = c("A", "C"))
  tx <- data.frame(person_id = c(11, 11), time = c(5, 2),
    group_origin = c("B", "A"), group_dest = c("A", "B"))
  unions <- data.frame(union_id = 1, wife_id = 11, husband_id = 29,
    start_time = 1, end_time = 3, end_reason = "dissolution")
  hammel_population(persons, unions, end_time = 6, transitions = tx)
}

test_that("directed age rates use origin exposure and retain zero-event directions", {
  pop <- transition_rate_population()
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 6)
  expect_false(any(out$group_origin == out$group_dest))
  expect_false(any(c("sex", "union_type") %in% names(out)))
  expect_equal(out$person_years[out$group_origin == "A"], c(3, 3))
  expect_equal(out$person_years[out$group_origin == "B"], c(3, 3))
  expect_equal(out$person_years[out$group_origin == "C"], c(4, 4))
  expect_equal(out$transitions[out$group_origin == "A" & out$group_dest == "B"], 1)
  expect_equal(out$transitions[out$group_origin == "B" & out$group_dest == "A"], 1)
  expect_equal(sum(out$transitions), 2)
  expect_equal(out$rate_transition[out$transitions > 0], c(1/3, 1/3))
  clipped <- calculate_transition_rates(pop, c(3, 6), age_breaks = c(0, Inf))
  expect_equal(clipped$person_years[clipped$group_origin == "A"], c(1, 1))
  expect_equal(clipped$person_years[clipped$group_origin == "B"], c(2, 2))
  expect_equal(clipped$person_years[clipped$group_origin == "C"], c(1, 1))
  expect_equal(sum(clipped$transitions), 1)
})

test_that("duration counts reset at group entry, retain initial age and ignore union changes", {
  out <- calculate_transition_rates(transition_rate_population(), c(0, 6), dur_breaks = c(0, 2))
  ab <- subset(out, group_origin == "A" & group_dest == "B")
  ba <- subset(out, group_origin == "B" & group_dest == "A")
  expect_equal(ab$dur_end, c(2, Inf))
  expect_equal(ab$person_years, c(3, 0))
  expect_equal(ab$transitions, c(0, 1))
  expect_equal(ab$rate_transition[1], 0)
  expect_true(is.na(ab$rate_transition[2]))
  expect_equal(ba$person_years, c(2, 1))
  expect_equal(ba$transitions, c(0, 1))
  c_group <- subset(out, group_origin == "C" & group_dest == "A")
  expect_equal(c_group$person_years, c(2, 2))
  late <- calculate_transition_rates(transition_rate_population(), c(3, 6), dur_breaks = c(0, 2, Inf))
  expect_equal(subset(late, group_origin == "B" & group_dest == "A")$person_years, c(1, 1))
})

test_that("sex and union-type components recover pooled counts and exposure", {
  pop <- transition_rate_population()
  original <- pop
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf),
    by_sex = TRUE, by_union_type = TRUE)
  ab <- subset(out, group_origin == "A" & group_dest == "B" & sex == "f")
  ba <- subset(out, group_origin == "B" & group_dest == "A" & sex == "f")
  expect_equal(ab$person_years, c(2, 1))
  expect_equal(ab$transitions, c(0, 1))
  expect_equal(ba$person_years, c(2, 1))
  expect_equal(ba$transitions, c(1, 0))
  expect_true(all(out$person_years[out$sex == "m" & out$group_origin == "A"] == 0))
  expect_true(all(is.na(out$rate_transition[out$person_years == 0])))
  pooled <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  for (j in seq_len(nrow(pooled))) {
    rows <- out$group_origin == pooled$group_origin[j] & out$group_dest == pooled$group_dest[j]
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
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf), by_union_type = TRUE)
  ab <- subset(out, group_origin == "A" & group_dest == "B")
  expect_equal(ab$union_type, c("unpartnered", "cohab", "marriage"))
  expect_equal(ab$person_years, c(2, 1, 0))
  expect_equal(ab$transitions, c(0, 0, 1))
  pop$transitions$time[1] <- 3 # B -> A exactly at dissolution
  pop$transitions$time[2] <- 1 # A -> B exactly at union formation
  out <- calculate_transition_rates(pop, c(0, 6), dur_breaks = c(0, Inf), by_union_type = TRUE)
  expect_equal(subset(out, group_origin == "A" & group_dest == "B" & union_type == "cohab")$transitions, 1)
  expect_equal(subset(out, group_origin == "B" & group_dest == "A" & union_type == "unpartnered")$transitions, 1)
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
  out <- calculate_transition_rates(pop, c(0, 6), age_breaks = c(0, Inf))
  expect_identical(out$group_origin, c(10L, 20L))
  expect_equal(out$person_years, c(6, 4))
  expect_equal(out$transitions, c(0, 0))
  expect_equal(out$rate_transition, c(0, 0))
  pop$persons <- pop$persons[0, ]; pop$unions <- pop$unions[0, ]
  pop$ancestry <- pop$ancestry[0, ]
  empty <- calculate_transition_rates(pop, c(0, 6), dur_breaks = c(0, Inf), by_sex = TRUE, by_union_type = TRUE)
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 0)
  expect_type(empty$group_origin, "integer")
  expect_type(empty$rate_transition, "double")
  one <- hammel_population(data.frame(person_id = 1, sex = "m", birth_time = 0,
    death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_, group = "A"), end_time = 6)
  expect_equal(nrow(calculate_transition_rates(one, c(0, 6), age_breaks = c(0, Inf))), 0)
})

test_that("transition counter rejects invalid axes, endpoints and flags", {
  pop <- transition_rate_population()
  expect_error(calculate_transition_rates(pop, c(0, 6)), "exactly one")
  expect_error(calculate_transition_rates(pop, c(0, 6), c(0, Inf), c(0, Inf)), "exactly one")
  expect_error(calculate_transition_rates(pop, c(0, 7), c(0, Inf)), "endpoint")
  expect_error(calculate_transition_rates(pop, c(0, 6), c(0, 0)), "increasing")
  expect_error(calculate_transition_rates(pop, c(0, 6), c(-1, Inf)), "nonnegative")
  expect_error(calculate_transition_rates(pop, c(0, 6), c(0, Inf), by_sex = NA), "by_sex")
  expect_error(calculate_transition_rates(pop, c(0, 6), c(0, Inf), by_union_type = 1), "by_union_type")
  pop$persons$group <- NULL; pop$transitions <- NULL
  expect_error(calculate_transition_rates(pop, c(0, 6), c(0, Inf)), "group")
})

test_that("simulated competing hazards recover directed rates with shared origin exposure", {
  model <- make_model(10000)
  model$population$persons$group <- rep("A", 10000)
  model <- set_group_transition(model, data.frame(group_origin = "A", group_dest = c("B", "C"),
    age_start = 0, age_end = Inf, rate = c(.2, .3)))
  pop <- simulate(model, 2, seed = 317)
  out <- calculate_transition_rates(pop, c(0, 2), age_breaks = c(0, Inf))
  a <- subset(out, group_origin == "A")
  expect_equal(a$person_years[1], a$person_years[2])
  targets <- c(B = .2, C = .3)[a$group_dest]
  expect_lt(max(abs(a$rate_transition - targets) / sqrt(targets / a$person_years)), 5)
  expect_equal(sum(out$transitions), nrow(pop$transitions))
  # No reproductive or death events: origin exposures sum to 20,000 years,
  # after removing the intentional denominator duplication by destination.
  expect_equal(sum(out$person_years) / 2, 20000)
})
