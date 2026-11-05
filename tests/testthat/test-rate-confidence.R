confidence_population <- function() {
  persons <- tibble::tibble(person_id = 1:3, sex = c("f", "m", "f"),
    birth_time = c(-20, -25, 1), death_time = c(NA, 7, NA),
    mother_id = c(NA, NA, 1), father_id = c(NA, NA, 2), group = c("A", "B", "A"))
  unions <- tibble::tibble(union_id = 1:3, wife_id = 1, husband_id = 2,
    start_time = c(0, 2, 3), end_time = c(2, 3, 6),
    end_reason = c("conversion", "dissolution", "dissolution"),
    union_type = c("cohab", "marriage", "cohab"),
    previous_cohab_id = c(NA, 1, NA))
  hammel_population(persons, unions, end_time = 8,
    transitions = tibble::tibble(person_id = 1, time = 4,
      group_origin = "B", group_dest = "A"), groups = if ("group" %in% names(persons)) "group")
}

confidence_calls <- function(population) {
  list(
    mortality = list(calculate_mortality_rates, list(population, c(0, 2, 8), c(0, 20, Inf), by = c("group", "union_type")), "deaths", "person_years"),
    fertility = list(calculate_fertility_rates, list(population, c(0, 2, 8), c(0, 20, Inf), parity_breaks = c(0, 1, Inf)), "births", "person_years"),
    union = list(calculate_union_rates, list(population, c(0, 2, 8), c(0, Inf), by = c("union_type", "previous_union")), "unions", "person_years"),
    marriage = list(calculate_marriage_rates, list(population, c(0, 2, 8), c(0, Inf), by = "group"), "marriages", "person_years"),
    dissolution = list(calculate_dissolution_rates, list(population, c(0, 2, 8), c(0, 2, Inf), by = "union_type"), "dissolutions", "union_years"),
    conversion = list(calculate_cohab_convert_rate, list(population, c(0, 2, 8), dur_breaks = c(0, 2, Inf)), "conversions", "union_years"),
    transition = list(calculate_transition_rates, list(population, c(0, 2, 8), age_breaks = c(0, Inf), by = "sex"), "transitions", "person_years")
  )
}

check_rate_intervals <- function(x, events, exposure, level) {
  expected <- matrix(NA_real_, nrow(x), 2)
  for (i in which(x[[exposure]] > 0))
    expected[i, ] <- stats::poisson.test(x[[events]][i], T = x[[exposure]][i],
      conf.level = level)$conf.int
  expect_equal(x$conf_low, expected[, 1])
  expect_equal(x$conf_high, expected[, 2])
}

test_that("all rate calculators append exact intervals without changing other output", {
  for (case in confidence_calls(confidence_population())) {
    baseline <- do.call(case[[1]], case[[2]])
    expect_identical(do.call(case[[1]], c(case[[2]], list(conf.level = NULL))), baseline)
    for (level in c(.95, .8)) {
      x <- do.call(case[[1]], c(case[[2]], list(conf.level = level)))
      expect_identical(x[names(baseline)], baseline)
      expect_identical(tail(names(x), 2), c("conf_low", "conf_high"))
      check_rate_intervals(x, case[[3]], case[[4]], level)
    }
  }
})

test_that("all calculators reject invalid confidence levels", {
  bad <- list(0, 1, -.1, 1.1, NA_real_, NaN, Inf, "0.95", TRUE,
    numeric(), c(.9, .95), matrix(.95))
  for (case in confidence_calls(confidence_population())) for (level in bad)
    expect_error(do.call(case[[1]], c(case[[2]], list(conf.level = level))), "conf.level")
})

test_that("zero events have positive upper limits and zero exposure has NA bounds", {
  persons <- tibble::tibble(person_id = 1, sex = "f", birth_time = 0,
    death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_)
  x <- calculate_mortality_rates(persons, c(0, 100), c(0, Inf), conf.level = .95)
  expect_equal(x$conf_low[x$sex == "f"], 0)
  expect_equal(x$conf_high[x$sex == "f"], -log(.025) / 100)
  expect_true(all(is.na(x$conf_high[x$person_years == 0])))
  # The conversion at the left endpoint has a count but no exposure in this band.
  x <- calculate_cohab_convert_rate(confidence_population(), c(2, 3),
    dur_breaks = c(0, 2, Inf), conf.level = .95)
  expect_true(any(x$conversions > 0 & x$union_years == 0))
  expect_true(all(is.na(x$conf_low)))
  expect_true(all(is.na(x$conf_high)))
})

test_that("empty grouped rate results retain numeric confidence columns", {
  pop <- confidence_population()
  pop$persons <- pop$persons[0, ]
  pop$unions <- pop$unions[0, ]
  pop <- hammel_population(pop$persons, pop$unions, end_time = 8, groups = if ("group" %in% names(pop$persons)) "group")
  for (case in confidence_calls(pop)) {
    group_by <- if (identical(case[[1]], calculate_transition_rates)) NULL else
      if (identical(case[[1]], calculate_dissolution_rates)) c("group_w", "group_h") else
      if (identical(case[[1]], calculate_cohab_convert_rate)) "group_w" else "group"
    case[[2]]$by <- unique(c(case[[2]]$by, group_by))
    x <- do.call(case[[1]], c(case[[2]], list(conf.level = .95)))
    expect_equal(nrow(x), 0)
    expect_identical(x$conf_low, numeric())
    expect_identical(x$conf_high, numeric())
  }
})

test_that("fertility intervals use observed births with pooled or unrealized output", {
  pop <- confidence_population()
  pop$simulation <- list(start_time = 0, unrealized_births = tibble::tibble(
    period_start = 0, period_end = 2, age_start = 20, age_end = 21,
    unrealized_births = 2.5))
  for (parity in list(NULL, c(0, 1, Inf))) {
    x <- calculate_fertility_rates(pop, c(0, 2, 8), c(0, Inf),
      parity_breaks = parity, unrealized = TRUE, conf.level = .95)
    check_rate_intervals(x, "births", "person_years", .95)
    expect_true("rate_fertility_adjusted" %in% names(x))
    if (is.null(parity)) {
      expect_equal(x$unrealized_births[1], 2.5)
      expect_gt(x$rate_fertility_adjusted[1], x$rate_fertility[1])
    }
  }
})
