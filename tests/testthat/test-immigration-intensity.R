# Intensity tables with counts or crude rates, and crude immigration rates
# (D194).

# A destination of `n` men aged 20 who neither die nor have children.
static_model <- function(n = 200, death_rate = 0) {
  make_model(n, death_rate = death_rate)
}

rate_template <- function() {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "m", "f", "m"), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
  )
  hammel_population(persons, end_time = 0)
}

test_that("rate intensities grow a closed population exponentially", {
  # With no births or deaths, arrivals at rate r per person-year make a
  # pure-birth process: E[N(T)] = N(0) exp(r T).
  model <- set_immigration(
    static_model(200),
    a = immigration_stream(
      data.frame(rate = 0.05), rate_template(),
      household = NULL
    )
  )
  sizes <- vapply(1:20, function(seed) {
    result <- run_simulation(model, 10, seed = seed)
    sum(is.na(result$persons$death_time))
  }, numeric(1))
  growth <- exp(0.05 * 10)
  expected <- 200 * growth
  sd_run <- sqrt(200 * growth * (growth - 1))
  expect_lt(abs(mean(sizes) - expected), 4 * sd_run / sqrt(20))
})

test_that("crude immigration rates recover a rate intensity", {
  # Deaths shrink the population, and immigration with it.
  model <- set_immigration(
    static_model(400, death_rate = 0.05),
    a = immigration_stream(
      data.frame(rate = 0.03), rate_template(),
      household = NULL
    )
  )
  result <- run_simulation(model, 20, seed = 1)
  rates <- calculate_immigration_rates(result, c(0, 20), conf_level = 0.95)
  expect_named(
    rates,
    c(
      "time_start", "time_end", "arrivals", "person_years",
      "rate_immigration", "dispersion", "conf_low", "conf_high"
    )
  )
  expect_gt(rates$arrivals, 100)
  expect_lt(
    abs(rates$rate_immigration - 0.03),
    4 * sqrt(rates$arrivals) / rates$person_years
  )
  expect_true(
    rates$conf_low <= rates$rate_immigration &&
      rates$rate_immigration <= rates$conf_high
  )
})

test_that("household rate streams count people", {
  # The worked-example households average 2 people when equally weighted.
  persons <- data.frame(
    person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-30, -32, -40),
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
  )
  unions <- data.frame(
    union_id = 1, wife_id = 1, husband_id = 2, start_time = -5,
    end_time = NA_real_, end_reason = NA_character_
  )
  template <- hammel_population(persons, unions = unions, end_time = 0)
  model <- set_immigration(
    static_model(300),
    a = immigration_stream(data.frame(rate = 0.04), template)
  )
  result <- run_simulation(model, 10, seed = 2)
  rates <- calculate_immigration_rates(result, c(0, 10))
  expect_lt(
    abs(rates$rate_immigration - 0.04),
    4 * sqrt(2 * rates$arrivals) / rates$person_years
  )
})

test_that("streams can mix counts and rates", {
  model <- set_immigration(
    static_model(300),
    counted = immigration_stream(
      per_year(20), rate_template(),
      household = NULL
    ),
    scaled = immigration_stream(
      data.frame(rate = 0.02), rate_template(),
      household = NULL
    )
  )
  result <- run_simulation(model, 10, seed = 3)
  rates <- calculate_immigration_rates(result, c(0, 10), streams = TRUE)
  expect_identical(rates$birthplace, c("counted", "scaled"))
  expect_lt(abs(rates$arrivals[1] - 200), 4 * sqrt(200))
  scaled <- rates[rates$birthplace == "scaled", ]
  expect_lt(
    abs(scaled$rate_immigration - 0.02),
    4 * sqrt(scaled$arrivals) / scaled$person_years
  )
  # Streams share the total person-years, so their rates sum to the total.
  total <- calculate_immigration_rates(result, c(0, 10))
  expect_equal(sum(rates$rate_immigration), total$rate_immigration)
  expect_equal(unique(rates$person_years), total$person_years)
})

test_that("calendar rate tables pause and resume", {
  intensity <- data.frame(
    time_start = c(0, 3), time_end = c(3, Inf), rate = c(0, 0.05)
  )
  model <- set_immigration(
    static_model(200),
    a = immigration_stream(intensity, rate_template(), household = NULL)
  )
  result <- run_simulation(model, 6, seed = 4)
  arrived <- result$persons$arrival_time[
    !is.na(result$persons$arrival_time)
  ]
  expect_gt(length(arrived), 0)
  expect_true(all(arrived >= 3))
  rates <- calculate_immigration_rates(result, c(0, 3, 6))
  expect_identical(rates$arrivals[1], 0)
})

test_that("crude immigration rates count person-years exactly", {
  model <- set_immigration(
    make_model(20, death_rate = 0.1),
    a = immigration_stream(per_year(5), rate_template(), household = NULL)
  )
  result <- run_simulation(model, 4, seed = 5)
  rates <- calculate_immigration_rates(result, c(0, 1, 4))
  persons <- result$persons
  lived <- function(lo, hi) {
    total <- 0
    for (i in seq_len(nrow(persons))) {
      start <- max(persons$birth_time[i], persons$arrival_time[i], na.rm = TRUE)
      stop <- min(persons$death_time[i], result$end_time, na.rm = TRUE)
      total <- total + max(0, min(stop, hi) - max(start, lo))
    }
    total
  }
  expect_equal(rates$person_years, c(lived(0, 1), lived(1, 4)))
  arrival <- persons$arrival_time
  expect_identical(
    rates$arrivals,
    c(
      sum(arrival >= 0 & arrival < 1, na.rm = TRUE),
      sum(arrival >= 1 & arrival < 4, na.rm = TRUE)
    ) + 0
  )
})

test_that("crude immigration rate arguments are validated", {
  plain <- run_simulation(static_model(10), 1, seed = 6)
  rates <- calculate_immigration_rates(plain, c(0, 1))
  expect_identical(rates$arrivals, 0)
  expect_hammel_error(
    calculate_immigration_rates(plain, c(0, 1), streams = TRUE),
    "groups", "calculate_immigration_rates"
  )
  expect_hammel_error(
    calculate_immigration_rates(plain, c(0, 1), by = "race"),
    "argument", "calculate_immigration_rates"
  )
  migrated <- run_simulation(
    set_immigration(
      static_model(10),
      a = immigration_stream(per_year(5), rate_template(), household = NULL)
    ),
    1,
    seed = 7
  )
  error <- expect_hammel_error(
    calculate_immigration_rates(migrated, c(0, 1), by = "birthplace"),
    "argument", "calculate_immigration_rates"
  )
  expect_identical(error$argument, "by")
  summary_text <- cli::cli_fmt(print(summary(immigration(
    a = immigration_stream(data.frame(rate = 0.01), rate_template())
  ))))
  expect_true(any(grepl("0.01 per person-year", summary_text, fixed = TRUE)))
})

test_that("immigration rates stratify arrivals and person-years together", {
  # Arrivals settle in the north; everyone moves south at rate 0.5.
  model <- static_model(100)
  model$population$persons$region <- "n"
  model <- model |>
    set_groups(region = group(categories = c("n", "s"))) |>
    set_group_transition(
      data.frame(
        origin = "n", destination = "s", age_start = 0, age_end = Inf,
        rate = 0.5
      ),
      group = "region"
    ) |>
    set_immigration(
      a = immigration_stream(
        per_year(30), rate_template(),
        assignment = list(
          region = migrant_assignment("fixed", fixed_value = "n")
        ),
        household = NULL
      )
    )
  result <- run_simulation(model, 4, seed = 8)
  by_region <- calculate_immigration_rates(result, c(0, 4), by = "region")
  total <- calculate_immigration_rates(result, c(0, 4))
  expect_identical(by_region$region, c("n", "s"))
  # Everyone arrives in the north, even those who have since moved south.
  expect_gt(sum(result$persons$region[!is.na(result$persons$arrival_time)] ==
                  "s"), 0)
  expect_identical(by_region$arrivals, c(total$arrivals, 0))
  expect_equal(sum(by_region$person_years), total$person_years)
  expect_equal(
    by_region$rate_immigration[1],
    total$arrivals / by_region$person_years[1]
  )
  # By sex and stream: strata partition arrivals and person-years.
  split <- calculate_immigration_rates(
    result, c(0, 2, 4),
    by = c("sex", "region"),
    streams = TRUE
  )
  expect_named(
    split,
    c(
      "time_start", "time_end", "sex", "region", "birthplace", "arrivals",
      "person_years", "rate_immigration"
    )
  )
  expect_identical(sum(split$arrivals), total$arrivals)
  women <- split[split$sex == "f", ]
  expect_identical(
    sum(women$arrivals),
    sum(result$persons$sex == "f" & !is.na(result$persons$arrival_time)) + 0
  )
})

# A template of two-person households: a couple and a mother with a child.
couple_template <- function() {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "m", "f", "m"),
    birth_time = c(-30, -32, -35, -5), death_time = NA_real_,
    mother_id = c(NA, NA, NA, 3), father_id = NA_integer_
  )
  unions <- data.frame(
    union_id = 1, wife_id = 1, husband_id = 2, start_time = -5,
    end_time = NA_real_, end_reason = NA_character_
  )
  hammel_population(persons, unions = unions, end_time = 0)
}

test_that("individual arrivals keep exact Poisson intervals", {
  model <- set_immigration(
    static_model(200),
    a = immigration_stream(
      data.frame(rate = 0.03), rate_template(),
      household = NULL
    )
  )
  result <- run_simulation(model, 5, seed = 9)
  rates <- calculate_immigration_rates(result, c(0, 5), conf_level = 0.9)
  expect_identical(rates$dispersion, 1)
  exact <- stats::poisson.test(
    rates$arrivals, rates$person_years,
    conf.level = 0.9
  )$conf.int
  expect_equal(c(rates$conf_low, rates$conf_high), as.numeric(exact))
})

test_that("household arrivals widen intervals by their dispersion", {
  model <- set_immigration(
    static_model(300),
    a = immigration_stream(data.frame(rate = 0.03), couple_template())
  )
  result <- run_simulation(model, 5, seed = 10)
  rates <- calculate_immigration_rates(result, c(0, 5), conf_level = 0.95)
  # Every household brings two people: dispersion 4 / 2.
  expect_equal(rates$dispersion, 2)
  plain <- stats::poisson.test(
    rates$arrivals, rates$person_years,
    conf.level = 0.95
  )$conf.int
  expect_gt(rates$conf_high - rates$conf_low, diff(plain))
  # Split by sex, each household has one member in each stratum.
  by_sex <- calculate_immigration_rates(
    result, c(0, 5),
    by = "sex",
    conf_level = 0.95
  )
  expect_equal(by_sex$dispersion, c(1, 1))
  # Without intervals, dispersion is not reported.
  expect_false(
    "dispersion" %in% names(calculate_immigration_rates(result, c(0, 5)))
  )
  # Across replications the adjusted 90% intervals cover the input rate.
  covered <- vapply(1:40, function(seed) {
    r <- calculate_immigration_rates(
      run_simulation(model, 5, seed = 100 + seed), c(0, 5),
      conf_level = 0.9
    )
    r$conf_low <= 0.03 && 0.03 <= r$conf_high
  }, logical(1))
  expect_gt(mean(covered), 0.75)
})
