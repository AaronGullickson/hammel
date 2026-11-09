test_that("Lexis counts use exact exposure and maternal ages", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  expect_equal(exposure$person_years[exposure$sex == "f"], c(1, 0.25))
  expect_equal(exposure$person_years[exposure$sex == "m"], c(0.75, 0))
  expect_equal(sum(exposure$person_years), 2, tolerance = 0)
  mx <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  expect_equal(mx$deaths[mx$sex == "f" & mx$age_start == 1], 1, tolerance = 0)
  expect_equal(mx$rate_mortality[mx$sex == "f" & mx$age_start == 1], 4, tolerance = 0)
  expect_true(is.na(mx$rate_mortality[mx$sex == "m" & mx$age_start == 1]))
  rate_fertility <- demographic_counts(population, c(0, 1), c(0, 1, 2))
  expect_identical(rate_fertility$births[rate_fertility$sex == "f"], c(0, 1, 0))
  expect_true(all(is.na(rate_fertility$births[rate_fertility$sex == "m"])))
  expect_true(all(is.na(rate_fertility$rate_fertility[rate_fertility$sex == "m"])))
  expect_equal(rate_fertility$rate_fertility[rate_fertility$sex == "f"][2], 4, tolerance = 0)
  expect_equal(attr(rate_fertility, "unassigned_births")$births, 1, tolerance = 0)
})

test_that("fractional periods and ages conserve exposure", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  split <- demographic_counts(population, c(0, 0.5, 1), c(0, 0.5, 1, Inf))
  expect_equal(sum(split$person_years), sum(exposure$person_years))
})

test_that("Lexis event boundaries are half open", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  p <- persons[1, ]
  p$birth_time <- 0
  p$death_time <- 1
  boundary <- demographic_counts(p, c(0, 1, 2), c(0, 1, Inf))
  expect_equal(sum(boundary$deaths), 1, tolerance = 0)
  expect_equal(boundary$deaths[boundary$sex == "f" & boundary$time_start == 1 & boundary$age_start ==
    1], 1, tolerance = 0)
  excluded <- demographic_counts(p, c(0, 1), c(0, 1, Inf))
  expect_equal(sum(excluded$deaths), 0, tolerance = 0)
  expect_equal(sum(excluded$person_years), 1, tolerance = 0)
  expect_error(demographic_counts(population, c(0, 2), c(0, Inf)))
})

test_that("duplicate identities in stacked histories are rejected", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  stacked <- rbind(persons, persons)
  expect_error(demographic_counts(stacked, c(0, 1), c(0, Inf)))
})

test_that("TFR uses band widths, propagates missing rates and groups periods", {
  rates <- data.frame(age_start = c(15, 20, 30), age_end = c(20, 30, 50), rate_fertility = c(0.1,
    0.2, 0.05))
  expect_equal(calculate_tfr(rates)$tfr, 3.5, tolerance = 0)
  bad <- rates
  bad$rate_fertility[2] <- NA_real_
  expect_identical(calculate_tfr(bad)$tfr, NA_real_)
  bad$rate_fertility[] <- NA_real_
  expect_identical(calculate_tfr(bad)$tfr, NA_real_)
  bad$rate_fertility[] <- 0
  expect_identical(calculate_tfr(bad)$tfr, 0)
  expect_error(calculate_tfr(rates[c(1, 3), ]))
  bad <- rates
  bad$age_end[3] <- Inf
  bad$rate_fertility[3] <- NA_real_
  expect_equal(calculate_tfr(bad)$tfr, 2.5, tolerance = 0)
  rates$time_start <- 0
  rates$time_end <- 1
  later <- rates
  later$time_start <- 1
  later$time_end <- 2
  later$rate_fertility <- later$rate_fertility * 2
  expect_identical(calculate_tfr(rbind(rates, later))$tfr, c(3.5, 7))
})

test_that("life tables use only rate or rate_mortality as input rates", {
  rates <- data.frame(age_start = c(0, 5), age_end = c(5, Inf), rate = 0.1)
  expected <- calculate_life_table(rates)
  for (column in c("rate", "rate_mortality")) {
    input <- rates
    names(input)[names(input) == "rate"] <- column
    expect_identical(calculate_life_table(input), expected)
    expect_error(
      calculate_life_table(input[FALSE, ]),
      "A rate table needs at least one age band.",
      fixed = TRUE
    )
    input$mx <- "not an input mortality rate"
    expect_identical(calculate_life_table(input), expected)
    input[[column]] <- "invalid"
    expect_error(calculate_life_table(input), paste0("rates\\$", column))
  }
  mx_only <- rates
  names(mx_only)[names(mx_only) == "rate"] <- "mx"
  expect_error(
    calculate_life_table(mx_only), "Supply `rate_mortality` or `rate`.",
    fixed = TRUE
  )
  expect_error(calculate_life_table(mx_only[FALSE, ]),
    "Supply `rate_mortality` or `rate`.", fixed = TRUE
  )
  both <- rates
  both$rate_mortality <- both$rate
  both$rate <- 0.2
  both$mx <- 0.5
  expect_identical(calculate_life_table(both), expected)
  expect_identical(expected$mx, rep(0.1, 2))
  expect_equal(expected$ex, rep(10, 2))
})

test_that("life tables reject missing stratum values", {
  rates <- data.frame(
    region = c("a", NA), age_start = c(0, 1), age_end = c(1, Inf), rate = 0.1
  )
  expect_error(
    calculate_life_table(rates),
    "`rates$region` must not contain missing values.",
    fixed = TRUE
  )
})

test_that("interval errors name the failing stratum", {
  labelled <- data.frame(
    sex = "f", age_start = c(0, 1), age_end = c(1, Inf),
    rate = 0.01, note = c("a", "b")
  )
  expect_error(
    calculate_life_table(labelled),
    "explicit open interval in stratum sex = f, note = a \\(columns other"
  )
  expect_error(
    calculate_life_table(labelled[2, c("age_start", "age_end", "rate")]),
    "explicit open interval\\.$"
  )
  overlap <- data.frame(
    age_start = c(15, 20), age_end = c(25, 30), rate = 0.1, region = "x"
  )
  expect_error(
    calculate_tfr(overlap),
    "nonoverlapping in stratum region = x \\("
  )
})

test_that("life tables follow piecewise constant hazards and handle limiting cases", {
  rates <- data.frame(
    age_start = c(0, 1, 5), age_end = c(1, 5, Inf), rate = 0.1
  )
  lt <- calculate_life_table(rates)
  expect_true(all(abs(lt$ex - 10) < 1e-12))
  expect_true(abs(lt$lx[3] - exp(-0.5)) < 1e-12)
  expect_true(abs(sum(lt$dx) - 1) < 1e-12)
  expect_true(abs(sum(lt$Lx) - 10) < 1e-12)
  scaled <- calculate_life_table(rates, radix = 1e+05)
  expect_equal(scaled$ex, lt$ex)
  expect_equal(scaled$lx, lt$lx * 1e+05)
  tiny <- rates
  tiny$rate <- 1e-12
  lt <- calculate_life_table(tiny)
  expect_true(all(abs(lt$ex/1e+12 - 1) < 1e-12))
  expect_true(abs(lt$ax[1] - 0.5) < 1e-10)
  zero <- rates
  zero$rate <- 0
  lt <- calculate_life_table(zero)
  expect_true(all(is.infinite(lt$ex)))
  expect_true(all(lt$qx == 0))
  expect_true(all(lt$lx == 1))
  expect_identical(lt$Lx[1:2], c(1, 4))
  missing <- rates
  missing$rate[2] <- NA_real_
  expect_true(all(is.na(calculate_life_table(missing)$ex)))
  closed <- rates
  closed$age_end[3] <- 100
  expect_error(calculate_life_table(closed))
  expect_error(calculate_life_table(rates[-1, ]))
  expect_error(calculate_life_table(rates, radix = 0))
  female <- rates
  female$sex <- "f"
  male <- rates
  male$sex <- "m"
  male$rate <- 0.2
  tables <- calculate_life_table(rbind(female, male))
  expect_true(all(abs(tables$ex[tables$sex == "f"] - 10) < 1e-12))
  expect_true(all(abs(tables$ex[tables$sex == "m"] - 5) < 1e-12))
  high <- data.frame(age_start = c(0, 1), age_end = c(1, Inf), rate = 40)
  expect_true(calculate_life_table(high)$lx[2] > 0)
  expect_true(all(abs(calculate_life_table(high)$ex - 1/40) < 1e-12))
})

test_that("empty populations have zero exposure and missing mortality rates", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  empty <- hammel_population(persons[FALSE, ], end_time = 1, groups = if ("group" %in% names(persons[FALSE, ])) "group")
  empty_rates <- demographic_counts(empty, c(0, 1), c(0, 1, Inf))
  expect_true(all(empty_rates$person_years == 0))
  expect_true(all(empty_rates$deaths == 0))
  expect_true(all(is.na(empty_rates$rate_mortality)))
})

test_that("Lexis summaries feed TFR and life tables directly", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  summary <- demographic_counts(population, c(0, 1), c(0, 1, 2, Inf))
  expect_equal(calculate_tfr(summary)$tfr, 4, tolerance = 0)
  expect_equal(nrow(calculate_life_table(summary)), nrow(summary), tolerance = 0)
  expect_identical(names(calculate_life_table(summary)), c("time_start", "time_end", "sex",
    "age_start", "age_end", "mx", "n", "ax", "qx", "lx", "lx_end", "dx", "Lx", "Tx",
    "ex"))
  expect_identical(calculate_tfr(summary), calculate_tfr(summary[is.finite(summary$age_end),
    ]))
})

test_that("finite age breaks append an open interval", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  expect_identical(demographic_counts(population, c(0, 1), c(0, 1)), demographic_counts(population,
    c(0, 1), c(0, 1, Inf)))
})

test_that("rate tables and typed processes give equivalent summaries", {
  input_fertility <- data.frame(age_start = c(15, 20, 30), age_end = c(20, 30, 50), rate = c(0.1,
    0.2, 0.05))
  expected <- input_fertility
  expected$rate_fertility <- expected$rate
  expected$rate <- NULL
  expect_identical(calculate_tfr(input_fertility), calculate_tfr(expected))
  expect_identical(calculate_tfr(fertility(input_fertility)), calculate_tfr(expected))
  input_mortality <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = c(0.1,
    0.2))
  expect_identical(calculate_life_table(mortality(input_mortality)), calculate_life_table(input_mortality))
  expect_error(calculate_tfr(mortality(input_mortality)))
  expect_error(calculate_life_table(fertility(input_fertility)))
})

test_that("population counts include births and exclude deaths at each instant", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  empty <- hammel_population(persons[FALSE, ], end_time = 1, groups = if ("group" %in% names(persons[FALSE, ])) "group")
  counts <- calculate_population_size(population, c(0, 0.25, 0.5, 0.75, 1))
  expect_identical(counts$population, c(1L, 2L, 3L, 2L, 2L))
  expect_identical(calculate_population_size(population, c(1, 0, 1))$population, c(2L,
    1L, 2L))
  expect_equal(calculate_population_size(population, -1)$population, 0L, tolerance = 0)
  expect_identical(calculate_population_size(persons, c(0, 1)), calculate_population_size(population,
    c(0, 1)))
  expect_equal(calculate_population_size(empty, 1)$population, 0L, tolerance = 0)
  expect_equal(nrow(calculate_population_size(population, numeric())), 0L, tolerance = 0)
  expect_error(calculate_population_size(population, 2))
  expect_error(calculate_population_size(population, NA_real_))
  expect_error(calculate_population_size(population, Inf))
})

test_that("all public demographic tables are tibbles including empty outputs", {
  persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-0.5, 0.25,
    0.5), death_time = c(0.75, NA, NA), mother_id = c(NA_integer_, NA_integer_, 1L),
    father_id = NA_integer_)
  population <- hammel_population(persons, end_time = 1, groups = if ("group" %in% names(persons)) "group")
  exposure <- demographic_counts(population, c(0, 1), c(0, 1, Inf))
  empty <- hammel_population(persons[FALSE, ], end_time = 1, groups = if ("group" %in% names(persons[FALSE, ])) "group")
  summary <- demographic_counts(population, c(0, 1), c(0, 1, 2, Inf))
  counts <- calculate_population_size(population, c(0, 0.25, 0.5, 0.75, 1))
  expect_true(tibble::is_tibble(counts))
  expect_true(tibble::is_tibble(summary))
  expect_true(tibble::is_tibble(attr(summary, "unassigned_births")))
  expect_true(tibble::is_tibble(calculate_tfr(summary)))
  expect_true(tibble::is_tibble(calculate_life_table(summary)))
  expect_true(tibble::is_tibble(calculate_population_size(empty, numeric())))
  expect_error(calculate_tfr(summary[FALSE, ]), "at least one age band")
  expect_error(
    calculate_life_table(summary[FALSE, ]),
    "at least one age band"
  )
})


test_that("TFR separates union statuses in input tables and fertility processes", {
  rates <- data.frame(union_status = c("unpartnered", "partnered", "unpartnered", "partnered"),
    age_start = c(15, 20, 20, 25), age_end = c(20, 25, 25, 30),
    rate = c(0.1, 0.3, 0.2, 0.4))
  expected <- tibble::tibble(union_status = c("unpartnered", "partnered"),
    age_start = c(15, 20), age_end = c(25, 30), tfr = c(1.5, 3.5))
  expect_equal(calculate_tfr(rates), expected)
  process_rates <- rates
  names(process_rates)[names(process_rates) == "union_status"] <- "union_type"
  expected_types <- expected
  names(expected_types)[names(expected_types) == "union_status"] <- "union_type"
  expect_equal(calculate_tfr(fertility(process_rates)), expected_types)
  later <- transform(rates, rate = rate * 2, time_start = 2005, time_end = 2010)
  first <- transform(rates, time_start = 2000, time_end = 2005)
  result <- calculate_tfr(rbind(first, later))
  expect_equal(result$tfr, c(1.5, 3.5, 3, 7))
  expect_equal(result$union_status, rep(c("unpartnered", "partnered"), 2))
  expect_equal(result$time_start, c(2000, 2000, 2005, 2005))
  expect_error(calculate_tfr(rates[FALSE, ]), "at least one age band")
  bad <- rates
  bad$age_start[3] <- 21
  expect_error(calculate_tfr(bad), "contiguous")
  rates$rate[3] <- NA_real_
  expect_equal(calculate_tfr(rates)$tfr, c(NA_real_, 3.5))
})

test_that("TFR warns only for known nonzero female open-age rates", {
  rates <- data.frame(age_start=c(15,20,50),age_end=c(20,50,Inf),rate=c(.1,.1,.01))
  expect_warning(x <- calculate_tfr(rates),"TFR.*underestimated")
  expect_equal(x$tfr,3.5)
  expect_equal(x$age_end,50)
  rates$rate[3] <- 0
  expect_warning(calculate_tfr(rates),NA)
  rates$rate[3] <- NA_real_
  expect_warning(calculate_tfr(rates),NA)
  rates$rate[3] <- .01
  rates$sex <- "m"
  female <- transform(rates, sex = "f", rate = c(.1, .1, 0))
  expect_warning(calculate_tfr(rbind(female, rates)), NA)
  expect_error(calculate_tfr(rates), "female age band")
  rates$sex <- "f"
  rates$union_type <- "partnered"
  expect_warning(calculate_tfr(rbind(rates,transform(rates,union_type="unpartnered"))),
    "TFR.*underestimated")
  only <- rates[3, ]
  expect_warning(
    expect_error(calculate_tfr(only), "finite end"),
    "TFR.*underestimated"
  )
})
