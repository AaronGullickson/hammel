test_that("bundled USA rates are tibbles accepted by their constructors", {
  datasets <- list(
    usa_2000_mortality, usa_2000_fertility, usa_2000_fertility_union,
    usa_2000_marriage, usa_2000_migration
  )
  for (rates in datasets) {
    expect_s3_class(rates, "tbl_df")
    expect_true(all(is.finite(rates$rate) & rates$rate >= 0))
  }
  expect_s3_class(mortality(usa_2000_mortality), "hammel_mortality")
  expect_s3_class(fertility(usa_2000_fertility), "hammel_process")
  expect_s3_class(fertility(usa_2000_fertility_union), "hammel_process")
  expect_s3_class(
    union_formation("rate", rates = usa_2000_marriage),
    "hammel_process"
  )
  expect_s3_class(
    group_transition(usa_2000_migration, group = "region"),
    "hammel_process"
  )
})

test_that("bundled USA mortality reproduces the source life expectancy", {
  e0 <- calculate_life_table(usa_2000_mortality)
  e0 <- e0[e0$age_start == 0, ]
  # Published e0, 1999-2001: men 74.13, women 79.47. The abridged constant
  # hazards are within 0.1 years.
  expect_equal(e0$ex[e0$sex == "m"], 74.13, tolerance = 0.1 / 74.13)
  expect_equal(e0$ex[e0$sex == "f"], 79.47, tolerance = 0.1 / 79.47)
  expect_equal(
    calculate_tfr(usa_2000_fertility)$tfr,
    sum(usa_2000_fertility$rate * 5)
  )
})

test_that("England/Wales example data have the documented table schemas", {
  inputs <- england_wales_1992_onward
  expect_named(inputs, c(
    "fertility", "marriage", "premarital_cohabitation", "cohabitation_outcomes",
    "marital_divorce_death_2005", "marital_divorce_death_2010"
  ))
  expect_identical(
    unname(vapply(inputs, nrow, integer(1))), c(18L, 672L, 373L, 3L, 52L, 60L)
  )
  for (table in inputs) {
    expect_s3_class(table, "tbl_df")
  }
  expect_named(
    inputs$fertility, c("union_type", "age_start", "age_end", "rate")
  )
  expect_named(inputs$marriage, c(
    "time_start", "time_end", "age_start", "age_end", "previous_union", "rate"
  ))
  expect_named(inputs$premarital_cohabitation, c(
    "year", "age_groups", "prop_cohabit"
  ))
  expect_named(inputs$cohabitation_outcomes, c(
    "cohort_start", "cohort_end", "marriage_5", "separation_5", "continuing_5",
    "marriage_10", "separation_10", "continuing_10"
  ))
  for (year in c("2005", "2010")) {
    expect_named(inputs[[paste0("marital_divorce_death_", year)]], c(
      "year", "anniversary", "cumulative_divorce", "cumulative_death"
    ))
  }
})

test_that("cleaned rate tables preserve the observations without tuning", {
  rates <- england_wales_1992_onward$fertility
  expect_setequal(rates$union_type, c("marriage", "cohab", "unpartnered"))
  expect_identical(
    sort(unique(rates$age_start)), c(16L, 20L, 25L, 30L, 35L, 40L)
  )
  expect_identical(sort(unique(rates$age_end)), c(20L, 25L, 30L, 35L, 40L, 45L))
  expect_equal(rates$rate[rates$union_type == "cohab" & rates$age_start == 16],
    0.192
  )
  expect_s3_class(fertility(rates), "hammel_process")
  rates <- england_wales_1992_onward$marriage
  expect_setequal(rates$time_start, 1992:2019)
  expect_true(all(rates$time_end == rates$time_start + 1))
  expect_setequal(rates$previous_union, c("none", "dissolved"))
  expect_equal(sum(is.na(rates$rate)), 48L)
  expect_true(all(rates$previous_union[is.na(rates$rate)] == "dissolved"))
  first <- rates$time_start == 1992 & rates$age_start == 16 &
    rates$previous_union == "none"
  expect_equal(rates$rate[first], 14981 / 1207600)
  early <- rates[rates$time_start == 1992, ]
  expect_s3_class(union_formation("rate", rates = early), "hammel_process")
})

test_that("observation tables preserve source rounding and missing values", {
  inputs <- england_wales_1992_onward
  outcomes <- inputs$cohabitation_outcomes
  expect_identical(outcomes$cohort_start, c(1990L, 1995L, 2000L))
  expect_identical(outcomes$cohort_end, c(1995L, 2000L, 2005L))
  expect_equal(
    outcomes$marriage_5 + outcomes$separation_5 + outcomes$continuing_5,
    c(1, 1, 1.001)
  )
  expect_true(all(is.na(outcomes[3, c(
    "marriage_10", "separation_10", "continuing_10"
  )])))
  expect_true("All ages" %in% inputs$premarital_cohabitation$age_groups)
  expect_true("Age not stated" %in%
    inputs$premarital_cohabitation$age_groups)
  first <- inputs$marital_divorce_death_2005
  expect_identical(tail(first$anniversary, 2), c(55L, 60L))
  expect_equal(first$cumulative_divorce[1], 0.00022)
  expect_equal(first$cumulative_death[1], 0.00259)
  expect_identical(inputs$marital_divorce_death_2010$anniversary, 1:60)
})
