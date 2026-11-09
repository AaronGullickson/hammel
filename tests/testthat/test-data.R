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
