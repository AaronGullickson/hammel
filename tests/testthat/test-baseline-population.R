test_that("baseline populations are reproducible and follow the specified distributions",
  {
    set.seed(42)
    population <- create_baseline_population(50000, end_time = 2000)
    set.seed(42)
    expect_identical(population, create_baseline_population(50000, end_time = 2000))
    validate_population(population)
    p <- population$persons
    ages <- population$end_time - p$birth_time
    expect_identical(p$person_id, seq_len(50000))
    expect_true(all(ages >= 0 & ages <= 70))
    expect_true(all(is.na(p$mother_id)))
    expect_true(all(is.na(p$father_id)))
    expect_true(all(is.na(p$death_time)))
    expect_equal(nrow(population$unions), 0L, tolerance = 0)
    expect_identical(set_population(hammel_model(), population)$population, population)
    expect_true(abs(mean(p$sex == "f") - 0.5) < 6 * sqrt(0.25/nrow(p)))
    expect_true(abs(mean(ages) - 35) < 6 * sqrt(70^2/(12 * nrow(p))))
  })

test_that("baseline generator supports newborns, empty populations and fractional ages",
  {
    newborns <- create_baseline_population(3, max_age = 0, end_time = 2000, female_probability = 1)
    expect_identical(newborns$persons$sex, rep("f", 3))
    expect_identical(newborns$persons$birth_time, rep(2000, 3))
    expect_true(all(create_baseline_population(3, female_probability = 0)$persons$sex ==
      "m"))
    empty <- create_baseline_population(0)
    validate_population(empty)
    expect_equal(nrow(empty$persons), 0L, tolerance = 0)
    expect_equal(nrow(empty$unions), 0L, tolerance = 0)
    fractional <- create_baseline_population(10, max_age = 0.5, end_time = -1.25)
    expect_true(all(fractional$persons$birth_time <= -1.25))
    expect_true(all(fractional$persons$birth_time >= -1.75))
  })

test_that("baseline generator rejects invalid arguments", {
  for (size in list(-1, 1.5, NA_real_, Inf, "10", c(1, 2), 2^31)) {
    expect_error(create_baseline_population(size))
  }
  for (age in list(-1, NA_real_, Inf, "70", c(1, 2))) {
    expect_error(create_baseline_population(2, max_age = age))
  }
  for (time in list(NULL, NA_real_, Inf, "2000", c(1, 2))) {
    expect_error(create_baseline_population(2, end_time = time))
  }
  for (probability in list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0, 1))) {
    expect_error(create_baseline_population(2, female_probability = probability))
  }
})

