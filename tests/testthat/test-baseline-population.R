test_that("baseline populations are reproducible and follow the specified distributions",
  {
    set.seed(42)
    population <- generate_uniform_population(50000, end_time = 2000)
    set.seed(42)
    expect_identical(population, generate_uniform_population(50000, end_time =
      2000))
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
    newborns <- generate_uniform_population(3, max_age = 0, end_time = 2000,
      female_probability = 1)
    expect_identical(newborns$persons$sex, rep("f", 3))
    expect_identical(newborns$persons$birth_time, rep(2000, 3))
    expect_true(all(generate_uniform_population(3, female_probability =
      0)$persons$sex ==
      "m"))
    empty <- generate_uniform_population(0)
    validate_population(empty)
    expect_equal(nrow(empty$persons), 0L, tolerance = 0)
    expect_equal(nrow(empty$unions), 0L, tolerance = 0)
    fractional <- generate_uniform_population(10, max_age = 0.5, end_time = -
      1.25)
    expect_true(all(fractional$persons$birth_time <= -1.25))
    expect_true(all(fractional$persons$birth_time >= -1.75))
  })

test_that("baseline generator rejects invalid arguments", {
  for (size in list(-1, 1.5, NA_real_, Inf, "10", c(1, 2), 2^31)) {
    expect_error(generate_uniform_population(size))
  }
  for (age in list(-1, NA_real_, Inf, "70", c(1, 2))) {
    expect_error(generate_uniform_population(2, max_age = age))
  }
  for (time in list(NULL, NA_real_, Inf, "2000", c(1, 2))) {
    expect_error(generate_uniform_population(2, end_time = time))
  }
  for (probability in list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0, 1))) {
    expect_error(generate_uniform_population(2, female_probability =
      probability))
  }
})


test_that("uniform ages respect positive fractional bounds", {
  set.seed(42)
  population <- generate_uniform_population(
    50000, min_age = 18.25, max_age = 65.75, end_time = 2000
  )
  ages <- population$end_time - population$persons$birth_time
  width <- 65.75 - 18.25
  expect_true(all(ages >= 18.25 & ages <= 65.75))
  expect_lt(abs(mean(ages) - 42), 6 * width / sqrt(12 * length(ages)))
  expect_lt(
    abs(var(ages) - width^2 / 12),
    6 * width^2 / sqrt(180 * length(ages))
  )
  expect_true(all(is.na(population$persons$death_time)))
  expect_equal(nrow(population$unions), 0L)
  expect_invisible(validate_population(population))
})

test_that("equal age bounds and empty populations are supported", {
  fixed <- generate_uniform_population(
    3, min_age = 25.5, max_age = 25.5, end_time = 2000
  )
  expect_identical(fixed$persons$birth_time, rep(1974.5, 3))
  empty <- generate_uniform_population(0, min_age = 18, max_age = 65)
  expect_equal(nrow(empty$persons), 0L)
  expect_s3_class(empty$persons, "tbl_df")
  expect_invisible(validate_population(empty))
})

test_that("zero minimum age preserves output and RNG state", {
  weights <- list(region = c(A = 1, B = 2))
  set.seed(42)
  positional <- generate_uniform_population(
    100, 50, 2000, female_probability = 0.6, groups = weights
  )
  rng <- .Random.seed
  set.seed(42)
  explicit <- generate_uniform_population(
    size = 100,
    max_age = 50,
    end_time = 2000,
    female_probability = 0.6,
    groups = weights,
    min_age = 0
  )
  expect_identical(explicit, positional)
  expect_identical(.Random.seed, rng)
})

test_that("invalid age bounds fail before consuming random numbers", {
  invalid <- list(
    NULL, numeric(), -1, NA_real_, NaN, Inf, -Inf, "18", TRUE,
    c(1, 2), 0 + 0i, matrix(0), array(0)
  )
  for (age in invalid) {
    set.seed(42)
    rng <- .Random.seed
    expect_error(generate_uniform_population(2, min_age = age), "min_age")
    expect_identical(.Random.seed, rng)
    expect_error(generate_uniform_population(2, max_age = age), "max_age")
    expect_identical(.Random.seed, rng)
  }
  expect_error(
    generate_uniform_population(0, min_age = 20, max_age = 19),
    "min_age.*must not exceed.*max_age"
  )
})


test_that("minimum age immediately follows the reference time", {
  expect_identical(
    names(formals(generate_uniform_population)),
    c("size", "max_age", "end_time", "min_age", "female_probability", "groups")
  )
  set.seed(42)
  positional <- generate_uniform_population(10, 50, 2000, 18, 1)
  set.seed(42)
  named <- generate_uniform_population(
    10, max_age = 50, end_time = 2000, min_age = 18, female_probability = 1
  )
  expect_identical(positional, named)
  expect_true(all(positional$persons$sex == "f"))
})
