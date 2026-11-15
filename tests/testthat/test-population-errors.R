two_persons <- function() {
  tibble::tibble(
    person_id = c(1L, 2L),
    sex = c("f", "m"),
    birth_time = c(-30, -32),
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_
  )
}

test_that("population errors name the function the user called", {
  persons <- two_persons()
  persons$person_id[2] <- 1L
  error <- expect_hammel_error(
    hammel_population(persons, end_time = 0),
    "population",
    "hammel_population"
  )
  expect_identical(error$rows, 1:2)
  expect_hammel_error(hammel_population(two_persons()), "argument",
    "hammel_population")
  population <- hammel_population(two_persons(), end_time = 0)
  population$persons$birth_time[2] <- 5
  expect_hammel_error(
    set_population(hammel_model(), population),
    "population",
    "set_population"
  )
  expect_hammel_error(
    validate_population(population),
    "population",
    "validate_population"
  )
  expect_hammel_error(
    set_population(hammel_model(), list()),
    "population",
    "set_population"
  )
})

test_that("population relationship errors report their rows", {
  persons <- two_persons()
  persons$mother_id[2] <- 9L
  error <- expect_hammel_error(
    hammel_population(persons, end_time = 0),
    "population",
    "hammel_population"
  )
  expect_identical(error$rows, 2L)
  expect_identical(error$argument, "persons")
  unions <- tibble::tibble(
    union_id = 1L,
    wife_id = 2L,
    husband_id = 1L,
    start_time = -5,
    end_time = NA_real_,
    end_reason = NA_character_
  )
  error <- expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "population",
    "hammel_population"
  )
  expect_identical(error$argument, "unions")
  unions$wife_id <- 1L
  unions$husband_id <- 2L
  unions$end_reason <- "dissolution"
  error <- expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "population",
    "hammel_population"
  )
  expect_match(conditionMessage(error), "end_reason")
})

test_that("generator arguments are checked before drawing", {
  expect_hammel_error(
    generate_uniform_population(-1),
    "argument",
    "generate_uniform_population"
  )
  error <- expect_hammel_error(
    generate_uniform_population(10, max_age = 20, min_age = 30),
    "argument",
    "generate_uniform_population"
  )
  expect_identical(error$argument, "min_age")
  error <- expect_hammel_error(
    generate_uniform_population(10, groups = list(race = c(A = -1, B = 1))),
    "argument",
    "generate_uniform_population"
  )
  expect_identical(error$argument, "groups$race")
})

test_that("population messages keep their layout", {
  persons <- two_persons()
  persons$death_time[1] <- -40
  expect_snapshot(hammel_population(persons, end_time = 0), error = TRUE)
})

test_that("union end reasons and conversion links are checked", {
  unions <- tibble::tibble(
    union_id = 1L,
    wife_id = 1L,
    husband_id = 2L,
    start_time = -5,
    end_time = -1,
    end_reason = 1
  )
  expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "table",
    "hammel_population"
  )
  unions$end_reason <- "divorce"
  expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "table",
    "hammel_population"
  )
  unions$end_reason <- "conversion"
  expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "population",
    "hammel_population"
  )
  unions$union_type <- "cohab"
  expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "population",
    "hammel_population"
  )
  unions <- rbind(
    unions,
    tibble::tibble(
      union_id = 2L, wife_id = 1L, husband_id = 2L, start_time = -0.5,
      end_time = NA_real_, end_reason = NA_character_, union_type = "marriage"
    )
  )
  unions$previous_cohab_id <- c(NA, 1L)
  error <- expect_hammel_error(
    hammel_population(two_persons(), unions, end_time = 0),
    "population",
    "hammel_population"
  )
  expect_identical(error$rows, 2L)
})
