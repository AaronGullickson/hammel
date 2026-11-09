test_that("model construction retains population and rate specifications", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  expect_true(inherits(configured, "hammel_model"))
  expect_true(is.null(empty$population))
  expect_identical(empty$processes, list())
  expect_identical(configured$population$persons, persons)
  expect_identical(configured$population$unions, unions)
  expect_identical(configured$population$end_time, 0)
  expect_identical(configured$processes$mortality$rates, mortality)
  expect_identical(configured$processes$fertility$rates, fertility)
  expect_identical(set_mortality(empty, mortality), set_process(empty, hammel::mortality(mortality)))
  expect_identical(set_fertility(empty, fertility), set_process(empty, hammel::fertility(fertility)))
})

test_that("replacing a process preserves the original model and other processes", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  replacement <- mortality
  replacement$rate <- 0.02
  updated <- set_mortality(configured, replacement)
  expect_identical(configured$processes$mortality$rates, mortality)
  expect_identical(updated$processes$mortality$rates, replacement)
  expect_identical(updated$processes$fertility$rates, fertility)
  expect_equal(length(updated$processes), 2L, tolerance = 0)
})

test_that("population restarts preserve histories, endpoints and processes", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  restarted <- set_population(configured, results)
  expect_true(inherits(results, "hammel_population"))
  expect_identical(restarted$population$persons, persons)
  expect_identical(restarted$population$unions, unions)
  expect_identical(restarted$population$end_time, 120)
  expect_identical(restarted$processes, configured$processes)
  expect_identical(configured$population$end_time, 0)
  unknown <- set_population(restarted, hammel_population(persons, unions, groups = if ("group" %in% names(persons)) "group"))
  expect_true("end_time" %in% names(unknown$population))
  expect_true(is.null(unknown$population$end_time))
})

test_that("population constructors normalize input subclasses to tibbles", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  custom_persons <- structure(persons, class = c("custom_persons", "data.frame"))
  expect_identical(set_population(empty, hammel_population(custom_persons, unions, groups = if ("group" %in% names(custom_persons)) "group"))$population$persons,
    tibble::as_tibble(custom_persons))
})

test_that("configuration rejects invalid objects and reference times", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  expect_error(hammel_population(NULL, unions, groups = if ("group" %in% names(NULL)) "group"))
  expect_error(hammel_population(persons, list(), groups = if ("group" %in% names(persons)) "group"))
  expect_error(set_population(empty, persons, list()))
  expect_error(set_population(empty, results, unions))
  expect_error(set_population(empty, results, start_time = 0))
  expect_error(set_process(list(), hammel::mortality(mortality)))
  expect_error(set_process(empty, mortality))
  expect_error(set_process(empty, NULL))
  for (month in list(NA_real_, Inf, c(0, 1), "120")) {
    expect_error(hammel_population(persons, unions, end_time = month, groups = if ("group" %in% names(persons)) "group"))
    expect_error(set_population(empty, persons, unions, start_time = month))
  }
})

test_that("omitted unions produce an empty union history", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  unmarried <- set_population(configured, hammel_population(persons, groups = if ("group" %in% names(persons)) "group"))
  expect_identical(unmarried, set_population(configured, hammel_population(persons, unions = NULL, groups = if ("group" %in% names(persons)) "group")))
  expect_identical(unmarried$population$persons, persons)
  expect_true(is.data.frame(unmarried$population$unions))
  expect_equal(nrow(unmarried$population$unions), 0L, tolerance = 0)
  expect_identical(unmarried$processes, configured$processes)
  expect_identical(configured$population$unions, unions)
})

test_that("fractional year endpoints are retained", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  expect_identical(set_population(empty, hammel_population(persons, end_time = 1.5, groups = if ("group" %in% names(persons)) "group"))$population$end_time,
    1.5)
})

test_that("population construction and attachment reject inconsistent histories", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  check_rejected <- function(p = persons, u = unions) {
    expect_error(hammel_population(p, u, groups = if ("group" %in% names(p)) "group"))
    invalid <- structure(list(persons = p, unions = u, end_time = NULL), class = "hammel_population")
    expect_error(set_population(empty, invalid))
  }
  check_rejected(persons[setdiff(names(persons), "sex")])
  p <- persons
  p$person_id[2] <- p$person_id[1]
  check_rejected(p)
  p <- persons
  p$person_id[1] <- 1.5
  check_rejected(p)
  p <- persons
  p$sex[1] <- 2
  check_rejected(p)
  p <- persons
  p$birth_time[1] <- Inf
  check_rejected(p)
  p <- persons
  p$death_time[1] <- -30
  check_rejected(p)
  p <- persons
  p$mother_id[1] <- 999L
  check_rejected(p)
  p <- persons
  p$father_id[1] <- 1L
  check_rejected(p)
  u <- data.frame(union_id = 1L, wife_id = 2L, husband_id = 1L, start_time = -1, end_time = NA_real_,
    end_reason = NA_character_)
  expect_identical(hammel_population(persons, u, end_time = 0, groups = if ("group" %in% names(persons)) "group")$unions, tibble::as_tibble(u))
  v <- u
  v$husband_id <- 999L
  check_rejected(u = v)
  v <- u
  v$husband_id <- 2L
  check_rejected(u = v)
  v <- u
  v$end_time <- -2
  v$end_reason <- "dissolution"
  check_rejected(u = v)
  v <- u
  v$end_reason <- "death_husband"
  check_rejected(u = v)
  v <- u
  v$start_time <- -30
  check_rejected(u = v)
  check_rejected(u = rbind(u, u))
  expect_error(hammel_population(persons, u, end_time = -2, groups = if ("group" %in% names(persons)) "group"))
  expect_error(set_population(empty, persons, u, start_time = -2))
})

test_that("custom columns and empty populations retain their structure", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  p <- persons
  p$custom_attribute <- c("a", "b")
  expect_identical(hammel_population(p, groups = if ("group" %in% names(p)) "group")$persons, p)
  expect_equal(nrow(hammel_population(persons[FALSE, ], groups = if ("group" %in% names(persons[FALSE, ])) "group")$persons), 0L, tolerance = 0)
})

test_that("population validation is read only and rejects malformed objects", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  checked <- withVisible(validate_population(results))
  expect_identical(checked$value, results)
  expect_true(!checked$visible)
  invalid <- results
  invalid$unions <- NULL
  expect_error(validate_population(invalid))
  expect_error(set_population(empty, invalid))
  invalid <- results
  invalid$end_time <- Inf
  expect_error(validate_population(invalid))
  invalid <- results
  invalid$persons$person_id[2] <- 1L
  expect_error(validate_population(invalid))
  expect_error(validate_population(persons))
  expect_identical(hammel_population(persons, data.frame(), groups = if ("group" %in% names(persons)) "group")$unions, unions)
})

test_that("population attachment preserves optional metadata", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  annotated <- results
  annotated$observations <- list(description = "preserve supplied components")
  stored <- set_population(empty, annotated)
  expect_identical(stored$population, annotated)
  expect_identical(names(stored), c("population", "processes"))
  expect_identical(set_population(empty, hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group"))$population,
    results)
  expect_error(set_population(empty, persons))
})

test_that("sex labels are explicit and public population tables are tibbles", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(persons, groups = if ("group" %in% names(persons)) "group")$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  check_rejected <- function(p = persons, u = unions) {
    expect_error(hammel_population(p, u, groups = if ("group" %in% names(p)) "group"))
    invalid <- structure(list(persons = p, unions = u, end_time = NULL), class = "hammel_population")
    expect_error(set_population(empty, invalid))
  }
  for (sex in list(c(0L, 1L), c("m", NA_character_), c("male", "female"))) {
    p <- persons
    p$sex <- sex
    check_rejected(p)
  }
  expect_true(tibble::is_tibble(results$persons))
  expect_true(tibble::is_tibble(results$unions))
  expect_true(tibble::is_tibble(hammel_population(as.data.frame(persons))$persons))
  expect_true(tibble::is_tibble(hammel_population(persons[FALSE, ], groups = if ("group" %in% names(persons[FALSE, ])) "group")$persons))
})


test_that("union wife and husband columns enforce female and male roles", {
  persons <- data.frame(person_id = c(40, 10), sex = c("m", "f"),
    birth_time = c(-30, -25), death_time = NA_real_,
    mother_id = NA_real_, father_id = NA_real_)
  unions <- data.frame(union_id = 3, wife_id = 10, husband_id = 40,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  population <- hammel_population(persons, unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")
  expect_equal(population$unions$wife_id, 10)
  expect_equal(population$unions$husband_id, 40)
  swapped <- transform(unions, wife_id = 40, husband_id = 10)
  expect_error(hammel_population(persons, swapped, end_time = 0, groups = if ("group" %in% names(persons)) "group"), "wife_id.*female")
  swapped$end_time <- -0.5
  swapped$end_reason <- "dissolution"
  expect_error(hammel_population(persons, swapped, end_time = 0, groups = if ("group" %in% names(persons)) "group"), "wife_id.*female")
  expect_identical(names(generate_uniform_population(0)$unions), names(unions))
})
