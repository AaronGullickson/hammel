test_that("model construction retains population and rate specifications", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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

})

test_that("population constructors normalize input subclasses to tibbles", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  custom_persons <- structure(persons, class = c("custom_persons", "data.frame"))
  expect_identical(set_population(empty, hammel_population(
    custom_persons,
    end_time = 0,
    unions,
    groups = if ("group" %in% names(custom_persons)) "group"
  ))$population$persons,
    tibble::as_tibble(custom_persons))
})

test_that("configuration rejects invalid objects and reference times", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  expect_error(hammel_population(
    NULL,
    end_time = 0,
    unions,
    groups = if ("group" %in% names(NULL)) "group"
  ))
  expect_error(hammel_population(
    persons,
    end_time = 0,
    list(),
    groups = if ("group" %in% names(persons)) "group"
  ))
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  unmarried <- set_population(configured, hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  ))
  expect_identical(unmarried, set_population(configured, hammel_population(
    persons,
    end_time = 0,
    unions = NULL,
    groups = if ("group" %in% names(persons)) "group"
  )))
  expect_identical(unmarried$population$persons, persons)
  expect_true(is.data.frame(unmarried$population$unions))
  expect_equal(nrow(unmarried$population$unions), 0L, tolerance = 0)
  expect_identical(unmarried$processes, configured$processes)
  expect_identical(configured$population$unions, unions)
})

test_that("fractional year endpoints are retained", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  check_rejected <- function(p = persons, u = unions) {
    expect_error(hammel_population(
      p,
      end_time = 0,
      u,
      groups = if ("group" %in% names(p)) "group"
    ))
    invalid <- structure(
      list(persons = p, unions = u, end_time = 0),
      class = "hammel_population"
    )
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  p <- persons
  p$custom_attribute <- c("a", "b")
  expect_identical(hammel_population(
    p,
    end_time = 0,
    groups = if ("group" %in% names(p)) "group"
  )$persons, p)
  expect_equal(nrow(hammel_population(
    persons[FALSE, ],
    end_time = 0,
    groups = if ("group" %in% names(persons[FALSE, ])) "group"
  )$persons), 0L, tolerance = 0)
})

test_that("population validation is read only and rejects malformed objects", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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
  expect_identical(hammel_population(
    persons,
    end_time = 0,
    data.frame(),
    groups = if ("group" %in% names(persons)) "group"
  )$unions, unions)
})

test_that("population attachment preserves optional metadata", {
  persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"), birth_time = c(-25,
    -23), death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
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
  unions <- hammel_population(
    persons,
    end_time = 0,
    groups = if ("group" %in% names(persons)) "group"
  )$unions
  mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
  fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)
  empty <- hammel_model()
  configured <- set_fertility(set_mortality(set_population(empty, hammel_population(persons,
    unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")), mortality), fertility)
  results <- hammel_population(persons, unions, end_time = 120, groups = if ("group" %in% names(persons)) "group")
  check_rejected <- function(p = persons, u = unions) {
    expect_error(hammel_population(
      p,
      end_time = 0,
      u,
      groups = if ("group" %in% names(p)) "group"
    ))
    invalid <- structure(
      list(persons = p, unions = u, end_time = 0),
      class = "hammel_population"
    )
    expect_error(set_population(empty, invalid))
  }
  for (sex in list(c(0L, 1L), c("m", NA_character_), c("male", "female"))) {
    p <- persons
    p$sex <- sex
    check_rejected(p)
  }
  expect_true(tibble::is_tibble(results$persons))
  expect_true(tibble::is_tibble(results$unions))
  expect_true(tibble::is_tibble(hammel_population(
    as.data.frame(persons),
    end_time = 0
  )$persons))
  expect_true(tibble::is_tibble(hammel_population(
    persons[FALSE, ],
    end_time = 0,
    groups = if ("group" %in% names(persons[FALSE, ])) "group"
  )$persons))
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
  expect_error(hammel_population(persons, swapped, end_time = 0, groups = if ("group" %in% names(persons)) "group"),
    "wife_id must refer to a woman", class = "hammel_error_population")
  swapped$end_time <- -0.5
  swapped$end_reason <- "dissolution"
  expect_error(hammel_population(persons, swapped, end_time = 0, groups = if ("group" %in% names(persons)) "group"),
    "wife_id must refer to a woman", class = "hammel_error_population")
  expect_identical(names(generate_uniform_population(0)$unions), names(unions))
})


test_that("population reference times are required and validated everywhere", {
  population <- generate_uniform_population(2, end_time = 0)
  persons <- population$persons
  expect_error(hammel_population(persons), "end_time.*required")
  bad_times <- list(
    NULL, NA_real_, NaN, Inf, -Inf, numeric(), c(0, 1), "0", TRUE,
    0 + 0i, matrix(0), array(0)
  )
  for (time in bad_times) {
    error <- expect_error(
      hammel_population(persons, end_time = time),
      class = "hammel_error_argument"
    )
    expect_identical(error$argument, "end_time")
    invalid <- population
    invalid["end_time"] <- list(time)
    expect_error(validate_population(invalid), "end_time")
    expect_error(set_population(hammel_model(), invalid), "end_time")
  }
  for (time in c(-0.5, 0, 0.5, 2000)) {
    persons$birth_time <- rep(time - 25, nrow(persons))
    valid <- hammel_population(persons, end_time = time)
    expect_identical(valid$end_time, time)
    expect_invisible(validate_population(valid))
  }
  expect_error(
    hammel_population(persons, end_time = 1900),
    "Recorded events must not occur after end_time",
    class = "hammel_error_population"
  )
})

test_that("model summaries combine components and report readiness", {
  empty <- summary(hammel_model())
  expect_s3_class(empty, "summary.hammel_model")
  expect_null(empty$population)
  expect_length(empty$processes, 0L)
  expect_identical(
    empty$missing,
    c("population", "mortality", "fertility", "union_formation")
  )
  text <- cli::cli_fmt(value <- print(empty))
  expect_identical(value, empty)
  expect_true(any(grepl("needs the following", text, fixed = TRUE)))
  expect_true(any(grepl("set_union_formation()", text, fixed = TRUE)))

  model <- hammel_model() |>
    set_population(generate_uniform_population(50)) |>
    set_union_formation("on_demand") |>
    set_mortality(usa_2000_mortality)
  partial <- summary(model)
  expect_s3_class(partial$population, "summary.hammel_population")
  expect_named(partial$processes, c("mortality", "union_formation"))
  expect_identical(partial$missing, "fertility")
  expect_error(run_simulation(model, years = 1), "fertility")

  ready <- summary(set_fertility(model, usa_2000_fertility))
  expect_named(
    ready$processes,
    c("mortality", "fertility", "union_formation")
  )
  expect_length(ready$missing, 0L)
  text <- cli::cli_fmt(print(ready))
  expect_true(any(grepl("Fertility process", text, fixed = TRUE)))
  expect_true(any(grepl("elements needed to run", text, fixed = TRUE)))
})

test_that("model summaries combine group declarations and transitions", {
  population <- generate_uniform_population(
    60,
    groups = list(ethnicity = c(A = 1, B = 1), region = c(n = 1, s = 1))
  )
  rates <- tibble::tribble(
    ~age_start, ~age_end, ~origin, ~destination, ~rate,
    0,          Inf,      "n",     "s",          0.05
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_groups(ethnicity = group(ancestry = TRUE), region = group()) |>
    set_group_transition(rates, group = "region")
  result <- summary(model)
  expect_named(result$groups, c("ethnicity", "region"))
  expect_identical(result$groups$ethnicity$group$categories, c("A", "B"))
  expect_true(result$groups$ethnicity$group$ancestry)
  expect_null(result$groups$ethnicity$transition)
  expect_s3_class(
    result$groups$region$transition,
    "summary.hammel_group_transition"
  )
  expect_named(result$population$group_counts, c("ethnicity", "region"))
  text <- cli::cli_fmt(print(result))
  groups <- grep("Groups", text, fixed = TRUE)
  expect_gt(groups, grep("Living population by age", text, fixed = TRUE))
  expect_true(any(grepl("Transition process: none", text, fixed = TRUE)))
  expect_s3_class(result$groups$region$group, "summary.hammel_group")
  expect_true(any(grepl("Inheritance: random_parent", text, fixed = TRUE)))
  expect_true(any(grepl("(1 of 2 possible): n", text, fixed = TRUE)))

  expect_null(summary(hammel_model())$groups)
})

test_that("model readiness requires declared transition dimensions", {
  population <- generate_uniform_population(
    60,
    groups = list(region = c(n = 1, s = 1, w = 1))
  )
  rates <- tibble::tribble(
    ~age_start, ~age_end, ~origin, ~destination, ~rate,
    0,          Inf,      "n",     "s",          0.05
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand") |>
    set_group_transition(rates)
  result <- summary(model)
  expect_length(result$missing, 0L)
  expect_equal(nrow(result$groups$region$transition$missing_transitions), 5L)
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl("(1 of 6 possible)", text, fixed = TRUE)))
  expect_no_error(run_simulation(model, years = 1, seed = 1))

  model$processes$group_transition$ethnicity <- group_transition(
    rates,
    group = "ethnicity"
  )
  result <- summary(model)
  expect_identical(result$missing, "transition:ethnicity")
  expect_null(result$groups$ethnicity$group)
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl("dimension ethnicity", text, fixed = TRUE)))
  expect_error(run_simulation(model, years = 1), "ethnicity")

  unassigned <- hammel_model() |> set_group_transition(rates)
  expect_identical(
    summary(unassigned)$missing,
    c(missing_model_elements(unassigned), "transition:.default")
  )
  text <- cli::cli_fmt(print(summary(unassigned)))
  expect_true(any(grepl("Unassigned dimension", text, fixed = TRUE)))
  expect_true(any(grepl("set_group_transition()", text, fixed = TRUE)))
})

test_that("population summaries match the simulated histories", {
  set.seed(7)
  population <- generate_uniform_population(300)
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand")
  result <- run_simulation(model, years = 10, seed = 11)
  persons <- result$persons
  alive <- is.na(persons$death_time)
  out <- summary(result)

  expect_s3_class(out, "summary.hammel_population")
  expect_identical(out$first_time, min(persons$birth_time))
  expect_identical(out$end_time, 10)
  expect_identical(out$persons_alive, sum(alive))
  expect_identical(out$n_births, nrow(persons))
  expect_identical(out$n_deaths, sum(!alive))
  expect_identical(out$n_unions, nrow(result$unions))
  expect_identical(
    out$n_dissolutions,
    sum(result$unions$end_reason %in% "dissolution")
  )
  expect_identical(
    is.na(out$n_births_u),
    is.null(result$unrealized_births)
  )
  expect_s3_class(out$age_counts, "tbl_df")
  expect_equal(sum(out$age_counts$male), sum(alive & persons$sex == "m"))
  expect_equal(sum(out$age_counts$female), sum(alive & persons$sex == "f"))
  expect_identical(out$simulation$start, 0)
  expect_identical(out$simulation$end, 10)
  expect_identical(out$simulation$pop_size_start, 300L)
  expect_identical(out$simulation$pop_size_end, sum(alive))
  expect_null(out$group_counts)

  text <- cli::cli_fmt(value <- print(out))
  expect_identical(value, out)
  expect_true(any(grepl("Latest simulation run", text, fixed = TRUE)))
  expect_true(any(grepl("Time period: 0-10", text, fixed = TRUE)))
  expect_true(any(grepl(
    paste("Persons alive at end time:", sum(alive)), text, fixed = TRUE
  )))
  expect_identical(
    any(grepl("Total unrealized births", text, fixed = TRUE)),
    !is.na(out$n_births_u)
  )

  start <- summary(population)
  expect_null(start$simulation)
  expect_true(is.na(start$n_births_u))
  text <- cli::cli_fmt(print(start))
  expect_false(any(grepl("Latest simulation run", text, fixed = TRUE)))
  expect_false(any(grepl("unrealized", text, fixed = TRUE)))
})

test_that("population summaries count declared categories without members", {
  population <- generate_uniform_population(
    40,
    groups = list(region = c(n = 1, s = 1))
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_groups(region = group(categories = c("n", "s", "w")))
  counts <- summary(model$population)$group_counts$region
  expect_s3_class(counts, "tbl_df")
  expect_identical(counts$category, c("n", "s", "w"))
  expect_equal(sum(counts$count), 40)
  expect_equal(counts$count[3], 0)
  expect_equal(sum(counts$percent), 100)
  text <- cli::cli_fmt(print(summary(model$population)))
  expect_true(any(grepl("w 0 (0.0%)", text, fixed = TRUE)))
})
