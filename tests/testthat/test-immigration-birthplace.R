# The automatic birthplace group of models with immigration (D190).

# A destination of 10 native couples who can have children, and an
# unpartnered template of women aged 25.
birthplace_model <- function(birth_rate = 0) {
  model <- make_model(20, birth_rate = birth_rate, birth_interval = 0)
  model$population$persons$sex <- rep(c("f", "m"), 10)
  model
}

birthplace_template <- function() {
  persons <- data.frame(
    person_id = 1:5, sex = "f", birth_time = -25, death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_
  )
  hammel_population(persons, end_time = 0)
}

test_that("birthplace is native for the born and the stream for arrivals", {
  model <- birthplace_model(birth_rate = 1) |>
    set_immigration(
      asia = immigration_stream(per_year(10), birthplace_template()),
      latam = immigration_stream(per_year(10), birthplace_template())
    )
  result <- run_simulation(model, 5, seed = 1)
  persons <- result$persons
  expect_identical(
    result$groups$birthplace$categories, c("native", "asia", "latam")
  )
  expect_true(result$groups$birthplace$automatic)
  expect_false(result$groups$birthplace$ancestry)
  expect_null(result$ancestry$birthplace)
  arrived <- !is.na(persons$arrival_time)
  expect_setequal(persons$birthplace[arrived], c("asia", "latam"))
  expect_true(all(persons$birthplace[!arrived] == "native"))
  # Immigrant mothers' children are native.
  kids <- persons[persons$mother_id %in% persons$person_id[arrived], ]
  expect_gt(nrow(kids), 0)
  expect_true(all(kids$birthplace == "native"))
  expect_false("stream" %in% names(persons))
})

test_that("birthplace ancestry is opt-in and passes to descendants", {
  model <- birthplace_model(birth_rate = 1) |>
    set_immigration(
      asia = immigration_stream(per_year(10), birthplace_template()),
      birthplace_ancestry = TRUE
    )
  result <- run_simulation(model, 5, seed = 2)
  shares <- result$ancestry$birthplace
  expect_named(shares, c("person_id", "native", "asia", "unknown"))
  persons <- result$persons
  row <- function(ids) shares[match(ids, shares$person_id), ]
  arrived <- persons[!is.na(persons$arrival_time), ]
  expect_true(all(row(arrived$person_id)$asia == 1))
  starting <- persons[1:20, ]
  expect_true(all(row(starting$person_id)$native == 1))
  # A newborn has the average of the parents' shares; an unknown father's
  # half is unknown.
  kids <- persons[!is.na(persons$mother_id) & persons$birth_time > 0, ]
  expect_true(any(kids$mother_id %in% arrived$person_id))
  mother <- row(kids$mother_id)
  father <- row(kids$father_id)
  child <- row(kids$person_id)
  known <- !is.na(kids$father_id)
  expect_equal(child$asia[known], (mother$asia[known] + father$asia[known]) / 2)
  expect_equal(child$asia[!known], mother$asia[!known] / 2)
  expect_equal(child$unknown[!known], rep(0.5, sum(!known)))
  expect_true(summary(immigration(
    a = immigration_stream(per_year(1), birthplace_template()),
    birthplace_ancestry = TRUE
  ))$birthplace_ancestry)
})

test_that("a template's own birthplace is ignored", {
  first <- run_simulation(
    set_immigration(
      birthplace_model(),
      asia = immigration_stream(per_year(10), birthplace_template())
    ),
    2, seed = 3
  )
  # Drawn from a template whose people are native or from asia.
  # Its entry rows are transitions that clones leave behind.
  expect_warning(
    stream <- immigration_stream(per_year(20), first, household = NULL),
    class = "hammel_warning_no_effect"
  )
  model <- set_immigration(birthplace_model(), europe = stream)
  result <- run_simulation(model, 2, seed = 4)
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  expect_gt(nrow(arrived), 0)
  expect_true(all(arrived$birthplace == "europe"))
  expect_identical(
    result$groups$birthplace$categories, c("native", "europe")
  )
})

test_that("birthplace persists through restarts and supplied values", {
  model <- set_immigration(
    birthplace_model(),
    asia = immigration_stream(per_year(10), birthplace_template())
  )
  first <- run_simulation(model, 2, seed = 5)
  # Without immigration, the dimension and its records continue.
  plain <- run_simulation(set_population(birthplace_model(), first), 2, 6)
  expect_identical(plain$persons$birthplace[seq_len(nrow(first$persons))],
                   first$persons$birthplace)
  expect_true(plain$groups$birthplace$automatic)
  # A new stream adds its category; earlier ones stay.
  again <- set_immigration(
    set_population(birthplace_model(), first),
    latam = immigration_stream(per_year(10), birthplace_template())
  )
  second <- run_simulation(again, 2, seed = 7)
  expect_identical(
    second$groups$birthplace$categories, c("native", "asia", "latam")
  )
  # A starting population may supply birthplaces.
  supplied <- birthplace_model()
  supplied$population$persons$birthplace <- rep(c("native", "asia"), 10)
  supplied <- set_immigration(
    supplied, latam = immigration_stream(per_year(5), birthplace_template())
  )
  result <- run_simulation(supplied, 1, seed = 8)
  expect_identical(
    result$persons$birthplace[1:20], rep(c("native", "asia"), 10)
  )
  expect_setequal(
    result$groups$birthplace$categories, c("native", "asia", "latam")
  )
})

test_that("rates can be stratified by birthplace", {
  # Only immigrants die.
  mortality <- data.frame(
    sex = rep(c("m", "f"), each = 2), age_start = 0, age_end = Inf,
    birthplace = c("native", "asia"), rate = c(0, 2)
  )
  model <- birthplace_model() |>
    set_immigration(
      asia = immigration_stream(per_year(20), birthplace_template())
    ) |>
    set_mortality(mortality)
  result <- run_simulation(model, 3, seed = 9)
  dead <- result$persons[!is.na(result$persons$death_time), ]
  expect_gt(nrow(dead), 0)
  expect_true(all(dead$birthplace == "asia"))
  rates <- calculate_mortality_rates(
    result, c(0, 3), c(0, Inf),
    by = "birthplace"
  )
  expect_identical(sum(rates$deaths[rates$birthplace == "native"]), 0)
})

test_that("birthplace never becomes the only dimension", {
  model <- birthplace_model()
  model$population$persons$race <- rep(c("A", "B"), each = 10)
  model <- model |>
    set_groups(race = group(categories = c("A", "B"))) |>
    set_immigration(
      asia = immigration_stream(
        per_year(10), birthplace_template(),
        assignment = list(race = migrant_assignment("fixed", fixed_value = "A"))
      )
    ) |>
    set_union_formation(
      "on_demand",
      matching = matching_model(
        ~age_difference, c(age_difference = 0),
        group_exogamy = c("A:B" = -1)
      ),
      cancel_birth_if_no_match = FALSE
    ) |>
    set_group_transition(
      data.frame(
        origin = "B", destination = "A", age_start = 0, age_end = Inf,
        rate = 0.5
      )
    )
  expect_identical(names(model$processes$group_transition), "race")
  result <- run_simulation(model, 2, seed = 10)
  expect_setequal(names(result$groups), c("race", "birthplace"))
  by_race <- calculate_population_size(result, c(0, 2), by = "race")
  expect_named(by_race, c("time", "race", "population"), ignore.order = TRUE)
  both <- calculate_population_size(
    result, c(0, 2),
    by = c("race", "birthplace")
  )
  expect_equal(
    sum(both$population[both$time == 2]),
    sum(by_race$population[by_race$time == 2])
  )
  rates <- calculate_transition_rates(result, c(0, 2), c(0, Inf))
  expect_true(all(rates$origin %in% c("A", "B")))
})

test_that("a transition attached before its groups skips birthplace", {
  # Without `group`, the transition resolves at run time to the only
  # declared dimension; the automatic birthplace doesn't count (D190).
  model <- birthplace_model() |>
    set_group_transition(
      data.frame(
        origin = "B", destination = "A", age_start = 0, age_end = Inf,
        rate = 0.5
      )
    )
  model$population$persons$race <- rep(c("A", "B"), each = 10)
  model <- model |>
    set_groups(race = group(categories = c("A", "B"))) |>
    set_immigration(
      asia = immigration_stream(
        per_year(10), birthplace_template(),
        assignment = list(race = migrant_assignment("fixed", fixed_value = "A"))
      )
    )
  result <- run_simulation(model, 2, seed = 10)
  expect_true(all(result$transitions$group[!is.na(
    result$transitions$origin
  )] == "race"))
})

test_that("birthplace can't have transitions", {
  model <- birthplace_model() |>
    set_immigration(
      asia = immigration_stream(per_year(10), birthplace_template()),
      europe = immigration_stream(per_year(10), birthplace_template())
    )
  for (destination in c("native", "europe")) {
    moving <- set_group_transition(
      model,
      data.frame(
        origin = "asia", destination = destination, age_start = 0,
        age_end = Inf, rate = 1
      ),
      group = "birthplace"
    )
    error <- expect_hammel_error(
      run_simulation(moving, 1, seed = 1),
      "groups",
      "run_simulation"
    )
    expect_identical(error$group, "birthplace")
  }
  # Without groups of its own, a transition without `group` can't fall back
  # to birthplace.
  implicit <- set_group_transition(
    model,
    data.frame(
      origin = "asia", destination = "native", age_start = 0, age_end = Inf,
      rate = 1
    )
  )
  expect_hammel_error(run_simulation(implicit, 1), "groups", "run_simulation")
})

test_that("birthplace names and settings are validated", {
  template <- birthplace_template()
  expect_hammel_error(
    immigration(native = immigration_stream(per_year(1), template)),
    "specification", "immigration"
  )
  expect_hammel_error(
    immigration(
      a = immigration_stream(per_year(1), template),
      birthplace_ancestry = 1
    ),
    "argument", "immigration"
  )
  stream <- immigration_stream(
    per_year(1), template,
    assignment = list(
      birthplace = migrant_assignment("fixed", fixed_value = "x")
    )
  )
  expect_hammel_error(
    run_simulation(set_immigration(birthplace_model(), a = stream), 1),
    "specification", "run_simulation"
  )
  model <- birthplace_model()
  model$population$persons$birthplace <- "native"
  model <- model |>
    set_groups(birthplace = group()) |>
    set_immigration(a = immigration_stream(per_year(1), template))
  expect_hammel_error(run_simulation(model, 1), "groups", "run_simulation")
})

test_that("redeclaring groups keeps birthplace", {
  model <- birthplace_model(birth_rate = 1)
  model$population$persons$race <- "A"
  model <- model |>
    set_groups(race = group(categories = "A", ancestry = TRUE)) |>
    set_immigration(
      asia = immigration_stream(
        per_year(10), birthplace_template(),
        assignment = list(race = migrant_assignment("fixed", fixed_value = "A"))
      ),
      birthplace_ancestry = TRUE
    )
  first <- run_simulation(model, 3, seed = 11)
  # A restart that redeclares the user's groups, without immigration.
  again <- set_population(birthplace_model(birth_rate = 1), first) |>
    set_groups(race = group(categories = "A", ancestry = TRUE))
  expect_true(again$population$groups$birthplace$automatic)
  second <- run_simulation(again, 2, seed = 12)
  expect_setequal(names(second$groups), c("race", "birthplace"))
  expect_setequal(names(second$ancestry), c("race", "birthplace"))
  old <- match(first$persons$person_id, second$ancestry$birthplace$person_id)
  expect_equal(
    second$ancestry$birthplace$asia[old], first$ancestry$birthplace$asia
  )
})
