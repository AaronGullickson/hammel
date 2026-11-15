grouped_model <- function(categories = c("A", "B")) {
  population <- generate_uniform_population(
    40,
    groups = list(race = stats::setNames(rep(1, length(categories)), categories))
  )
  hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand")
}

test_that("group declaration errors name the function the user called", {
  expect_hammel_error(group(ancestry = NA), "argument", "group")
  expect_hammel_error(
    group(categories = c("A", "A")),
    "groups",
    "group"
  )
  expect_hammel_error(group(inheritance = "bogus"), "specification", "group")
  error <- expect_hammel_error(
    set_groups(grouped_model(), region = group()),
    "groups",
    "set_groups"
  )
  expect_identical(error$group, "region")
  error <- expect_hammel_error(
    set_groups(grouped_model(), sex = group()),
    "groups",
    "set_groups"
  )
  expect_identical(error$group, "sex")
  error <- expect_hammel_error(
    set_groups(grouped_model(), race = group(categories = "A")),
    "groups",
    "set_groups"
  )
  expect_identical(error$categories, "B")
  model <- set_groups(grouped_model(), race = group(inheritance = "ancestry"))
  error <- expect_hammel_error(
    run_simulation(model, 1),
    "specification",
    "run_simulation"
  )
  expect_identical(error$group, "race")
})

test_that("inheritance errors name group_inheritance or run_simulation", {
  expect_hammel_error(
    group_inheritance("hierarchy"),
    "specification",
    "group_inheritance"
  )
  expect_hammel_error(
    group_inheritance("mother", passing = 0.2),
    "specification",
    "group_inheritance"
  )
  model <- set_groups(
    grouped_model(),
    race = group(inheritance = group_inheritance("fixed", fixed_value = "C"))
  )
  error <- expect_hammel_error(run_simulation(model, 1), "groups",
    "run_simulation")
  expect_identical(error$categories, "C")
  rule <- group_inheritance(
    "custom",
    fun = function(mother_group, ...) {
      "Z"
    }
  )
  model <- set_groups(
    set_fertility(grouped_model(), data.frame(age_start = 0, age_end = Inf,
      rate = 5)),
    race = group(inheritance = rule)
  )
  error <- expect_hammel_error(
    run_simulation(model, 1, seed = 1),
    "runtime",
    "run_simulation"
  )
  expect_identical(error$categories, "Z")
})

test_that("matching errors name the constructor or run_simulation", {
  error <- expect_hammel_error(
    matching_model(~age_difference, c(age_differenc = 1)),
    "specification",
    "matching_model"
  )
  expect_identical(error$argument, "coefficients")
  expect_hammel_error(
    age_gap_matching_model(0, -1),
    "argument",
    "age_gap_matching_model"
  )
  expect_hammel_error(
    matching_model(~1, numeric(), group_exogamy = c("A:A" = NA)),
    "specification",
    "matching_model"
  )
  expect_hammel_error(
    union_formation("on_demand", matching = "close"),
    "specification",
    "union_formation"
  )
  matching <- matching_model(
    ~ female_race,
    c(female_raceB = 1, female_raceC = 1)
  )
  model <- set_union_formation(grouped_model(), "on_demand",
    matching = matching)
  error <- expect_hammel_error(
    run_simulation(model, 0),
    "specification",
    "run_simulation"
  )
  expect_match(conditionMessage(error), "In the matching model.",
    fixed = TRUE)
  matching <- matching_model(~ same_region, c(same_region = 1))
  model <- set_union_formation(grouped_model(), "on_demand",
    matching = matching)
  error <- expect_hammel_error(
    run_simulation(model, 0),
    "groups",
    "run_simulation"
  )
  expect_identical(error$group, "region")
})

test_that("simulation setup errors name run_simulation", {
  model <- grouped_model()
  expect_hammel_error(run_simulation(model, -1), "argument", "run_simulation")
  expect_hammel_error(
    run_simulation(model, 1, seed = 1.5),
    "argument",
    "run_simulation"
  )
  expect_hammel_error(
    run_simulation(model, 1, progress = NA),
    "argument",
    "run_simulation"
  )
  error <- expect_hammel_error(
    run_simulation(hammel_model(), 1),
    "specification",
    "run_simulation"
  )
  expect_identical(
    error$missing,
    c("population", "mortality", "fertility", "union_formation")
  )
  rates <- usa_2000_mortality
  rates$age_end[rates$age_end == Inf] <- 110
  model <- set_mortality(grouped_model(), rates)
  expect_hammel_error(run_simulation(model, 1), "coverage", "run_simulation")
})

test_that("group and matching messages keep their layout", {
  expect_snapshot(
    set_groups(grouped_model(), race = group(categories = "A")),
    error = TRUE
  )
  matching <- matching_model(~ female_race, c(female_raceC = 1))
  model <- set_union_formation(grouped_model(), "on_demand",
    matching = matching)
  expect_snapshot(run_simulation(model, 0), error = TRUE)
  expect_snapshot(run_simulation(hammel_model(), 1), error = TRUE)
})
