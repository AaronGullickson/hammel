grouped_model <- function(categories = c("A", "B")) {
  population <- generate_uniform_population(
    40,
    groups = list(
      race = stats::setNames(rep(1, length(categories)), categories)
    )
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

test_that("matching formulas that can't be evaluated are classified", {
  error <- expect_hammel_error(
    matching_model(
      ~ nonexistent(age_difference),
      c("nonexistent(age_difference)" = 1)
    ),
    "specification",
    "matching_model"
  )
  expect_identical(error$argument, "formula")
  expect_match(conditionMessage(error$parent), "nonexistent")
})

test_that("remaining inheritance and matching checks name their functions", {
  expect_hammel_error(
    group_inheritance("mother", order = c("A", "B")),
    "specification",
    "group_inheritance"
  )
  expect_hammel_error(
    group_inheritance("mother", fun = function(...) "A"),
    "specification",
    "group_inheritance"
  )
  expect_hammel_error(
    matching_model("age_difference", c(age_difference = 1)),
    "specification",
    "matching_model"
  )
  expect_hammel_error(
    matching_model(~1, c(a = 1)),
    "specification",
    "matching_model"
  )
  expect_hammel_error(
    matching_model(~ I(1 / age_difference), c("I(1/age_difference)" = 1)),
    "specification",
    "matching_model"
  )
  expect_hammel_error(
    matching_model(~1, numeric(), group_exogamy = c("A: " = 1)),
    "specification",
    "matching_model"
  )
  model <- matching_model(~ I(age_difference^2),
    c("I(age_difference^2)" = 1e308))
  expect_error(
    calculate_matching_utilities(model, 1e10),
    class = "hammel_error_runtime"
  )
})

test_that("multi-dimension matching needs two categories for female_ terms", {
  population <- generate_uniform_population(
    40,
    groups = list(race = c(A = 1, B = 1), region = c(N = 1))
  )
  matching <- matching_model(
    ~ female_region + same_race,
    c(same_race = 1)
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand", matching = matching)
  expect_hammel_error(run_simulation(model, 0), "groups", "run_simulation")
})

test_that("simulation checks union types against cohabitation settings", {
  population <- generate_uniform_population(20)
  women <- population$persons$person_id[population$persons$sex == "f"]
  men <- population$persons$person_id[population$persons$sex == "m"]
  unions <- tibble::tibble(
    union_id = 1L, wife_id = women[1], husband_id = men[1], start_time = -1,
    end_time = NA_real_, end_reason = NA_character_, union_type = "cohab"
  )
  population <- hammel_population(population$persons, unions, end_time = 0)
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("rate", data.frame(age_start = 15, age_end = Inf,
      rate = 0.1))
  expect_hammel_error(
    run_simulation(model, 1),
    "specification",
    "run_simulation"
  )
  typed <- rbind(
    transform(usa_2000_fertility, union_type = "unpartnered"),
    transform(usa_2000_fertility, union_type = "cohab"),
    transform(usa_2000_fertility, union_type = "marriage")
  )
  model <- hammel_model() |>
    set_population(generate_uniform_population(20)) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(typed) |>
    set_union_formation("rate", data.frame(age_start = 15, age_end = Inf,
      rate = 0.1))
  expect_hammel_error(
    run_simulation(model, 1),
    "specification",
    "run_simulation"
  )
  model <- set_population(model, generate_uniform_population(5,
    end_time = 1e300))
  expect_hammel_error(run_simulation(model, 1), "argument", "run_simulation")
})
