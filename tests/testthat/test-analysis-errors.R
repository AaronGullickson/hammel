simulated_population <- function() {
  withr::local_seed(11)
  model <- hammel_model() |>
    set_population(generate_uniform_population(30)) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand")
  run_simulation(model, 2, seed = 3)
}

test_that("rate calculators name themselves through shared helpers", {
  population <- simulated_population()
  expect_hammel_error(
    calculate_mortality_rates(population, c(0, 3), c(0, Inf)),
    "argument",
    "calculate_mortality_rates"
  )
  error <- expect_hammel_error(
    calculate_fertility_rates(population, c(0, 2), c(0, Inf), by = "sex"),
    "argument",
    "calculate_fertility_rates"
  )
  expect_identical(error$variables, "sex")
  expect_hammel_error(
    calculate_marriage_rates(population, c(0, 3), c(0, Inf)),
    "argument",
    "calculate_marriage_rates"
  )
  expect_hammel_error(
    calculate_union_rates(population, c(0, 2), c(0, Inf), conf_level = 2),
    "argument",
    "calculate_union_rates"
  )
  expect_hammel_error(
    calculate_dissolution_rates(population, c(0, 2), c(0, Inf),
      by = "group_w"),
    "groups",
    "calculate_dissolution_rates"
  )
  expect_hammel_error(
    calculate_population_size(population$persons, 0),
    "population",
    "calculate_population_size"
  )
  expect_hammel_error(
    calculate_population_size(population, 0, by = NA),
    "argument",
    "calculate_population_size"
  )
  expect_hammel_error(count_kin(population, 5), "argument", "count_kin")
  expect_hammel_error(extract_kin(population, 0), "table", "extract_kin")
})

test_that("summary measures name themselves for per-stratum problems", {
  rates <- data.frame(
    age_start = c(15, 20, 15, 22),
    age_end = c(20, 25, 20, 25),
    rate = 0.1,
    region = c("a", "a", "b", "b")
  )
  error <- expect_hammel_error(calculate_tfr(rates), "rates", "calculate_tfr")
  expect_match(conditionMessage(error), "region = b", fixed = TRUE)
  rates <- data.frame(age_start = c(0, 1), age_end = c(1, 50), rate = 0.1)
  expect_hammel_error(
    calculate_life_table(rates),
    "coverage",
    "calculate_life_table"
  )
  expect_hammel_error(
    calculate_life_expectancy(usa_2000_mortality, age = 3),
    "argument",
    "calculate_life_expectancy"
  )
  expect_warning(
    calculate_tfr(
      data.frame(age_start = c(15, 50), age_end = c(50, Inf), rate = 0.1)
    ),
    class = "hammel_warning_estimate"
  )
})

test_that("stratification errors name stratify_rates", {
  rates <- data.frame(age_start = 0, age_end = Inf, rate = 0.1)
  expect_hammel_error(
    stratify_rates(rates, list(race = c(A = -1))),
    "argument",
    "stratify_rates"
  )
  expect_hammel_error(
    stratify_rates(transform(rates, race = "A"), list(race = c(A = 1))),
    "table",
    "stratify_rates"
  )
})

test_that("database errors name the download function", {
  expect_hammel_error(
    get_mortality_hmd("USA", 2000.5, "user", "pass"),
    "argument",
    "get_mortality_hmd"
  )
  expect_hammel_error(
    get_mortality_hmd("USA", 2000, "user", "pass", shape = "2x2"),
    "argument",
    "get_mortality_hmd"
  )
  error <- expect_hammel_error(
    get_fertility_hfd("USA", 2000, "", "pass"),
    "argument",
    "get_fertility_hfd"
  )
  expect_identical(error$argument, "username")
  skip_if_not_installed("HMDHFDplus")
  with_database_reader(
    "readHMDweb",
    function(...) {
      stop("HTTP 401")
    },
    {
      error <- expect_hammel_error(
        get_mortality_hmd("USA", 2000, "user", "pass"),
        "database",
        "get_mortality_hmd"
      )
      expect_identical(error$database, "HMD")
      expect_match(conditionMessage(error), "Couldn't download", fixed = TRUE)
      expect_s3_class(error$parent, "error")
    }
  )
  with_database_reader(
    "readHFDweb",
    function(...) {
      data.frame(Year = 1999, Age = "15", ASFR = 0.1)
    },
    {
      expect_hammel_error(
        get_fertility_hfd("USA", 2000, "user", "pass"),
        "database",
        "get_fertility_hfd"
      )
    }
  )
})

test_that("analysis and database messages keep their layout", {
  rates <- data.frame(age_start = c(15, 22), age_end = c(20, 25), rate = 0.1)
  expect_snapshot(calculate_tfr(rates), error = TRUE)
  skip_if_not_installed("HMDHFDplus")
  with_database_reader(
    "readHMDweb",
    function(...) {
      stop("HTTP 401")
    },
    {
      expect_snapshot(
        get_mortality_hmd("USA", 2000, "user", "pass"),
        error = TRUE
      )
    }
  )
})

test_that("history augmentation wrappers name themselves", {
  expect_hammel_error(
    extract_unpartnered_births(NULL),
    "population",
    "extract_unpartnered_births"
  )
  expect_hammel_error(augment_persons(NULL), "population", "augment_persons")
})
