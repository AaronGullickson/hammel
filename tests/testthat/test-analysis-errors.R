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

test_that("summary measures check their rate tables", {
  rates <- data.frame(age_start = c(15, 30), age_end = c(30, 50), rate = 0.1)
  expect_hammel_error(
    calculate_tfr(rates[c("age_start", "age_end")]),
    "table",
    "calculate_tfr"
  )
  bad <- rates
  bad$age_end <- as.character(bad$age_end)
  expect_hammel_error(calculate_tfr(bad), "table", "calculate_tfr")
  bad <- rates
  bad$rate[1] <- -0.1
  expect_hammel_error(calculate_tfr(bad), "table", "calculate_tfr")
  bad <- rates
  bad$age_end[1] <- 10
  expect_hammel_error(calculate_tfr(bad), "rates", "calculate_tfr")
  parity <- data.frame(
    age_start = 15,
    age_end = 50,
    parity_start = c(0, 1),
    parity_end = c(1, Inf),
    rate = 0.1
  )
  expect_hammel_error(
    calculate_patfr(parity[setdiff(names(parity), "rate")]),
    "table",
    "calculate_patfr"
  )
  parity$age_end <- 10
  expect_hammel_error(calculate_patfr(parity), "rates", "calculate_patfr")
  mortality <- data.frame(age_start = c(0, 1), age_end = c(1, Inf),
    rate = c(0.01, -0.01))
  expect_hammel_error(
    calculate_life_table(mortality),
    "table",
    "calculate_life_table"
  )
})

test_that("database arguments are checked before any download", {
  expect_hammel_error(
    get_fertility_hfd("USA", 2000, "user", "pass", update = "2020"),
    "argument",
    "get_fertility_hfd"
  )
  error <- expect_hammel_error(
    get_mortality_hmd(c("USA", "SWE"), 2000, "user", "pass"),
    "argument",
    "get_mortality_hmd"
  )
  expect_identical(error$argument, "country")
})

test_that("unexpected HMD downloads are database errors", {
  five <- data.frame(
    Year = 2000,
    Age = c("0", "1-4", "5-9", "12", "15+"),
    Male = 0.01,
    Female = 0.01
  )
  expect_error(
    mortality_from_hmd(five, 2000, "5x1"),
    "five-year ages",
    class = "hammel_error_database"
  )
  mx <- data.frame(
    Year = 2000,
    Age = c("0", "1", "2+"),
    Male = c("a", "b", "c"),
    Female = 0.01
  )
  expect_error(
    mortality_from_hmd(mx, 2000),
    "must be numeric",
    class = "hammel_error_database"
  )
  x <- list(
    mx = data.frame(Year = 2000, Age = c("0", "1", "2+"), Male = 0.01,
      Female = 0.01),
    deaths = data.frame(Year = 2000, Age = c("0", "1", "2+"),
      Male = c(1, 1e308, 0), Female = 1),
    exposures = data.frame(Year = 2000, Age = c("0", "1", "2+"),
      Male = c(10, 1e-300, 1e-300), Female = 10)
  )
  expect_error(
    mortality_from_hmd(x$mx, 2000, "1x1", NULL, 1, x$deaths, x$exposures),
    "not finite",
    class = "hammel_error_database"
  )
})
