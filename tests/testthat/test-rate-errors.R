mortality_input <- function() {
  tibble::tibble(
    sex = c("m", "f"),
    age_start = 0,
    age_end = Inf,
    rate = 0.01
  )
}

fertility_input <- function() {
  tibble::tibble(age_start = c(15, 30), age_end = c(30, 50), rate = 0.1)
}

dissolution_input <- function() {
  tibble::tibble(dur_start = c(0, 5), dur_end = c(5, Inf), rate = 0.05)
}

test_that("constructors and setters name themselves", {
  bad <- mortality_input()[c("sex", "age_start", "rate")]
  model <- hammel_model()
  expect_hammel_error(mortality(bad), "table", "mortality")
  expect_hammel_error(set_mortality(model, bad), "table", "set_mortality")
  expect_hammel_error(
    model |> set_mortality(bad),
    "table",
    "set_mortality"
  )
  bad <- fertility_input()
  bad$rate[2] <- -1
  expect_hammel_error(fertility(bad), "table", "fertility")
  expect_hammel_error(set_fertility(model, bad), "table", "set_fertility")
  expect_hammel_error(
    set_fertility(model, fertility_input(), ffab = 2),
    "argument",
    "set_fertility"
  )
  expect_hammel_error(
    fertility(fertility_input(), birth_interval = -1),
    "argument",
    "fertility"
  )
  bad <- dissolution_input()
  bad$dur_end[1] <- 10
  expect_hammel_error(union_dissolution(bad), "rates", "union_dissolution")
  expect_hammel_error(
    set_union_dissolution(model, bad),
    "rates",
    "set_union_dissolution"
  )
  process <- mortality(mortality_input())
  process$rates$rate <- NA_real_
  expect_hammel_error(set_process(model, process), "table", "set_process")
  expect_hammel_error(
    set_process(model, list()),
    "specification",
    "set_process"
  )
  expect_hammel_error(set_mortality(list(), mortality_input()), "specification",
    "set_mortality")
})

test_that("pipelines name the setter that failed, not run_simulation", {
  population <- generate_uniform_population(20)
  bad <- mortality_input()[c("sex", "age_start", "rate")]
  expect_hammel_error(
    hammel_model() |>
      set_population(population) |>
      set_mortality(bad) |>
      run_simulation(1),
    "table",
    "set_mortality"
  )
})

test_that("row problems report the input's rows", {
  rates <- tibble::tibble(
    sex = rep(c("m", "f"), 2),
    region = rep(c("n", "s"), each = 2),
    age_start = 0,
    age_end = Inf,
    rate = c(0.01, 0.01, 0.01, NA),
    time_start = 2000,
    time_end = Inf
  )
  error <- expect_hammel_error(mortality(rates), "table", "mortality")
  expect_identical(error$rows, 4L)
  rates$rate[4] <- 0.01
  rates$age_end[3] <- 0
  error <- expect_hammel_error(mortality(rates), "rates", "mortality")
  expect_identical(error$rows, 3L)
})

test_that("schedule problems name their period and stratum", {
  rates <- tibble::tibble(
    sex = c("m", "f", "m", "m", "f"),
    region = c("n", "n", "s", "s", "s"),
    age_start = c(0, 0, 0, 10, 0),
    age_end = c(Inf, Inf, 20, Inf, Inf),
    rate = 0.01,
    time_start = 2000,
    time_end = Inf
  )
  error <- expect_hammel_error(mortality(rates), "rates", "mortality")
  message <- conditionMessage(error)
  expect_match(message, "Bands [0, 20) and [10, Inf) overlap", fixed = TRUE)
  expect_match(message, 'for sex "m"', fixed = TRUE)
  expect_match(message, "In group stratum region = s.", fixed = TRUE)
  expect_match(message, "In calendar period [2000, Inf).", fixed = TRUE)

  rates <- mortality_input()[1, ]
  error <- expect_hammel_error(mortality(rates), "rates", "mortality")
  expect_match(conditionMessage(error), 'No rates for sex "f"', fixed = TRUE)
})

test_that("parity bands must partition parities in every schedule", {
  rates <- tibble::tibble(
    age_start = 15,
    age_end = 50,
    parity_start = c(0, 2),
    parity_end = c(1, Inf),
    rate = 0.1
  )
  error <- expect_hammel_error(fertility(rates), "rates", "fertility")
  expect_match(conditionMessage(error), "[0, 1) and [2, Inf)", fixed = TRUE)
  rates$parity_start[2] <- 1.5
  error <- expect_hammel_error(fertility(rates), "rates", "fertility")
  expect_identical(error$rows, 2L)
})

test_that("coverage errors name run_simulation and the uncovered time", {
  rates <- mortality_input()
  rates$time_start <- 0
  rates$time_end <- 1
  model <- hammel_model() |>
    set_population(generate_uniform_population(10)) |>
    set_mortality(rates) |>
    set_fertility(fertility_input()) |>
    set_union_formation("on_demand")
  error <- expect_hammel_error(
    run_simulation(model, 2),
    "coverage",
    "run_simulation"
  )
  expect_identical(error$time, 1)
  expect_match(conditionMessage(error), "Mortality rates", fixed = TRUE)
})

test_that("model schedules report the row with an invalid model", {
  schedule <- tibble::tibble(
    time_start = c(0, 10),
    time_end = c(10, Inf),
    model = list(weight_model(~partnered, c(partnered = 0)), "x")
  )
  error <- expect_hammel_error(
    mortality(mortality_input(), weights = schedule),
    "specification",
    "mortality"
  )
  expect_identical(error$rows, 2L)
})

test_that("rate messages keep their layout", {
  rates <- tibble::tibble(
    sex = c("m", "f", "m"),
    age_start = c(0, 0, 10),
    age_end = c(20, Inf, Inf),
    rate = 0.01,
    time_start = 2000,
    time_end = Inf
  )
  expect_snapshot(mortality(rates), error = TRUE)
  rates <- mortality_input()
  rates$time_start <- 0
  rates$time_end <- 1
  model <- hammel_model() |>
    set_population(generate_uniform_population(10)) |>
    set_mortality(rates) |>
    set_fertility(fertility_input()) |>
    set_union_formation("on_demand")
  expect_snapshot(run_simulation(model, 2), error = TRUE)
})
