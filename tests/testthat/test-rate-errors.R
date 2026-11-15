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

search_input <- function() {
  tibble::tibble(age_start = c(15, 30), age_end = c(30, Inf), rate = 0.1)
}

test_that("union formation errors name the constructor or setter", {
  model <- hammel_model()
  expect_hammel_error(
    union_formation("sometimes"),
    "argument",
    "union_formation"
  )
  expect_hammel_error(
    set_union_formation(model, "rate"),
    "specification",
    "set_union_formation"
  )
  expect_hammel_error(
    set_union_formation(model, "on_demand", rates = search_input()),
    "specification",
    "set_union_formation"
  )
  expect_hammel_error(
    union_formation("on_demand", min_male_age = -1),
    "argument",
    "union_formation"
  )
  expect_hammel_error(
    union_formation("on_demand", diagnostics = NA),
    "argument",
    "union_formation"
  )
  expect_hammel_error(
    union_formation("on_demand", kinship_generations = 1.5),
    "argument",
    "union_formation"
  )
  expect_hammel_error(
    union_formation(
      "on_demand",
      weights = weight_model(~partnered, c(partnered = 0))
    ),
    "specification",
    "union_formation"
  )
  process <- union_formation("rate", search_input())
  process$rates$rate[2] <- -1
  expect_hammel_error(set_process(model, process), "table", "set_process")
})

test_that("search rows report the input's rows and clocks", {
  rates <- tibble::tibble(
    age_start = c(15, NA, 30),
    age_end = c(30, NA, Inf),
    dur_start = c(NA, 0, 5),
    dur_end = c(NA, 5, Inf),
    previous_union = c("none", "dissolved", "dissolved"),
    rate = 0.1
  )
  error <- expect_hammel_error(
    union_formation("rate", rates),
    "rates",
    "union_formation"
  )
  expect_identical(error$rows, 3L)
  rates$age_start[3] <- rates$age_end[3] <- NA
  rates$previous_union[2] <- "none"
  error <- expect_hammel_error(
    union_formation("rate", rates),
    "rates",
    "union_formation"
  )
  expect_identical(error$rows, 2L)
})

test_that("cohabitation tables are named in their errors", {
  cohab <- tibble::tibble(age_start = 15, age_end = Inf, prob = 1.5)
  error <- expect_hammel_error(
    union_formation("rate", search_input(), cohab_prob = cohab),
    "table",
    "union_formation"
  )
  expect_identical(error$argument, "cohab_prob$prob")
  cohab$prob <- 0.5
  cohab$age_end <- 10
  error <- expect_hammel_error(
    union_formation("rate", search_input(), cohab_prob = cohab),
    "rates",
    "union_formation"
  )
  expect_identical(error$argument, "cohab_prob")
  convert <- tibble::tibble(
    dur_start = c(0, 3),
    dur_end = c(5, Inf),
    rate = 0.1
  )
  cohab$age_end <- Inf
  error <- expect_hammel_error(
    set_union_formation(
      hammel_model(), "rate", search_input(),
      cohab_prob = cohab,
      union_convert_rate = convert
    ),
    "rates",
    "set_union_formation"
  )
  expect_identical(error$argument, "union_convert_rate")
  expect_match(
    conditionMessage(error),
    "Bands [0, 5) and [3, Inf) overlap",
    fixed = TRUE
  )
  convert$age_start <- c(15, NA)
  convert$age_end <- c(Inf, NA)
  convert$dur_start[1] <- convert$dur_end[1] <- NA
  error <- expect_hammel_error(
    union_formation(
      "rate", search_input(),
      cohab_prob = cohab,
      union_convert_rate = convert
    ),
    "rates",
    "union_formation"
  )
  expect_match(conditionMessage(error), "one clock", fixed = TRUE)
  expect_hammel_error(
    union_formation("on_demand", union_convert_rate = convert),
    "specification",
    "union_formation"
  )
})

test_that("union formation messages keep their layout", {
  rates <- tibble::tibble(
    age_start = 15,
    age_end = Inf,
    previous_union = c("none", "separated"),
    rate = 0.1
  )
  expect_snapshot(union_formation("rate", rates), error = TRUE)
})

transition_input <- function() {
  tibble::tibble(
    origin = c("A", "B"),
    destination = c("B", "A"),
    age_start = 0,
    age_end = Inf,
    rate = 0.1
  )
}

test_that("transition errors name the constructor or setter", {
  bad <- transition_input()
  bad$destination[1] <- "A"
  error <- expect_hammel_error(
    group_transition(bad),
    "rates",
    "group_transition"
  )
  expect_identical(error$rows, 1L)
  expect_hammel_error(
    set_group_transition(hammel_model(), bad),
    "rates",
    "set_group_transition"
  )
  expect_hammel_error(
    group_transition(transition_input(), household = "yes"),
    "specification",
    "group_transition"
  )
  expect_hammel_error(household(-1), "argument", "household")
  bad <- transition_input()
  bad$dur_start <- 0
  error <- expect_hammel_error(
    group_transition(bad),
    "rates",
    "group_transition"
  )
  expect_match(conditionMessage(error), "age bands or duration bands")
})

test_that("transition schedules report overlapping bands by stratum", {
  rates <- tibble::tibble(
    origin = "A",
    destination = "B",
    age_start = c(0, 10),
    age_end = c(20, Inf),
    rate = 0.1
  )
  error <- expect_hammel_error(
    group_transition(rates),
    "rates",
    "group_transition"
  )
  expect_match(
    conditionMessage(error),
    "overlap for origin = A, destination = B",
    fixed = TRUE
  )
})

test_that("weight model errors name the function that received them", {
  expect_hammel_error(weight_model(), "argument", "weight_model")
  expect_hammel_error(
    weight_model(~ parity - 1, c(parity = 1)),
    "specification",
    "weight_model"
  )
  weights <- weight_model(~parity, c(paritx = 1))
  error <- expect_hammel_error(
    mortality(mortality_input(), weights = weights),
    "specification",
    "mortality"
  )
  expect_identical(error$argument, "coefficients")
  weights <- weight_model(~same_group, c(same_group = 1))
  error <- expect_hammel_error(
    set_mortality(hammel_model(), mortality_input(), weights = weights),
    "specification",
    "set_mortality"
  )
  expect_identical(error$variables, "same_group")
  expect_warning(
    mortality(mortality_input(), weights = weight_model(~sex, c(sexm = 1))),
    class = "hammel_warning_no_effect"
  )
})

test_that("weight and transition messages keep their layout", {
  weights <- weight_model(~parity, c(paritx = 1))
  expect_snapshot(
    mortality(mortality_input(), weights = weights),
    error = TRUE
  )
  rates <- transition_input()
  rates$dur_start <- 0
  expect_snapshot(group_transition(rates), error = TRUE)
})

test_that("tables are checked before they are converted", {
  rates <- data.frame(
    dur_start = 0,
    dur_end = Inf,
    rate = 0.1,
    rate = 0.2,
    check.names = FALSE
  )
  error <- expect_hammel_error(
    union_dissolution(rates),
    "table",
    "union_dissolution"
  )
  expect_identical(error$columns, "rate")
  rates <- data.frame(
    origin = "A",
    destination = "B",
    age_start = 0,
    age_end = Inf,
    rate = 1,
    rate = 2,
    check.names = FALSE
  )
  error <- expect_hammel_error(
    set_group_transition(hammel_model(), rates),
    "table",
    "set_group_transition"
  )
  expect_identical(error$columns, "rate")
  expect_hammel_error(group_transition(list()), "table", "group_transition")
})

test_that("remaining rate-table checks name their functions", {
  rates <- mortality_input()
  rates$age_end <- as.character(rates$age_end)
  expect_hammel_error(mortality(rates), "table", "mortality")
  rates <- fertility_input()
  rates$sex <- "f"
  expect_hammel_error(fertility(rates), "rates", "fertility")
  rates <- data.frame(
    age_start = 15,
    age_end = 50,
    parity_start = c(0, 3e9),
    parity_end = c(3e9, Inf),
    rate = 0.1
  )
  expect_hammel_error(fertility(rates), "rates", "fertility")
  search <- search_input()[0, ]
  expect_hammel_error(
    union_formation("rate", search),
    "rates",
    "union_formation"
  )
  search <- search_input()
  search$sex <- "f"
  expect_hammel_error(
    union_formation("rate", search),
    "rates",
    "union_formation"
  )
  cohab <- tibble::tibble(age_start = 15, age_end = Inf, prob = 0.5)
  expect_hammel_error(
    union_formation("rate", search_input(), cohab_prob = cohab,
      union_convert_rate = data.frame(dur_start = numeric(),
        dur_end = numeric(), rate = numeric())),
    "rates",
    "union_formation"
  )
  process <- mortality(mortality_input())
  process$name <- NA_character_
  expect_hammel_error(
    set_process(hammel_model(), process),
    "specification",
    "set_process"
  )
})

test_that("remaining transition checks name their functions", {
  expect_hammel_error(
    group_transition(transition_input()[0, ]),
    "rates",
    "group_transition"
  )
  rates <- transition_input()
  rates$union_status <- "single"
  expect_hammel_error(group_transition(rates), "rates", "group_transition")
  rates <- transition_input()[c("origin", "destination", "rate")]
  expect_hammel_error(group_transition(rates), "rates", "group_transition")
  rates <- transition_input()
  rates$dur_start <- 0
  rates$dur_end <- Inf
  rates$age_start <- rates$age_end <- NULL
  expect_hammel_error(
    group_transition(rates, household = household()),
    "specification",
    "group_transition"
  )
  rates <- rbind(
    transform(transition_input(), union_type = "partnered"),
    transform(transition_input(), union_type = "cohab")
  )
  expect_hammel_error(group_transition(rates), "rates", "group_transition")
})

test_that("weight and matching models with no columns take no coefficients", {
  expect_hammel_error(
    mortality(mortality_input(), weights = weight_model(~1, c(a = 1))),
    "specification",
    "mortality"
  )
})

test_that("model schedules need a list of models", {
  schedule <- data.frame(time_start = 0, time_end = Inf, model = "x")
  expect_hammel_error(
    union_formation("on_demand", matching = schedule),
    "specification",
    "union_formation"
  )
})
