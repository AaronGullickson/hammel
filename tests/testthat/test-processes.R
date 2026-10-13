test_that("process constructors and setters preserve rates and replace independently", {
  rates <- data.frame(sex = c("m", "m", "f", "f"), age_start = c(0, 1, 0, 1), age_end = c(1,
    Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
  births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
  death_process <- mortality(rates)
  birth_process <- fertility(births)
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  expect_true(inherits(death_process, "hammel_process"))
  expect_true(inherits(death_process, "hammel_mortality"))
  expect_true(inherits(birth_process, "hammel_fertility"))
  expect_identical(death_process$rates, tibble::as_tibble(rates))
  expect_identical(birth_process$rate_unit, "person_year")
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  expect_identical(model, set_fertility(set_mortality(hammel_model(), rates), births))
  replacement <- rates
  replacement$rate <- 0.03
  updated <- set_mortality(model, replacement)
  expect_identical(model$processes$mortality, death_process)
  expect_identical(updated$processes$fertility, birth_process)
  expect_identical(updated$processes$mortality$rates, tibble::as_tibble(replacement))
})

test_that("rate validation preserves custom columns and input order", {
  rates <- data.frame(sex = c("m", "m", "f", "f"), age_start = c(0, 1, 0, 1), age_end = c(1,
    Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
  births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
  death_process <- mortality(rates)
  birth_process <- fertility(births)
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  custom <- rates[c(4, 2, 1, 3), ]
  custom$source <- "example"
  class(custom) <- c("custom_rates", "data.frame")
  expect_identical(mortality(custom)$rates, tibble::as_tibble(custom))
})

test_that("invalid rates and edited process specifications are rejected", {
  rates <- data.frame(sex = c("m", "m", "f", "f"), age_start = c(0, 1, 0, 1), age_end = c(1,
    Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
  births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
  death_process <- mortality(rates)
  birth_process <- fertility(births)
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  expect_error(mortality(rates[rates$sex == "f", ]))
  expect_error(mortality(rates[FALSE, ]))
  expect_error(fertility(births[setdiff(names(births), "rate")]))
  for (value in list(-1, Inf, NA_real_, "0.1")) {
    bad <- births
    bad$rate <- value
    expect_error(fertility(bad))
  }
  for (value in list(-1, Inf, NA_real_)) {
    bad <- births
    bad$age_start[1] <- value
    expect_error(fertility(bad))
  }
  bad <- births
  bad$age_end[1] <- 15
  expect_error(fertility(bad))
  bad <- births
  bad$age_end[1] <- 21
  expect_error(fertility(bad))
  expect_error(fertility(rbind(births, births[1, ])))
  expect_error(fertility(data.frame(age_start = c(0, 1, 3), age_end = c(10, 2, 4), rate = 0.1)))
  bad <- rates
  bad$sex <- c(0, 0, 1, 1)
  expect_error(mortality(bad))
  bad <- rates
  bad$time_start <- 2000
  expect_error(mortality(bad))
  bad <- births
  bad$time_end <- 2100
  expect_error(fertility(bad))
  edited <- death_process
  edited$rates$rate[1] <- -1
  expect_error(set_process(model, edited))
  edited <- death_process
  edited$rate_unit <- "month"
  expect_error(set_process(model, edited))
  expect_error(set_process(model, rates))
})

test_that("gaps are retained and public rate tables are tibbles", {
  rates <- data.frame(sex = c("m", "m", "f", "f"), age_start = c(0, 1, 0, 1), age_end = c(1,
    Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
  births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
  death_process <- mortality(rates)
  birth_process <- fertility(births)
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  gap <- births
  gap$age_start[2] <- 25
  expect_identical(fertility(gap)$rates, tibble::as_tibble(gap))
  expect_true(tibble::is_tibble(death_process$rates))
  expect_true(tibble::is_tibble(birth_process$rates))
  edited <- death_process
  edited$rates <- as.data.frame(edited$rates)
  expect_true(tibble::is_tibble(set_process(model, edited)$processes$mortality$rates))
})

test_that("fertility spacing is stored and revalidated on attachment", {
  rates <- data.frame(sex = c("m", "m", "f", "f"), age_start = c(0, 1, 0, 1), age_end = c(1,
    Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
  births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
  death_process <- mortality(rates)
  birth_process <- fertility(births)
  model <- set_process(set_process(hammel_model(), death_process), birth_process)
  expect_equal(birth_process$birth_interval, 0.75, tolerance = 0)
  spaced <- fertility(births, birth_interval = 1.5)
  expect_identical(set_process(model, spaced), set_fertility(model, births, birth_interval = 1.5))
  expect_equal(fertility(births, birth_interval = 0)$birth_interval, 0, tolerance = 0)
  for (value in list(NULL, NA_real_, Inf, -1, 0.1, "1", c(1, 2))) {
    expect_error(fertility(births, birth_interval = value))
    expect_error(set_fertility(model, births, birth_interval = value))
    edited <- spaced
    edited$birth_interval <- value
    expect_error(set_process(model, edited))
  }
})

