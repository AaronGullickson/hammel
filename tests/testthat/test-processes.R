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
  expect_equal(fertility(births, birth_interval = 0.1)$birth_interval, 0.1)
  for (value in list(NULL, NA_real_, Inf, -1, "1", c(1, 2), matrix(1))) {
    expect_error(fertility(births, birth_interval = value))
    expect_error(set_fertility(model, births, birth_interval = value))
    edited <- spaced
    edited$birth_interval <- value
    expect_error(set_process(model, edited))
  }
})


test_that("fertility summaries flag positive open-ended age bands", {
  rates <- tibble::tibble(age_start = 15, age_end = Inf, rate = 0.1)
  expect_no_warning(process <- fertility(rates))
  expect_no_warning(result <- summary(process))
  expect_s3_class(result, "summary.hammel_fertility")
  expect_true(result$open_band_non_zero)
  expect_identical(process$rates, rates)
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl("open-ended age band", text, fixed = TRUE)))

  rates$rate <- 0
  expect_false(summary(fertility(rates))$open_band_non_zero)
  rates$rate <- 0.1
  rates$age_end <- 50
  rates$parity_start <- 0
  rates$parity_end <- Inf
  rates$time_start <- 2000
  rates$time_end <- Inf
  result <- summary(fertility(rates))
  expect_false(result$open_band_non_zero)
  text <- cli::cli_fmt(print(result))
  expect_false(any(grepl("open-ended age band", text, fixed = TRUE)))
  rates$age_end <- Inf
  rates$race <- "A"
  expect_true(summary(fertility(rates))$open_band_non_zero)
})

test_that("dissolution summaries show bands by union type", {
  rates <- tibble::tribble(
    ~union_type, ~dur_start, ~dur_end, ~rate,
    "cohab",     0,          5,        0.08,
    "cohab",     5,          Inf,      0.02,
    "marriage",  0,          10,       0.02,
    "marriage",  10,         Inf,      0.01
  )
  result <- summary(union_dissolution(rates))
  expect_s3_class(result, "summary.hammel_union_dissolution")
  expect_identical(result$rate_strata$stratified_by, "union_type")
  expect_true(result$rate_strata$bands_by_stratum)
  expect_null(result$weights)
  text <- cli::cli_fmt(value <- print(result))
  expect_identical(value, result)
  expect_true(any(grepl("Union dissolution process", text, fixed = TRUE)))
  expect_true(any(grepl(
    "Stratified by: union_type (cohab, marriage)", text, fixed = TRUE
  )))
  expect_true(any(grepl("union_type = marriage", text, fixed = TRUE)))
  expect_true(any(grepl("[10, Inf)", text, fixed = TRUE)))
  expect_true(any(grepl("Weight model: none", text, fixed = TRUE)))
})

test_that("mortality and fertility summaries report their settings", {
  result <- summary(mortality(usa_2000_mortality))
  expect_identical(result$rate_strata$stratified_by, "sex")
  text <- cli::cli_fmt(value <- print(result))
  expect_identical(value, result)
  expect_true(any(grepl("Mortality process", text, fixed = TRUE)))

  result <- summary(fertility(usa_2000_fertility, birth_interval = 1.5))
  expect_identical(result$birth_interval, 1.5)
  expect_identical(result$ffab, 0.4886)
  expect_false(result$open_band_non_zero)
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl(
    "Minimum birth interval: 1.5 years", text, fixed = TRUE
  )))
  expect_true(any(grepl(
    "Fraction female at birth: 0.4886", text, fixed = TRUE
  )))
})

test_that("union formation summaries report cohabitation and conversion", {
  cohab_prob <- tibble::tibble(
    age_start = c(15, 30),
    age_end = c(30, Inf),
    prob = c(0.6, 0.3)
  )
  convert <- tibble::tibble(dur_start = 0, dur_end = Inf, rate = 0.1)
  process <- union_formation(
    "rate",
    usa_2000_marriage,
    cohab_prob = cohab_prob,
    union_convert_rate = convert,
    weights_conversion = weight_model(fun = function(members, bin) 0),
    cohab_as_previous_union = TRUE,
    min_male_age = 18,
    kinship_generations = 1
  )
  result <- summary(process)
  expect_identical(result$method, "rate")
  expect_identical(result$min_male_age, 18)
  expect_identical(result$kinship_generations, 1)
  expect_true(result$cohab_as_previous_union)
  expect_identical(result$cohabitation$cohab_prob_range, c(0.3, 0.6))
  expect_s3_class(
    result$cohabitation$conversion_weights,
    "summary.hammel_weight_model"
  )
  text <- cli::cli_fmt(value <- print(result))
  expect_identical(value, result)
  expect_true(any(grepl("Cohabitation: enabled", text, fixed = TRUE)))
  expect_true(any(grepl("Probability range: 0.3-0.6", text, fixed = TRUE)))
  expect_true(any(grepl("Counts as previous union: yes", text, fixed = TRUE)))
  expect_true(any(grepl("Conversion to marriage rates", text, fixed = TRUE)))
  expect_true(any(grepl("Weight model: custom function", text, fixed = TRUE)))
  expect_true(any(grepl("Minimum male age: 18", text, fixed = TRUE)))

  process <- union_formation(
    "rate",
    usa_2000_marriage,
    cohab_prob = cohab_prob
  )
  result <- summary(process)
  expect_null(result$cohabitation$convert_rate_strata)
  expect_null(summary(union_formation("on_demand"))$rate_strata)
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl(
    "Conversion to marriage: not allowed", text, fixed = TRUE
  )))
  text <- cli::cli_fmt(print(summary(union_formation("on_demand"))))
  expect_true(any(grepl("Cohabitation: disabled", text, fixed = TRUE)))
})
