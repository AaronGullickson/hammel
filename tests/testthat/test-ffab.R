test_that("fertility stores and validates the female-birth probability", {
  rates <- data.frame(age_start = 15, age_end = 50, rate = 0.2)
  expect_identical(fertility(rates)$ffab, 0.4886)
  for (probability in c(0, 0.25, 0.5, 1)) {
    process <- fertility(rates, ffab = probability)
    expect_identical(process$ffab, probability)
    model <- set_fertility(hammel_model(), rates, ffab = probability)
    expect_identical(model$processes$fertility, process)
  }
  for (bad in list(NULL, numeric(), c(0.2, 0.8), NA_real_, NaN, Inf,
                  -Inf, -0.01, 1.01, TRUE, "0.5", matrix(0.5))) {
    expect_error(fertility(rates, ffab = bad), "ffab")
    expect_error(set_fertility(hammel_model(), rates, ffab = bad), "ffab")
    process <- fertility(rates)
    process["ffab"] <- list(bad)
    expect_error(set_process(hammel_model(), process), "ffab")
    model <- make_model(0)
    model$processes$fertility <- process
    expect_error(simulate(model, 0, seed = 1), "ffab")
  }
})

test_that("default and older fertility specifications preserve seeded runs", {
  model <- make_model(40, sex = "f", birth_rate = 1)
  expected <- simulate(model, 2, seed = 735)
  explicit <- set_fertility(model, model$processes$fertility$rates, ffab = 0.4886)
  expect_identical(simulate(explicit, 2, seed = 735), expected)
  older <- model
  older$processes$fertility$ffab <- NULL
  expect_identical(set_process(older, older$processes$fertility), model)
  expect_identical(simulate(older, 2, seed = 735), expected)
  expect_identical(expected$simulation$ffab, 0.4886)
  expect_identical(simulate(set_fertility(model, model$processes$fertility$rates,
    ffab = 0.5), 0)$simulation$ffab, 0.5)
})

test_that("ffab endpoints control newborn sex in both initiation modes", {
  for (mode in c("on_demand", "rate")) {
    model <- make_model(100, birth_rate = 2, birth_interval = 0)
    model$population$persons$sex <- rep(c("f", "m"), 50)
    if (mode == "rate") {
      rates <- model$processes$fertility$rates
      model <- set_fertility(model, rbind(transform(rates, union_type = "partnered"),
        transform(rates, union_type = "unpartnered")), birth_interval = 0)
      model <- set_union_formation(model, "rate", rates = data.frame(
        age_start = 15, age_end = 50, rate = 1), cancel_birth_if_no_match = FALSE)
    }
    for (probability in c(0, 1)) {
      model$processes$fertility$ffab <- probability
      set.seed(123)
      rng <- .Random.seed
      result <- simulate(model, 1, seed = 736)
      expect_identical(.Random.seed, rng)
      children <- result$persons$birth_time >= 0
      expect_gt(sum(children), 0)
      expect_true(all(result$persons$sex[children] == if (probability == 1) "f" else "m"))
      expect_identical(result$persons$sex[!children], model$population$persons$sex)
      expect_identical(result$simulation$ffab, probability)
      expect_identical(simulate(model, 1, seed = 736), result)
    }
  }
})

test_that("newborn sex follows a nondefault Bernoulli probability", {
  model <- make_model(2000, sex = "f", birth_rate = 1, birth_interval = 0)
  model$processes$fertility$ffab <- 0.25
  result <- simulate(model, 1, seed = 737)
  children <- result$persons$birth_time >= 0
  n <- sum(children)
  expect_gt(n, 1000)
  # Six binomial standard errors detects gross miswiring without a tight
  # stochastic acceptance threshold for demographic validation.
  expect_lt(abs(mean(result$persons$sex[children] == "f") - 0.25),
    6 * sqrt(0.25 * 0.75 / n))
})

test_that("configured sex is available to built-in and custom inheritance", {
  for (probability in c(0, 1)) {
    for (method in c("same_sex", "opposite_sex", "custom")) {
      inheritance <- if (method == "custom") group_inheritance("custom",
        fun = function(child_sex, ...) if (child_sex == "f") "A" else "B") else method
      model <- make_model(2, birth_rate = 20, birth_interval = 0)
      model$population$persons$sex <- c("f", "m")
      model$population$persons$group <- c("A", "B")
      model$population$groups <- list(group = group())
      model <- model |>
        set_fertility(data.frame(age_start = 15, age_end = 50, rate = 20),
          birth_interval = 0, ffab = probability) |> set_inheritance(inheritance)
      result <- simulate(model, 0.3, seed = 654)
      children <- result$persons$birth_time >= 0
      expect_true(any(children))
      expected_group <- if (probability == 1) "A" else "B"
      if (method == "opposite_sex") expected_group <- if (probability == 1) "B" else "A"
      expect_true(all(result$persons$group[children] == expected_group))
    }
  }
})

test_that("changing ffab affects only subsequent births and independent runs", {
  model <- make_model(20, sex = "f", birth_rate = 20, birth_interval = 0)
  model$processes$fertility$ffab <- 1
  first <- simulate(model, 0.3, seed = 738)
  continued_model <- set_population(model, first)
  continued_model$processes$fertility$ffab <- 0
  continued <- simulate(continued_model, 0.3, seed = 739)
  new <- continued$persons$birth_time >= first$end_time
  expect_true(any(new))
  expect_true(all(continued$persons$sex[new] == "m"))
  expect_identical(continued$persons$sex[!new], first$persons$sex)
  expect_identical(continued$simulation$ffab, 0)
  expect_identical(model$processes$fertility$ffab, 1)
  expect_identical(simulate(model, 0.3, seed = 738), first)
})
