prop_hazard_fertility_rates <- function() {
  data.frame(age_start = 15, age_end = 50, rate = 0.1)
}

prop_hazard_mortality_rates <- function() {
  data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
}

test_that("processes retain group multipliers and setters accept them", {
  factors <- c(A = 1.1, B = 0.9)
  rates <- prop_hazard_mortality_rates()
  attr(rates, "source") <- list(example = "synthetic")
  process <- mortality(rates, group_prop_hazard = factors)
  expect_identical(process$group_prop_hazard, factors)
  expect_identical(process$rates, tibble::as_tibble(rates))
  expect_identical(process$source, attr(rates, "source"))
  model <- set_mortality(hammel_model(), rates, group_prop_hazard = factors)
  expect_identical(model$processes$mortality, process)
  births <- prop_hazard_fertility_rates()
  expect_identical(fertility(births, group_prop_hazard = factors)$group_prop_hazard, factors)
  model <- set_fertility(model, births, group_prop_hazard = factors)
  expect_identical(model$processes$fertility$group_prop_hazard, factors)
  searches <- data.frame(age_start = 15, age_end = Inf, rate = 0.1)
  process <- union_formation("rate", rates = searches, group_prop_hazard = factors)
  expect_identical(process$group_prop_hazard, factors)
  model <- set_union_formation(model, "rate", rates = searches, group_prop_hazard = factors)
  expect_identical(model$processes$union_formation, process)
  expect_identical(model$processes$mortality$group_prop_hazard, factors)
})

test_that("group multipliers validate names, shape, and numerical values", {
  for (bad in list(numeric(), c(1, 2), c(A = NA_real_), c(A = NaN),
                   c(A = Inf), c(A = -0.1), c(A = 1, A = 2),
                   setNames(1, NA_character_), setNames(1, ""),
                   setNames(1, "  "), c(A = "1"), list(A = 1), matrix(1))) {
    expect_error(mortality(prop_hazard_mortality_rates(), group_prop_hazard = bad), "group_prop_hazard")
    expect_error(fertility(prop_hazard_fertility_rates(), group_prop_hazard = bad), "group_prop_hazard")
    expect_error(union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
                                 group_prop_hazard = bad), "group_prop_hazard")
  }
  expect_no_error(mortality(prop_hazard_mortality_rates(), group_prop_hazard = c(A = 0)))
  expect_no_error(fertility(prop_hazard_fertility_rates(), group_prop_hazard = c(A = 2L)))
})

test_that("edited process multiplier vectors are revalidated", {
  model <- make_rate_model()
  for (name in c("mortality", "fertility", "union_formation")) {
    process <- model$processes[[name]]
    process$group_prop_hazard <- c(A = -1)
    expect_error(set_process(model, process), "group_prop_hazard")
    edited <- model
    edited$processes[[name]] <- process
    expect_error(simulate(edited, 0), "group_prop_hazard")
  }
  expect_error(mortality(prop_hazard_mortality_rates(),
                         group_prop_hazard = c(A = 1+1i)), "group_prop_hazard")
})

prop_hazard_model <- function() {
  model <- make_rate_model(6, search_rate = 0.4, birth_interval = 0.75)
  model$population$persons$sex <- rep(c("f", "m"), 3)
  model$population$persons$group <- c("A", "B", "C", "A", "B", "C")
  model$processes$fertility$rates$rate <- c(0.2, 0.1)
  model
}

# Independently construct explicit tables for comparison with execution scaling.
explicit_prop_hazard_model <- function(model) {
  groups <- unique(as.character(model$population$persons$group))
  for (name in c("mortality", "fertility", "union_formation")) {
    process <- model$processes[[name]]
    factors <- process$group_prop_hazard
    if (is.null(factors)) next
    tables <- lapply(groups, function(group) {
      rates <- process$rates
      if ("group" %in% names(rates)) rates <- rates[as.character(rates$group) == group, , drop = FALSE]
      if (is.list(factors)) {
        factor <- vapply(rates$union_type, function(status) {
          vector <- factors[[status]]
          if (group %in% names(vector)) unname(vector[group]) else 1
        }, numeric(1), USE.NAMES = FALSE)
      } else {
        factor <- if (group %in% names(factors)) unname(factors[group]) else 1
      }
      rates$rate <- rates$rate * factor
      rates$group <- group
      rates
    })
    process$rates <- do.call(rbind, tables)
    process$group_prop_hazard <- NULL
    model$processes[[name]] <- process
  }
  model
}

test_that("multipliers equal explicit schedules for each process and together", {
  for (stratified in c(FALSE, TRUE)) {
    model <- prop_hazard_model()
    if (stratified) {
      for (name in c("mortality", "fertility", "union_formation")) {
        rates <- model$processes[[name]]$rates
        model$processes[[name]]$rates <- rbind(transform(rates, group = "C", rate = rate * 1.2),
          transform(rates, group = "B", rate = rate * 0.8), transform(rates, group = "A"))
      }
    }
    for (selected in list("mortality", "fertility", "union_formation",
                         c("mortality", "fertility", "union_formation"))) {
      modified <- model
      for (name in selected) modified$processes[[name]]$group_prop_hazard <- c(A = 1.1, B = 0.9, absent = 3)
      original <- modified
      explicit <- explicit_prop_hazard_model(modified)
      expect_identical(simulate(modified, 20, seed = 654), simulate(explicit, 20, seed = 654))
      expect_identical(modified, original)
    }
  }
})

test_that("unit and absent multipliers leave RNG consumption and histories unchanged", {
  original <- prop_hazard_model()
  for (factors in list(c(A = 1, B = 1, C = 1), c(absent = 2))) {
    modified <- original
    for (name in c("mortality", "fertility", "union_formation"))
      modified$processes[[name]]$group_prop_hazard <- factors
    set.seed(123)
    state <- .Random.seed
    expect_identical(simulate(modified, 5, seed = 654), simulate(original, 5, seed = 654))
    expect_identical(.Random.seed, state)
    set.seed(321)
    a <- simulate(original, 5)
    after <- .Random.seed
    set.seed(321)
    b <- simulate(modified, 5)
    expect_identical(a, b)
    expect_identical(.Random.seed, after)
  }
})

test_that("fertility scaling precedes cell accounting, including zero and high rates", {
  model <- make_model(2, sex = "f", birth_rate = 2, birth_interval = 0.75)
  model$population$persons$group <- c("A", "B")
  # Zero suppresses group A; B is scaled to .2.
  model$processes$fertility$group_prop_hazard <- c(A = 0, B = 0.1)
  result <- simulate(model, 20, seed = 654)
  expect_gt(result$simulation$births, 0)
  # Woman 1 (factor 0) never gives birth; daughters of woman 2 may, later.
  expect_false(any(result$persons$mother_id %in% 1))
  expect_gt(sum(result$persons$mother_id %in% 2), 0)
  expect_identical(result, simulate(explicit_prop_hazard_model(model), 20, seed = 654))
  model$processes$fertility$rates$rate <- 0.2
  model$processes$fertility$group_prop_hazard <- c(B = 10)
  expect_no_error(simulate(model, 1, seed = 654))
})

test_that("zero mortality and search multipliers suppress only the selected risk", {
  model <- make_model(2, death_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$processes$mortality$group_prop_hazard <- c(A = 0)
  result <- simulate(model, 1, seed = 654)
  expect_true(is.na(result$persons$death_time[1]))
  expect_false(is.na(result$persons$death_time[2]))
  model <- make_rate_model(2, search_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$processes$union_formation$group_prop_hazard <- c(A = 0)
  expect_equal(nrow(simulate(model, 1, seed = 654)$unions), 0)
  model$processes$union_formation$group_prop_hazard <- c(B = 0)
  expect_equal(nrow(simulate(model, 1, seed = 654)$unions), 1)
})

test_that("multipliers without population groups warn and are ignored", {
  original <- prop_hazard_model()
  original$population$persons$group <- NULL
  for (duration in c(0, 5)) {
    modified <- original
    for (name in c("mortality", "fertility", "union_formation"))
      modified$processes[[name]]$group_prop_hazard <- c(A = 0, B = 10)
    warnings <- character()
    result <- withCallingHandlers(simulate(modified, duration, seed = 654), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_length(warnings, 3)
    expect_true(all(grepl("group_prop_hazard.*no groups found.*ignoring", warnings)))
    expect_identical(result, simulate(original, duration, seed = 654))
    expect_false("group" %in% names(result$persons))
  }
  # Grouped-rate alphabetical fallback remains separate from ignoring factors.
  original$processes$mortality$rates <- rbind(
    transform(original$processes$mortality$rates, group = "B", rate = 1e8),
    transform(original$processes$mortality$rates, group = "A"))
  original$processes$mortality$group_prop_hazard <- c(A = 0)
  warnings <- character()
  withCallingHandlers(simulate(original, 0), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  expect_length(warnings, 2)
  expect_match(warnings[1], "using group 'A'")
  expect_match(warnings[2], "ignoring multipliers")
})

test_that("multipliers support numeric codes, empty populations, and restarts", {
  model <- prop_hazard_model()
  model$population$persons$group <- match(model$population$persons$group, c("A", "B", "C"))
  for (name in c("mortality", "fertility", "union_formation"))
    model$processes[[name]]$group_prop_hazard <- c("1" = 1.1, "2" = 0.9)
  result <- simulate(model, 5, seed = 654)
  explicit <- explicit_prop_hazard_model(model)
  expect_identical(result, simulate(explicit, 5, seed = 654))
  expect_identical(simulate(set_population(model, result), 5, seed = 123),
                   simulate(set_population(explicit, result), 5, seed = 123))
  model$population$persons <- model$population$persons[FALSE, ]
  expect_no_warning(empty <- simulate(model, 1))
  expect_equal(nrow(empty$persons), 0)
  expect_type(empty$persons$group, "integer")
})

test_that("multipliers cannot supply missing schedules or overflow rates", {
  model <- prop_hazard_model()
  model$processes$mortality$rates$group <- "A"
  model$processes$mortality$group_prop_hazard <- c(B = 1, C = 1)
  expect_error(simulate(model, 0), "missing schedules")
  model <- prop_hazard_model()
  model$processes$mortality$rates$rate <- 1e308
  model$processes$mortality$group_prop_hazard <- c(A = 10)
  expect_error(simulate(model, 0), "must remain finite")
  expect_error(union_formation("on_demand", group_prop_hazard = c(A = 1)), "rate-based")
  edited <- union_formation("on_demand")
  edited$group_prop_hazard <- c(A = 1)
  expect_error(set_process(hammel_model(), edited), "rate-based")
})

test_that("fertility status lists match explicit conditional group schedules", {
  for (stratified in c(FALSE, TRUE)) {
    for (factors in list(list(partnered = c(A = 1.1, B = 0.9), unpartnered = c(A = 1.05)),
                        list(partnered = c(A = 1.1)), list(unpartnered = c(B = 0)))) {
      model <- prop_hazard_model()
      if (stratified) {
        rates <- model$processes$fertility$rates
        model$processes$fertility$rates <- rbind(transform(rates, group = "C", rate = rate * 1.2),
          transform(rates, group = "B", rate = rate * 0.8), transform(rates, group = "A"))
      }
      model$processes$fertility$group_prop_hazard <- factors
      explicit <- explicit_prop_hazard_model(model)
      expect_identical(simulate(model, 20, seed = 654), simulate(explicit, 20, seed = 654))
    }
  }
  model <- prop_hazard_model()
  factors <- list(partnered = c(A = 1.1, B = 0.9), unpartnered = c(A = 1.05))
  rates <- model$processes$fertility$rates
  process <- fertility(rates, group_prop_hazard = factors)
  expect_identical(process$group_prop_hazard, factors)
  model <- set_fertility(model, rates, group_prop_hazard = factors)
  expect_identical(model$processes$fertility$group_prop_hazard, factors)
  explicit <- explicit_prop_hazard_model(model)
  expect_equal(explicit$processes$fertility$rates$rate,
               c(0.22, 0.105, 0.18, 0.1, 0.2, 0.1)) # groups A, B, C
})

test_that("a shared fertility vector equals the same vector for both statuses", {
  model <- prop_hazard_model()
  factors <- c(A = 1.1, B = 0.9)
  model$processes$fertility$group_prop_hazard <- factors
  shared <- simulate(model, 20, seed = 654)
  model$processes$fertility$group_prop_hazard <- list(partnered = factors, unpartnered = factors)
  expect_identical(simulate(model, 20, seed = 654), shared)
  model$processes$fertility$group_prop_hazard <- list(unpartnered = c(A = 1))
  original <- prop_hazard_model()
  expect_identical(simulate(model, 20, seed = 654), simulate(original, 20, seed = 654))
})

test_that("fertility multiplier lists validate status names and nested vectors", {
  rates <- prop_hazard_model()$processes$fertility$rates
  for (bad in list(list(), list(c(A = 1)), list(overall = c(A = 1)),
                   list(partnered = c(A = 1), partnered = c(B = 1)),
                   list(partnered = NULL), list(unpartnered = c(A = NA_real_)),
                   list(partnered = list(A = 1)), list(unpartnered = numeric()))) {
    expect_error(fertility(rates, group_prop_hazard = bad), "group_prop_hazard")
  }
  factors <- list(partnered = c(A = 1.1))
  expect_error(fertility(prop_hazard_fertility_rates(), group_prop_hazard = factors), "overall fertility")
  expect_error(mortality(prop_hazard_mortality_rates(), group_prop_hazard = factors), "group_prop_hazard")
  expect_error(union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0.1),
                               group_prop_hazard = factors), "group_prop_hazard")
  model <- prop_hazard_model()
  model$processes$fertility$group_prop_hazard <- list(partnered = c(A = -1))
  expect_error(set_process(model, model$processes$fertility), "group_prop_hazard")
  expect_error(simulate(model, 0), "group_prop_hazard")
  model$processes$fertility$group_prop_hazard <- list(partnered = c(A = 10))
  expect_no_error(simulate(model, 1, seed = 654))
})

test_that("fertility status lists without population groups warn once and are ignored", {
  original <- prop_hazard_model()
  original$population$persons$group <- NULL
  model <- original
  model$processes$fertility$group_prop_hazard <- list(partnered = c(A = 10), unpartnered = c(B = 0))
  for (duration in c(0, 5)) {
    warnings <- character()
    result <- withCallingHandlers(simulate(model, duration, seed = 654), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    expect_length(warnings, 1)
    expect_match(warnings, "fertility.*ignoring multipliers")
    expect_identical(result, simulate(original, duration, seed = 654))
  }
})
