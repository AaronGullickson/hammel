group_status_rates <- function(partnered = 0, unpartnered = 0) {
  data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(partnered, unpartnered))
}

# Group labels deliberately differ from row order and public person IDs.
group_model <- function() {
  model <- make_rate_model(4, search_rate = 0)
  model$population$persons$sex <- c("f", "m", "f", "m")
  model$population$persons$group <- c("B", "A", "A", "B")
  model
}

group_rates <- function(rates, a = 0, b = 0) {
  rbind(transform(rates, group = "B", rate = b),
        transform(rates, group = "A", rate = a))
}

test_that("group labels validate and survive population construction", {
  p <- group_model()$population$persons
  for (labels in list(c("A", "B", "A", "B"), c(20, 10, 20, 10))) {
    p$group <- labels
    expect_identical(hammel_population(p, end_time = 0)$persons$group, labels)
  }
  for (labels in list(c("A", NA, "B", "A"), c(1, NA, 2, 1),
                     c("A", "", "B", "A"), rep(TRUE, 4))) {
    p$group <- labels
    expect_error(hammel_population(p, end_time = 0), "group")
  }
})

test_that("group schedules validate overlap and required strata independently", {
  model <- group_model()
  mortality_rates <- group_rates(model$processes$mortality$rates)
  expect_s3_class(mortality(mortality_rates)$rates, "tbl_df")
  expect_error(mortality(rbind(mortality_rates, mortality_rates[1, ])), "overlap")
  expect_error(mortality(mortality_rates[-1, ]), "both")
  births <- group_rates(group_status_rates())
  expect_s3_class(fertility(births), "hammel_fertility")
  expect_error(fertility(births[-1, ]), "union_status")
  searches <- group_rates(model$processes$union_formation$rates)
  expect_s3_class(set_union_formation(model, "rate", rates = searches), "hammel_model")
  expect_error(set_union_formation(model, "rate", rates = rbind(searches, searches[1, ])), "overlap")
  searches$group[1] <- NA_character_
  expect_error(set_union_formation(model, "rate", rates = searches), "group")
  mortality_rates$group[1] <- NA_character_
  expect_error(mortality(mortality_rates), "group")
  births$group[1] <- NA_character_
  expect_error(fertility(births), "group")
})

test_that("shared and identical grouped schedules preserve fixed-seed histories", {
  model <- group_model()
  model <- set_fertility(model, group_status_rates(0.3, 0.2), 0.75)
  model <- set_union_formation(model, "rate", rates = data.frame(
    age_start = 15, age_end = Inf, rate = 1), sample_size = 0)
  model <- set_mortality(model, data.frame(sex = c("m", "f"),
    age_start = 0, age_end = Inf, rate = 0.01))
  original <- model
  original$population$persons$group <- NULL
  ungrouped <- simulate(original, 5, seed = 123)
  shared <- simulate(model, 5, seed = 123)
  for (name in c("mortality", "fertility", "union_formation")) {
    rates <- model$processes[[name]]$rates
    model$processes[[name]]$rates <- rbind(transform(rates, group = "B"),
                                        transform(rates, group = "A"))
  }
  grouped <- simulate(model, 5, seed = 123)
  expect_identical(grouped, shared)
  shared$persons$group <- NULL
  expect_identical(shared, ungrouped)
})

test_that("group-specific mortality and newborn mortality use own group", {
  model <- group_model()
  model <- set_mortality(model, group_rates(model$processes$mortality$rates, b = 1e8))
  result <- simulate(model, 1, seed = 12)
  expect_true(all(is.na(result$persons$death_time[result$persons$group == "A"])))
  expect_true(all(!is.na(result$persons$death_time[result$persons$group == "B"])))

  model <- set_fertility(group_model(), group_status_rates(0, 20), 0)
  # Only newborn B children face mortality; adults are outside this band.
  rates <- data.frame(sex = rep(c("m", "f"), each = 2),
    age_start = c(0, 1), age_end = c(1, Inf), rate = 0)
  rates <- rbind(transform(rates, group = "A"), transform(rates, group = "B"))
  rates$rate[rates$group == "B" & rates$age_start == 0] <- 1e8
  model <- set_mortality(model, rates)
  result <- simulate(model, 0.5, seed = 12)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 0)
  expect_true(any(children$group == "B"))
  expect_true(all(!is.na(children$death_time[children$group == "B"])))
  expect_true(all(is.na(children$death_time[children$group == "A"])))
})

test_that("fertility selects maternal groups and inheritance persists over restarts", {
  model <- group_model()
  model$population$persons$group <- c(20, 10, 10, 20)
  rates <- rbind(transform(group_status_rates(0, 20), group = 20),
                 transform(group_status_rates(), group = 10))
  model <- set_fertility(model, rates, 0)
  result <- simulate(model, 0.5, seed = 12)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 0)
  expect_true(all(children$mother_id == 1))
  expect_true(all(children$group == 20))
  expect_type(result$persons$group, "double")
  continued <- simulate(set_population(model, result), 0.5, seed = 13)
  expect_identical(head(continued$persons$group, nrow(result$persons)), result$persons$group)
  kids <- continued$persons[continued$persons$birth_time >= 0, ]
  expect_equal(kids$group, continued$persons$group[match(kids$mother_id, continued$persons$person_id)])
})

test_that("group search hazards allow intergroup unions", {
  model <- make_rate_model(2, search_rate = 0)
  model$population$persons$group <- c("B", "A")
  model <- set_union_formation(model, "rate", sample_size = 0,
    rates = group_rates(model$processes$union_formation$rates, b = 1e8))
  result <- simulate(model, 0.1, seed = 12)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$wife_id, 1)
  expect_equal(result$unions$husband_id, 2)
  model$processes$union_formation$rates$rate <- c(0, 1e8)
  expect_equal(nrow(simulate(model, 0.1, seed = 12)$unions), 0)
})

test_that("ungrouped populations use alphabetical first rates and warn per process", {
  model <- make_rate_model(2, search_rate = 0)
  for (name in c("mortality", "fertility", "union_formation")) {
    model$processes[[name]]$rates <- group_rates(model$processes[[name]]$rates, b = 1e8)
  }
  warnings <- character()
  result <- withCallingHandlers(simulate(model, 1, seed = 12), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 3)
  expect_true(all(grepl("no groups found.*'A'", warnings)))
  expect_equal(result$simulation$births, 0)
  expect_equal(result$simulation$deaths, 0)
  expect_equal(nrow(result$unions), 0)
  expect_false("group" %in% names(result$persons))
})

test_that("missing group schedules fail and mortality coverage is per group", {
  for (name in c("mortality", "fertility", "union_formation")) {
    model <- group_model()
    model$processes[[name]]$rates$group <- "A"
    expect_error(simulate(model, 0), "missing schedules.*B")
  }
  model <- group_model()
  rates <- group_rates(model$processes$mortality$rates)
  rates$age_start[rates$group == "B"] <- 1
  model <- set_mortality(model, rates)
  expect_error(simulate(model, 0), "cover")
})

test_that("group fertility compensation stays within each group and status", {
  model <- group_model()
  model <- set_fertility(model, group_rates(group_status_rates(), a = 0.2, b = 0.4), 0.75)
  expect_no_warning(simulate(model, 0.1, seed = 12))
  model$processes$fertility$rates$rate[1] <- 2
  expect_error(simulate(model, 0.1), "midpoint eligibility")
})

test_that("empty grouped populations remain valid and preserve group type", {
  model <- group_model()
  model$population$persons <- model$population$persons[FALSE, ]
  model <- set_fertility(model, group_rates(group_status_rates()))
  result <- simulate(model, 1, seed = 12)
  expect_equal(nrow(result$persons), 0)
  expect_identical(result$persons$group, character())
  expect_s3_class(result$persons, "tbl_df")
})

test_that("baseline group weights are normalized, sampled and reproducible", {
  set.seed(135)
  p <- create_baseline_population(10000, group_weights = c(A = 2, B = 1, C = 0))
  expect_s3_class(p$persons, "tbl_df")
  expect_true(all(p$persons$group %in% c("A", "B")))
  # Five binomial standard deviations; tests sampling rather than exact quotas.
  expect_lt(abs(mean(p$persons$group == "A") - 2/3), 5 * sqrt((2/3) * (1/3) / 10000))
  set.seed(135)
  expect_identical(p, create_baseline_population(10000, group_weights = c(A = 2, B = 1, C = 0)))
  set.seed(135)
  scaled <- create_baseline_population(10000, group_weights = c(A = 20, B = 10, C = 0))
  expect_identical(scaled, p)
  set.seed(135)
  ungrouped <- create_baseline_population(10000)
  p$persons$group <- NULL
  expect_identical(p, ungrouped)
  expect_identical(create_baseline_population(0, group_weights = c(A = 1))$persons$group,
                   character())
  expect_true(all(create_baseline_population(10, group_weights = c(A = 0, B = 1))$persons$group == "B"))
  expect_no_warning(create_baseline_population(10, group_weights = c(A = 1e308, B = 1e308)))
})

test_that("invalid baseline group weights fail without consuming RNG", {
  set.seed(321)
  state <- .Random.seed
  for (weights in list(numeric(), c(1, 2), c(A = 0, B = 0), c(A = -1, B = 2),
                       c(A = NA_real_), c(A = Inf), c(A = 1, A = 2),
                       setNames(1, ""), setNames(1, NA_character_),
                       list(A = 1), c(A = "1"), matrix(1, 1))) {
    expect_error(create_baseline_population(10, group_weights = weights), "group_weights")
    expect_identical(.Random.seed, state)
  }
})

test_that("on-demand intergroup unions pass maternal identity to newborns", {
  model <- make_model(2, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m")
  model$population$persons$group <- c("B", "A")
  result <- simulate(model, 0.2, seed = 12)
  expect_equal(nrow(result$unions), 1)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 0)
  expect_true(all(children$group == "B"))
  expect_true(all(children$mother_id == 1 & children$father_id == 2))
})

test_that("group and union status jointly select fertility schedules", {
  model <- group_model()
  model$population$unions <- tibble::tibble(union_id = 1,
    wife_id = 1, husband_id = 2, start_time = -1,
    end_time = NA_real_, end_reason = NA_character_)
  rates <- rbind(transform(group_status_rates(20, 0), group = "B"),
                 transform(group_status_rates(0, 0), group = "A"))
  model <- set_fertility(model, rates, 0)
  result <- simulate(model, 0.5, seed = 12)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 0)
  expect_true(all(children$mother_id == 1))
})

test_that("group-specific survival agrees with constant-hazard probabilities", {
  model <- make_model(10000)
  model$population$persons$group <- rep(c("A", "B"), each = 5000)
  model <- set_mortality(model, group_rates(model$processes$mortality$rates, a = 0.1, b = 0.3))
  result <- simulate(model, 2, seed = 654)
  for (group in c("A", "B")) {
    hazard <- if (group == "A") 0.1 else 0.3
    probability <- 1 - exp(-2 * hazard)
    observed <- mean(!is.na(result$persons$death_time[result$persons$group == group]))
    expect_lt(abs(observed - probability), 5 * sqrt(probability * (1-probability) / 5000))
  }
})

test_that("numeric codes use alphabetical fallback and empty runs warn", {
  model <- make_model(0)
  rates <- rbind(transform(model$processes$mortality$rates, group = 2),
                 transform(model$processes$mortality$rates, group = 10))
  model <- set_mortality(model, rates)
  expect_warning(result <- simulate(model, 0), "using group '10'")
  expect_equal(nrow(result$persons), 0)
  expect_false("group" %in% names(result$persons))
})
