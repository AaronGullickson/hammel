group_status_rates <- function(partnered = 0, unpartnered = 0) {
  data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(partnered, unpartnered))
}

# Group labels deliberately differ from row order and public person IDs.
group_model <- function() {
  model <- make_rate_model(4, search_rate = 0)
  model$population$persons$sex <- c("f", "m", "f", "m")
  model$population$persons$group <- c("B", "A", "A", "B")
  model$population$groups <- list(group = group())
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
    expect_identical(hammel_population(p, end_time = 0, groups = if ("group" %in% names(p)) "group")$persons$group, labels)
  }
  for (labels in list(c("A", NA, "B", "A"), c(1, NA, 2, 1),
                     c("A", "", "B", "A"), rep(TRUE, 4))) {
    p$group <- labels
    expect_error(hammel_population(p, end_time = 0, groups = if ("group" %in% names(p)) "group"), "group")
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
  expect_error(fertility(births[-1, ]), "union_type")
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

test_that("shared and identical grouped schedules agree in distribution", {
  # Cells stratify by group only when the rates do, so shared schedules pool
  # groups while explicitly identical grouped schedules split them. The two
  # agree in distribution, not seed by seed.
  model <- group_model()
  model <- set_inheritance(set_fertility(model, group_status_rates(0.3, 0.2), 0.75), "mother")
  model <- set_union_formation(model, "rate", rates = data.frame(
    age_start = 15, age_end = Inf, rate = 1),
    matching = create_age_gap_matching_model(0, 5))
  model <- set_mortality(model, data.frame(sex = c("m", "f"),
    age_start = 0, age_end = Inf, rate = 0.01))
  original <- model
  original$population$persons$group <- NULL
  original$population$groups <- NULL
  original$groups <- NULL
  ungrouped <- simulate(original, 5, seed = 123)
  shared <- simulate(model, 5, seed = 123)
  original_grouped <- model
  for (name in c("mortality", "fertility", "union_formation")) {
    rates <- model$processes[[name]]$rates
    model$processes[[name]]$rates <- rbind(transform(rates, group = "B"),
                                        transform(rates, group = "A"))
  }
  expect_equivalent_unions(original_grouped, model, 5)
  expect_s3_class(ungrouped, "hammel_population")
  expect_s3_class(shared, "hammel_population")
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
  model$population$groups <- list(group = group())
  rates <- rbind(transform(group_status_rates(0, 20), group = 20),
                 transform(group_status_rates(), group = 10))
  model <- set_inheritance(set_fertility(model, rates, 0), "mother")
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
  model$population$groups <- list(group = group())
  model <- set_union_formation(model, "rate",
    matching = create_age_gap_matching_model(0, 5),
    rates = group_rates(model$processes$union_formation$rates, b = 1e8))
  result <- simulate(model, 0.1, seed = 12)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$wife_id, 1)
  expect_equal(result$unions$husband_id, 2)
  model$processes$union_formation$rates$rate <- c(0, 1e8)
  expect_equal(nrow(simulate(model, 0.1, seed = 12)$unions), 0)
})

test_that("group-stratified rates require a declared group dimension", {
  model <- make_rate_model(2, search_rate = 0)
  for (name in c("mortality", "fertility", "union_formation")) {
    edited <- model
    edited$processes[[name]]$rates <- group_rates(edited$processes[[name]]$rates, b = 1e8)
    expect_error(simulate(edited, 1, seed = 12), "not a declared group dimension")
  }
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

test_that("group fertility cells stay within each group and status", {
  model <- group_model()
  model <- set_fertility(model, group_rates(group_status_rates(), a = 0.2, b = 0.4), 0.75)
  expect_no_warning(simulate(model, 0.1, seed = 12))
  # High rates need no feasibility check under cell accounting.
  model$processes$fertility$rates$rate[1] <- 2
  expect_no_error(simulate(model, 0.1, seed = 12))
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
  p <- create_baseline_population(10000, groups = list(group = c(A = 2, B = 1, C = 0)))
  expect_s3_class(p$persons, "tbl_df")
  expect_true(all(p$persons$group %in% c("A", "B")))
  # Five binomial standard deviations; tests sampling rather than exact quotas.
  expect_lt(abs(mean(p$persons$group == "A") - 2/3), 5 * sqrt((2/3) * (1/3) / 10000))
  set.seed(135)
  expect_identical(p, create_baseline_population(10000, groups = list(group = c(A = 2, B = 1, C = 0))))
  set.seed(135)
  scaled <- create_baseline_population(10000, groups = list(group = c(A = 20, B = 10, C = 0)))
  expect_identical(scaled, p)
  set.seed(135)
  ungrouped <- create_baseline_population(10000)
  p$persons$group <- NULL
  p$groups <- NULL
  p$ancestry <- NULL
  # Grouped populations now also carry an empty transition history.
  p$transitions <- NULL
  expect_identical(p, ungrouped)
  expect_identical(create_baseline_population(0, groups = list(group = c(A = 1)))$persons$group,
                   character())
  expect_true(all(create_baseline_population(10, groups = list(group = c(A = 0, B = 1)))$persons$group == "B"))
  expect_no_warning(create_baseline_population(10, groups = list(group = c(A = 1e308, B = 1e308))))
})

test_that("invalid baseline group weights fail without consuming RNG", {
  set.seed(321)
  state <- .Random.seed
  for (weights in list(numeric(), c(1, 2), c(A = 0, B = 0), c(A = -1, B = 2),
                       c(A = NA_real_), c(A = Inf), c(A = 1, A = 2),
                       setNames(1, ""), setNames(1, NA_character_),
                       list(A = 1), c(A = "1"), matrix(1, 1))) {
    expect_error(create_baseline_population(10, groups = list(group = weights)), "groups")
    expect_identical(.Random.seed, state)
  }
})

test_that("on-demand intergroup unions support explicit maternal inheritance", {
  model <- make_model(2, birth_rate = 20, birth_interval = 0)
  model$population$persons$sex <- c("f", "m")
  model$population$persons$group <- c("B", "A")
  model$population$groups <- list(group = group())
  model <- set_inheritance(model, create_group_inheritance("mother"))
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
  model$population$groups <- list(group = group())
  model <- set_mortality(model, group_rates(model$processes$mortality$rates, a = 0.1, b = 0.3))
  result <- simulate(model, 2, seed = 654)
  for (group in c("A", "B")) {
    hazard <- if (group == "A") 0.1 else 0.3
    probability <- 1 - exp(-2 * hazard)
    observed <- mean(!is.na(result$persons$death_time[result$persons$group == group]))
    expect_lt(abs(observed - probability), 5 * sqrt(probability * (1-probability) / 5000))
  }
})

test_that("numeric group rates without a declaration are rejected", {
  model <- make_model(0)
  rates <- rbind(transform(model$processes$mortality$rates, group = 2),
                 transform(model$processes$mortality$rates, group = 10))
  model <- set_mortality(model, rates)
  expect_error(simulate(model, 0), "set_groups")
})

