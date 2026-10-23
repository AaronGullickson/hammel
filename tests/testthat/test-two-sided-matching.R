# Independent fixed-matrix reference; one uniform draw per retained match.
reference_two_sided_market <- function(utility, ego = 1L, draws = 1L, scorer = NULL) {
  if (!is.null(scorer)) utility <- t(vapply(seq_len(nrow(utility)), function(i)
    scorer(seq_len(nrow(utility)), seq_len(ncol(utility)), TRUE, i), numeric(ncol(utility))))
  normalize <- function(x) {
    if (all(!is.finite(x))) return(rep(0, length(x)))
    w <- exp(x - max(x)); w / sum(w)
  }
  a <- normalize(utility[ego, ])
  b <- vapply(seq_len(ncol(utility)), function(j) normalize(utility[, j])[ego], numeric(1))
  w <- a / nrow(utility) + b / ncol(utility)
  choices <- vapply(seq_len(draws), function(i)
    which(runif(1) * sum(w) < cumsum(w))[1], integer(1))
  list(choices = choices, steps = rep(1L, draws), chooser = rep(-1L, draws),
    scored_pairs = rep(as.double(length(utility)), draws))
}

test_that("fixed-pool kernel agrees with independently normalized matrix probabilities", {
  for (utility in list(matrix(c(0, 1, -2, 4, 0, -1, 2, -3, 1), 3),
                       matrix(c(0, -Inf, 2, -Inf, 0, 1, -2, -Inf, 0), 3),
                       matrix(c(1000, -1000), 1))) {
    set.seed(361)
    expected <- reference_two_sided_market(utility, draws = 1000L)
    expected_rng <- .Random.seed
    set.seed(361)
    observed <- inspect_two_sided_matrix_cpp(utility, 1L, 1000L)
    expect_identical(observed, expected)
    expect_identical(.Random.seed, expected_rng)
  }
})

test_that("integrated choices match analytical probabilities and repeated rejection draws", {
  weights <- matrix(c(1, 2, 3, 4), 2)
  expected <- 49 / 148
  set.seed(362)
  result <- inspect_two_sided_matrix_cpp(log(weights), 1L, 12000L)
  expect_lt(abs(mean(result$choices == 1L) - expected),
    6 * sqrt(expected * (1 - expected) / 12000))
  expect_true(all(result$steps == 1L))
  # Independent literal implementation: keep both pools fixed until ego appears.
  rejection <- replicate(12000, {
    repeat {
      female <- runif(1) < 0.5
      initiator <- sample.int(2, 1)
      if (female) {
        partner <- sample.int(2, 1, prob = weights[initiator, ])
        if (initiator == 1) break
      } else {
        woman <- sample.int(2, 1, prob = weights[, initiator])
        if (woman == 1) { partner <- initiator; break }
      }
    }
    partner
  })
  expect_lt(abs(mean(rejection == 1) - expected),
    6 * sqrt(expected * (1 - expected) / 12000))
  expect_error(inspect_two_sided_matrix_cpp(matrix(NA_real_, 1), 1L, 1L), "finite")
  expect_error(inspect_two_sided_matrix_cpp(matrix(0, 1), 0L, 1L), "ego")
  expect_true(is.na(inspect_two_sided_matrix_cpp(matrix(-Inf, 1), 1L, 1L)$choices))
})

test_that("native and custom scoring retain female-first orientation in male choices", {
  fb <- c(-20, -30, -40); mb <- c(-22, -28, -38)
  fg <- c(0L, 1L, 0L); mg <- c(1L, 0L, 1L)
  effects <- c("A:B" = 2, "B:A" = -1)
  native <- create_matching_model(~ I(age_difference^2), c("I(age_difference^2)" = -0.02),
    group_exogamy = effects)
  fallback <- create_matching_model(~ I((age_difference + 0)^2),
    c("I((age_difference + 0)^2)" = -0.02), group_exogamy = effects)
  models <- list(native, fallback,
    create_matching_model(~ female_group + male_group + same_group,
      c(female_groupB = 1, male_groupB = -2, same_group = 0.5),
      group_exogamy = effects, group_levels = c("A", "B")))
  outputs <- list()
  for (k in seq_along(models)) {
    model <- models[[k]]
    utility <- do.call(rbind, lapply(seq_along(fb), function(i)
      calculate_matching_utilities(model, fb[i] - mb, c("A", "B")[fg[i] + 1L],
        c("A", "B")[mg + 1L])))
    set.seed(363)
    expected <- inspect_two_sided_matrix_cpp(utility, 1L, 1000L)
    set.seed(363)
    outputs[[k]] <- inspect_two_sided_matching_cpp(fb, mb, fg, mg,
      prepare_matching_model(model, c("A", "B"), TRUE), 1000L)
    expect_identical(outputs[[k]], expected)
  }
  expect_identical(outputs[[1]], outputs[[2]])
})

test_that("matching method is optional, validated and backwards compatible", {
  expect_identical(union_formation("on_demand")$matching_method, "female_choice")
  expect_identical(set_union_formation(hammel_model(), "on_demand",
    matching_method = "two_sided")$processes$union_formation$matching_method, "two_sided")
  for (bad in list(NULL, NA, TRUE, "both", c("female_choice", "two_sided"))) {
    expect_error(union_formation("on_demand", matching_method = bad), "matching_method")
    model <- make_model(0)
    model$processes$union_formation["matching_method"] <- list(bad)
    expect_error(simulate(model, 0), "matching_method")
  }
  model <- make_model(30, sex = "f", birth_rate = 1)
  expected <- simulate(model, 1, seed = 364)
  model$processes$union_formation$matching_method <- NULL
  expect_identical(simulate(model, 1, seed = 364), expected)
})

make_two_sided_model <- function(mode = "on_demand", extra_women = 5L) {
  model <- make_model(3L + extra_women, birth_rate = 100, birth_interval = 0)
  model$population$persons$sex <- c("f", "m", "m", rep("f", extra_women))
  model$population$persons$birth_time <- c(-20, -25, -30, rep(-40, extra_women))
  # Only ego can actually have births/searches during this short run.
  model <- set_fertility(model, data.frame(age_start = 19, age_end = 21, rate = 100), birth_interval = 0)
  set_union_formation(model, mode, matching = create_age_gap_matching_model(0, 5),
    sample_size = 0, matching_method = "two_sided", diagnostics = TRUE,
    rates = if (mode == "rate") data.frame(age_start = 19, age_end = 21, rate = 100) else NULL)
}

test_that("only ego's actual union changes state in either initiation mode", {
  for (mode in c("on_demand", "rate")) {
    model <- make_two_sided_model(mode)
    set.seed(365); before <- .Random.seed
    result <- suppressWarnings(simulate(model, 0.2, seed = 366))
    expect_identical(.Random.seed, before)
    expect_identical(result, suppressWarnings(simulate(model, 0.2, seed = 366)))
    expect_equal(nrow(result$unions), 1)
    expect_equal(result$unions$wife_id, 1)
    expect_true(result$unions$husband_id %in% 2:3)
    expect_equal(result$persons[4:8, ], model$population$persons[4:8, ])
    expect_identical(result$unions$.n_eligible, 2L)
    chosen <- result$unions$husband_id
    expect_identical(result$unions$.utility_rank, if (chosen == 2) 1L else 2L)
    model$processes$union_formation$diagnostics <- FALSE
    plain <- suppressWarnings(simulate(model, 0.2, seed = 366))
    expect_identical(plain$persons, result$persons)
    expect_identical(plain$unions, result$unions[setdiff(names(result$unions), c(".n_eligible", ".utility_rank"))])
  }
})

test_that("shortage reduction, age cutoffs, and ego exclusion retain eligible partners", {
  for (extra in 0:1) {
    model <- make_two_sided_model(extra_women = extra)
    result <- simulate(model, 0.2, seed = 367)
    expect_equal(nrow(result$unions), 1)
    expect_true(result$unions$husband_id %in% 2:3)
    # Diagnostics retain the original screened count even after reduction.
    expect_identical(result$unions$.n_eligible, 2L)
  }
  model <- make_two_sided_model(extra_women = 1L)
  model$population$persons$birth_time[4] <- -14.8
  # The alternate ages into the pool during the run. Zero birth rate keeps
  # births out of the accounting; only ego has scheduled searches.
  model <- set_fertility(model, data.frame(age_start = 19, age_end = 21, rate = 0))
  model$processes$union_formation$initiation <- "rate"
  model$processes$union_formation$rates <- tibble::tibble(age_start = 19, age_end = 21, rate = 1)
  result <- suppressWarnings(simulate(model, 1, seed = 368))
  expect_true(all(result$unions$wife_id == 1))
  expect_true(all(result$unions$husband_id %in% 2:3))
})

test_that("hard group exclusions can strand an alternate but never pair ego illegally", {
  model <- make_two_sided_model(extra_women = 1L)
  model$population$persons$group <- c("B", "B", "B", "A")
  model$processes$union_formation$matching <- create_age_gap_matching_model(0, 5,
    prevent_exogamy = TRUE)
  result <- simulate(model, 0.2, seed = 369)
  expect_equal(nrow(result$unions), 1)
  expect_equal(result$unions$wife_id, 1)
  expect_true(result$unions$husband_id %in% 2:3)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] %in% 2:3))
})

test_that("nonmarital fathers share two-sided selection and respect overrides", {
  model <- make_two_sided_model("rate", extra_women = 2L)
  model$population$persons$group <- c("A", "A", "B", "B", "A")
  model$processes$union_formation$rates$rate <- 0
  model$processes$union_formation$matching <- create_age_gap_matching_model(0, 5,
    prevent_exogamy = TRUE)
  # A father-specific restriction leaves only the A father.
  model$processes$union_formation$father_matching <- create_age_gap_matching_model(0, 5,
    prevent_exogamy = "A:B")
  result <- suppressWarnings(simulate(model, 0.3, seed = 370))
  expect_equal(nrow(result$unions), 0)
  children <- result$persons$birth_time >= 0
  expect_true(any(children))
  expect_true(all(result$persons$father_id[children] == 2))
  model$processes$union_formation$father_matching <- create_age_gap_matching_model(0, 5,
    group_exogamy = c("A:B" = 1000))
  result <- suppressWarnings(simulate(model, 0.3, seed = 370))
  expect_equal(nrow(result$unions), 0)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] %in% 2:3))
})


test_that("R scoring evaluates each fixed female row once", {
  fb <- c(-20, -30, -40); mb <- c(-22, -28, -38)
  model <- create_matching_model(~ I((age_difference - mean(age_difference))^2),
    c("I((age_difference - mean(age_difference))^2)" = -0.02))
  prepared <- prepare_matching_model(model, two_sided = TRUE)
  evaluate <- prepared$evaluate
  observed_lengths <- integer()
  prepared$evaluate <- function(gaps, ...) {
    observed_lengths <<- c(observed_lengths, length(gaps))
    evaluate(gaps, ...)
  }
  scorer <- function(women, men, female, initiator) {
    gaps <- if (female) fb[women[initiator]] - mb[men] else fb[women] - mb[men[initiator]]
    calculate_matching_utilities(model, gaps)
  }
  set.seed(371)
  expected <- reference_two_sided_market(matrix(0, 3, 3), draws = 100L, scorer = scorer)
  set.seed(371)
  observed <- inspect_two_sided_matching_cpp(fb, mb, rep(0L, 3), rep(0L, 3), prepared, 100L)
  expect_identical(observed, expected)
  expect_true(all(observed_lengths %in% 1:3))
  expect_true(all(observed_lengths == 3))
  expect_length(observed_lengths, 300)
})

test_that("materialized and streamed native matrices give identical fixed-pool outcomes", {
  prepared <- prepare_matching_model(create_age_gap_matching_model(3, 5), two_sided = TRUE)
  set.seed(372)
  lazy <- inspect_two_sided_matching_cpp(c(-20, -30, -40), c(-22, -28, -38),
    rep(0L, 3), rep(0L, 3), prepared, 1000L)
  state <- .Random.seed
  set.seed(372)
  eager <- inspect_two_sided_matching_cpp(c(-20, -30, -40), c(-22, -28, -38),
    rep(0L, 3), rep(0L, 3), prepared, 1000L, precompute = TRUE)
  expect_identical(eager, lazy)
  expect_identical(.Random.seed, state)
})

test_that("two-sided selection retains ego kin screening, failure rules and independent state", {
  model <- make_two_sided_model(extra_women = 1L)
  model$population$persons$father_id[1] <- 2
  model$population$persons$birth_time[2] <- -50
  original <- model
  first <- simulate(model, 0.2, seed = 373)
  expect_equal(first$unions$husband_id, 3)
  expect_identical(first$unions$.n_eligible, 1L)
  expect_identical(model, original)
  continued <- simulate(set_population(model, first), 0.2, seed = 374)
  expect_equal(nrow(continued$unions), 1)
  expect_identical(simulate(model, 0.2, seed = 373), first)
  # Remove all eligible men without giving any actual unions to alternates.
  model$population$persons$sex[2:3] <- "f"
  empty <- simulate(model, 0.2, seed = 373)
  expect_equal(nrow(empty$unions), 0)
  expect_equal(empty$simulation$births, 0)
  model$processes$union_formation$cancel_birth_if_no_match <- FALSE
  allowed <- simulate(model, 0.2, seed = 373)
  expect_gt(allowed$simulation$births, 0)
  expect_true(all(is.na(allowed$persons$father_id[allowed$persons$birth_time >= 0])))
})

test_that("alternate eligibility ages in and removes actual deaths before searches", {
  original_prepare <- prepare_matching_model
  lengths <- integer()
  local_mocked_bindings(prepare_matching_model = function(...) {
    prepared <- original_prepare(...)
    evaluate <- prepared$evaluate
    prepared$evaluate <- function(gaps, ...) {
      lengths <<- c(lengths, length(gaps))
      evaluate(gaps, ...)
    }
    prepared
  }, .package = "hammel")
  model <- make_two_sided_model(extra_women = 1L)
  model$processes$union_formation$diagnostics <- FALSE
  # The observer wraps the R evaluator without changing utilities or RNG.
  model$processes$union_formation$matching <- create_matching_model(
    ~ abs(age_difference), c("abs(age_difference)" = 0))
  model <- set_fertility(model, data.frame(age_start = 20.3, age_end = 21, rate = 100),
    birth_interval = 0)
  model$population$persons$birth_time[4] <- -14
  result <- simulate(model, 0.5, seed = 375)
  expect_equal(nrow(result$unions), 1)
  expect_true(all(lengths == 1L))
  model$population$persons$birth_time[4] <- -14.8
  lengths <- integer()
  aged <- simulate(model, 0.5, seed = 375)
  expect_equal(nrow(aged$unions), 1)
  expect_identical(lengths[1], 2L)
  # Restart before the threshold must reconstruct the same eligibility pool.
  checkpoint <- simulate(model, 0.2, seed = 375)
  lengths <- integer()
  aged_restart <- simulate(set_population(model, checkpoint), 0.3, seed = 376)
  expect_equal(nrow(aged_restart$unions), 1)
  expect_identical(lengths[1], 2L)

  model$population$persons$group <- c("A", "A", "A", "B")
  model$population$persons$birth_time[4] <- -40
  model <- set_mortality(model, data.frame(sex = rep(c("m", "f"), 2),
    age_start = 0, age_end = Inf, group = rep(c("A", "B"), each = 2),
    rate = c(0, 0, 0, 1e6)))
  lengths <- integer()
  died <- simulate(model, 0.5, seed = 377)
  expect_lt(died$persons$death_time[4], died$unions$start_time[1])
  expect_true(all(lengths == 1L))
})

test_that("widowhood restores the female pool without altering prior unions", {
  model <- make_two_sided_model(extra_women = 1L)
  model$population$persons$group <- c("A", "D", "A", "A")
  model$population$unions <- tibble::tibble(union_id = 10, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  model <- set_fertility(model, data.frame(age_start = 20.3, age_end = 21, rate = 100),
    birth_interval = 0)
  model <- set_mortality(model, data.frame(sex = rep(c("m", "f"), 2),
    age_start = 0, age_end = Inf, group = rep(c("A", "D"), each = 2),
    rate = c(0, 0, 1e6, 0)))
  result <- simulate(model, 0.5, seed = 378)
  expect_equal(nrow(result$unions), 2)
  expect_equal(result$unions$wife_id, c(1, 1))
  expect_equal(result$unions$husband_id, c(2, 3))
  expect_equal(result$unions$end_reason, c("death_husband", NA_character_))
  expect_equal(result$unions$start_time[1], -1)
})

test_that("bulk native columns agree with materialized matrices across restrictions and scales", {
  models <- list(
    create_matching_model(~ age_difference, c(age_difference = 0.1)),
    create_age_gap_matching_model(3, 5,
      group_exogamy = c("A:B" = 1000, "B:A" = -1000)),
    create_age_gap_matching_model(0, 5, prevent_exogamy = "A:B"),
    create_age_gap_matching_model(0, 5, prevent_exogamy = TRUE))
  for (model in models) {
    prepared <- prepare_matching_model(model, c("A", "B"), TRUE)
    set.seed(821)
    fb <- -runif(17, 15, 80); mb <- -runif(13, 15, 80)
    fg <- sample(0:1, 17, TRUE); mg <- sample(0:1, 13, TRUE)
    set.seed(822)
    bulk <- inspect_two_sided_matching_cpp(fb, mb, fg, mg, prepared, 300L)
    state <- .Random.seed
    set.seed(822)
    matrix <- inspect_two_sided_matching_cpp(fb, mb, fg, mg, prepared, 300L, TRUE)
    expect_identical(bulk, matrix)
    expect_identical(.Random.seed, state)
  }
  prepared <- prepare_matching_model(create_age_gap_matching_model(0, 5), two_sided = TRUE)
  expect_error(inspect_two_sided_matching_cpp(-1e200, 1e200, 0L, 0L, prepared, 1L), "finite")
})
