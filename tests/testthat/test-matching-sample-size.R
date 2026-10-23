test_that("matching sample size defaults and validation belong to matching models", {
  expect_equal(create_matching_model(~ 1, numeric())$sample_size, 1000)
  expect_equal(create_age_gap_matching_model(0, 5)$sample_size, 1000)
  expect_equal(union_formation("on_demand")$matching$sample_size, 1000)
  for (size in c(0, 1, 50, .Machine$integer.max)) {
    matching <- create_matching_model(~ 1, numeric(), sample_size = size)
    expect_equal(create_age_gap_matching_model(0, 5, sample_size = size)$sample_size, size)
    expect_identical(union_formation("on_demand", matching)$matching, matching)
  }
  for (size in list(NULL, -1, 0.5, Inf, NA_real_, NaN, "100", TRUE,
                   c(1, 2), matrix(1), 1 + 1i, .Machine$integer.max + 1)) {
    expect_error(create_matching_model(~ 1, numeric(), sample_size = size), "sample_size")
    expect_error(create_age_gap_matching_model(0, 5, sample_size = size), "sample_size")
    model <- make_model(0)
    model$processes$union_formation$matching$sample_size <- size
    expect_error(set_process(hammel_model(), model$processes$union_formation), "sample_size")
    expect_error(simulate(model, 0), "sample_size")
    model <- make_rate_model(0)
    model$processes$union_formation$father_matching <- create_matching_model(~ 1, numeric())
    model$processes$union_formation$father_matching$sample_size <- size
    expect_error(simulate(model, 0), "sample_size")
  }
  expect_error(union_formation("on_demand", sample_size = 250), "unused argument")
  expect_error(set_union_formation(hammel_model(), "on_demand", sample_size = 250), "unused argument")
  legacy <- union_formation("on_demand")
  legacy$sample_size <- 1000
  expect_error(set_process(hammel_model(), legacy), "matching model")
})

test_that("the default samples 1000 men in both matching methods", {
  model <- make_model(2201, birth_rate = 100, birth_interval = 0)
  model$population$persons$sex <- c("f", rep("m", 1100), rep("f", 1100))
  model$population$persons$birth_time <- c(-20, rep(-25, 1100), rep(-60, 1100))
  for (method in c("female_choice", "weighted_choice")) {
    model <- set_union_formation(model, "on_demand", diagnostics = TRUE,
      matching = create_age_gap_matching_model(0, 5, matching_method = method))
    result <- simulate(model, 0.2, seed = 52)
    expect_equal(nrow(result$unions), 1)
    expect_identical(result$unions$.n_eligible, 1000L)
  }
})

test_that("union and father models use independent sample sizes with NULL inheritance", {
  for (method in c("female_choice", "weighted_choice")) {
    union_sizes <- father_sizes <- integer()
    record_union <- function(gap) { union_sizes <<- c(union_sizes, length(gap)); gap }
    record_father <- function(gap) { father_sizes <<- c(father_sizes, length(gap)); gap }
    general <- create_matching_model(~ record_union(age_difference),
      c("record_union(age_difference)" = 0), sample_size = 3, matching_method = method)
    fathers <- create_matching_model(~ record_father(age_difference),
      c("record_father(age_difference)" = 0), sample_size = 1, matching_method = method)
    model <- make_model(6)
    model$population$persons$sex <- c("f", "m", "m", "m", "f", "f")
    model$population$persons$birth_time <- c(-20, -20, -25, -30, -60, -60)
    # Female-only newborns keep the male encounter pool fixed at three adults.
    model <- set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
      age_start = 19, age_end = 21, rate = 100), birth_interval = 0, ffab = 1)
    model <- set_union_formation(model, "rate", general, father_matching = fathers,
      rates = data.frame(age_start = 20.2, age_end = 21, rate = 1000),
      diagnostics = TRUE)
    # Constructor and execution validation evaluate formulas on a 101-row grid;
    # keep those validation calls distinct from actual candidate-set sizes.
    result <- simulate(model, 0.4, seed = 52)
    expect_equal(nrow(result$unions), 1)
    expect_identical(result$unions$.n_eligible, 3L)
    expect_true(3L %in% union_sizes)
    expect_true(1L %in% father_sizes)
    expect_true(all(union_sizes %in% c(3L, 101L)))
    expect_true(all(father_sizes %in% c(1L, 101L)))
    children <- result$persons[result$persons$birth_time >= 0, ]
    expect_true(any(children$birth_time < result$unions$start_time))
    expect_equal(model$processes$union_formation$matching$sample_size, 3)
    expect_equal(model$processes$union_formation$father_matching$sample_size, 1)
    model$processes$union_formation$father_matching <- NULL
    inherited <- simulate(model, 0.4, seed = 52)
    model$processes$union_formation$father_matching <- general
    expect_identical(simulate(model, 0.4, seed = 52), inherited)
  }
})
