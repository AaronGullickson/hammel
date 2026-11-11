test_that("Hammel execution leaves the stats simulate generic available", {
  withr::local_seed(42)
  expect_identical(simulate, stats::simulate)
  fitted <- stats::lm(mpg ~ wt, data = datasets::mtcars)
  responses <- simulate(fitted, nsim = 2, seed = 42)
  expect_identical(dim(responses), c(32L, 2L))
})

test_that("simulation handles empty populations and zero risks", {
  expect_error(run_simulation(hammel_model(), 1))
  expect_equal(run_simulation(make_model(), 0.1)$end_time, 0.1, tolerance = 0)
  expect_equal(run_simulation(make_model(birth_interval = 0.1),
    0.1)$simulation$birth_interval, 0.1)
  model <- make_model()
  zero <- run_simulation(model, 1, seed = 42)
  expect_equal(count_run_events(zero, "births"), 0, tolerance = 0)
  expect_equal(count_run_events(zero, "deaths"), 0, tolerance = 0)
  expect_true(all(is.na(zero$persons$death_time)))
  expect_equal(zero$end_time, 1, tolerance = 0)
  expect_identical(model$population$end_time, 0)
  expect_identical(run_simulation(model, 0)$persons, model$population$persons)
  empty <- run_simulation(make_model(0), 1, seed = 42)
  expect_equal(nrow(empty$persons), 0, tolerance = 0)
  expect_equal(empty$end_time, 1, tolerance = 0)
})

test_that("mortality survival is analytical and seeded runs preserve RNG and histories",
  {
    n <- 50000
    model <- make_model(n, death_rate = 0.2)
    set.seed(99)
    before <- .Random.seed
    result <- run_simulation(model, 2, seed = 42)
    expect_identical(.Random.seed, before)
    expect_identical(result, run_simulation(model, 2, seed = 42))
    expected <- exp(-0.4)
    expect_true(abs(mean(is.na(result$persons$death_time)) - expected) < 6 * sqrt(expected *
      (1 - expected)/n))
    expect_true(all(result$persons$death_time < 2, na.rm = TRUE))
    restart <- set_fertility(set_mortality(set_population(hammel_model(), result), model$processes$mortality$rates),
      model$processes$fertility$rates) |>
      set_process(model$processes$union_formation)
    later <- run_simulation(restart, 1, seed = 43)
    dead <- !is.na(result$persons$death_time)
    expect_identical(later$persons$death_time[dead], result$persons$death_time[dead])
    expect_equal(later$end_time, 3, tolerance = 0)
  })

test_that("mortality times follow the process clock's exponential draws", {
  # All 100 men share one cell. The mortality clock spends an Exp(1) threshold
  # at the total intensity rate * members; each death draws two uniforms (cell,
  # then member) before the next threshold, and the dead man leaves the cell.
  set.seed(42)
  thresholds <- numeric(100)
  thresholds[1] <- rexp(1)
  for (k in 2:100) {
    runif(2)
    thresholds[k] <- rexp(1)
  }
  expected <- cumsum(thresholds / (1e+08 * (100:1)))
  exact <- run_simulation(make_model(100, death_rate = 1e+08), 1/12, seed = 42)
  expect_equal(sort(exact$persons$death_time), expected, tolerance = 1e-14)
  expect_true(all(exact$persons$death_time > 0))
})

test_that("events at the run endpoint are excluded", {
  set.seed(42)
  endpoint <- rexp(1)
  expect_true(is.na(run_simulation(make_model(1, death_rate = 1), endpoint,
    seed = 42)$persons$death_time))
})

test_that("birth precedes a later death within a former monthly bin", {
  # The woman's fertility-cell clock draws its threshold before her death
  # threshold; with one member the cell intensity equals her rate.
  set.seed(42)
  draws <- rexp(2)
  ordered <- run_simulation(make_model(1, sex = "f", death_rate = draws[2]/0.06,
    birth_rate = draws[1]/0.02,
    birth_interval = 0.05), 0.03, seed = 42) # spacing leaves room for one birth
  expect_equal(count_run_events(ordered, "births"), 1, tolerance = 0)
  expect_true(abs(ordered$persons$birth_time[2] - 0.02) < 1e-14)
})

test_that("births preserve maternal IDs and allow initially unknown fathers", {
  model <- make_model(1000, sex = "f", birth_rate = 0.5)
  model$population$persons$custom <- seq_len(1000)
  result <- run_simulation(model, 1, seed = 42)
  new <- result$persons$person_id > 1000
  expect_true(any(new))
  expect_true(all(result$persons$mother_id[new] %in% seq_len(1000)))
  expect_true(all(is.na(result$persons$custom[new])))
  expect_identical(result$persons$custom[!new], seq_len(1000))
  expect_true(all(result$unions$start_time >= 0))
  expect_true(is.na(result$persons$father_id[which(new)[1]]))
})

test_that("birth accounting, father availability and spacing survive restarts", {
  expect_valid_spacing <- function(times, mother_ids) {
    intervals <- lapply(split(times, mother_ids), function(x) diff(sort(x)))
    expect_gt(sum(lengths(intervals)), 0L)
    invalid <- vapply(intervals, function(x) any(x < 0.75 - 1e-09), logical(1))
    expect_false(any(invalid), info = paste(
      "Birth spacing violated for mother IDs:", paste(names(intervals)[invalid], collapse = ", ")))
  }
  model <- make_model(2000, birth_rate = 0.2)
  model$population$persons$sex[seq_len(1000)] <- "f"
  result <- run_simulation(model, 5, seed = 42)
  new <- result$persons$person_id > 2000
  expect_true(!anyDuplicated(result$persons$person_id))
  expect_true(all(result$persons$father_id[new] %in% result$persons$person_id[result$persons$sex ==
    "m"]))
  father_rows <- match(result$persons$father_id[new], result$persons$person_id)
  child_times <- result$persons$birth_time[new]
  expect_true(all(result$persons$birth_time[father_rows] <= child_times))
  expect_true(all(is.na(result$persons$death_time[father_rows]) | result$persons$death_time[father_rows] >=
    child_times))
  expect_valid_spacing(result$persons$birth_time[new], result$persons$mother_id[new])
  validate_population(result)
  restarted <- set_fertility(set_mortality(set_population(hammel_model(), result), model$processes$mortality$rates),
    model$processes$fertility$rates) |>
      set_process(model$processes$union_formation)
  continued <- run_simulation(restarted, 1, seed = 43)
  children <- !is.na(continued$persons$mother_id)
  expect_valid_spacing(continued$persons$birth_time[children],
    continued$persons$mother_id[children])
})

test_that("adjusted fertility follows the continuous renewal expectation", {
  model <- make_model(10000, sex = "f", birth_rate = 0.2)
  result <- run_simulation(model, 20, seed = 123)
  eligible_hazard <- 0.2/(1 - 0.75 * 0.2)
  expected_rate <- 1/(0.75 + 1/eligible_hazard)
  observed <- sum(result$persons$mother_id <= 10000 & result$persons$birth_time >= 5, na.rm = TRUE)/(10000 *
    15)
  expect_true(abs(observed - expected_rate) < 0.005)
  expect_equal(expected_rate, 0.2, tolerance = 1e-14)
})

test_that("death closes existing unions without forming new unions", {
  model <- make_model(2, death_rate = 1e+08)
  model$population$persons$sex[2] <- "f"
  model$population$unions <- data.frame(union_id = 1L, wife_id = 2L, husband_id = 1L,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  result <- run_simulation(model, 1/12, seed = 42)
  expect_equal(nrow(result$unions), 1L, tolerance = 0)
  expect_equal(result$unions$end_time, min(result$persons$death_time), tolerance = 0)
  expected_reason <- if (result$persons$death_time[1] < result$persons$death_time[2])
    "death_husband" else "death_wife"
  expect_equal(result$unions$end_reason, expected_reason)
})

test_that("fractional run durations and tibble outputs are supported", {
  result <- run_simulation(make_model(), 0.123, seed = 42)
  expect_equal(result$end_time, 0.123, tolerance = 0)
  expect_identical(
    names(result$simulation),
    c("start_time", "end_time", "birth_interval", "ffab")
  )
  expect_s3_class(result$unrealized_births, "tbl_df")
  expect_equal(result$unrealized_since, 0)
  expect_error(run_simulation(make_model(), -1))
  expect_error(run_simulation(make_model(), 1, step_duration = 1))
  expect_true(tibble::is_tibble(result$persons))
  expect_true(tibble::is_tibble(result$unions))
  expect_true(tibble::is_tibble(run_simulation(make_model(0), 1)$persons))
})

test_that("nondefault process spacing is enforced and revalidated on restart", {
  model <- make_model(500, sex = "f", birth_rate = 0.3, birth_interval = 1.13)
  result <- run_simulation(model, 5, seed = 42)
  continued <- run_simulation(set_population(model, result), 3, seed = 43)
  expect_equal(result$simulation$birth_interval, 1.13, tolerance = 0)
  expect_equal(continued$simulation$birth_interval, 1.13, tolerance = 0)
  expect_true(count_run_events(continued, "births") > 0)
  children <- continued$persons[!is.na(continued$persons$mother_id), ]
  intervals <- unlist(lapply(split(children$birth_time, children$mother_id), function(times) diff(sort(times))))
  expect_true(length(intervals) > 0)
  expect_true(all(intervals >= 1.13 - 1e-09))
  edited <- model
  edited$processes$fertility$birth_interval <- NULL
  expect_error(run_simulation(edited, 0))
  edited$processes$fertility$birth_interval <- -1
  expect_error(run_simulation(edited, 1))
})



test_that("restarts keep settings but discard legacy run counts", {
  model <- make_model(2)
  result <- run_simulation(model, 1, seed = 42)
  counts <- c("births", "deaths", "dissolutions", "conversions", "transitions")
  for (name in counts) {
    result$simulation[[name]] <- 10
  }
  continued <- run_simulation(set_population(model, result), 0)
  expect_false(any(counts %in% names(continued$simulation)))
  expect_identical(continued$simulation$birth_interval, 0.75)
  expect_identical(continued$simulation$ffab, result$simulation$ffab)
  expect_identical(continued$persons, result$persons)
  expect_identical(continued$unions, result$unions)
  expect_identical(continued$unrealized_births, result$unrealized_births)
  expect_identical(continued$unrealized_since, result$unrealized_since)
})
