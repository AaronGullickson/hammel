# Flat probability calendars, using the existing age/group/previous strata.
cohab_calendar <- function(values = c(0, 1), times = c(0, .5, Inf)) {
  data.frame(age_start = 0, age_end = Inf, prob = values,
    time_start = head(times, -1), time_end = tail(times, -1))
}
cohab_calendar_model <- function(prob = cohab_calendar(), n = 2, search = 1) {
  m <- make_rate_model(n, search_rate = search)
  if (n) m$population$persons$sex <- rep(c("f", "m"), length.out = n)
  set_union_formation(m, "rate", rates = data.frame(age_start = 0, age_end = Inf, rate = search),
    cohab_prob = prob, union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = 0),
    kinship_generations = 0)
}

test_that("cohabitation calendars validate every active probability slice", {
  p <- cohab_calendar()
  attr(p, "source") <- list(example = "synthetic probability schedule")
  m <- cohab_calendar_model(p)
  expect_s3_class(m$processes$union_formation$cohab_prob, "tbl_df")
  expect_identical(attr(m$processes$union_formation$cohab_prob, "source"), attr(p, "source"))
  for (bad in list(p[setdiff(names(p), "time_end")],
      transform(p, time_start = NA_real_), transform(p, time_end = time_start),
      transform(p, time_start = Inf), transform(p, time_end = NaN),
      transform(p, prob = c(0, 1.1)), transform(p, prob = c(0, NA_real_)),
      transform(p, rate = .1), transform(p, dur_start = 0, dur_end = Inf),
      transform(p, sex = "f"), rbind(p, p[1, ])))
    expect_error(cohab_calendar_model(bad))
  edited <- m
  edited$processes$union_formation$cohab_prob$prob[2] <- -1
  expect_error(simulate(edited, 0), "prob")
  # Nonoverlap applies to simultaneous age/stratum risks, not calendar alone.
  p <- rbind(transform(p, age_end = 20), transform(p, age_start = 20))
  expect_s3_class(cohab_calendar_model(p), "hammel_model")
  p$age_start <- 0
  expect_error(cohab_calendar_model(p), "overlap")
})

test_that("probability coverage follows rate calendar and group/category rules", {
  m <- cohab_calendar_model(cohab_calendar(times = c(.1, .5, 1)))
  expect_error(simulate(m, 0), "cohabitation_probability.*calendar coverage")
  m$population$end_time <- .1
  expect_s3_class(simulate(m, .9, seed = 3), "hammel_population")
  expect_error(simulate(m, 1, seed = 3), "calendar coverage")
  p <- cohab_calendar(); p$time_start[2] <- .7
  expect_error(simulate(cohab_calendar_model(p), 1), "calendar coverage")
  expect_s3_class(simulate(cohab_calendar_model(p), .4, seed = 3), "hammel_population")
  m <- cohab_calendar_model(rbind(transform(cohab_calendar(), group = "A"),
    transform(cohab_calendar(), group = "B")))
  m$population$persons$group <- c("A", "B")
  p <- m$processes$union_formation$cohab_prob
  m$processes$union_formation$cohab_prob <- p[!(p$group == "B" & p$time_start == .5), ]
  expect_error(simulate(m, 1), "cohabitation probability.*missing schedules.*B")
  p <- rbind(transform(cohab_calendar(), previous_union = "none"),
    transform(cohab_calendar(), previous_union = "separated"),
    transform(cohab_calendar(), previous_union = "widowed"))
  expect_s3_class(cohab_calendar_model(p), "hammel_model")
  expect_error(cohab_calendar_model(p[-nrow(p), ]), "all three")
})

test_that("type probabilities change at fractional dates while searches recover their hazard", {
  n <- 8000
  m <- cohab_calendar_model(n = 2 * n, search = .7)
  out <- simulate(m, 1.3, seed = 911)
  u <- out$unions
  expect_gt(sum(u$start_time < .5), 0)
  expect_gt(sum(u$start_time >= .5), 0)
  expect_true(all(u$union_type[u$start_time < .5] == "marriage"))
  expect_true(all(u$union_type[u$start_time >= .5] == "cohab"))
  expected <- 1 - exp(-.7 * 1.3)
  expect_lt(abs(nrow(u) / n - expected), 5 * sqrt(expected * (1 - expected) / n))
  expect_equal(out$simulation$conversions, 0)
  # Existing marriages remain marriage after the probability change.
  expect_true(all(is.na(u$end_time)))
})

test_that("non-extreme calendar probabilities recover type shares in each period", {
  m <- cohab_calendar_model(cohab_calendar(c(.2, .8)), n = 16000, search = .7)
  out <- simulate(m, 1.3, seed = 554)
  for (late in c(FALSE, TRUE)) {
    u <- out$unions[(out$unions$start_time >= .5) == late, ]
    p <- if (late) .8 else .2
    expect_gt(nrow(u), 100)
    expect_lt(abs(mean(u$union_type == "cohab") - p), 5 * sqrt(p * (1 - p) / nrow(u)))
  }
})

test_that("probability changes precede an exact tied search without rescheduling it", {
  seed <- 783
  set.seed(seed)
  # Female: first threshold is death, second is search. No positive fertility.
  boundary <- rexp(2)[2]
  m <- cohab_calendar_model(cohab_calendar(times = c(0, boundary, Inf)))
  early <- simulate(m, boundary, seed = seed)
  expect_equal(nrow(early$unions), 0)
  out <- simulate(m, boundary + .01, seed = seed)
  expect_equal(nrow(out$unions), 1)
  expect_equal(out$unions$start_time, boundary, tolerance = 0)
  expect_equal(out$unions$union_type, "cohab")
})

test_that("probability-only changes preserve existing union clocks and histories", {
  m <- cohab_calendar_model(data.frame(age_start = 0, age_end = Inf, prob = 1),
    n = 200, search = 0)
  m$population$unions <- tibble::tibble(union_id = 1:100,
    wife_id = seq(1, 199, 2), husband_id = seq(2, 200, 2), start_time = 0,
    end_time = NA_real_, end_reason = NA_character_, union_type = "cohab")
  m$processes$union_formation$union_convert_rate$rate <- .2
  m <- set_union_dissolution(m, data.frame(dur_start = 0, dur_end = Inf, rate = .3))
  reference <- simulate(m, 2, seed = 839)
  m$processes$union_formation$cohab_prob <- cohab_calendar(c(0, 1))
  out <- simulate(m, 2, seed = 839)
  expect_identical(out, reference)
})

test_that("probability changes refresh weighted matching contributions immediately", {
  n <- 400
  p <- rbind(transform(cohab_calendar(c(.9, .2), c(0, 1e-6, Inf)), group = "A"),
    transform(cohab_calendar(c(.1, .8), c(0, 1e-6, Inf)), group = "B"))
  m <- cohab_calendar_model(p, n = n, search = .4)
  m$population$persons$sex <- rep(c("f", "m"), each = n / 2)
  m$population$persons$group <- rep(c("A", "B"), n / 2)
  m$processes$union_formation$diagnostics <- TRUE
  reference <- m
  reference$processes$union_formation$cohab_prob <- data.frame(
    age_start = 0, age_end = Inf, group = c("A", "B"), prob = c(.2, .8))
  expected <- simulate(reference, .5, seed = 921)
  expect_true(all(expected$unions$start_time > 1e-6))
  expect_gt(nrow(expected$unions), 10)
  withr::local_options(hammel.check_cells = TRUE)
  actual <- simulate(m, .5, seed = 921)
  expect_equal(actual$unions, expected$unions, tolerance = 1e-10)
  expect_equal(actual$persons, expected$persons)
  expect_true(all(c(".log_weight_A", ".log_weight_B") %in% names(actual$unions)))
})

test_that("active age and previous-union strata apply after a calendar change", {
  p <- rbind(transform(cohab_calendar(c(0, 1)), previous_union = "none"),
    transform(cohab_calendar(c(1, 0)), previous_union = "separated"),
    transform(cohab_calendar(c(0, 0)), previous_union = "widowed"))
  m <- cohab_calendar_model(p, n = 4, search = 1e5)
  m$population$unions <- tibble::tibble(union_id = 1L, wife_id = 1L, husband_id = 2L,
    start_time = -2, end_time = -1, end_reason = "dissolution")
  # Both searches start after .5, when the women reach age 20.5.
  m$processes$union_formation$rates$age_start <- 20.5
  out <- simulate(m, .6, seed = 45)
  u <- out$unions[out$unions$start_time >= 0, ]
  expect_equal(nrow(u), 2)
  expect_equal(u$union_type[match(c(1, 3), u$wife_id)], c("marriage", "cohab"))
  # An uncovered age selects probability zero even in the active calendar slice.
  m <- cohab_calendar_model(transform(cohab_calendar(c(1, 1)), age_start = 50), search = 1e5)
  expect_equal(simulate(m, .6, seed = 1)$unions$union_type, "marriage")
})

test_that("probability calendars handle zero/empty runs, restarts and RNG isolation", {
  m <- cohab_calendar_model(n = 100)
  set.seed(731); before <- .Random.seed
  out <- simulate(m, 1.2, seed = 92)
  expect_identical(.Random.seed, before)
  expect_identical(simulate(m, 1.2, seed = 92), out)
  expect_equal(nrow(simulate(m, 0, seed = 92)$unions), 0)
  expect_equal(nrow(simulate(cohab_calendar_model(n = 0), 1.2, seed = 92)$persons), 0)
  early <- simulate(m, .5, seed = 33)
  restarted <- simulate(set_population(m, early), .5, seed = 34)
  old <- nrow(early$unions)
  expect_equal(restarted$unions[seq_len(old), names(early$unions)], early$unions)
  added <- restarted$unions[restarted$unions$start_time >= .5, ]
  expect_gt(nrow(added), 0)
  expect_true(all(added$union_type == "cohab"))
})

test_that("shared and asynchronous rate/probability boundaries apply the current slices", {
  n <- 3000
  m <- cohab_calendar_model(cohab_calendar(c(0, 1, 0), c(0, .3, .7, Inf)),
    n = 2 * n, search = 0)
  m$processes$union_formation$rates <- data.frame(age_start = 0, age_end = Inf,
    rate = c(0, .8), time_start = c(0, .3), time_end = c(.3, Inf))
  out <- simulate(m, 1.2, seed = 593)
  u <- out$unions
  expect_true(all(u$start_time >= .3))
  expect_gt(sum(u$start_time < .7), 100)
  expect_gt(sum(u$start_time >= .7), 100)
  expect_true(all(u$union_type == ifelse(u$start_time < .7, "cohab", "marriage")))
  expected <- 1 - exp(-.8 * .9)
  expect_lt(abs(nrow(u) / n - expected), 5 * sqrt(expected * (1 - expected) / n))
})
