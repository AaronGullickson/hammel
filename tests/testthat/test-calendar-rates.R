# Synthetic calendar schedules: independent of database preparation.
calendar_epochs <- function(rates, values, times = c(0, .5, Inf)) {
  do.call(rbind, lapply(seq_along(values), function(i) {
    x <- rates
    x$rate <- values[i]
    x$time_start <- times[i]
    x$time_end <- times[i + 1L]
    x
  }))
}

calendar_pairs <- function(n = 2000, cohab = FALSE) {
  m <- make_rate_model(2 * n, search_rate = 0)
  m$population$persons$sex <- rep(c("f", "m"), each = n)
  u <- data.frame(union_id = seq_len(n), wife_id = seq_len(n),
    husband_id = n + seq_len(n), start_time = -3, end_time = NA_real_,
    end_reason = NA_character_)
  if (cohab) {
    u$union_type <- "cohab"
    m <- set_union_formation(m, "rate", rates = m$processes$union_formation$rates,
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = 0))
  }
  m$population <- hammel_population(m$population$persons, u, end_time = 0, groups = if ("group" %in% names(m$population$persons)) "group")
  m
}

test_that("calendar bounds validate and retain source tables", {
  r <- calendar_epochs(data.frame(sex = c("m", "f"), age_start = 0,
    age_end = Inf, rate = 0), c(.1, .2))
  attr(r, "source") <- list(source = "synthetic")
  p <- mortality(r)
  expect_identical(p$source, attr(r, "source"))
  expect_equal(p$rates$time_start, r$time_start)
  for (bad in list(r[setdiff(names(r), "time_end")],
      transform(r, time_start = NA_real_), transform(r, time_start = Inf),
      transform(r, time_end = time_start), transform(r, time_end = NA_real_),
      transform(r, time_start = as.character(time_start)), rbind(r, r[1, ])))
    expect_error(mortality(bad))
  future <- r
  future$rate[nrow(future)] <- -1
  expect_error(mortality(future), "nonnegative")
  expect_s3_class(mortality(transform(r, time_start = time_start - 2000,
    time_end = time_end - 2000)), "hammel_process")
})

test_that("overlap checks combine calendar and demographic bands", {
  r <- data.frame(sex = rep(c("m", "f"), 2), age_start = c(0, 0, 20, 20),
    age_end = c(20, 20, Inf, Inf), rate = .1,
    time_start = c(0, 0, .2, .2), time_end = Inf)
  expect_s3_class(mortality(r), "hammel_process")
  # Simultaneous rows for disjoint ages are valid; overlapping ages are not.
  expect_error(mortality(transform(r, age_start = 0)), "overlap")
  expect_error(run_simulation(set_mortality(make_model(1), r), 1), "cover")
})

test_that("run coverage, zero years and endpoint exclusion use fractional times", {
  r <- calendar_epochs(data.frame(sex = c("m", "f"), age_start = 0,
    age_end = Inf, rate = 0), c(0, 1e5), c(2015.5, 2019.8, 2021))
  m <- make_model(3)
  m$population$end_time <- 2015.5
  m <- set_mortality(m, r)
  expect_equal(run_simulation(m, 4.3, seed = 9)$simulation$deaths, 0)
  expect_equal(run_simulation(m, 4.31, seed = 9)$simulation$deaths, 3)
  expect_equal(run_simulation(m, 0, seed = 9)$end_time, 2015.5)
  expect_error(run_simulation(m, 6), "calendar coverage gap")
  gap <- r
  gap$time_start[3:4] <- 2020
  expect_error(run_simulation(set_mortality(m, gap), 5),
    "calendar coverage gap")
  # A gap outside this run is allowed.
  expect_equal(run_simulation(set_mortality(m, gap), 1, seed =
    1)$simulation$deaths, 0)
  m$population$end_time <- 2021
  expect_error(run_simulation(m, 0), "calendar coverage gap")
})

test_that("mortality integrates fractional calendar and age boundaries", {
  n <- 16000
  r <- calendar_epochs(data.frame(sex = c("m", "f"), age_start = 0,
    age_end = Inf, rate = 0), c(.1, .4), c(0, .4, Inf))
  m <- set_mortality(make_model(n), r)
  out <- run_simulation(m, 1.3, seed = 912)
  expected <- exp(-.1 * .4 - .4 * .9)
  expect_lt(abs(mean(is.na(out$persons$death_time)) - expected),
    5 * sqrt(expected * (1 - expected) / n))
  expect_identical(run_simulation(m, 1.3, seed = 912), out)
  set.seed(45)
  before <- .Random.seed
  run_simulation(m, .6, seed = 21)
  expect_identical(.Random.seed, before)
  # Sampling a candidate beyond the boundary must not execute an old hazard.
  zero_after <- calendar_epochs(m$processes$mortality$rates[1:2,
    c("sex", "age_start", "age_end", "rate")], c(.4, 0))
  out <- run_simulation(set_mortality(m, zero_after), 2.5, seed = 92)
  expect_true(all(out$persons$death_time[!is.na(out$persons$death_time)] < .5))
  expect_equal(run_simulation(make_model(2), .75, seed = 2)$simulation$deaths,
    0)
})

test_that("fertility changes retain spacing, weights and count unrealized births", {
  # Women in spacing until 1.9 under rates .2 then .4 from .5 (calendar epochs).
  n <- 2000
  m <- make_model(2, sex = "f", age = 30, birth_interval = 2) |>
    set_population(spacing_population(n))
  r <- calendar_epochs(data.frame(age_start = 15, age_end = 50, rate = 0), c(.2, .4))
  m <- set_fertility(m, r, birth_interval = 2, hetfert = het_fert(cv = 0, alpha = 1))
  m$population$persons$hetfert <- c(rep(2, n), rep(NA, n))
  out <- run_simulation(m, 1, seed = 81)
  u <- out$unrealized_births
  expect_count(u$unrealized_births, (.2 * .5 + .4 * .5) * n)
  # Counts sit on the one-year grid from the run start, across the rate change.
  expect_equal(u$period_start, 0)
  expect_equal(u$period_end, 1)
  expect_equal(out$persons$hetfert[1], 2)
  expect_equal(out$simulation$births, 0)
  expect_equal(out$persons$mother_id, m$population$persons$mother_id)
  restarted <- run_simulation(set_population(m, out), .5, seed = 82)
  # The restart keeps the first run's records and appends its own.
  record <- restarted$unrealized_births
  expect_identical(record[seq_len(nrow(u)), ], u)
  expect_equal(restarted$unrealized_since, 0)
  later <- record[record$period_start >= 1, ]
  expect_count(later$unrealized_births, .4 * .5 * n)
  # Adjusted rates now cover periods before the restart.
  rates <- calculate_fertility_rates(
    restarted, c(0, 1, 1.5), c(15, 50), unrealized = TRUE
  )
  by_period <- tapply(rates$unrealized_births, rates$time_start, sum)
  expect_equal(as.vector(by_period), c(
    sum(record$unrealized_births[record$period_start < 1]),
    sum(later$unrealized_births)
  ))
})

test_that("cell birth counts recover rates on either side of a calendar change", {
  n <- 4000
  m <- make_model(n, sex = "f", birth_interval = .2)
  m <- set_fertility(m, calendar_epochs(data.frame(age_start = 15, age_end = 50,
    rate = 0), c(.2, .6)), birth_interval = .2)
  withr::local_options(hammel.check_cells = TRUE)
  out <- run_simulation(m, 1.25, seed = 333)
  children <- out$persons[out$persons$mother_id %in% seq_len(n), ]
  counts <- c(sum(children$birth_time < .5), sum(children$birth_time >= .5))
  means <- n * c(.2 * .5, .6 * .75)
  expect_true(all(abs(counts - means) < 5 * sqrt(means)))
  expect_true(all(vapply(split(children$birth_time, children$mother_id),
    function(times) all(diff(sort(times)) >= .2 - 1e-12), logical(1))))
})

test_that("duration search clocks retain their origin across rate boundaries", {
  m <- calendar_pairs(1)
  m$population$unions$end_time <- -.2
  m$population$unions$end_reason <- "dissolution"
  r <- data.frame(previous_union = c("none", "separated", "widowed"),
    age_start = c(0, NA, 0), age_end = c(Inf, NA, Inf),
    dur_start = c(NA, 1, NA), dur_end = c(NA, Inf, NA), rate = 0)
  r <- calendar_epochs(r, c(0, 1e5))
  m <- set_union_formation(m, "rate", rates = r)
  out <- run_simulation(m, 1, seed = 122)
  expect_equal(nrow(out$unions), 2)
  expect_gte(out$unions$start_time[2], .8)
  expect_lt(out$unions$start_time[2], .81)
  bad <- r
  bad$age_start[5] <- 0
  bad$age_end[5] <- Inf
  bad$dur_start[5] <- bad$dur_end[5] <- NA_real_
  expect_error(set_union_formation(m, "rate", rates = bad), "same age or duration clock")
})

test_that("search intensity changes recover analytical formation probability", {
  n <- 5000
  m <- make_rate_model(2 * n)
  m$population$persons$sex <- rep(c("f", "m"), each = n)
  m <- set_union_formation(m, "rate", rates = calendar_epochs(
    data.frame(age_start = 0, age_end = Inf, rate = 0), c(.2, .8)))
  out <- run_simulation(m, 1.25, seed = 551)
  p <- 1 - exp(-.2 * .5 - .8 * .75)
  expect_lt(abs(nrow(out$unions) / n - p), 5 * sqrt(p * (1 - p) / n))
})

test_that("dissolution and conversion retain attained episode duration", {
  for (conversion in c(FALSE, TRUE)) {
    m <- calendar_pairs(1, cohab = conversion)
    # Existing episode starts at -3: attained duration reaches 3.75 at .75.
    r <- calendar_epochs(data.frame(dur_start = 3.75, dur_end = Inf, rate = 0),
      c(0, 1e5))
    if (conversion) {
      m$processes$union_formation$union_convert_rate <- r
    } else m <- set_union_dissolution(m, r)
    withr::local_options(hammel.check_cells = TRUE)
    out <- run_simulation(m, 1, seed = 15)
    expect_gte(out$unions$end_time[1], .75)
    expect_lt(out$unions$end_time[1], .76)
    expect_equal(out$unions$end_reason[1], if (conversion) "conversion" else "dissolution")
    if (conversion) {
      expect_equal(out$unions$wife_id, c(1, 1))
      expect_equal(out$unions$start_time[2], out$unions$end_time[1])
    }
  }
})

test_that("dissolution and conversion recover integrated calendar hazards", {
  n <- 5000
  for (conversion in c(FALSE, TRUE)) {
    m <- calendar_pairs(n, cohab = conversion)
    r <- calendar_epochs(data.frame(dur_start = 0, dur_end = Inf, rate = 0), c(.2, .6))
    if (conversion) m$processes$union_formation$union_convert_rate <- r
    else m <- set_union_dissolution(m, r)
    out <- run_simulation(m, 1.25, seed = 719)
    count <- if (conversion) out$simulation$conversions else out$simulation$dissolutions
    p <- 1 - exp(-.2 * .5 - .6 * .75)
    expect_lt(abs(count / n - p), 5 * sqrt(p * (1 - p) / n))
  }
})

test_that("directed transitions compete after a fractional rate boundary", {
  n <- 12000
  m <- make_model(n)
  m$population$persons$group <- "A"
  m$population$groups <- list(group = group(ancestry = TRUE))
  r <- data.frame(origin = "A", destination = c("B", "C"),
    age_start = 0, age_end = Inf, rate = 0)
  r <- calendar_epochs(r, c(0, 1), c(0, .3, Inf))
  r$rate[3:4] <- c(.2, .3)
  out <- run_simulation(set_group_transition(m, r), 1.3, seed = 91)
  expect_true(all(out$transitions$time >= .3))
  expected <- c(exp(-.5), .4 * (1 - exp(-.5)), .6 * (1 - exp(-.5)))
  actual <- as.numeric(table(factor(out$persons$group, levels = c("A", "B", "C")))) / n
  expect_true(all(abs(actual - expected) < 5 * sqrt(expected * (1 - expected) / n)))
  expect_equal(out$ancestry$group$A, rep(1, n))
})

test_that("grouped schedules require future group and pair coverage", {
  m <- make_model(2)
  m$population$persons$group <- c("A", "B")
  m$population$groups <- list(group = group())
  r <- calendar_epochs(data.frame(group = rep(c("A", "B"), each = 2),
    sex = rep(c("m", "f"), 2), age_start = 0, age_end = Inf, rate = 0), c(0, 0))
  r <- r[!(r$time_start == .5 & r$group == "B"), ]
  expect_error(run_simulation(set_mortality(m, r), 1), "missing schedules.*B")
  d <- calendar_epochs(directional_pairs(data.frame(group_pair = c("A:A", "A:B", "B:B"),
    dur_start = 0, dur_end = Inf, rate = 0)), c(0, 0))
  d <- d[!(d$time_start == .5 & d$group_w == "A" & d$group_h == "B"), ]
  expect_error(run_simulation(set_union_dissolution(m, d), 1),
    "missing schedules.*group_w=A, group_h=B")
})

test_that("probability calendar columns are accepted and empty populations cross changes", {
  m <- make_rate_model()
  prob <- data.frame(age_start = 0, age_end = Inf, prob = .5,
    time_start = 0, time_end = Inf)
  expect_s3_class(set_union_formation(m, "rate", rates = m$processes$union_formation$rates,
    cohab_prob = prob), "hammel_model")
  m <- make_model(0)
  m <- set_mortality(m, calendar_epochs(m$processes$mortality$rates, c(.2, .3)))
  expect_equal(nrow(run_simulation(m, 1.3, seed = 1)$persons), 0)
})

test_that("calendar changes precede exact person and cell clock ties", {
  for (birth in c(FALSE, TRUE)) {
    set.seed(682)
    boundary <- rexp(1) / 10
    m <- make_model(1, sex = if (birth) "f" else "m", birth_interval = 0)
    if (birth) m <- set_fertility(m, calendar_epochs(m$processes$fertility$rates,
      c(10, 0), c(0, boundary, Inf)), birth_interval = 0)
    else m <- set_mortality(m, calendar_epochs(m$processes$mortality$rates,
      c(10, 0), c(0, boundary, Inf)))
    out <- run_simulation(m, boundary + .1, seed = 682)
    expect_equal(out$simulation$births, 0)
    expect_equal(out$simulation$deaths, 0)
  }
})

test_that("group-duration transitions preserve latest-entry origin at boundaries", {
  m <- make_model(1)
  m$population$persons$group <- "A"
  m$population$groups <- list(group = group())
  r <- data.frame(origin = "A", destination = "B", dur_start = 20.8,
    dur_end = Inf, rate = 0)
  m <- set_group_transition(m, calendar_epochs(r, c(0, 1e5)))
  out <- run_simulation(m, 1, seed = 44)
  expect_gte(out$transitions$time, .8)
  expect_lt(out$transitions$time, .81)
})

test_that("asynchronous changes retain parity/type/group accounting and calibration", {
  n <- 400
  m <- make_rate_model(n)
  m$population$persons$sex <- rep(c("f", "m"), each = n / 2)
  m$population$persons$group <- rep(c("A", "B"), n / 2)
  m$population$groups <- list(group = group())
  f <- expand.grid(group = c("A", "B"), union_type = c("unpartnered", "partnered"),
    parity_start = c(0, 1), stringsAsFactors = FALSE)
  f$parity_end <- ifelse(f$parity_start == 0, 1, Inf)
  f$age_start <- 15
  f$age_end <- 50
  f$rate <- 0
  m <- set_fertility(m, calendar_epochs(f, c(.3, .6), c(0, .4, Inf)),
    birth_interval = .3, hetfert = het_fert(.5, .5))
  m <- set_union_formation(m, "rate", rates = calendar_epochs(
    data.frame(age_start = 15, age_end = Inf, rate = 0), c(.5, .9), c(0, .7, Inf)),
    diagnostics = TRUE, cancel_birth_if_no_match = FALSE)
  m <- set_mortality(m, calendar_epochs(m$processes$mortality$rates,
    c(0, .1), c(0, .9, Inf)))
  m <- set_union_dissolution(m, calendar_epochs(data.frame(dur_start = 0,
    dur_end = Inf, rate = 0), c(.2, .6), c(0, 1.1, Inf)))
  m <- set_group_transition(m, calendar_epochs(data.frame(origin = c("A", "B"),
    destination = c("B", "A"), age_start = 0, age_end = Inf, rate = 0),
    c(.2, .6), c(0, .2, Inf)))
  withr::local_options(hammel.check_cells = TRUE)
  out <- run_simulation(m, 1.5, seed = 13)
  expect_gt(out$simulation$births, 0)
  expect_gt(out$simulation$transitions, 0)
  expect_gt(out$simulation$dissolutions, 0)
  expect_true(all(is.finite(out$persons$hetfert[out$persons$sex == "f"])))
  expect_true(all(c(".log_weight_A", ".log_weight_B") %in% names(out$unions)))
  expect_identical(run_simulation(m, 1.5, seed = 13), out)
})

test_that("fertility can change parity partitions and pooled/typed union rates", {
  m <- calendar_pairs(1, cohab = TRUE)
  # Force no births before .5, then allow parity zero cohabitors to have one.
  before <- data.frame(age_start = 15, age_end = 50, parity_start = 0,
    parity_end = Inf, union_type = c("unpartnered", "partnered"), rate = 0,
    time_start = 0, time_end = .5)
  after <- expand.grid(parity_start = c(0, 1),
    union_type = c("unpartnered", "cohab", "marriage"), stringsAsFactors = FALSE)
  after$parity_end <- ifelse(after$parity_start == 0, 1, Inf)
  after$age_start <- 15
  after$age_end <- 50
  after$rate <- ifelse(after$union_type == "cohab" & after$parity_start == 0, 1e5, 0)
  after$time_start <- .5
  after$time_end <- Inf
  m <- set_fertility(m, rbind(before, after), birth_interval = .1)
  withr::local_options(hammel.check_cells = TRUE)
  out <- run_simulation(m, .7, seed = 124)
  expect_equal(out$simulation$births, 1)
  expect_gte(tail(out$persons$birth_time, 1), .5)
  expect_equal(tail(out$persons$mother_id, 1), 1)
})

test_that("rate changes recover rates in place and when band layouts change", {
  n <- 20000
  check <- function(rates) {
    out <- run_simulation(set_mortality(make_model(n, age = 30), rates), 2,
      seed = 71)
    d <- out$persons$death_time
    years <- function(a, b) sum(pmax(0, pmin(ifelse(is.na(d), 2, d), b) - a))
    deaths <- function(a, b) sum(!is.na(d) & d >= a & d < b)
    expect_count(deaths(0, 1), .05 * years(0, 1))
    expect_count(deaths(1, 2), .2 * years(1, 2))
  }
  withr::local_options(hammel.check_cells = TRUE)
  # Same bands in both periods: rates change in place.
  same <- data.frame(sex = rep(c("m", "f"), 2), age_start = 0, age_end = Inf,
    rate = c(.05, .05, .2, .2), time_start = rep(c(0, 1), each = 2), time_end = rep(c(1, Inf), each = 2))
  check(same)
  # Different band edges after the change: processes are rebuilt.
  split <- rbind(same[1:2, ],
    data.frame(sex = rep(c("m", "f"), 2), age_start = rep(c(0, 31.5), each = 2),
      age_end = rep(c(31.5, Inf), each = 2), rate = .2, time_start = 1, time_end = Inf))
  check(split)
})
