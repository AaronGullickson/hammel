dissolution_model <- function(n = 1L, duration = 0, rate = 0.1) {
  model <- make_rate_model(2L * n, search_rate = 0)
  p <- model$population$persons
  p$sex <- rep(c("f", "m"), n)
  model <- set_population(model, hammel_population(p, unions = data.frame(
    union_id = seq_len(n), wife_id = seq.int(1L, 2L*n, 2L),
    husband_id = seq.int(2L, 2L*n, 2L), start_time = -duration,
    end_time = NA_real_, end_reason = NA_character_), end_time = 0, groups = if ("group" %in% names(p)) "group"))
  set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = rate))
}

test_that("dissolution constructors validate duration schedules and retain provenance", {
  rates <- data.frame(dur_start = c(5, 0), dur_end = c(Inf, 5), rate = c(.1, 0), note = c("late", "early"))
  attr(rates, "source") <- list(description = "synthetic")
  p <- union_dissolution(rates)
  expect_s3_class(p, "hammel_process")
  expect_s3_class(p$rates, "tbl_df")
  expect_identical(p$rates$dur_start, rates$dur_start)
  expect_identical(p$source, attr(rates, "source"))
  expect_identical(set_process(hammel_model(), p), set_union_dissolution(hammel_model(), rates))
  for (column in c("sex", "age_start", "age_end", "union_status")) {
    bad <- rates; bad[[column]] <- 0
    expect_error(union_dissolution(bad), "duration bands|cannot be stratified")
  }
  for (column in c("time_start", "time_end")) {
    bad <- rates; bad[[column]] <- 0
    expect_error(union_dissolution(bad), "missing columns")
  }
  for (value in list(-1, NA_real_, Inf, "0", matrix(0, 2, 1))) {
    bad <- rates; bad$dur_start <- value
    expect_error(union_dissolution(bad))
  }
  for (value in list(-1, NA_real_, Inf, "0")) {
    bad <- rates; bad$rate <- value
    expect_error(union_dissolution(bad))
  }
  expect_error(union_dissolution(rates[FALSE, ]))
  expect_error(union_dissolution(transform(rates, dur_end = dur_start)))
  expect_error(union_dissolution(transform(rates[c("dur_start", "dur_end", "rate")], dur_end = Inf)), "overlap")
  # A custom column may be a group column, so it separates strata at construction;
  # at simulation only declared group dimensions stratify rates.
  noted <- union_dissolution(transform(rates, dur_end = Inf))
  model <- dissolution_model()
  model$processes$union_dissolution <- noted
  expect_error(simulate(model, 0), "overlapping bands within a stratum")
  model <- dissolution_model()
  model$processes$union_dissolution$rates$rate <- -1
  expect_error(simulate(model, 0), "duration bands")
  expect_error(set_process(model, model$processes$union_dissolution), "duration bands")
})

test_that("one annual dissolution hazard per couple matches analytical survival", {
  n <- 8000L
  model <- dissolution_model(n, duration = 8, rate = .2)
  set.seed(123); original_rng <- .Random.seed
  result <- simulate(model, 2, seed = 52)
  expected <- exp(-.4)
  survived <- is.na(result$unions$end_time)
  expect_lt(abs(mean(survived) - expected), 6 * sqrt(expected * (1 - expected) / n))
  expect_identical(.Random.seed, original_rng)
  expect_identical(result, simulate(model, 2, seed = 52))
  expect_equal(result$simulation$dissolutions, sum(!survived))
  expect_true(all(result$unions$end_reason[!survived] == "dissolution"))
  expect_true(all(result$unions$end_time[!survived] >= 0 & result$unions$end_time[!survived] < 2))
  expect_true(all(is.na(result$persons$death_time)))
  expect_equal(nrow(result$unions), n)
})

test_that("duration bands use union duration, gaps and attained restart duration", {
  model <- dissolution_model(100, duration = 0)
  model <- set_union_dissolution(model, data.frame(dur_start = 5, dur_end = Inf, rate = 1e4))
  early <- simulate(model, 4, seed = 52)
  expect_true(all(is.na(early$unions$end_time)))
  later <- simulate(set_population(model, early), 2, seed = 53)
  expect_true(all(later$unions$end_time >= 5 & later$unions$end_time < 6))
  expect_equal(later$simulation$dissolutions, 100)
  # Old input unions use duration six immediately, independently of age.
  model$population$unions$start_time <- -6
  old <- simulate(model, 1, seed = 52)
  expect_true(all(old$unions$end_time < 1))
  model$population$persons$birth_time <- -40
  expect_identical(simulate(model, 1, seed = 52)$unions, old$unions)
})

test_that("birth ties retain the pending dissolution and death wins same-person ties", {
  # Process clocks draw thresholds in refresh order: fertility, mortality (only
  # when positive), then dissolution. With one unit per process each due time is
  # its threshold over the rate. Ties sit at .5, away from integer-age cells.
  seed <- 71
  set.seed(seed); draws <- rexp(3)
  model <- dissolution_model(rate = draws[2] / .5)
  model <- set_fertility(model, data.frame(age_start = 15, age_end = 50,
    union_type = c("partnered", "unpartnered"), rate = c(draws[1] / .5, 0)), birth_interval = 0)
  tied <- simulate(model, 1.5, seed = seed)
  kids <- tied$persons[tied$persons$birth_time >= 0, ]
  expect_equal(nrow(kids), 1)
  expect_equal(kids$birth_time, .5, tolerance = 0)
  expect_equal(kids$father_id, 2)
  expect_equal(tied$unions$end_time, .5, tolerance = 0)
  expect_identical(tied$unions$end_reason, "dissolution")
  counts <- demographic_counts(tied, c(0, .5, 1.5), c(0, Inf), by = "union_type")
  maternal <- subset(counts, sex == "f" & time_start == .5)
  expect_equal(maternal$births[maternal$union_type == "partnered"], 1)
  expect_equal(maternal$births[maternal$union_type == "unpartnered"], 0)
  expect_equal(maternal$person_years[maternal$union_type == "partnered"], 0)
  # Observation endpoint is exclusive even for an exact tied event.
  endpoint <- simulate(model, .5, seed = seed)
  expect_equal(endpoint$simulation$births, 0)
  expect_equal(endpoint$simulation$dissolutions, 0)
  # Her positive death rate now draws second, so dissolution draws third.
  model <- set_mortality(model, data.frame(sex = c("f", "m"),
    age_start = 0, age_end = Inf, rate = c(draws[2] / .5, 0)))
  model$processes$union_dissolution$rates$rate <- draws[3] / .5
  dead <- simulate(model, 1.5, seed = seed)
  expect_equal(dead$simulation$births, 0)
  expect_equal(dead$simulation$dissolutions, 0)
  expect_identical(dead$unions$end_reason, "death_wife")
})

test_that("widowhood invalidates pending dissolution in both initiation modes", {
  # Zero rates draw no thresholds: the dissolution clock draws first (due at
  # 1), then the mortality clock when his positive-rate cell fills (due at .5).
  set.seed(71); draws <- rexp(2)
  for (mode in c("rate", "on_demand")) {
    model <- dissolution_model(rate = draws[1])
    if (mode == "on_demand") {
      model <- set_union_formation(model, mode)
      model <- set_fertility(model, data.frame(age_start = 15, age_end = 50, rate = 0))
    }
    model <- set_mortality(model, data.frame(sex = c("f", "m"),
      age_start = 0, age_end = Inf, rate = c(0, draws[2] * 2)))
    result <- simulate(model, 1.5, seed = 71)
    expect_equal(result$persons$death_time[2], .5, tolerance = 0)
    expect_identical(result$unions$end_reason, "death_husband")
    expect_equal(result$simulation$dissolutions, 0)
  }
})

test_that("births do not reset the union dissolution clock", {
  # A negligible positive rate gives the quiet model a cell clock too, so both
  # draw the initial dissolution after the same cell/death draws.
  model <- set_fertility(dissolution_model(), data.frame(age_start = 15, age_end = 50,
    union_type = c("partnered", "unpartnered"), rate = c(1e-12, 0)), birth_interval = 0)
  busy <- set_fertility(model, data.frame(age_start = 15, age_end = 50,
    union_type = c("partnered", "unpartnered"), rate = c(20, 0)), birth_interval = 0)
  quiet <- simulate(model, 10, seed = 52)
  result <- simulate(busy, 10, seed = 52)
  expect_gt(result$simulation$births, 0)
  expect_identical(result$unions, quiet$unions)
  births <- result$persons$birth_time[result$persons$birth_time >= 0]
  expect_true(all(births <= result$unions$end_time))
})

test_that("dissolution restores matching availability and repeat episodes start anew", {
  for (mode in c("on_demand", "rate")) {
    model <- make_model(2, birth_rate = 5, birth_interval = .1)
    model$population$persons$sex <- c("f", "m")
    if (mode == "rate") model <- make_rate_model(2, search_rate = 50, birth_rate = 5, birth_interval = .1)
    model <- set_union_dissolution(model, data.frame(dur_start = .1, dur_end = Inf, rate = 1e4))
    result <- simulate(model, 2, seed = 52)
    expect_gt(nrow(result$unions), 1)
    expect_true(all(result$unions$wife_id == 1 & result$unions$husband_id == 2))
    ended <- !is.na(result$unions$end_time)
    expect_true(all(result$unions$end_time[ended] - result$unions$start_time[ended] >= .1))
    expect_equal(result$simulation$dissolutions, sum(ended))
    expect_no_error(validate_population(result))
    births <- result$persons$birth_time[result$persons$mother_id %in% 1]
    expect_true(all(diff(births) >= .1 - 1e-12))
  }
})

test_that("dissolution uses directional group strata with future coverage", {
  model <- dissolution_model(4)
  model$population$persons$group <- c("A", "A", "A", "B", "B", "A", "B", "B")
  model$population$groups <- list(group = group())
  model <- set_population(model, model$population)
  # Symmetric rates: cross-group unions in either direction dissolve.
  rates <- directional_pairs(data.frame(group_pair = c("A:A", "A:B", "B:B"),
    dur_start = 0, dur_end = Inf, rate = c(0, 1e4, 0)))
  model <- set_union_dissolution(model, rates)
  result <- simulate(model, 1, seed = 52)
  expect_true(all(is.na(result$unions$end_time[c(1, 4)])))
  expect_true(all(result$unions$end_time[c(2, 3)] < 1))
  # Swapping wives' and husbands' identities keeps their (symmetric) risk.
  reversed <- model
  reversed$population$persons$group <- c("A", "A", "B", "A", "A", "B", "B", "B")
  reversed$population$groups <- list(group = group())
  # Directions are separate cells, so the draws differ but the same unions dissolve.
  expect_identical(!is.na(simulate(reversed, 1, seed = 52)$unions$end_time), !is.na(result$unions$end_time))
  # Directional rates: only one direction is at risk.
  persons <- result$persons
  wife <- persons$group[match(result$unions$wife_id, persons$person_id)]
  husband <- persons$group[match(result$unions$husband_id, persons$person_id)]
  one_way <- set_union_dissolution(model, transform(rates,
    rate = ifelse(group_w == "A" & group_h == "B", 1e4, 0)))
  out <- simulate(one_way, 1, seed = 52)
  expect_identical(!is.na(out$unions$end_time), wife == "A" & husband == "B")
  # Proportional strata from stratify_rates equal the explicit table.
  scaled <- set_union_dissolution(model, stratify_rates(
    data.frame(dur_start = 0, dur_end = Inf, rate = .2),
    list(group_w = c(A = 2, B = 1), group_h = c(A = 1, B = 3))))
  explicit <- set_union_dissolution(model, data.frame(group_w = c("A", "B", "A", "B"),
    group_h = c("A", "A", "B", "B"), dur_start = 0, dur_end = Inf, rate = .2 * c(2, 1, 6, 3)))
  expect_equal(histories(simulate(scaled, 5, seed = 52)), histories(simulate(explicit, 5, seed = 52)),
               tolerance = 1e-12)
  # Every combination must be covered, including future identities.
  partial <- set_union_dissolution(model, data.frame(group_w = "A", group_h = "A",
    dur_start = 0, dur_end = Inf, rate = 0))
  expect_error(simulate(partial, 0), "missing schedules.*group_w=B")
  model <- dissolution_model()
  model$population$persons$group <- "A"
  model$population$groups <- list(group = group())
  model <- set_inheritance(model, group_inheritance("custom", fun = function(...) "AB"), extra = "AB")
  model <- set_union_dissolution(model, data.frame(group_w = "A", group_h = "A",
    dur_start = 0, dur_end = Inf, rate = 0))
  expect_error(simulate(model, 0), "missing schedules.*AB")
})

test_that("optional process, zero-duration runs, empty populations and fallbacks work", {
  model <- dissolution_model()
  expect_equal(simulate(model, 0)$simulation$dissolutions, 0)
  model$processes$union_dissolution <- NULL
  absent <- simulate(model, 2, seed = 52)
  expect_null(absent$simulation$dissolutions)
  expect_true(all(is.na(absent$unions$end_time)))
  empty <- make_model(0)
  empty <- set_union_dissolution(empty, data.frame(dur_start = 0, dur_end = Inf, rate = .1))
  out <- simulate(empty, 2, seed = 52)
  expect_equal(out$simulation$dissolutions, 0)
  expect_s3_class(out$unions, "tbl_df")
  fallback <- set_union_dissolution(model, data.frame(group_pair = c("B:B", "A:A"), dur_start = 0,
    dur_end = Inf, rate = c(1e4, 0)))
  expect_error(simulate(fallback, 1, seed = 52), "not a declared group dimension")
})

test_that("dissolution histories conserve partnered exposure and survive later deaths", {
  model <- dissolution_model(100, rate = 2)
  result <- simulate(model, 1, seed = 52)
  counts <- demographic_counts(result, c(0, 1), c(0, Inf), by = "union_type")
  partnered <- subset(counts, union_type == "partnered")
  unpartnered <- subset(counts, union_type == "unpartnered")
  expected <- sum(ifelse(is.na(result$unions$end_time), 1, result$unions$end_time))
  expect_equal(partnered$person_years, rep(expected, 2))
  expect_equal(partnered$person_years + unpartnered$person_years, rep(100, 2))
  marriages <- calculate_union_rates(result, c(0, 1), c(0, Inf))
  expect_equal(sum(marriages$unions), 200) # Input formations at the window start.
  # Follow-up deaths must not rewrite already dissolved records or count them again.
  continuation <- set_population(model, result)
  continuation <- set_mortality(continuation,
    data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1e4))
  continued <- simulate(continuation, 1, seed = 53)
  ended <- !is.na(result$unions$end_time)
  expect_identical(continued$unions[ended, ], result$unions[ended, ])
})


test_that("directional group columns are validated", {
  rates <- data.frame(group_w = c("A", "B"), group_h = c("B", "A"), dur_start = c(0, 1),
                      dur_end = c(1, Inf), rate = c(.1, .2))
  expect_no_error(union_dissolution(rates))
  expect_error(union_dissolution(transform(rates, group_w = "A", group_h = "B", dur_start = 0)), "overlap")
  # Unordered pairs and person-level group columns are not union strata.
  model <- dissolution_model(4)
  model$population$persons$group <- c("A", "A", "A", "B", "B", "A", "B", "B")
  model$population$groups <- list(group = group())
  shared <- data.frame(dur_start = 0, dur_end = Inf, rate = .1)
  for (column in c("group_pair", "group")) {
    bad <- shared; bad[[column]] <- "A"
    expect_error(simulate(set_union_dissolution(model, bad), 0), "_w`")
  }
  expect_error(simulate(set_mortality(model, transform(model$processes$mortality$rates, group_w = "A")), 0),
               "own name")
  expect_error(stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = 1e308),
    list(group_w = c(A = 2))), "remain finite")
})
