# Weight models (milestone 8 step 2): selection within rate cells in proportion
# to exp(eta), without changing any cell's input intensity.

# Couples aged 70 (wives 1..n, husbands 2n+1..3n) plus n unpartnered people of
# each sex, with zero fertility and on-demand formation.
weighted_mortality_model <- function(n, groups = NULL) {
  persons <- data.frame(person_id = 1:(4 * n), sex = rep(c("f", "m"), each = 2 * n),
    birth_time = -70, death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  if (!is.null(groups)) persons$group <- rep(groups, length.out = 4 * n)
  unions <- data.frame(union_id = 1:n, wife_id = 1:n, husband_id = 2 * n + 1:n,
    start_time = -10, end_time = NA_real_, end_reason = NA_character_)
  hammel_model() |>
    set_population(hammel_population(persons, unions = unions, end_time = 0,
      groups = if ("group" %in% names(persons)) list(group = group(ancestry = TRUE)))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .05)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0)
}
mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .05)

test_that("weight models validate formulas, coefficients and variables", {
  expect_s3_class(weight_model(~ partnered, c(partnered = 1)), "hammel_weight_model")
  expect_s3_class(weight_model(fun = function(members, bin) 0), "hammel_weight_model")
  expect_error(weight_model(), "exactly one")
  expect_error(weight_model(~ parity, c(parity = 1), fun = function(m, b) 0), "exactly one")
  expect_error(weight_model(~ I(parity > 1), c(x = 1)), "custom function")
  expect_error(weight_model(~ log(parity), c(x = 1)), "custom function")
  expect_error(weight_model(~ parity - 1, c(parity = 1)), "intercept")
  expect_error(weight_model(~ parity, c(parity = Inf)), "finite")
  model <- weighted_mortality_model(2)
  # Names that could be group variables are checked against the declaration
  # at simulation; others fail when the model is set.
  education <- set_mortality(model, mortality_rates, weights = weight_model(~ education, c(education = 1)))
  expect_error(run_simulation(education, 0),
    "not available for this process: education")
  expect_error(set_mortality(model, mortality_rates,
    weights = weight_model(~ parity, c(paritx = 1))), "missing parity; unknown paritx")
  expect_error(set_mortality(model, mortality_rates,
    weights = weight_model(~ same_group, c(same_group = 1))), "not available")
  expect_error(set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .1),
    weights = weight_model(~ partnered, c(partnered = 1))), "not available")
  expect_error(set_mortality(model, mortality_rates, weights = list()), "weight_model")
  # Interactions follow R's model-matrix coding and names.
  w <- weight_model(~ union_type * parity,
    c(union_typecohab = 1, union_typemarriage = 2, parity = .1,
      "union_typecohab:parity" = .2, "union_typemarriage:parity" = .3))
  expect_silent(set_mortality(model, mortality_rates, weights = w))
})

test_that("terms the input rates already stratify by warn that they cancel", {
  model <- weighted_mortality_model(2)
  expect_warning(set_mortality(model, mortality_rates,
    weights = weight_model(~ sex + partnered, c(sexm = 1, partnered = 1))), "`sex` has no effect")
  status <- data.frame(union_type = c("partnered", "unpartnered"), age_start = 15, age_end = 50, rate = .1)
  expect_warning(set_fertility(model, status,
    weights = weight_model(~ partnered, c(partnered = 1))), "`partnered` has no effect")
  expect_warning(set_union_dissolution(model, data.frame(union_type = c("marriage", "cohab"),
    dur_start = 0, dur_end = Inf, rate = .1),
    weights = weight_model(~ union_type, c(union_typemarriage = 1))), "`union_type` has no effect")
  expect_silent(set_mortality(model, mortality_rates,
    weights = weight_model(~ sex:partnered, c("sexf:partnered" = 1, "sexm:partnered" = 1))))
})

test_that("mortality weights give within-cell hazard ratios and keep the cell rate", {
  n <- 20000
  model <- set_mortality(weighted_mortality_model(n), mortality_rates,
    weights = weight_model(~ partnered, c(partnered = log(0.5))))
  out <- run_simulation(model, 1, seed = 11)
  p <- out$persons
  d <- ifelse(is.na(p$death_time), 1, p$death_time)
  partner <- c(2 * n + 1:n, rep(NA, n), 1:n, rep(NA, n))
  # Exposure in the initial status ends at the first death in the couple.
  status_end <- ifelse(is.na(partner), d, pmin(d, d[partner]))
  died <- !is.na(p$death_time) & p$death_time <= status_end
  part <- !is.na(partner)
  rate <- function(sel) sum(died[sel]) / sum(status_end[sel])
  expect_lt(abs(log(rate(part) / rate(!part)) - log(0.5)), 4 * sqrt(1 / sum(died[part]) + 1 / sum(died[!part])))
  # The cell intensity is unchanged: total deaths match the input rate.
  deaths <- sum(!is.na(p$death_time))
  expected <- .05 * sum(d)
  expect_lt(abs(deaths - expected), 4 * sqrt(expected))
})

test_that("zero coefficients and constant custom weights reproduce unweighted runs", {
  model <- weighted_mortality_model(500)
  plain <- run_simulation(model, 2, seed = 5)
  zero <- set_mortality(model, mortality_rates,
    weights = weight_model(~ partnered, c(partnered = 0)))
  expect_identical(run_simulation(zero, 2, seed = 5), plain)
  constant <- set_mortality(model, mortality_rates,
    weights = weight_model(fun = function(members, bin) rep(3, nrow(members))))
  expect_identical(run_simulation(constant, 2, seed = 5), plain)
})

test_that("custom functions receive members and the cell, without exact times", {
  model <- weighted_mortality_model(200)
  seen <- NULL
  f <- function(members, bin) {
    if (is.null(seen)) seen <<- list(members = members, bin = bin)
    log(ifelse(members$partnered == 1, .5, 1))
  }
  out <- run_simulation(set_mortality(model, mortality_rates, weights =
    weight_model(fun = f)), 1, seed = 3)
  expect_s3_class(seen$members, "tbl_df")
  expect_true(all(c("person_id", "sex", "union_type", "partnered", "parity",
    "previous_union") %in% names(seen$members)))
  # Without declared groups there are no group columns.
  expect_false(any(c("group", "partner_group") %in% names(seen$members)))
  expect_false(any(c("birth_time", "age", "time") %in% names(seen$members)))
  expect_identical(seen$bin$process, "mortality")
  expect_true(all(c("sex", "age_start", "age_end", "time_start", "time_end") %in% names(seen$bin)))
  # The same model as a formula agrees in distribution. Formulas store weights
  # and pick by rejection sampling; custom functions scan the chosen cell.
  formula <- set_mortality(model, mortality_rates,
    weights = weight_model(~ partnered, c(partnered = log(.5))))
  share <- function(x) {
    dead <- x$persons[!is.na(x$persons$death_time), ]
    mean(dead$person_id %in% c(seq_len(200), 400 + seq_len(200)))   # partnered share
  }
  runs <- 1:15
  custom_share <- vapply(runs, function(s) {
    share(run_simulation(set_mortality(
      model, mortality_rates,
      weights = weight_model(fun = f)), 1, seed = s))
  }, 0)
  formula_share <- vapply(runs, function(s) {
    share(run_simulation(formula, 1,
      seed = s))
  }, 0)
  expect_lt(abs(mean(custom_share) - mean(formula_share)),
            4 * sqrt(var(custom_share) / 15 + var(formula_share) / 15) + 1e-9)
  bad_length <- weight_model(fun = function(members, bin) 0)
  expect_error(run_simulation(set_mortality(model, mortality_rates, weights =
    bad_length), 1, seed = 3),
    "returned 1 values")
  infinite <- weight_model(fun = function(members, bin) rep(Inf, nrow(members)))
  expect_error(run_simulation(set_mortality(model, mortality_rates, weights =
    infinite), 1, seed = 3), "infinite")
  # Draws inside the function share the simulation stream reproducibly.
  noisy <- set_mortality(model, mortality_rates,
    weights = weight_model(fun = function(members, bin) rnorm(nrow(members))))
  expect_identical(run_simulation(noisy, 1, seed = 8), run_simulation(noisy, 1,
    seed = 8))
})

test_that("missing variables contribute zero with one warning", {
  model <- weighted_mortality_model(200, groups = c("A", "B"))
  w <- weight_model(~ partner_group + partnered,
    c(partner_groupB = 1, partnered = 1))
  expect_warning(run_simulation(set_mortality(model, mortality_rates, weights =
    w), 1, seed = 2),
    "Missing values in the mortality weight formula")
  na_fun <- weight_model(fun = function(members, bin) rep(NA_real_, nrow(members)))
  model <- weighted_mortality_model(200)
  expect_warning(out <- run_simulation(set_mortality(model, mortality_rates,
    weights = na_fun), 1, seed = 2),
    "function result")
  expect_identical(out, run_simulation(model, 1, seed = 2))
  # Coefficients name the dimension's categories (the first is the reference).
  partial <- weight_model(~ group, c(groupC = 1))
  grouped <- weighted_mortality_model(20, groups = c("A", "B"))
  expect_error(run_simulation(set_mortality(grouped, mortality_rates, weights =
    partial), 1),
    "groupB; missing groupB; unknown groupC")
})

# Women aged 25 (one fertility cell for a year), with birth interval zero.
weighted_fertility_model <- function(n, groups = NULL) {
  persons <- data.frame(person_id = 1:(2 * n), sex = rep(c("f", "m"), each = n),
    birth_time = -25.2, death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  if (!is.null(groups)) persons$group <- rep(groups, length.out = 2 * n)
  hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0, cancel_birth_if_no_match = FALSE)
}

test_that("fertility weights select mothers at selection and at refresh under calibration", {
  rates <- data.frame(age_start = 15, age_end = 50, rate = .5)
  check <- function(model, n) {
    out <- run_simulation(model, .7, seed = 21)
    p <- out$persons
    births <- table(factor(p$group[match(p$mother_id[p$birth_time >= 0], p$person_id)],
                           levels = c("A", "B")))
    expected <- .5 * .7 * n
    expect_lt(abs(sum(births) - expected), 4 * sqrt(expected))
    expect_lt(abs(log(births[["B"]] / births[["A"]]) - log(2)),
              4 * sqrt(1 / births[["A"]] + 1 / births[["B"]]))
  }
  w <- weight_model(~ group, c(groupB = log(2)))
  n <- 4000
  # Two groups activate weighted matching calibration: weights stored at refresh.
  check(set_fertility(weighted_fertility_model(n, c("A", "B")), rates, birth_interval = 0, weights = w), n)
  # Calibration off (matching restricted to the woman's own group): selection time.
  restricted <- weighted_fertility_model(n, c("A", "B"))
  restricted <- set_union_formation(restricted, "on_demand", kinship_generations = 0,
    cancel_birth_if_no_match = FALSE,
    matching = matching_model(~ age_difference, c(age_difference = 0), group_exogamy = NA))
  check(set_fertility(restricted, rates, birth_interval = 0, weights = w), n)
  # A fertility weight model combines with hetfert, and zero coefficients keep it.
  model <- set_fertility(weighted_fertility_model(300, c("A", "B")), rates, hetfert = het_fert(.5, 1))
  zero <- set_fertility(model, rates, hetfert = het_fert(.5, 1),
    weights = weight_model(~ group, c(groupB = 0)))
  expect_identical(run_simulation(zero, 1, seed = 4), run_simulation(model, 1,
    seed = 4))
  expect_warning(run_simulation(set_fertility(model, rates, weights =
    weight_model(
    fun = function(members, bin) rep(0, nrow(members)))), .1, seed = 4), "every refresh")
})

# Marriages with an age gap of 0 or 5 years, no mortality or fertility.
weighted_union_model <- function(n, type = "marriage") {
  persons <- data.frame(person_id = 1:(2 * n), sex = rep(c("f", "m"), each = n),
    birth_time = c(rep(-30, n), rep(c(-30, -35), length.out = n)), death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_)
  unions <- data.frame(union_id = 1:n, wife_id = 1:n, husband_id = n + 1:n, start_time = -1,
    end_time = NA_real_, end_reason = NA_character_, union_type = type)
  hammel_model() |>
    set_population(hammel_population(persons, unions = unions, end_time = 0,
      groups = if ("group" %in% names(persons)) list(group = group(ancestry = TRUE)))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(union_type = c("partnered", "unpartnered"), age_start = 15,
      age_end = 50, rate = 0), birth_interval = 0)
}
gap_ratio <- function(unions, n, event) {
  gap <- rep(c(0, 5), length.out = n)[seq_len(n)]
  first <- unions[seq_len(n), ]
  ended <- !is.na(first$end_reason) & first$end_reason == event
  years <- ifelse(is.na(first$end_time), 1, first$end_time)
  rate <- function(g) sum(ended[gap == g]) / sum(years[gap == g])
  c(log_ratio = log(rate(5) / rate(0)), se = sqrt(1 / sum(ended[gap == 5]) + 1 / sum(ended[gap == 0])))
}

test_that("dissolution and conversion weights act on unions", {
  n <- 6000
  w <- weight_model(~ age_difference, c(age_difference = log(1.2)))
  model <- weighted_union_model(n) |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      kinship_generations = 0) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .3), weights = w)
  out <- suppressWarnings(run_simulation(model, 1, seed = 31))
  r <- gap_ratio(out$unions, n, "dissolution")
  expect_lt(abs(r[["log_ratio"]] - 5 * log(1.2)), 4 * r[["se"]])
  cohabs <- weighted_union_model(n, "cohab") |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3),
      kinship_generations = 0, weights_conversion = w)
  out <- suppressWarnings(run_simulation(cohabs, 1, seed = 32))
  r <- gap_ratio(out$unions, n, "conversion")
  expect_lt(abs(r[["log_ratio"]] - 5 * log(1.2)), 4 * r[["se"]])
})

test_that("transition weights can use person variables and the destination", {
  n <- 4000
  persons <- data.frame(person_id = 1:n, sex = rep(c("f", "m"), length.out = n), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_, group = "A")
  model <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0)
  rates <- data.frame(origin = "A", destination = c("B", "C"), age_start = 0, age_end = Inf, rate = .2)
  w <- weight_model(~ sex:destination,
    c("sexf:destinationA" = 0, "sexm:destinationA" = 0, "sexf:destinationB" = log(2),
      "sexm:destinationB" = 0, "sexf:destinationC" = 0, "sexm:destinationC" = 0))
  out <- run_simulation(set_group_transition(model, rates, weights = w), 1,
    seed = 41)
  tx <- out$transitions
  sex <- out$persons$sex[match(tx$person_id, out$persons$person_id)]
  b <- table(factor(sex[tx$destination == "B"], levels = c("f", "m")))
  c_moves <- table(factor(sex[tx$destination == "C"], levels = c("f", "m")))
  # Women are twice as likely as men to be chosen for a move to B, not to C.
  expect_lt(abs(log(b[["f"]] / b[["m"]]) - log(2)), 4 * sqrt(1 / b[["f"]] + 1 / b[["m"]]))
  expect_lt(abs(log(c_moves[["f"]] / c_moves[["m"]])), 4 * sqrt(1 / c_moves[["f"]] + 1 / c_moves[["m"]]))
  expect_warning(group_transition(rates, weights = weight_model(~ destination,
    c(destinationB = 1, destinationC = 1))), "no effect")
})

# Unpartnered women aged 22 in two groups, with no mortality or fertility.
weighted_search_model <- function(n) {
  persons <- data.frame(person_id = 1:(2 * n), sex = rep(c("f", "m"), each = n), birth_time = -22,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    group = rep(c("A", "B"), length.out = 2 * n))
  hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(union_type = c("partnered", "unpartnered"), age_start = 15,
      age_end = 50, rate = 0), birth_interval = 0)
}
first_union_ratio <- function(out) {
  u <- out$unions
  first <- tapply(u$start_time, u$wife_id, min)
  wives <- out$persons$person_id[out$persons$sex == "f"]
  start <- first[as.character(wives)]
  formed <- !is.na(start)
  start[!formed] <- out$simulation$end_time - out$simulation$start_time
  group <- out$persons$group[match(wives, out$persons$person_id)]
  events <- tapply(formed, group, sum)
  years <- tapply(start, group, sum)
  c(log_ratio = log((events[["B"]] / years[["B"]]) / (events[["A"]] / years[["A"]])),
    se = sqrt(1 / events[["A"]] + 1 / events[["B"]]))
}

test_that("search weights select searchers with and without matching calibration", {
  n <- 6000
  rates <- data.frame(age_start = 15, age_end = Inf, rate = .3)
  w <- weight_model(~ group, c(groupB = log(2)))
  # Two groups activate calibration: weights stored, calibration totals weighted.
  withr::local_options(hammel.check_cells = TRUE)
  calibrated <- set_union_formation(weighted_search_model(n), "rate", rates = rates,
    kinship_generations = 0, weights = w)
  r <- first_union_ratio(suppressWarnings(run_simulation(calibrated, 1, seed =
    3)))
  expect_lt(abs(r[["log_ratio"]] - log(2)), 4 * r[["se"]])
  # Matching within groups only: no calibration, weights at selection.
  within <- matching_model(~ age_difference, c(age_difference = 0), group_exogamy = NA)
  restricted <- set_union_formation(weighted_search_model(n), "rate", rates = rates,
    matching = within, kinship_generations = 0, weights = w)
  r <- first_union_ratio(suppressWarnings(run_simulation(restricted, 1, seed =
    4)))
  expect_lt(abs(r[["log_ratio"]] - log(2)), 4 * r[["se"]])
  withr::local_options(hammel.check_cells = FALSE)
  plain <- set_union_formation(weighted_search_model(500), "rate", rates = rates, kinship_generations = 0)
  zero <- set_union_formation(weighted_search_model(500), "rate", rates = rates, kinship_generations = 0,
    weights = weight_model(~ group, c(groupB = 0)))
  expect_identical(suppressWarnings(run_simulation(zero, 1, seed = 5)),
    suppressWarnings(run_simulation(plain, 1, seed = 5)))
})

test_that("search weights require rate mode and warn about cancelled terms and slow paths", {
  model <- weighted_search_model(50)
  w <- weight_model(~ parity, c(parity = 1))
  expect_error(set_union_formation(model, "on_demand", weights = w), "rate-based")
  rates <- data.frame(age_start = 15, age_end = Inf, rate = .3)
  expect_warning(set_union_formation(model, "rate", rates = rates,
    weights = weight_model(~ partnered, c(partnered = 1))), "`partnered` has no effect")
  slow <- set_union_formation(model, "rate", rates = rates, kinship_generations = 0,
    weights = weight_model(fun = function(members, bin) rep(0, nrow(members))))
  expect_warning(run_simulation(slow, .2, seed = 1),
    "every refresh of every searching woman")
})

test_that("weight schedules switch models at calendar boundaries without changing clocks", {
  n <- 20000
  half <- weight_model(~ partnered, c(partnered = log(.5)))
  off <- weight_model(~ partnered, c(partnered = 0))
  schedule <- function(a, b) tibble::tibble(time_start = c(0, 1), time_end = c(1, Inf), model = list(a, b))
  model <- weighted_mortality_model(n)
  # Identical models in both intervals reproduce the static model exactly.
  small <- weighted_mortality_model(300)
  expect_identical(run_simulation(set_mortality(small, mortality_rates,
    weights = schedule(half, half)), 2, seed = 4),
                   run_simulation(set_mortality(small, mortality_rates,
                     weights = half), 2, seed = 4))
  out <- run_simulation(set_mortality(model, mortality_rates, weights =
    schedule(half, off)), 2, seed = 4)
  p <- out$persons
  partner <- c(2 * n + 1:n, rep(NA, n), 1:n, rep(NA, n))
  part <- !is.na(partner)
  d <- ifelse(is.na(p$death_time), 2, p$death_time)
  status_end <- ifelse(is.na(partner), d, pmin(d, d[partner]))
  log_ratio <- function(a, b) {
    years <- pmax(0, pmin(status_end, b) - a)
    died <- !is.na(p$death_time) & p$death_time >= a & p$death_time < b & p$death_time <= status_end
    c(log(sum(died[part]) / sum(years[part]) / (sum(died[!part]) / sum(years[!part]))),
      sqrt(1 / sum(died[part]) + 1 / sum(died[!part])))
  }
  before <- log_ratio(0, 1)
  after <- log_ratio(1, 2)
  expect_lt(abs(before[1] - log(.5)), 4 * before[2])
  expect_lt(abs(after[1]), 4 * after[2])
  # Fertility weights stored under matching calibration also switch.
  rates <- data.frame(age_start = 15, age_end = 50, rate = .5)
  b2 <- weight_model(~ group, c(groupB = log(2)))
  b1 <- weight_model(~ group, c(groupB = 0))
  fert <- set_fertility(weighted_fertility_model(4000, c("A", "B")), rates, birth_interval = 0,
    weights = tibble::tibble(time_start = c(0, .35), time_end = c(.35, Inf), model = list(b2, b1)))
  out <- run_simulation(fert, .7, seed = 21)
  kids <- out$persons[out$persons$birth_time >= 0, ]
  group <- out$persons$group[match(kids$mother_id, out$persons$person_id)]
  early <- table(factor(group[kids$birth_time < .35], levels = c("A", "B")))
  late <- table(factor(group[kids$birth_time >= .35], levels = c("A", "B")))
  expect_lt(abs(log(early[["B"]] / early[["A"]]) - log(2)), 4 * sqrt(sum(1 / early)))
  expect_lt(abs(log(late[["B"]] / late[["A"]])), 4 * sqrt(sum(1 / late)))
  # Custom-function schedules switch functions.
  fun_schedule <- tibble::tibble(time_start = c(0, 1), time_end = c(1, Inf), model = list(
    weight_model(fun = function(members, bin) rep(0, nrow(members))),
    weight_model(fun = function(members, bin) log(ifelse(members$partnered == 1, .5, 1)))))
  expect_s3_class(run_simulation(set_mortality(small, mortality_rates, weights =
    fun_schedule), 2, seed = 4),
                  "hammel_population")
})

test_that("weight schedules are validated like model schedules", {
  model <- weighted_mortality_model(20)
  half <- weight_model(~ partnered, c(partnered = log(.5)))
  custom <- weight_model(fun = function(members, bin) 0)
  expect_error(set_mortality(model, mortality_rates, weights = tibble::tibble(
    time_start = c(0, 1), time_end = c(1, Inf), model = list(half, custom))), "all be formulas")
  expect_error(set_mortality(model, mortality_rates, weights = tibble::tibble(
    time_start = c(0, .5), time_end = c(1, Inf), model = list(half, half))), "overlap")
  expect_error(set_mortality(model, mortality_rates, weights = tibble::tibble(
    time_start = 0, time_end = Inf, model = list(weight_model(~ parity, c(paritx = 1))))),
    "weights\\$model\\[\\[1\\]\\]")
  gap <- set_mortality(model, mortality_rates, weights = tibble::tibble(time_start = 1, time_end = Inf,
    model = list(half)))
  expect_error(run_simulation(gap, 2, seed = 1), "coverage gap")
  expect_s3_class(set_mortality(model, mortality_rates, weights = data.frame(time_start = 0, time_end = Inf,
    model = I(list(half))))$processes$mortality$weights, "tbl_df")
})

test_that("weight switches refresh stored weights when rates also change in place", {
  # Regression: a weight-model switch at the same instant as an in-place rate
  # update must still recompute every stored weight. An identical-rate calendar
  # split then reproduces the static-rate run exactly.
  model <- weighted_mortality_model(2000)
  w <- tibble::tibble(time_start = c(0, 1), time_end = c(1, Inf), model = list(
    weight_model(~ partnered, c(partnered = log(.1))),
    weight_model(~ partnered, c(partnered = 0))))
  split <- data.frame(sex = rep(c("m", "f"), 2), age_start = 0, age_end = Inf, rate = .05,
    time_start = rep(c(0, 1), each = 2), time_end = rep(c(1, Inf), each = 2))
  expect_identical(run_simulation(set_mortality(model, split, weights = w), 2,
    seed = 9)$persons,
                   run_simulation(set_mortality(model, mortality_rates,
                     weights = w), 2, seed = 9)$persons)
})

test_that("rejection sampling recovers after an extreme-weight member leaves", {
  # Regression: after the only high-weight member dies, the cell's weight bound
  # is stale; repeated rejections reset it to the exact maximum.
  persons <- data.frame(person_id = 1:4, sex = "m", birth_time = -30, death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_, group = c("B", "A", "A", "A"))
  model <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0,
                                     groups = list(group = group(categories = c("A", "B"))))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1),
      weights = weight_model(~ group, c(groupB = 40))) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand")
  withr::local_options(hammel.check_cells = TRUE)
  out <- run_simulation(model, 20, seed = 1)
  expect_equal(count_run_events(out, "deaths"), 4)
  # The extreme-weight member is chosen first.
  expect_equal(which.min(out$persons$death_time), 1L)
})

test_that("remaining validation rules for weight models", {
  expect_error(weight_model(fun = "not a function"), "must be a function")
  expect_error(weight_model(y ~ parity, c(parity = 1)), "one-sided")
  expect_error(weight_model(~ parity + offset(parity), c(parity = 1)), "offsets")
  model <- weighted_mortality_model(10)
  grouped_formula <- weight_model(~ group, c(groupB = 1))
  expect_error(run_simulation(set_mortality(model, mortality_rates, weights =
    grouped_formula), 1, seed = 1),
               "not available for this process: group")
  # Group levels are the dimension's categories; coefficients must name them.
  grouped <- weighted_mortality_model(200, groups = c("A", "B"))
  extra <- weight_model(~ group, c(groupB = log(2), groupC = log(5)))
  expect_error(run_simulation(set_mortality(grouped, mortality_rates, weights =
    extra), 1, seed = 1),
               "unknown groupC")
  declared <- set_groups(grouped, group = group(categories = c("A", "B", "C")))
  expect_s3_class(run_simulation(set_mortality(declared, mortality_rates,
    weights = extra), 1, seed = 1),
                  "hammel_population")
  # Ancestry variables need a tracked dimension and a known category.
  missing <- weight_model(~ ancestry_group_C, c(ancestry_group_C = 1))
  expect_error(run_simulation(set_mortality(grouped, mortality_rates, weights =
    missing), 1, seed = 1),
               "not available for this process: ancestry_group_C")
})

# Couples in two groups: half the unions are exogamous. No mortality or fertility.
grouped_union_model <- function(n, type = "marriage") {
  persons <- data.frame(person_id = 1:(2 * n), sex = rep(c("f", "m"), each = n),
    birth_time = c(rep(-30, n), rep(c(-30, -35), length.out = n)), death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_,
    group = c(rep(c("A", "B"), length.out = n), rep(c("A", "A", "B", "B"), length.out = n)))
  unions <- data.frame(union_id = 1:n, wife_id = 1:n, husband_id = n + 1:n, start_time = -1,
    end_time = NA_real_, end_reason = NA_character_, union_type = type)
  hammel_model() |>
    set_population(hammel_population(persons, unions = unions, end_time = 0,
      groups = if ("group" %in% names(persons)) list(group = group(ancestry = TRUE)))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(union_type = c("partnered", "unpartnered"), age_start = 15,
      age_end = 50, rate = 0), birth_interval = 0) |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3),
      kinship_generations = 0)
}

test_that("union weight formulas use both partners' variables", {
  n <- 8000
  model <- grouped_union_model(n)
  # Exogamous unions are three times as likely to dissolve within the cell.
  exo <- weight_model(~ same_group, c(same_group = -log(3)))
  out <- suppressWarnings(run_simulation(set_union_dissolution(model,
    data.frame(dur_start = 0, dur_end = Inf, rate = .2), weights = exo), 1, seed = 61))
  p <- out$persons
  first <- out$unions[seq_len(n), ]
  same <- p$group[match(first$wife_id, p$person_id)] == p$group[match(first$husband_id, p$person_id)]
  ended <- first$end_reason %in% "dissolution"
  years <- ifelse(is.na(first$end_time), 1, first$end_time)
  ratio <- (sum(ended[!same]) / sum(years[!same])) / (sum(ended[same]) / sum(years[same]))
  expect_lt(abs(log(ratio) - log(3)), 4 * sqrt(1 / sum(ended[!same]) + 1 / sum(ended[same])))
  # Every union variable is accepted and evaluated, for dissolution and conversion.
  all_terms <- weight_model(
    ~ group_w + group_h + parity_w + previous_union_w + previous_union_h + same_group +
      age_difference + union_type + ancestry_group_A_w + ancestry_group_A_h + ancestry_group_unknown_w,
    c(group_wB = .1, group_hB = -.1, parity_w = .05, previous_union_wseparated = .1,
      previous_union_wwidowed = .1, previous_union_hseparated = -.1, previous_union_hwidowed = -.1,
      same_group = -.2, age_difference = .02, union_typemarriage = .1, ancestry_group_A_w = .3,
      ancestry_group_A_h = -.3, ancestry_group_unknown_w = .1))
  small <- grouped_union_model(400, type = "cohab")
  both <- set_union_dissolution(small, data.frame(dur_start = 0, dur_end = Inf, rate = .2),
    weights = all_terms)
  # union_type is constant for conversion, so its term cancels with a warning.
  expect_warning(both <- set_union_formation(both, "rate",
    rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
    union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3),
    kinship_generations = 0, weights_conversion = all_terms), "union_type")
  withr::local_options(hammel.check_cells = TRUE)
  out <- suppressWarnings(run_simulation(both, 2, seed = 62))
  expect_gt(sum(out$unions$end_reason %in% "dissolution"), 0)
  expect_gt(sum(out$unions$end_reason %in% "conversion"), 0)
})

test_that("person formulas use parity, previous unions and ancestry", {
  model <- weighted_mortality_model(300, groups = c("A", "B"))
  w <- weight_model(~ parity + previous_union + ancestry_group_A + ancestry_group_unknown + partner_group,
    c(parity = .1, previous_unionseparated = .2, previous_unionwidowed = .3, ancestry_group_A = .5,
      ancestry_group_unknown = .1, partner_groupB = .2))
  withr::local_options(hammel.check_cells = TRUE)
  out <- suppressWarnings(run_simulation(set_mortality(model, mortality_rates,
    weights = w), 2, seed = 63))
  expect_gt(count_run_events(out, "deaths"), 0)
  # Group B has no group-A ancestry: a strong ancestry effect shifts deaths to A.
  strong <- weight_model(~ ancestry_group_A, c(ancestry_group_A = log(4)))
  big <- weighted_mortality_model(5000, groups = c("A", "B"))
  out <- run_simulation(set_mortality(big, mortality_rates, weights = strong),
    1, seed = 64)
  dead <- out$persons$group[!is.na(out$persons$death_time)]
  expect_gt(sum(dead == "A") / sum(dead == "B"), 2.5)
})

test_that("custom functions receive members and cells for every process", {
  seen <- new.env()
  record <- function(name) {
    function(members, bin) {
      if (is.null(seen[[name]])) {
        seen[[name]] <- list(members = members, bin = bin)
      }
      rep(0, nrow(members))
    }
  }
  # Dissolution and conversion: union members and their cells.
  model <- grouped_union_model(200, type = "cohab")
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .3),
    weights = weight_model(fun = record("dissolution")))
  model <- set_union_formation(model, "rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
    union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3),
    kinship_generations = 0, weights_conversion = weight_model(fun = record("conversion")))
  suppressWarnings(run_simulation(model, 2, seed = 65))
  for (name in c("dissolution", "conversion")) {
    expect_true(all(c("union_id", "wife_id", "husband_id", "group_w", "group_h", "parity_w",
      "previous_union_w", "previous_union_h", "same_group", "age_difference", "union_type",
      "ancestry_group_A_w", "ancestry_group_unknown_h") %in% names(seen[[name]]$members)))
    expect_true(all(c("dur_start", "dur_end", "time_start", "time_end") %in% names(seen[[name]]$bin)) ||
                all(c("age_start", "age_end") %in% names(seen[[name]]$bin)))
  }
  expect_true("union_type" %in% names(seen$dissolution$bin))
  # Search without calibration (one group), with calendar rates and previous-union schedules.
  rates <- data.frame(previous_union = rep(c("none", "separated", "widowed"), 2), age_start = 15,
    age_end = Inf, rate = .5, time_start = rep(c(0, 1), each = 3), time_end = rep(c(1, Inf), each = 3))
  search <- make_rate_model(200) |>
    set_union_formation("rate", rates = rates, kinship_generations = 0,
      weights = weight_model(fun = record("search")))
  search$population$persons$sex <- rep(c("f", "m"), 100)
  suppressWarnings(run_simulation(search, 2, seed = 66))
  expect_identical(seen$search$bin$process, "search")
  expect_true(seen$search$bin$previous_union %in% c("none", "separated", "widowed"))
  expect_equal(c(seen$search$bin$time_start, seen$search$bin$time_end), c(0, 1))
  # Transitions: the destination appears as a member column and in the cell.
  persons <- data.frame(person_id = 1:200, sex = rep(c("f", "m"), 100), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_, group = "A")
  moves <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand") |>
    set_group_transition(data.frame(origin = "A", destination = "B", sex = c("f", "m"),
      age_start = 0, age_end = Inf, rate = .5), weights = weight_model(fun = record("transition")))
  run_simulation(moves, 1, seed = 67)
  expect_true(all(seen$transition$members$destination == "B"))
  expect_identical(c(seen$transition$bin$origin, seen$transition$bin$destination), c("A", "B"))
  expect_true(seen$transition$bin$sex %in% c("f", "m"))
  # Fertility without calibration: custom weights at selection.
  fert <- set_fertility(weighted_fertility_model(300), data.frame(age_start = 15, age_end = 50, rate = .5),
    birth_interval = 0,
    weights = weight_model(fun = function(members, bin) {
      if (is.null(seen$fertility)) seen$fertility <- list(members = members, bin = bin)
      log(seq_len(nrow(members)))
    }))
  expect_gt(count_run_events(run_simulation(fert, 1, seed = 68), "births"), 0)
  expect_identical(seen$fertility$bin$process, "fertility")
  expect_true(all(c("parity_start", "parity_end", "age_start", "age_end") %in% names(seen$fertility$bin)))
})
