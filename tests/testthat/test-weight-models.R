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
    set_population(hammel_population(persons, unions = unions, end_time = 0)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .05)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0)
}
mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .05)

test_that("weight models validate formulas, coefficients and variables", {
  expect_s3_class(create_weight_model(~ partnered, c(partnered = 1)), "hammel_weight_model")
  expect_s3_class(create_weight_model(fun = function(members, bin) 0), "hammel_weight_model")
  expect_error(create_weight_model(), "exactly one")
  expect_error(create_weight_model(~ parity, c(parity = 1), fun = function(m, b) 0), "exactly one")
  expect_error(create_weight_model(~ I(parity > 1), c(x = 1)), "custom function")
  expect_error(create_weight_model(~ log(parity), c(x = 1)), "custom function")
  expect_error(create_weight_model(~ parity - 1, c(parity = 1)), "intercept")
  expect_error(create_weight_model(~ education, c(education = 1)), "Unknown weight variables")
  expect_error(create_weight_model(~ parity, c(parity = Inf)), "finite")
  expect_error(create_weight_model(~ group, c(groupB = 1)), "group_levels")
  model <- weighted_mortality_model(2)
  expect_error(set_mortality(model, mortality_rates,
    weights = create_weight_model(~ parity, c(paritx = 1))), "missing parity; unknown paritx")
  expect_error(set_mortality(model, mortality_rates,
    weights = create_weight_model(~ same_group, c(same_group = 1))), "not available")
  expect_error(set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .1),
    weights = create_weight_model(~ partnered, c(partnered = 1))), "not available")
  expect_error(set_mortality(model, mortality_rates, weights = list()), "create_weight_model")
  # Interactions follow R's model-matrix coding and names.
  w <- create_weight_model(~ union_type * parity,
    c(union_typecohab = 1, union_typemarriage = 2, parity = .1,
      "union_typecohab:parity" = .2, "union_typemarriage:parity" = .3))
  expect_silent(set_mortality(model, mortality_rates, weights = w))
})

test_that("terms the input rates already stratify by warn that they cancel", {
  model <- weighted_mortality_model(2)
  expect_warning(set_mortality(model, mortality_rates,
    weights = create_weight_model(~ sex + partnered, c(sexm = 1, partnered = 1))), "`sex` has no effect")
  status <- data.frame(union_type = c("partnered", "unpartnered"), age_start = 15, age_end = 50, rate = .1)
  expect_warning(set_fertility(model, status,
    weights = create_weight_model(~ partnered, c(partnered = 1))), "`partnered` has no effect")
  expect_warning(set_union_dissolution(model, data.frame(union_type = c("marriage", "cohab"),
    dur_start = 0, dur_end = Inf, rate = .1),
    weights = create_weight_model(~ union_type, c(union_typemarriage = 1))), "`union_type` has no effect")
  expect_silent(set_mortality(model, mortality_rates,
    weights = create_weight_model(~ sex:partnered, c("sexf:partnered" = 1, "sexm:partnered" = 1))))
})

test_that("mortality weights give within-cell hazard ratios and keep the cell rate", {
  n <- 20000
  model <- set_mortality(weighted_mortality_model(n), mortality_rates,
    weights = create_weight_model(~ partnered, c(partnered = log(0.5))))
  out <- simulate(model, 1, seed = 11)
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
  plain <- simulate(model, 2, seed = 5)
  zero <- set_mortality(model, mortality_rates,
    weights = create_weight_model(~ partnered, c(partnered = 0)))
  expect_identical(simulate(zero, 2, seed = 5), plain)
  constant <- set_mortality(model, mortality_rates,
    weights = create_weight_model(fun = function(members, bin) rep(3, nrow(members))))
  expect_identical(simulate(constant, 2, seed = 5), plain)
})

test_that("custom functions receive members and the cell, without exact times", {
  model <- weighted_mortality_model(200)
  seen <- NULL
  f <- function(members, bin) {
    if (is.null(seen)) seen <<- list(members = members, bin = bin)
    log(ifelse(members$partnered == 1, .5, 1))
  }
  out <- simulate(set_mortality(model, mortality_rates, weights = create_weight_model(fun = f)), 1, seed = 3)
  expect_s3_class(seen$members, "tbl_df")
  expect_true(all(c("person_id", "sex", "group", "union_type", "partnered", "parity",
    "previous_union", "partner_group") %in% names(seen$members)))
  expect_false(any(c("birth_time", "age", "time") %in% names(seen$members)))
  expect_identical(seen$bin$process, "mortality")
  expect_true(all(c("sex", "age_start", "age_end", "time_start", "time_end") %in% names(seen$bin)))
  # The same model written as a formula gives the same run.
  formula <- set_mortality(model, mortality_rates,
    weights = create_weight_model(~ partnered, c(partnered = log(.5))))
  expect_identical(simulate(formula, 1, seed = 3), out)
  bad_length <- create_weight_model(fun = function(members, bin) 0)
  expect_error(simulate(set_mortality(model, mortality_rates, weights = bad_length), 1, seed = 3),
    "returned 1 values")
  infinite <- create_weight_model(fun = function(members, bin) rep(Inf, nrow(members)))
  expect_error(simulate(set_mortality(model, mortality_rates, weights = infinite), 1, seed = 3), "infinite")
  # Draws inside the function share the simulation stream reproducibly.
  noisy <- set_mortality(model, mortality_rates,
    weights = create_weight_model(fun = function(members, bin) rnorm(nrow(members))))
  expect_identical(simulate(noisy, 1, seed = 8), simulate(noisy, 1, seed = 8))
})

test_that("missing variables contribute zero with one warning", {
  model <- weighted_mortality_model(200, groups = c("A", "B"))
  w <- create_weight_model(~ partner_group + partnered,
    c(partner_groupB = 1, partnered = 1), group_levels = c("A", "B"))
  expect_warning(simulate(set_mortality(model, mortality_rates, weights = w), 1, seed = 2),
    "Missing values in the mortality weight formula")
  na_fun <- create_weight_model(fun = function(members, bin) rep(NA_real_, nrow(members)))
  model <- weighted_mortality_model(200)
  expect_warning(out <- simulate(set_mortality(model, mortality_rates, weights = na_fun), 1, seed = 2),
    "function result")
  expect_identical(out, simulate(model, 1, seed = 2))
  # Group levels must cover the population.
  partial <- create_weight_model(~ group, c(groupC = 1), group_levels = c("A", "C"))
  grouped <- weighted_mortality_model(20, groups = c("A", "B"))
  expect_error(simulate(set_mortality(grouped, mortality_rates, weights = partial), 1),
    "include every population group")
})

# Women aged 25 (one fertility cell for a year), with birth interval zero.
weighted_fertility_model <- function(n, groups = NULL) {
  persons <- data.frame(person_id = 1:(2 * n), sex = rep(c("f", "m"), each = n),
    birth_time = -25.2, death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  if (!is.null(groups)) persons$group <- rep(groups, length.out = 2 * n)
  hammel_model() |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0, cancel_birth_if_no_match = FALSE)
}

test_that("fertility weights select mothers at selection and at refresh under calibration", {
  rates <- data.frame(age_start = 15, age_end = 50, rate = .5)
  check <- function(model, n) {
    out <- simulate(model, .7, seed = 21)
    p <- out$persons
    births <- table(factor(p$group[match(p$mother_id[p$birth_time >= 0], p$person_id)],
                           levels = c("A", "B")))
    expected <- .5 * .7 * n
    expect_lt(abs(sum(births) - expected), 4 * sqrt(expected))
    expect_lt(abs(log(births[["B"]] / births[["A"]]) - log(2)),
              4 * sqrt(1 / births[["A"]] + 1 / births[["B"]]))
  }
  w <- create_weight_model(~ group, c(groupB = log(2)), group_levels = c("A", "B"))
  n <- 4000
  # Two groups activate weighted matching calibration: weights stored at refresh.
  check(set_fertility(weighted_fertility_model(n, c("A", "B")), rates, birth_interval = 0, weights = w), n)
  # Calibration off (matching restricted to the woman's own group): selection time.
  restricted <- weighted_fertility_model(n, c("A", "B"))
  restricted <- set_union_formation(restricted, "on_demand", kinship_generations = 0,
    cancel_birth_if_no_match = FALSE,
    matching = create_matching_model(~ age_difference, c(age_difference = 0), prevent_exogamy = TRUE))
  check(set_fertility(restricted, rates, birth_interval = 0, weights = w), n)
  # A fertility weight model combines with hetfert, and zero coefficients keep it.
  model <- set_fertility(weighted_fertility_model(300, c("A", "B")), rates, hetfert = het_fert(.5, 1))
  zero <- set_fertility(model, rates, hetfert = het_fert(.5, 1),
    weights = create_weight_model(~ group, c(groupB = 0), group_levels = c("A", "B")))
  expect_identical(simulate(zero, 1, seed = 4), simulate(model, 1, seed = 4))
  expect_warning(simulate(set_fertility(model, rates, weights = create_weight_model(
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
    set_population(hammel_population(persons, unions = unions, end_time = 0)) |>
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
  w <- create_weight_model(~ age_difference, c(age_difference = log(1.2)))
  model <- weighted_union_model(n) |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      kinship_generations = 0) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .3), weights = w)
  out <- suppressWarnings(simulate(model, 1, seed = 31))
  r <- gap_ratio(out$unions, n, "dissolution")
  expect_lt(abs(r[["log_ratio"]] - 5 * log(1.2)), 4 * r[["se"]])
  cohabs <- weighted_union_model(n, "cohab") |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 1),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3),
      kinship_generations = 0, weights_conversion = w)
  out <- suppressWarnings(simulate(cohabs, 1, seed = 32))
  r <- gap_ratio(out$unions, n, "conversion")
  expect_lt(abs(r[["log_ratio"]] - 5 * log(1.2)), 4 * r[["se"]])
})

test_that("transition weights can use person variables and the destination", {
  n <- 4000
  persons <- data.frame(person_id = 1:n, sex = rep(c("f", "m"), length.out = n), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_, group = "A")
  model <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand", kinship_generations = 0)
  rates <- data.frame(group_origin = "A", group_dest = c("B", "C"), age_start = 0, age_end = Inf, rate = .2)
  w <- create_weight_model(~ sex:group_dest,
    c("sexf:group_destA" = 0, "sexm:group_destA" = 0, "sexf:group_destB" = log(2),
      "sexm:group_destB" = 0, "sexf:group_destC" = 0, "sexm:group_destC" = 0),
    group_levels = c("A", "B", "C"))
  out <- simulate(set_group_transition(model, rates, weights = w), 1, seed = 41)
  tx <- out$transitions
  sex <- out$persons$sex[match(tx$person_id, out$persons$person_id)]
  b <- table(factor(sex[tx$group_dest == "B"], levels = c("f", "m")))
  c_moves <- table(factor(sex[tx$group_dest == "C"], levels = c("f", "m")))
  # Women are twice as likely as men to be chosen for a move to B, not to C.
  expect_lt(abs(log(b[["f"]] / b[["m"]]) - log(2)), 4 * sqrt(1 / b[["f"]] + 1 / b[["m"]]))
  expect_lt(abs(log(c_moves[["f"]] / c_moves[["m"]])), 4 * sqrt(1 / c_moves[["f"]] + 1 / c_moves[["m"]]))
  expect_warning(group_transition(rates, weights = create_weight_model(~ group_dest,
    c(group_destB = 1, group_destC = 1), group_levels = c("A", "B", "C"))), "no effect")
})
