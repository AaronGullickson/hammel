# weighted_choice (D066): female choice with male-group weights calibrated at
# each search so expected male group flows meet the mirrored target.

prepare_weighted <- function(model, groups) {
  prepare_matching_model(as_internal_matching(model, groups), groups)
}

calibrate <- function(women, men, matching, rates, groups, now = 0, cutoff = 15,
                      width = 0.25, fertility = FALSE) {
  inspect_weighted_calibration_cpp(women$born, women$group, women$fertile,
    men$born, men$group, rates, prepare_weighted(matching, groups), length(groups),
    now, cutoff, width, fertility)
}

# Independent individual-level reference (two groups): exact ages, no bins.
reference_weights <- function(women, men, matching, rates, groups, now = 0,
                              cutoff = 15, fertility = FALSE) {
  women <- women[now - women$born >= cutoff, ]
  men <- men[now - men$born >= cutoff, ]
  hazard <- function(age) {
    i <- which(age >= rates$age_start & age < rates$age_end)
    if (length(i)) rates$rate[i[1]] else 0
  }
  lambda <- vapply(now - women$born, hazard, numeric(1)) * if (fertility) women$fertile else 1
  intensity <- tapply(lambda, factor(women$group, 0:1), sum)
  women_n <- tabulate(women$group + 1L, 2)
  men_n <- tabulate(men$group + 1L, 2)
  mirror <- men_n * intensity / women_n
  mirror <- mirror / sum(mirror)
  female_share <- intensity / sum(intensity)
  d <- vapply(groups, function(g) if (g %in% names(matching$delta)) matching$delta[[g]] else 0, numeric(1))
  target <- mirror * exp(d - max(d))
  target <- target / sum(target)
  U <- vapply(seq_len(nrow(men)), function(j) calculate_matching_utilities(matching,
    women$born - men$born[j], groups[women$group + 1L],
    rep(groups[men$group[j] + 1L], nrow(women))), numeric(nrow(women)))
  W <- exp(U - max(U))
  S <- cbind(rowSums(W[, men$group == 0, drop = FALSE]), rowSums(W[, men$group == 1, drop = FALSE]))
  flow_b <- function(log_b) {
    P <- S * rep(c(1, exp(log_b)), each = nrow(S))
    sum(lambda * P[, 2] / rowSums(P)) / sum(lambda)
  }
  root <- uniroot(function(x) flow_b(x) - target[2], c(-30, 30), tol = 1e-13)$root
  list(log_weights = c(0, root), target = as.numeric(target))
}

group_only <- function(women_share, men_share, effect, delta = NULL) {
  women <- data.frame(born = -30, group = rep(0:1, round(1000 * women_share)), fertile = TRUE)
  men <- data.frame(born = -30, group = rep(0:1, round(1000 * men_share)))
  matching <- matching_model(~ 1, numeric(), group_exogamy = c("A:B" = effect),
                                    delta = delta)
  list(women = women, men = men, matching = matching,
       rates = data.frame(age_start = 15, age_end = 50, rate = 0.1))
}

test_that("calibration reproduces the analytical group-only examples", {
  # Women and men 75/25, A:B = -3: symmetric table, b_B = 1.6830 (spec example).
  x <- group_only(c(0.75, 0.25), c(0.75, 0.25), -3)
  result <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_false(result$infeasible)
  expect_equal(result$log_weights[1], 0)
  expect_equal(exp(result$log_weights[2]), 1.682979, tolerance = 1e-6)
  expect_equal(result$flow, c(0.75, 0.25), tolerance = 1e-9)
  x <- group_only(c(.75, .25), c(.70, .30), -3)
  result <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_equal(exp(result$log_weights[2]), 3.875051, tolerance = 1e-6)
  expect_equal(result$target, c(.7, .3), tolerance = 1e-12)
  # Delta can offset a 70/30 male versus 75/25 female composition imbalance.
  x <- group_only(c(.75, .25), c(.70, .30), -3, c(B = log(7/9)))
  balanced <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_equal(balanced$target, c(.75, .25), tolerance = 1e-12)
  expect_equal(exp(balanced$log_weights[2]), 1.308984, tolerance = 1e-6)
  # Independently derived target for AB/BA = 3 at log cross-product -6.
  x <- group_only(c(.75, .25), c(.75, .25), -3, c(B = .12274682836))
  result <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_equal(result$target[2], .27371293659, tolerance = 1e-10)
  w <- exp(result$log_weights[2]); t <- exp(-3)
  AB <- .75 * .25 * w * t / (.75 + .25 * w * t)
  BA <- .25 * .75 * t / (.75 * t + .25 * w)
  expect_equal(AB / BA, 3, tolerance = 1e-7)
  expect_equal(log(AB * BA / ((.75 - AB) * (.25 - BA))), -6)
  x$matching$delta <- c(A = 2, B = 2 + .12274682836)
  shifted <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_equal(shifted$target, result$target, tolerance = 1e-12)
  expect_equal(shifted$log_weights, result$log_weights, tolerance = 1e-9)
  x$matching$delta <- c(B = -.12274682836)
  negative <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_lt(negative$target[2], .25)
  x$matching$delta <- c(B = 1000)
  extreme <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_true(all(is.finite(extreme$target)))
  expect_equal(sum(extreme$target), 1)

})

test_that("binned calibration agrees with an exact individual-level reference", {
  set.seed(901)
  women <- data.frame(born = -runif(400, 10, 80), group = sample(0:1, 400, TRUE, c(0.7, 0.3)),
                      fertile = runif(400) < 0.8)
  men <- data.frame(born = -runif(350, 10, 85), group = sample(0:1, 350, TRUE, c(0.6, 0.4)))
  rates <- data.frame(age_start = c(15, 30), age_end = c(30, 50), rate = c(0.2, 0.05))
  matching <- age_gap_matching_model(2, 5, group_exogamy = c("A:B" = -1.5))
  for (fertility in c(FALSE, TRUE)) {
    expected <- reference_weights(women, men, matching, rates, c("A", "B"), fertility = fertility)
    fine <- calibrate(women, men, matching, rates, c("A", "B"), width = 0.001, fertility = fertility)
    expect_false(fine$infeasible)
    expect_equal(fine$log_weights, expected$log_weights, tolerance = 2e-3)
    expect_equal(fine$target, expected$target, tolerance = 2e-3)
    default <- calibrate(women, men, matching, rates, c("A", "B"), fertility = fertility)
    expect_equal(default$log_weights, expected$log_weights, tolerance = 0.05)
    expect_equal(default$flow, default$target, tolerance = 1e-9)
  }
  # Custom R formulas use the same grid through the formula evaluator.
  custom <- matching_model(~ abs(age_difference - 2),
    c("abs(age_difference - 2)" = -0.3), group_exogamy = c("A:B" = -1))
  expected <- reference_weights(women, men, custom, rates, c("A", "B"))
  fine <- calibrate(women, men, custom, rates, c("A", "B"), width = 0.001)
  expect_equal(fine$log_weights, expected$log_weights, tolerance = 2e-3)
})

test_that("incremental pool sums match a full rebuild after many pool changes", {
  set.seed(909)
  women <- data.frame(born = -runif(300, 10, 60), group = sample(0:1, 300, TRUE, c(0.7, 0.3)))
  men <- data.frame(born = -runif(400, 5, 80), group = sample(0:1, 400, TRUE, c(0.6, 0.4)))
  # Men leave and rejoin the pool while time advances 30 years, so younger
  # cohorts also cross the age cutoff part-way through.
  n_ops <- 4000
  times <- sort(runif(n_ops, 0, 30))
  present <- rep(TRUE, nrow(men))
  who <- sample.int(nrow(men), n_ops, TRUE)
  deltas <- numeric(n_ops)
  for (i in seq_len(n_ops)) {
    deltas[i] <- if (present[who[i]]) -1 else 1
    present[who[i]] <- !present[who[i]]
  }
  matching <- age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -2))
  result <- inspect_weighted_incremental_cpp(women$born, women$group, men$born, men$group,
    times, who - 1L, deltas, data.frame(age_start = 15, age_end = 50, rate = 0.1),
    prepare_weighted(matching, c("A", "B")), 2L, 15)
  expect_lt(result$max_sum_difference, 1e-9)
  expect_equal(result$incremental, result$rebuilt, tolerance = 1e-8)
})

test_that("three-group calibration meets targets and handles zero and unattainable targets", {
  set.seed(902)
  women <- data.frame(born = -runif(600, 15, 45), group = sample(0:2, 600, TRUE, c(5, 3, 1)),
                      fertile = TRUE)
  men <- data.frame(born = -runif(600, 15, 50), group = sample(0:2, 600, TRUE, c(4, 4, 1)))
  rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
  matching <- age_gap_matching_model(0, 5,
    group_exogamy = c("A:B" = -3, "B:C" = -.75), delta = c(B = .3, C = -.2))
  result <- calibrate(women, men, matching, rates, c("A", "B", "C"))
  expect_false(result$infeasible)
  expect_equal(result$flow, result$target, tolerance = 1e-9)
  expect_equal(sum(result$target), 1)
  expected <- tabulate(men$group + 1L, 3) * exp(c(0, .3, -.2))
  expect_equal(result$target, expected / sum(expected), tolerance = 1e-12)
  # Group-stratified search hazards also enter mirrored baseline targets.
  stratified <- data.frame(.hammel_group = 0:2, age_start = 15, age_end = 50,
                            rate = c(.1, .2, .05))
  result_rates <- calibrate(women, men, matching, stratified, c("A", "B", "C"))
  expected <- expected * stratified$rate
  expect_equal(result_rates$target, expected / sum(expected), tolerance = 1e-12)
  # Men of a group with no eligible women still receive the population-wide
  # per-capita rate; with all exogamy prevented they cannot be reached.
  women$group <- 0L
  blocked <- matching_model(~ 1, numeric(), group_exogamy = NA)
  result <- calibrate(women, men, blocked, rates, c("A", "B", "C"))
  expect_true(result$infeasible)
  expect_equal(result$log_weights[2:3], c(50, 50))
  # Zero baseline hazards remain zero even with positive target adjustments.
  grouped_rates <- data.frame(.hammel_group = 0:2, age_start = 15,
                              age_end = 50, rate = c(.1, 0, 0))
  women_zero <- women
  women_zero$group <- rep(0:2, length.out = nrow(women_zero))
  zero <- matching_model(~ 1, numeric(), group_exogamy = NA,
                                 delta = c(B = 10, C = 10))
  result_zero <- calibrate(women_zero, men, zero, grouped_rates, c("A", "B", "C"))
  expect_equal(result_zero$target, c(1, 0, 0))
  # No women at risk leaves unit weights.
  idle <- calibrate(women, men, matching, data.frame(age_start = 15, age_end = 50, rate = 0),
                    c("A", "B", "C"))
  expect_equal(idle$log_weights, c(0, 0, 0))
})

weighted_population <- function(n = 600, weights = NULL, seed = 903) {
  set.seed(seed)
  create_baseline_population(n, max_age = 60, end_time = 0,
    groups = if (!is.null(weights)) list(group = weights))
}

weighted_model <- function(population, initiation = "on_demand",
                           matching = age_gap_matching_model(0, 5),
                           father_matching = NULL, ...) {
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01))
  identity(if (initiation == "on_demand") {
    model |>
      set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.12)) |>
      set_union_formation("on_demand", matching, father_matching = father_matching, ...)
  } else {
    model |>
      set_fertility(data.frame(age_start = 15, age_end = 45, rate = c(0.15, 0.05),
                               union_type = c("partnered", "unpartnered"))) |>
      set_union_formation("rate", matching, rates = data.frame(age_start = 15, age_end = 50, rate = 0.15),
                          father_matching = father_matching, ...)
  })
}

test_that("one group or no groups use unit weights without changing RNG isolation", {
  for (weights in list(NULL, c(A = 1))) {
    population <- weighted_population(400, weights)
    for (initiation in c("on_demand", "rate")) {
      model <- weighted_model(population, initiation)
      ordinary <- run_simulation(model, 15, seed = 904)
      # A common target adjustment cannot change choices in a one-group model.
      if (!is.null(weights)) model$processes$union_formation$matching$delta <- c(A = 2)
      set.seed(1); before <- .Random.seed
      result <- run_simulation(model, 15, seed = 904)
      expect_identical(.Random.seed, before)
      expect_identical(result, ordinary)
    }
  }
})

test_that("grouped weighted choice runs in both modes, restarts, and is reproducible", {
  population <- weighted_population(600, c(A = 3, B = 1))
  matching <- age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -2), delta = c(B = .2))
  for (initiation in c("on_demand", "rate")) {
    model <- weighted_model(population, initiation, matching = matching)
    first <- run_simulation(model, 20, seed = 905)
    expect_identical(first, run_simulation(model, 20, seed = 905))
    expect_gt(nrow(first$unions), nrow(population$unions))
    restart <- run_simulation(set_population(model, first), 5, seed = 906)
    expect_equal(restart$end_time, first$end_time + 5)
    expect_identical(restart$persons[seq_len(nrow(first$persons)), c("person_id", "birth_time")],
                     first$persons[c("person_id", "birth_time")])
    # Diagnostics consume no random draws and do not change histories.
    diagnosed <- run_simulation(weighted_model(population, initiation,
                                         matching = matching, diagnostics = TRUE), 20, seed = 905)
    expect_identical(diagnosed$persons, first$persons)
    expect_identical(diagnosed$unions[names(first$unions)], first$unions)
    expect_true(all(c(".log_weight_A", ".log_weight_B") %in% names(diagnosed$unions)))
    expect_false(any(c(".n_eligible", ".utility_rank") %in% names(diagnosed$unions)))
  }
  # Rate mode: unpartnered births select fathers through their own calibration.
  father <- age_gap_matching_model(3, 4, group_exogamy = c("A:B" = -1), delta = c(B = -.3))
  result <- run_simulation(weighted_model(population, "rate", matching =
    matching,
                                    father_matching = father), 20, seed = 907)
  persons <- result$persons
  born <- persons[persons$birth_time >= 0, ]
  in_union <- mapply(function(mother, time) any(result$unions$wife_id == mother &
    result$unions$start_time <= time & (is.na(result$unions$end_time) | result$unions$end_time > time)),
    born$mother_id, born$birth_time)
  expect_true(any(!in_union))
  expect_true(all(!is.na(born$father_id)))
})

test_that("father models retain separate delta and inherit the general model when absent", {
  population <- weighted_population(600, c(A = 3, B = 1))
  matching <- age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -2), delta = c(B = .3))
  inherited <- weighted_model(population, "rate", matching = matching)
  explicit <- weighted_model(population, "rate", matching = matching, father_matching = matching)
  expect_identical(run_simulation(inherited, 20, seed = 910), run_simulation(
    explicit, 20, seed = 910))
  father <- matching
  father$delta <- c(B = -.6)
  different <- weighted_model(population, "rate", matching = matching, father_matching = father)
  result <- run_simulation(different, 20, seed = 910)
  expect_identical(result, run_simulation(different, 20, seed = 910))
  expect_false(identical(result$persons, run_simulation(inherited, 20, seed =
    910)$persons))
})

test_that("weighted choice warns once when weights must be approximated", {
  population <- weighted_population(300, c(A = 1))
  persons <- population$persons
  persons$group[persons$sex == "m"][1:40] <- "B"
  population <- hammel_population(persons, population$unions, end_time = 0, groups = if ("group" %in% names(persons)) "group")
  # A partial restriction keeps calibration active (full prevention skips it).
  matching <- matching_model(~ age_difference, c(age_difference = 0),
                                    group_exogamy = c("A:B" = NA))
  model <- weighted_model(population, matching = matching)
  warnings <- character()
  withCallingHandlers(run_simulation(model, 10, seed = 908),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_length(grep("could not meet male-group targets", warnings), 1)
})

test_that("removed method arguments are rejected and delta inputs validated", {
  expect_false("matching_method" %in% names(matching_model(~ 1, numeric())))
  expect_error(matching_model(~ 1, numeric(), matching_method = "female_choice"), "unused argument")
  expect_error(age_gap_matching_model(0, 5, matching_method = "weighted_choice"), "unused argument")
  expect_error(union_formation("on_demand", matching_method = "female_choice"), "unused argument")
  expect_error(set_union_formation(hammel_model(), "on_demand", matching_method = "female_choice"), "unused argument")
  for (bad in list(1, NA_real_, c(B = Inf), c(B = NA_real_), c(B = "1"),
                   c(B = 1, B = 2), setNames(1, ""), matrix(1)))
    expect_error(age_gap_matching_model(0, 5, delta = bad), "delta")
  expect_null(matching_model(~ 1, numeric())$delta)
  expect_error(matching_model(~ 1, numeric(), imbalance_absorption = 1), "unused argument")
  edited <- age_gap_matching_model(0, 5)
  edited$delta <- c(B = Inf)
  expect_error(union_formation("on_demand", edited), "delta")
  expect_error(prepare_matching_model(matching_model(~ 1, numeric(), delta = c(X = 1)),
                                      c("A", "B")), "delta")
  expect_error(prepare_matching_model(matching_model(~ 1, numeric(), delta = c(A = 1))), "delta")

})

test_that("male-group formula terms warn about absorbed utilities", {
  expect_warning(matching_model(~ male_group, c(male_groupB = 1)),
                 "male_group` are absorbed")
  expect_silent(age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -3)))
  expect_silent(matching_model(~ same_group, c(same_group = 1)))
})

test_that("strict endogamy has zero group log weights even with delta adjustments", {
  population <- weighted_population(600, c(A = 3, B = 1))
  matching <- age_gap_matching_model(0, 5, group_exogamy = NA)
  for (initiation in c("on_demand", "rate")) {
    baseline <- run_simulation(weighted_model(population, initiation, matching =
      matching,
                                        diagnostics = TRUE), 15, seed = 912)
    adjusted <- matching
    adjusted$delta <- c(B = 2)
    expect_no_warning(result <- run_simulation(weighted_model(population,
      initiation,
      matching = adjusted, diagnostics = TRUE), 15, seed = 912))
    expect_identical(result, baseline)
    expect_true(all(as.matrix(result$unions[c(".log_weight_A", ".log_weight_B")]) == 0))
    unions <- get_detailed_unions(result)
    expect_true(all(unions$group_w == unions$group_h))
  }
})

test_that("weighted choice handles identities created during the run", {
  # Men of a group with no counted women use the population-wide female rate.
  women <- data.frame(born = rep(-30, 400), group = rep(0:1, c(300, 100)), fertile = TRUE)
  men <- data.frame(born = rep(-30, 450), group = rep(0:2, c(300, 100, 50)))
  result <- calibrate(women, men, age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1)),
                      data.frame(age_start = 15, age_end = 50, rate = 0.1), c("A", "B", "AB"))
  expect_false(result$infeasible)
  expect_equal(result$target[3], 50 / 450, tolerance = 1e-12)
  expect_equal(result$flow, result$target, tolerance = 1e-9)
  # Mixed couples have AB children, who later enter the pools and form unions.
  rule <- group_inheritance("custom", fun = function(mother_group, father_group, ...) {
    if (is.na(father_group) || mother_group == father_group) mother_group else "AB"
  })
  population <- weighted_population(800, c(A = 3, B = 1))
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.15)) |> set_inheritance(rule, extra = "AB") |>
    set_union_formation("on_demand",
      age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1, "A:AB" = -0.5, "B:AB" = -0.5)),
      diagnostics = TRUE)
  expect_no_warning(result <- run_simulation(model, 40, seed = 913))
  expect_identical(result, run_simulation(model, 40, seed = 913))
  expect_true(any(result$persons$group == "AB"))
  unions <- get_detailed_unions(result)
  expect_true(any(unions$group_w == "AB" | unions$group_h == "AB"))
  expect_true(all(c(".log_weight_A", ".log_weight_B", ".log_weight_AB") %in% names(result$unions)))
  restart <- run_simulation(set_population(model, result), 5, seed = 914)
  expect_equal(restart$end_time, result$end_time + 5)
})

test_that("diagnostics record the group log weights applied at each union", {
  population <- weighted_population(600, c(A = 3, B = 1))
  matching <- age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -2))
  model <- weighted_model(population, matching = matching, diagnostics = TRUE)
  first <- run_simulation(model, 10, seed = 915)
  columns <- c(".log_weight_A", ".log_weight_B")
  expect_true(all(columns %in% names(first$unions)))
  weights <- as.matrix(first$unions[columns])
  expect_true(nrow(weights) > 0 && all(is.finite(weights)))
  # One group is the reference (log weight 0); only differences between groups matter.
  expect_true(any(weights[, 1] != weights[, 2]))
  expect_true(all(apply(weights, 1, function(w) any(w == 0))))
  # Diagnostics add columns without changing histories.
  plain <- run_simulation(weighted_model(population, matching = matching), 10,
    seed = 915)
  expect_identical(first$persons, plain$persons)
  expect_identical(first$unions[names(plain$unions)], plain$unions)
  # Restarts keep earlier rows and fill only new unions.
  second <- run_simulation(set_population(model, first), 5, seed = 916)
  expect_identical(second$unions[seq_len(nrow(first$unions)), columns], first$unions[columns])
  expect_true(all(is.finite(as.matrix(second$unions[-seq_len(nrow(first$unions)), columns]))))
  # Zero-duration runs retain diagnostic columns; no groups means no columns.
  expect_true(all(columns %in% names(run_simulation(model, 0)$unions)))
  ungrouped <- run_simulation(weighted_model(weighted_population(300),
                                       diagnostics = TRUE), 5, seed = 917)
  expect_false(any(grepl("^\\.log_weight_", names(ungrouped$unions))))
})

test_that("unmet targets hidden by near-zero weights are reported as approximate", {
  # B women are present but inactive, so B men's target is zero; C women can
  # only partner B men (A:C prohibited, no C men), so that target cannot hold.
  women <- data.frame(born = c(rep(-30, 60), rep(-60, 20), rep(-30, 20)),
                      group = rep(0:2, c(60, 20, 20)), fertile = TRUE)
  men <- data.frame(born = -30, group = rep(0:1, c(60, 20)))
  rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
  matching <- age_gap_matching_model(0, 5, group_exogamy = c("A:C" = NA))
  result <- calibrate(women, men, matching, rates, c("A", "B", "C"))
  expect_equal(result$target, c(1, 0, 0))
  expect_true(result$infeasible)
  expect_equal(result$flow, c(0.75, 0.25, 0), tolerance = 1e-9)
})

test_that("calibration uses actual hazards inside narrow age bands", {
  # A band much narrower than a bin: actual age, not bin-average hazard, applies.
  x <- group_only(c(0.75, 0.25), c(0.75, 0.25), -3)
  x$women$born <- -30.01
  narrow <- data.frame(age_start = c(30, 30.05), age_end = c(30.05, 60), rate = c(1, 0))
  result <- calibrate(x$women, x$men, x$matching, narrow, c("A", "B"))
  expect_false(result$infeasible)
  expect_equal(exp(result$log_weights[2]), 1.682979, tolerance = 1e-6)
  # Different exact ages in the same calibration bin can have different risks.
  x$women$born[x$women$group == 1] <- -30.1
  different <- calibrate(x$women, x$men, x$matching, narrow, c("A", "B"))
  expect_equal(different$target, c(1, 0))
  expect_equal(different$flow, c(1, 0), tolerance = 1e-9)
  # One-year bands put a band edge inside many bins; agreement with the exact
  # individual-level reference stays close at the default width.
  set.seed(918)
  women <- data.frame(born = -runif(500, 12, 60), group = sample(0:1, 500, TRUE, c(0.7, 0.3)),
                      fertile = runif(500) < 0.8)
  men <- data.frame(born = -runif(450, 12, 70), group = sample(0:1, 450, TRUE, c(0.6, 0.4)))
  yearly <- data.frame(age_start = 15:49, age_end = 16:50, rate = runif(35, 0, 0.3))
  matching <- age_gap_matching_model(1, 5, group_exogamy = c("A:B" = -2))
  for (fertility in c(FALSE, TRUE)) {
    expected <- reference_weights(women, men, matching, yearly, c("A", "B"), fertility = fertility)
    default <- calibrate(women, men, matching, yearly, c("A", "B"), fertility = fertility)
    expect_equal(default$log_weights, expected$log_weights, tolerance = 0.02)
  }
})
