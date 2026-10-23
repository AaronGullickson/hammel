# weighted_choice (D066): female choice with male-group weights calibrated at
# each search so expected male group flows meet the mirrored target.

prepare_weighted <- function(model, groups) {
  prepare_matching_model(model, groups, weighted = TRUE)
}

calibrate <- function(women, men, matching, rates, groups, now = 0, cutoff = 15,
                      width = 0.25, fertility = FALSE) {
  inspect_weighted_calibration_cpp(women$born, women$group, women$fertile,
    men$born, men$group, rates, prepare_weighted(matching, groups), length(groups),
    now, cutoff, width, fertility)
}

# Independent individual-level reference (two groups): exact ages, no bins.
reference_weights <- function(women, men, matching, rates, groups, now = 0,
                              cutoff = 15, fertility = FALSE, absorption = 1) {
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
  target <- female_share + absorption * (mirror - female_share)
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

group_only <- function(women_share, men_share, effect, absorption = 1) {
  women <- data.frame(born = -30, group = rep(0:1, round(1000 * women_share)), fertile = TRUE)
  men <- data.frame(born = -30, group = rep(0:1, round(1000 * men_share)))
  matching <- create_matching_model(~ 1, numeric(), group_exogamy = c("A::B" = effect),
                                    imbalance_absorption = absorption)
  list(women = women, men = men, matching = matching,
       rates = data.frame(age_start = 15, age_end = 50, rate = 0.1))
}

test_that("calibration reproduces the analytical group-only examples", {
  # Women and men 75/25, A::B = -3: symmetric table, b_B = 1.6830 (spec example).
  x <- group_only(c(0.75, 0.25), c(0.75, 0.25), -3)
  result <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
  expect_false(result$infeasible)
  expect_equal(result$log_weights[1], 0)
  expect_equal(exp(result$log_weights[2]), 1.682979, tolerance = 1e-6)
  expect_equal(result$flow, c(0.75, 0.25), tolerance = 1e-9)
  # Men 70/30: full absorption b_B = 3.875; zero absorption targets female shares.
  for (case in list(list(1, 3.875051, c(0.70, 0.30)), list(0, 1.308984, c(0.75, 0.25)),
                    list(0.5, 2.400323, c(0.725, 0.275)))) {
    x <- group_only(c(0.75, 0.25), c(0.70, 0.30), -3, absorption = case[[1]])
    result <- calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))
    expect_equal(exp(result$log_weights[2]), case[[2]], tolerance = 1e-6)
    expect_equal(result$target, case[[3]], tolerance = 1e-12)
    expect_equal(result$flow, case[[3]], tolerance = 1e-9)
  }
  # Balanced composition is symmetric at every absorption value (weights equal).
  for (absorption in c(0, 0.5, 1)) {
    x <- group_only(c(0.5, 0.5), c(0.5, 0.5), -3, absorption)
    expect_equal(calibrate(x$women, x$men, x$matching, x$rates, c("A", "B"))$log_weights,
                 c(0, 0), tolerance = 1e-9)
  }
})

test_that("binned calibration agrees with an exact individual-level reference", {
  set.seed(901)
  women <- data.frame(born = -runif(400, 10, 80), group = sample(0:1, 400, TRUE, c(0.7, 0.3)),
                      fertile = runif(400) < 0.8)
  men <- data.frame(born = -runif(350, 10, 85), group = sample(0:1, 350, TRUE, c(0.6, 0.4)))
  rates <- data.frame(age_start = c(15, 30), age_end = c(30, 50), rate = c(0.2, 0.05))
  matching <- create_age_gap_matching_model(2, 5, group_exogamy = c("A:B" = -2, "B:A" = -1))
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
  custom <- create_matching_model(~ abs(age_difference - 2),
    c("abs(age_difference - 2)" = -0.3), group_exogamy = c("A::B" = -1))
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
  matching <- create_age_gap_matching_model(0, 5, group_exogamy = c("A::B" = -2))
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
  matching <- create_age_gap_matching_model(0, 5,
    group_exogamy = c("A::B" = -3, "B:C" = -1, "C:B" = -0.5), imbalance_absorption = 0.5)
  result <- calibrate(women, men, matching, rates, c("A", "B", "C"))
  expect_false(result$infeasible)
  expect_equal(result$flow, result$target, tolerance = 1e-9)
  expect_equal(sum(result$target), 1)
  # Men of a group with no eligible women still receive the population-wide
  # per-capita rate; with all exogamy prevented they cannot be reached.
  women$group <- 0L
  blocked <- create_matching_model(~ 1, numeric(), prevent_exogamy = TRUE)
  result <- calibrate(women, men, blocked, rates, c("A", "B", "C"))
  expect_true(result$infeasible)
  expect_equal(result$log_weights[2:3], c(50, 50))
  # With zero absorption those groups have zero target and effectively zero weight.
  zero <- create_matching_model(~ 1, numeric(), prevent_exogamy = TRUE, imbalance_absorption = 0)
  result <- calibrate(women, men, zero, rates, c("A", "B", "C"))
  expect_false(result$infeasible)
  expect_equal(result$target, c(1, 0, 0))
  expect_true(all(result$log_weights[2:3] <= -1000))
  # No women at risk leaves unit weights.
  idle <- calibrate(women, men, matching, data.frame(age_start = 15, age_end = 50, rate = 0),
                    c("A", "B", "C"))
  expect_equal(idle$log_weights, c(0, 0, 0))
})

weighted_population <- function(n = 600, weights = NULL, seed = 903) {
  set.seed(seed)
  create_baseline_population(n, max_age = 60, end_time = 0, group_weights = weights)
}

weighted_model <- function(population, method, initiation = "on_demand", ...) {
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01))
  if (initiation == "on_demand") {
    model |>
      set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.12)) |>
      set_union_formation("on_demand", matching_method = method, ...)
  } else {
    model |>
      set_fertility(data.frame(age_start = 15, age_end = 45, rate = c(0.15, 0.05),
                               union_status = c("partnered", "unpartnered"))) |>
      set_union_formation("rate", rates = data.frame(age_start = 15, age_end = 50, rate = 0.15),
                          matching_method = method, ...)
  }
}

test_that("one group or no groups reproduce explicit female choice exactly", {
  for (weights in list(NULL, c(A = 1))) {
    population <- weighted_population(400, weights)
    for (initiation in c("on_demand", "rate")) {
      female <- simulate(weighted_model(population, "female_choice", initiation), 15, seed = 904)
      set.seed(1); before <- .Random.seed
      weighted <- simulate(weighted_model(population, "weighted_choice", initiation), 15, seed = 904)
      expect_identical(.Random.seed, before)
      expect_identical(weighted, female)
    }
  }
})

test_that("grouped weighted choice runs in both modes, restarts, and is reproducible", {
  population <- weighted_population(600, c(A = 3, B = 1))
  matching <- create_age_gap_matching_model(0, 5, group_exogamy = c("A::B" = -2))
  for (initiation in c("on_demand", "rate")) {
    model <- weighted_model(population, "weighted_choice", initiation, matching = matching)
    first <- simulate(model, 20, seed = 905)
    expect_identical(first, simulate(model, 20, seed = 905))
    expect_gt(nrow(first$unions), nrow(population$unions))
    expect_false(identical(first$unions,
      simulate(weighted_model(population, "female_choice", initiation, matching = matching),
               20, seed = 905)$unions))
    restart <- simulate(set_population(model, first), 5, seed = 906)
    expect_equal(restart$end_time, first$end_time + 5)
    expect_identical(restart$persons[seq_len(nrow(first$persons)), c("person_id", "birth_time")],
                     first$persons[c("person_id", "birth_time")])
    # Diagnostics consume no random draws and do not change histories.
    diagnosed <- simulate(weighted_model(population, "weighted_choice", initiation,
                                         matching = matching, diagnostics = TRUE), 20, seed = 905)
    expect_identical(diagnosed$persons, first$persons)
    expect_identical(diagnosed$unions[names(first$unions)], first$unions)
    added <- !is.na(diagnosed$unions$.n_eligible)
    expect_true(any(added))
    expect_true(all(diagnosed$unions$.utility_rank[added] >= 1))
  }
  # Rate mode: unpartnered births select fathers through their own calibration.
  father <- create_age_gap_matching_model(3, 4, group_exogamy = c("A::B" = -1))
  result <- simulate(weighted_model(population, "weighted_choice", "rate", matching = matching,
                                    father_matching = father), 20, seed = 907)
  persons <- result$persons
  born <- persons[persons$birth_time >= 0, ]
  in_union <- mapply(function(mother, time) any(result$unions$wife_id == mother &
    result$unions$start_time <= time & (is.na(result$unions$end_time) | result$unions$end_time > time)),
    born$mother_id, born$birth_time)
  expect_true(any(!in_union))
  expect_true(all(!is.na(born$father_id)))
})

test_that("weighted choice warns once when weights must be approximated", {
  population <- weighted_population(300, c(A = 1))
  persons <- population$persons
  persons$group[persons$sex == "m"][1:40] <- "B"
  population <- hammel_population(persons, population$unions, end_time = 0)
  matching <- create_matching_model(~ age_difference, c(age_difference = 0),
                                    prevent_exogamy = TRUE)
  model <- weighted_model(population, "weighted_choice", matching = matching)
  warnings <- character()
  withCallingHandlers(simulate(model, 10, seed = 908),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_length(grep("weighted_choice found no exact", warnings), 1)
})

test_that("matching method and absorption inputs are validated", {
  expect_identical(union_formation("on_demand", matching_method = "weighted_choice")$matching_method,
                   "weighted_choice")
  expect_error(union_formation("on_demand", matching_method = "weighted"), "weighted_choice")
  for (bad in list(-0.1, 1.1, NA_real_, c(0.5, 0.5), "1", Inf))
    expect_error(create_age_gap_matching_model(0, 5, imbalance_absorption = bad),
                 "imbalance_absorption")
  expect_identical(create_matching_model(~ 1, numeric())$imbalance_absorption, 1)
  # Older matching objects without the field use the default.
  legacy <- create_age_gap_matching_model(0, 5)
  legacy$imbalance_absorption <- NULL
  expect_identical(prepare_matching_model(legacy, weighted = TRUE)$imbalance_absorption, 1)
  edited <- create_age_gap_matching_model(0, 5)
  edited$imbalance_absorption <- 2
  expect_error(union_formation("on_demand", edited), "imbalance_absorption")
})
