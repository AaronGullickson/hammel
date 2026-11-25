# Sampled matching weights and reconstructed choice sets (D200).

choice_population <- function(n = 4000, seed = 33) {
  set.seed(seed)
  persons <- data.frame(
    person_id = seq_len(n),
    sex = rep(c("f", "m"), n / 2),
    birth_time = -runif(n, 0, 40),
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_,
    race = sample(c("A", "B"), n, replace = TRUE, prob = c(0.7, 0.3)),
    heritage = sample(c("irish", "nw", "se"), n, replace = TRUE)
  )
  hammel_population(
    persons,
    end_time = 0,
    groups = list(race = group(), heritage = group(ancestry = TRUE))
  )
}

choice_model <- function(population, matching) {
  hammel_model() |>
    set_population(population) |>
    set_mortality(
      data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
    ) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.12)) |>
    set_union_formation("on_demand", matching = matching)
}

ancestry_matching <- function() {
  age_gap_matching_model(
    3, 5,
    group_exogamy = list(race = c("A:B" = -2)),
    ancestry_endogamy = list(heritage = c(irish = 1, nw = 0.5, se = 1.5))
  )
}

# Conditional logit by maximum likelihood: one stratum per choice set.
fit_clogit <- function(x, chosen, strata, offset) {
  strata <- match(strata, unique(strata))
  loglik <- function(beta) {
    eta <- drop(x %*% beta) + offset
    top <- tapply(eta, strata, max)
    log_sum <- log(tapply(exp(eta - top[strata]), strata, sum)) + top
    -sum(eta[chosen] - log_sum[strata[chosen]])
  }
  fit <- stats::optim(
    numeric(ncol(x)), loglik, method = "BFGS", hessian = TRUE,
    control = list(maxit = 500, reltol = 1e-12)
  )
  list(
    coef = stats::setNames(fit$par, colnames(x)),
    se = sqrt(diag(solve(fit$hessian)))
  )
}

test_that("union_weight_sample is validated", {
  model <- choice_model(choice_population(400), ancestry_matching())
  for (bad in list(-0.1, 1.5, NA_real_, c(0.1, 0.2), "a")) {
    expect_hammel_error(
      run_simulation(model, 1, seed = 1, union_weight_sample = bad),
      "argument",
      "run_simulation"
    )
  }
})

test_that("recording weights leaves the simulation unchanged", {
  model <- choice_model(choice_population(2000), ancestry_matching())
  plain <- run_simulation(model, 20, seed = 5)
  kept <- run_simulation(model, 20, seed = 5, union_weight_sample = 1)
  part <- run_simulation(model, 20, seed = 5, union_weight_sample = 0.3)
  expect_null(plain$matching_weights)
  for (result in list(kept, part)) {
    expect_identical(result$persons, plain$persons)
    expect_identical(result$unions, plain$unions)
    expect_identical(result$ancestry, plain$ancestry)
  }

  # Every union formed in the run is kept at 1; a fraction keeps a subset,
  # chosen by union ID.
  formed <- plain$unions$union_id[plain$unions$start_time > 0]
  expect_setequal(kept$matching_weights$unions$union_id, formed)
  ids <- part$matching_weights$unions$union_id
  expect_true(all(ids %in% formed))
  expect_true(abs(length(ids) / length(formed) - 0.3) < 0.05)
  # Another history keeps the same union IDs.
  again <- run_simulation(model, 20, seed = 9, union_weight_sample = 0.3)
  formed_again <- again$unions$union_id[again$unions$start_time > 0]
  common <- intersect(formed, formed_again)
  expect_identical(
    common %in% again$matching_weights$unions$union_id,
    common %in% ids
  )

  # The record has typed weights for the chosen husbands' types.
  record <- kept$matching_weights
  expect_named(
    record$types,
    c("type_id", "race", "heritage_irish", "heritage_nw", "heritage_se",
      "heritage_unknown")
  )
  expect_false(anyNA(record$unions$calibration_id))
  expect_true(all(record$weights$calibration_id %in%
    record$unions$calibration_id))

  # A later run without recording drops the record.
  later <- run_simulation(set_population(model, kept), 1, seed = 2)
  expect_null(later$matching_weights)
})

test_that("models without calibrated weights record zero offsets", {
  population <- generate_uniform_population(1000)
  model <- choice_model(population, age_gap_matching_model(3, 5))
  result <- run_simulation(model, 15, seed = 3, union_weight_sample = 1)
  expect_true(all(is.na(result$matching_weights$unions$calibration_id)))
  expect_equal(nrow(result$matching_weights$weights), 0L)
  set.seed(4)
  sets <- extract_choice_sets(result, 50, 10)
  expect_true(all(sets$log_weight == 0))
})

test_that("choice sets hold the husband and eligible alternatives", {
  model <- choice_model(choice_population(2000), ancestry_matching())
  result <- run_simulation(model, 20, seed = 5, union_weight_sample = 0.5)
  set.seed(6)
  sets <- extract_choice_sets(result, 200, 20)
  expect_s3_class(sets, "tbl_df")
  expect_named(
    sets,
    c("choice_id", "union_type", "time", "chosen", "wife_id", "man_id",
      "birth_time_w", "birth_time_h", "race_w", "race_h",
      "heritage_irish_w", "heritage_irish_h", "heritage_nw_w",
      "heritage_nw_h", "heritage_se_w", "heritage_se_h",
      "heritage_unknown_w", "heritage_unknown_h", "log_weight")
  )
  expect_equal(length(unique(sets$choice_id)), 200L)
  per_set <- split(sets, sets$choice_id)
  for (set in per_set) {
    expect_equal(sum(set$chosen), 1L)
    expect_lte(nrow(set), 21L)
    expect_false(anyDuplicated(set$man_id) > 0)
  }
  # The chosen man is the union's husband; every candidate is eligible.
  u <- result$unions[match(sets$choice_id, result$unions$union_id), ]
  chosen <- sets[sets$chosen, ]
  expect_equal(chosen$man_id, u$husband_id[sets$chosen])
  expect_equal(sets$time, u$start_time)
  p <- result$persons[match(sets$man_id, result$persons$person_id), ]
  expect_true(all(p$sex == "m"))
  expect_true(all(sets$time - p$birth_time >= 15))
  expect_true(all(is.na(p$death_time) | p$death_time > sets$time))
  expect_false(any(is.na(sets$log_weight)))
  expect_true(all(sets$log_weight[sets$chosen] > -999))
  # Shares are rounded profiles that sum to 1.
  share_w <- rowSums(sets[grep("^heritage_.*_w$", names(sets))])
  expect_equal(share_w, rep(1, nrow(sets)))

  only <- extract_choice_sets(result, 10, 5, union_types = "cohabitation")
  expect_equal(nrow(only), 0L)
})

test_that("extract_choice_sets validates its inputs", {
  model <- choice_model(choice_population(400), ancestry_matching())
  plain <- run_simulation(model, 5, seed = 1)
  expect_hammel_error(
    extract_choice_sets(plain),
    "population",
    "extract_choice_sets"
  )
  kept <- run_simulation(model, 5, seed = 1, union_weight_sample = 1)
  expect_hammel_error(
    extract_choice_sets(kept, unions = 0),
    "argument",
    "extract_choice_sets"
  )
  expect_hammel_error(
    extract_choice_sets(kept, alternatives = 2.5),
    "argument",
    "extract_choice_sets"
  )
  expect_hammel_error(
    extract_choice_sets(kept, union_types = "dating"),
    "table",
    "extract_choice_sets"
  )
  expect_hammel_error(
    extract_choice_sets(list()),
    "population",
    "extract_choice_sets"
  )
})

test_that("a conditional logit on choice sets recovers the matching inputs", {
  model <- choice_model(choice_population(), ancestry_matching())
  result <- run_simulation(model, 40, seed = 5, union_weight_sample = 1)
  set.seed(9)
  sets <- extract_choice_sets(result, 3000, 30)
  gap <- sets$birth_time_w - sets$birth_time_h
  distance <- function(k) {
    -abs(sets[[paste0("heritage_", k, "_w")]] -
      sets[[paste0("heritage_", k, "_h")]])
  }
  x <- cbind(
    age = gap,
    age2 = gap^2,
    ab = as.numeric(sets$race_w != sets$race_h),
    irish = distance("irish"),
    nw = distance("nw"),
    se = distance("se")
  )
  fit <- fit_clogit(x, sets$chosen, sets$choice_id, sets$log_weight)
  coefficients <- age_gap_matching_model(3, 5)$coefficients
  truth <- c(
    age = coefficients[[1]],
    age2 = coefficients[[2]],
    ab = -2,
    irish = 1,
    nw = 0.5,
    se = 1.5
  )
  expect_true(all(abs(fit$coef - truth) < 4 * fit$se))
})
