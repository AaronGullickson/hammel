# Ancestry endogamy in matching models (D199): interface and checks.

endogamy_model <- function(
  matching = age_gap_matching_model(0, 5),
  groups = list(group = group(ancestry = TRUE)),
  father_matching = NULL
) {
  persons <- data.frame(
    person_id = 1:4,
    sex = c("f", "m", "f", "m"),
    birth_time = -20,
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_,
    group = c("A", "A", "B", "B"),
    region = c("N", "S", "N", "S")
  )
  population <- hammel_population(persons, end_time = 0, groups = groups)
  make_model(0, birth_rate = 1, birth_interval = 0) |>
    set_population(population) |>
    set_union_formation(
      "on_demand",
      matching = matching,
      father_matching = father_matching,
      cancel_birth_if_no_match = FALSE
    )
}

test_that("matching models store ancestry endogamy and its resolution", {
  plain <- matching_model(~ age_difference, c(age_difference = 0))
  expect_null(plain$ancestry_endogamy)
  expect_equal(plain$ancestry_resolution, 1 / 8)
  weights <- c(irish = 2, nw_europe = 2)
  model <- matching_model(
    ~ 1,
    numeric(),
    ancestry_endogamy = weights,
    ancestry_resolution = 1 / 16
  )
  expect_identical(model$ancestry_endogamy, weights)
  expect_equal(model$ancestry_resolution, 1 / 16)
  gap <- age_gap_matching_model(
    3, 5,
    ancestry_endogamy = list(birthplace = weights),
    ancestry_resolution = 0.1
  )
  expect_identical(gap$ancestry_endogamy, list(birthplace = weights))
  expect_equal(gap$ancestry_resolution, 0.1)
  # Negative weights and `unknown` are allowed.
  expect_silent(matching_model(
    ~ 1, numeric(), ancestry_endogamy = c(A = -1, unknown = 0.5)
  ))
})

test_that("malformed ancestry endogamy weights are specification errors", {
  bad <- list(
    c(2, 2),
    c(A = NA),
    c(A = Inf),
    c(A = 1, A = 2),
    stats::setNames(1, ""),
    numeric(),
    c(A = "1"),
    list(c(A = 1)),
    list(group = c(A = 1), group = c(B = 1)),
    list(group = c(1, 2)),
    list()
  )
  for (x in bad) {
    error <- expect_hammel_error(
      matching_model(~ 1, numeric(), ancestry_endogamy = x),
      "specification",
      "matching_model"
    )
    expect_identical(error$argument, "ancestry_endogamy")
  }
})

test_that("ancestry resolution must be one over a whole number", {
  for (r in list(1, 1 / 3, 1 / 16, 0.1, 0.5)) {
    expect_silent(matching_model(
      ~ 1, numeric(), ancestry_endogamy = c(A = 1), ancestry_resolution = r
    ))
  }
  for (r in list(0, -0.125, 0.3, 2, NA_real_, Inf, c(0.1, 0.2), "1/8")) {
    error <- expect_hammel_error(
      matching_model(~ 1, numeric(), ancestry_resolution = r),
      "argument",
      "matching_model"
    )
    expect_identical(error$argument, "ancestry_resolution")
  }
})

test_that("summaries report ancestry endogamy and resolution", {
  plain <- matching_summary_lines(summary(age_gap_matching_model(0, 5)))
  expect_false(any(grepl("^Ancestry", names(plain))))
  bare <- matching_summary_lines(summary(matching_model(
    ~ 1, numeric(), ancestry_endogamy = c(irish = 2, nw_europe = 1.5)
  )))
  expect_identical(bare[["Ancestry endogamy"]], "irish = 2, nw_europe = 1.5")
  expect_identical(bare[["Ancestry resolution"]], "1/8")
  listed <- matching_summary_lines(summary(matching_model(
    ~ 1,
    numeric(),
    ancestry_endogamy = list(birthplace = c(irish = 2)),
    ancestry_resolution = 0.1
  )))
  expect_identical(listed[["Ancestry endogamy (birthplace)"]], "irish = 2")
  expect_identical(listed[["Ancestry resolution"]], "1/10")
})

test_that("ancestry endogamy is checked against tracked dimensions", {
  term <- function(x) matching_model(~ 1, numeric(), ancestry_endogamy = x)
  # Valid terms run, including in a father model.
  for (x in list(c(A = 1), c(unknown = 1), list(group = c(B = 2, A = 1)))) {
    expect_no_error(run_simulation(endogamy_model(term(x)), 1, seed = 1))
  }
  expect_no_error(run_simulation(
    endogamy_model(father_matching = term(c(A = 1))), 1, seed = 1
  ))
  # Unknown categories and dimensions without ancestry are group errors.
  error <- expect_hammel_error(
    run_simulation(endogamy_model(term(c(C = 1))), 1, seed = 1),
    "groups",
    "run_simulation"
  )
  expect_identical(error$group, "group")
  error <- expect_hammel_error(
    run_simulation(endogamy_model(term(list(region = c(N = 1)))), 1, seed = 1),
    "groups",
    "run_simulation"
  )
  expect_identical(error$group, "region")
  expect_hammel_error(
    run_simulation(
      endogamy_model(term(c(A = 1)), groups = list(group = group())), 1,
      seed = 1
    ),
    "groups",
    "run_simulation"
  )
  # A bare vector needs a single dimension that tracks ancestry.
  both <- list(group = group(ancestry = TRUE), region = group(ancestry = TRUE))
  expect_hammel_error(
    run_simulation(endogamy_model(term(c(A = 1)), groups = both), 1, seed = 1),
    "groups",
    "run_simulation"
  )
  expect_no_error(run_simulation(
    endogamy_model(term(list(region = c(N = 1))), groups = both), 1, seed = 1
  ))
})

test_that("models with the term must share dimensions and resolution", {
  both <- list(group = group(ancestry = TRUE), region = group(ancestry = TRUE))
  general <- matching_model(
    ~ 1, numeric(), ancestry_endogamy = list(group = c(A = 1))
  )
  # Different weights are allowed.
  same <- matching_model(
    ~ 1, numeric(), ancestry_endogamy = list(group = c(B = 3))
  )
  expect_no_error(run_simulation(
    endogamy_model(general, both, father_matching = same), 1, seed = 1
  ))
  other_dim <- matching_model(
    ~ 1, numeric(), ancestry_endogamy = list(region = c(N = 1))
  )
  finer <- matching_model(
    ~ 1,
    numeric(),
    ancestry_endogamy = list(group = c(A = 1)),
    ancestry_resolution = 1 / 16
  )
  for (father in list(other_dim, finer)) {
    error <- expect_hammel_error(
      run_simulation(
        endogamy_model(general, both, father_matching = father), 1, seed = 1
      ),
      "specification",
      "run_simulation"
    )
    expect_identical(error$argument, "ancestry_endogamy")
  }
  # A model without the term may join a run that uses it.
  expect_no_error(run_simulation(
    endogamy_model(
      general, both,
      father_matching = age_gap_matching_model(0, 5)
    ), 1, seed = 1
  ))
})

# ENGINE ------------------------------------------------------------------

test_that("profiles round jointly to whole units, ties to earlier columns", {
  shares <- rbind(
    c(0.375, 0.625, 0),
    c(0.35, 0.35, 0.30),
    c(0.3, 0.7, 0),
    c(0.04, 0.96, 0),
    c(0.5, 0, 0.5),
    c(1 / 3, 1 / 3, 1 / 3)
  )
  eighths <- inspect_round_profile_cpp(shares, 8L)
  expect_equal(
    eighths,
    rbind(c(3, 5, 0), c(3, 3, 2), c(2, 6, 0), c(0, 8, 0), c(4, 0, 4),
      c(3, 3, 2))
  )
  expect_true(all(rowSums(eighths) == 8))
  tenths <- inspect_round_profile_cpp(shares[1:2, ], 10L)
  expect_equal(tenths, rbind(c(4, 6, 0), c(4, 3, 3)))
})

# Shares of irish, nw_europe and unknown, one row per person.
profiles <- list(
  I = c(1, 0, 0),
  N = c(0, 1, 0),
  H = c(0.5, 0.5, 0)
)
layout <- list(offsets = 0L, sizes = 3L, units = 8L)

prepared_endogamy <- function(lambda, model = matching_model(~ 1, numeric()),
                              groups = NULL) {
  prepared <- prepare_matching_model(
    as_internal_matching(model, groups), groups
  )
  prepared$ancestry_lambda <- lambda
  prepared
}

typed_calibration <- function(women, men, prepared, groups = 1L, now = 0) {
  inspect_weighted_calibration_cpp(
    women$born, women$group, rep(TRUE, nrow(women)), men$born, men$group,
    data.frame(age_start = 15, age_end = 50, rate = 0.1), prepared, groups,
    now, 15, 0.25, FALSE,
    ancestry_layout = layout,
    female_shares = do.call(rbind, profiles[women$profile]),
    male_shares = do.call(rbind, profiles[men$profile])
  )
}

profile_people <- function(shares, n = 1000) {
  counts <- round(n * shares)
  data.frame(
    born = -30,
    group = 0L,
    profile = rep(names(shares), counts)
  )
}

# Log weight of each profile's male type, relative to the first in `order`.
profile_weights <- function(result, men, order) {
  type <- result$male_type[match(order, men$profile)]
  w <- result$log_weights[type + 1L]
  stats::setNames(w - w[1], order)
}

test_that("pure profiles reproduce the group exogamy calibration", {
  # I-N pairs score -(1.5 + 1.5) = -3, as group_exogamy "A:B" = -3, so the
  # 75/25 example gives b_N = 1.682979 (weighted-choice specification).
  people <- profile_people(c(I = 0.75, N = 0.25))
  result <- typed_calibration(
    people, people, prepared_endogamy(c(1.5, 1.5, 0))
  )
  expect_false(result$infeasible)
  w <- profile_weights(result, people, c("I", "N"))
  expect_equal(exp(w[["N"]]), 1.682979, tolerance = 1e-6)
  expect_equal(result$flow, result$target, tolerance = 1e-9)
})

test_that("three-profile calibration gives the specification's table", {
  shares <- c(I = 0.5, N = 0.3, H = 0.2)
  people <- profile_people(shares)
  result <- typed_calibration(
    people, people, prepared_endogamy(c(1.5, 1.5, 0))
  )
  w <- exp(profile_weights(result, people, names(shares)))
  expect_equal(unname(w), c(1, 1.2556, 1.2767), tolerance = 1e-4)
  # The expected table is symmetric, not only balanced in its margins.
  U <- -1.5 * outer(
    names(shares), names(shares),
    Vectorize(function(a, b) sum(abs(profiles[[a]] - profiles[[b]])))
  )
  P <- exp(U) * rep(w * shares, each = 3)
  table <- shares * P / rowSums(P)
  expect_equal(table, t(table), tolerance = 1e-8)
  expect_equal(table[1, 3], 0.04948, tolerance = 1e-4)
})

test_that("typed calibration meets exact person-level targets", {
  # Ages, two categories with an exogamy effect, and three profiles. People
  # are born at bin centres, so binned and exact ages agree.
  set.seed(31)
  n <- 120
  person <- function(n) {
    data.frame(
      born = -(sample(60:160, n, replace = TRUE) + 0.5) / 4,
      group = sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3)),
      profile = sample(names(profiles), n, replace = TRUE,
        prob = c(0.5, 0.3, 0.2))
    )
  }
  women <- person(n)
  men <- person(n)
  groups <- c("A", "B")
  model <- age_gap_matching_model(2, 6, group_exogamy = c("A:B" = -1))
  lambda <- c(1.2, 0.8, 0)
  result <- typed_calibration(
    women, men, prepared_endogamy(lambda, model, groups), groups = 2L
  )
  expect_false(result$infeasible)
  # Reference: exact utilities for every pair (all are past the cutoff and
  # inside the rate band) and each type's mirrored target.
  internal <- as_internal_matching(model, groups)
  share <- function(x) do.call(rbind, profiles[x$profile])
  wf <- share(women)
  wm <- share(men)
  U <- vapply(seq_len(nrow(men)), function(j) {
    calculate_matching_utilities(
      internal, women$born - men$born[j], groups[women$group + 1L],
      rep(groups[men$group[j] + 1L], nrow(women))
    ) - drop(abs(wf - rep(wm[j, ], each = nrow(women))) %*% lambda)
  }, numeric(nrow(women)))
  types <- length(result$type_group)
  f_type <- factor(result$female_type, 0:(types - 1L))
  m_type <- factor(result$male_type, 0:(types - 1L))
  women_n <- as.numeric(table(f_type))
  men_n <- as.numeric(table(m_type))
  mirror <- men_n * ifelse(women_n > 0, 0.1, 0.1)
  target <- mirror / sum(mirror)
  b <- exp(result$log_weights)[as.integer(m_type)]
  W <- exp(U) * rep(b, each = nrow(women))
  P <- W / rowSums(W)
  flow <- vapply(
    levels(m_type), function(t) sum(P[, m_type == t]), numeric(1)
  ) / nrow(women)
  expect_equal(unname(flow), target, tolerance = 1e-8)
  expect_equal(result$target, target, tolerance = 1e-12)
})

test_that("exact draws add the ancestry term to each candidate", {
  set.seed(32)
  n <- 40
  men <- data.frame(
    born = -runif(n, 16, 40),
    group = sample(0:1, n, replace = TRUE),
    profile = sample(names(profiles), n, replace = TRUE)
  )
  groups <- c("A", "B")
  model <- age_gap_matching_model(3, 5, group_exogamy = c("A:B" = -0.5))
  lambda <- c(2, 1, 0)
  prepared <- prepared_endogamy(lambda, model, groups)
  shares <- do.call(rbind, profiles[men$profile])
  female <- profiles$H
  log_weights <- c(0, 0.4)
  result <- inspect_exact_draw_cpp(
    men$born, men$group, logical(n), -25, 0L, log_weights, prepared, 2L, 0,
    15, 60000,
    ancestry_layout = layout, male_shares = shares, female_shares = female
  )
  internal <- as_internal_matching(model, groups)
  score <- log_weights[men$group + 1L] +
    calculate_matching_utilities(
      internal, -25 - men$born, rep("A", n), groups[men$group + 1L]
    ) -
    drop(abs(shares - rep(female, each = n)) %*% lambda)
  expected <- exp(score - max(score))
  expected <- expected / sum(expected)
  expect_equal(result$probability, unname(expected), tolerance = 1e-12)
  expect_equal(result$approximate, 0)
  observed <- result$counts / sum(result$counts)
  expect_lt(max(abs(observed - expected)), 0.01)
})

test_that("simulations with ancestry endogamy are reproducible and favor it", {
  # Two equal groups of pure Irish and pure NW European ancestry, one matching
  # category. A strong term makes mixed unions rare.
  set.seed(33)
  n <- 800
  persons <- data.frame(
    person_id = seq_len(n),
    sex = rep(c("f", "m"), n / 2),
    birth_time = -runif(n, 15, 30),
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_,
    heritage = rep(c("irish", "irish", "nw", "nw"), n / 4)
  )
  population <- hammel_population(
    persons, end_time = 0, groups = list(heritage = group(ancestry = TRUE))
  )
  run <- function(lambda) {
    model <- hammel_model() |>
      set_population(population) |>
      set_mortality(
        data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf,
          rate = 0.005)
      ) |>
      set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.15)) |>
      set_union_formation(
        "on_demand",
        matching = age_gap_matching_model(
          2, 5,
          ancestry_endogamy = c(irish = lambda, nw = lambda)
        )
      )
    run_simulation(model, 10, seed = 34)
  }
  mixed_share <- function(result) {
    u <- result$unions[result$unions$start_time >= 0, ]
    a <- result$ancestry$heritage
    wife <- a$irish[match(u$wife_id, a$person_id)]
    husband <- a$irish[match(u$husband_id, a$person_id)]
    mean(wife != husband)
  }
  strong <- run(3)
  expect_identical(run(3), strong)
  expect_gt(nrow(strong$unions), 100)
  expect_lt(mixed_share(strong), 0.1)
  expect_gt(mixed_share(run(0)), 0.3)
  # Restarts rebuild profiles from the recorded ancestry.
  restart <- run_simulation(
    set_population(
      hammel_model() |>
        set_mortality(
          data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf,
            rate = 0.005)
        ) |>
        set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.15)) |>
        set_union_formation(
          "on_demand",
          matching = age_gap_matching_model(
            2, 5, ancestry_endogamy = c(irish = 3, nw = 3)
          )
        ),
      strong
    ),
    5, seed = 35
  )
  expect_gt(nrow(restart$unions), nrow(strong$unions))
})

test_that("rate-mode searches calibrate by type with consistent cells", {
  set.seed(36)
  n <- 600
  persons <- data.frame(
    person_id = seq_len(n),
    sex = rep(c("f", "m"), n / 2),
    birth_time = -runif(n, 15, 35),
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_,
    race = sample(c("A", "B"), n, replace = TRUE),
    heritage = sample(c("irish", "nw"), n, replace = TRUE)
  )
  population <- hammel_population(
    persons,
    end_time = 0,
    groups = list(race = group(), heritage = group(ancestry = TRUE))
  )
  fertility <- data.frame(
    union_type = c("partnered", "unpartnered"),
    age_start = 15,
    age_end = 45,
    rate = c(0.2, 0.05)
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(
      data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
    ) |>
    set_fertility(fertility) |>
    set_union_formation(
      "rate",
      rates = data.frame(age_start = 15, age_end = Inf, rate = 0.3),
      matching = age_gap_matching_model(
        2, 5,
        group_exogamy = list(race = c("A:B" = -1)),
        ancestry_endogamy = list(heritage = c(irish = 1, nw = 1))
      ),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 0.4),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf,
        rate = 0.2),
      cancel_birth_if_no_match = FALSE
    ) |>
    set_union_dissolution(
      data.frame(dur_start = 0, dur_end = Inf, rate = 0.05)
    )
  withr::local_options(hammel.check_cells = TRUE)
  expect_no_warning(result <- run_simulation(model, 15, seed = 37))
  expect_gt(nrow(result$unions), 200)
  expect_identical(run_simulation(model, 15, seed = 37), result)
})

test_that("typed weights are held for a month unless set otherwise", {
  set.seed(38)
  n <- 600
  persons <- data.frame(
    person_id = seq_len(n),
    sex = rep(c("f", "m"), n / 2),
    birth_time = -runif(n, 15, 30),
    death_time = NA_real_,
    mother_id = NA_integer_,
    father_id = NA_integer_,
    heritage = sample(c("irish", "nw", "other"), n, replace = TRUE)
  )
  population <- hammel_population(
    persons, end_time = 0, groups = list(heritage = group(ancestry = TRUE))
  )
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(
      data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf,
        rate = 0.005)
    ) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = 0.15)) |>
    set_union_formation(
      "on_demand",
      matching = age_gap_matching_model(
        2, 5, ancestry_endogamy = c(irish = 1, nw = 1, other = 1)
      )
    )
  monthly <- run_simulation(model, 10, seed = 39)
  withr::with_options(
    list(hammel.ancestry_recalibration = 1 / 12),
    expect_identical(run_simulation(model, 10, seed = 39), monthly)
  )
  # Recalibrating at every search is a different, equally valid history.
  every <- withr::with_options(
    list(hammel.ancestry_recalibration = 0),
    run_simulation(model, 10, seed = 39)
  )
  expect_false(identical(every$unions, monthly$unions))
  expect_gt(nrow(every$unions), 100)
})
