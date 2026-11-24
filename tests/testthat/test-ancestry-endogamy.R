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
  # Valid terms reach the guard that awaits the engine (stage 2).
  for (x in list(c(A = 1), c(unknown = 1), list(group = c(B = 2, A = 1)))) {
    error <- expect_hammel_error(
      run_simulation(endogamy_model(term(x)), 1, seed = 1),
      "specification",
      "run_simulation"
    )
    expect_identical(error$argument, "ancestry_endogamy")
  }
  error <- expect_hammel_error(
    run_simulation(
      endogamy_model(father_matching = term(c(A = 1))), 1, seed = 1
    ),
    "specification",
    "run_simulation"
  )
  expect_identical(error$argument, "ancestry_endogamy")
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
  expect_hammel_error(
    run_simulation(
      endogamy_model(term(list(region = c(N = 1))), groups = both), 1,
      seed = 1
    ),
    "specification",
    "run_simulation"
  )
  # Models without the term run as before.
  expect_no_error(run_simulation(endogamy_model(), 1, seed = 1))
})
