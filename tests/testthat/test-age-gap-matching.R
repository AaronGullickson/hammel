test_that("age-gap preferences give symmetric half-peak selection weights", {
  for (preferred in c(-3, 0, 3)) {
    matching <- age_gap_matching_model(preferred, 5)
    expect_s3_class(matching, "hammel_matching_model")
    probabilities <- calculate_matching_probabilities(matching,
      preferred + c(-10, -5, 0, 5, 10))
    expect_equal(probabilities / probabilities[3], c(1/16, 1/2, 1, 1/2, 1/16))
    expect_equal(prepare_matching_model(matching)$quadratic,
      c(2 * preferred * log(2) / 25, -log(2) / 25))
  }
})

test_that("age-gap helper rejects invalid or unrepresentable preferences", {
  for (bad in list(NA_real_, Inf, NaN, numeric(), c(1, 2), "3", TRUE, matrix(3))) {
    expect_error(age_gap_matching_model(bad, 5), "preferred_age_gap")
    expect_error(age_gap_matching_model(3, bad), "half_weight_distance")
  }
  for (bad in c(0, -1))
    expect_error(age_gap_matching_model(3, bad), "positive")
  for (bad in c(1e-300, 1e300))
    expect_error(age_gap_matching_model(3, bad), class = "hammel_error_argument")
  expect_error(
    age_gap_matching_model(1e308, 0.1),
    class = "hammel_error_argument"
  )
})
test_that("union constructors default to zero-gap five-year half-weight matching", {
  expected <- age_gap_matching_model(0, 5)
  rates <- data.frame(age_start = 18, age_end = Inf, rate = 0.1)
  for (mode in c("on_demand", "rate")) {
    schedule <- if (mode == "rate") rates else NULL
    process <- union_formation(mode, rates = schedule)
    expect_equal(process$matching$coefficients, expected$coefficients)
    expect_equal(calculate_matching_probabilities(process$matching, c(-5, 0, 5)),
      c(0.25, 0.5, 0.25))
    model <- set_union_formation(hammel_model(), mode, rates = schedule)
    expect_equal(model$processes$union_formation$matching$coefficients,
      expected$coefficients)
    custom <- age_gap_matching_model(3, 7)
    expect_identical(union_formation(mode, custom, rates = schedule)$matching, custom)
    expect_identical(set_union_formation(hammel_model(), mode, custom,
      rates = schedule)$processes$union_formation$matching, custom)
  }
})


test_that("matching summaries recover age-gap parameters and settings", {
  result <- summary(
    age_gap_matching_model(
      preferred_age_gap = 3,
      half_weight_distance = 5,
      group_exogamy = list(race = c("W:B" = -3, "W:O" = NA), region = NA),
      delta = list(race = c(W = 0.1, B = 0))
    )
  )
  expect_s3_class(result, "summary.hammel_matching_model")
  expect_equal(result$age_gap$preferred_age_gap, 3)
  expect_equal(result$age_gap$half_weight_distance, 5)
  text <- cli::cli_fmt(value <- print(result))
  expect_identical(value, result)
  expect_true(any(grepl("Preferred age gap: 3 years", text, fixed = TRUE)))
  expect_true(any(grepl(
    "Group exogamy (race): W:B = -3, W:O = prohibited", text, fixed = TRUE
  )))
  expect_true(any(grepl(
    "Group exogamy (region): all pairs prohibited", text, fixed = TRUE
  )))
  expect_true(any(grepl("Delta (race): W = 0.1, B = 0", text, fixed = TRUE)))

  linear <- summary(
    matching_model(
      ~ age_difference + same_race,
      c(age_difference = 0.1, same_race = 0.69),
      group_exogamy = c("A:B" = -1)
    )
  )
  expect_null(linear$age_gap)
  lines <- matching_summary_lines(linear)
  expect_identical(lines[["Group exogamy"]], "A:B = -1")
  expect_false("Delta" %in% names(lines))

  table <- data.frame(
    race_x = "W", region_x = "S", race_y = "B", region_y = "S",
    log_odds = NA
  )
  lines <- matching_summary_lines(
    summary(matching_model(~1, numeric(), group_exogamy = table))
  )
  expect_identical(lines[["Coefficients"]], "none")
  expect_identical(
    lines[["Group exogamy (race|region)"]],
    "W|S:B|S = prohibited"
  )

  long <- data.frame(
    race_x = "W", race_y = LETTERS[1:12], log_odds = -seq_len(12)
  )
  text <- format_exogamy_table(long, "race")
  expect_match(text, "W:J = -10, and 2 more", fixed = TRUE)
})

test_that("union formation summaries nest matching models", {
  process <- union_formation(
    "rate",
    usa_2000_marriage,
    cohab_prob = tibble::tibble(age_start = 15, age_end = Inf, prob = 0.4),
    cohab_matching = age_gap_matching_model(1, 4)
  )
  result <- summary(process)
  expect_s3_class(result$matching, "summary.hammel_matching_model")
  expect_null(result$father_matching)
  expect_null(result$search_weights)
  expect_s3_class(
    result$cohabitation$matching,
    "summary.hammel_matching_model"
  )
  text <- cli::cli_fmt(print(result))
  expect_true(any(grepl("Father matching model: same as", text, fixed = TRUE)))
  expect_true(any(grepl("Preferred age gap: 1 years", text, fixed = TRUE)))
  expect_true(any(grepl("Search weight model: none", text, fixed = TRUE)))
  expect_true(result$cancel_birth_if_no_match)
  expect_true(any(grepl(
    "Births without a matched father: cancelled", text, fixed = TRUE
  )))

  process <- union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  text <- cli::cli_fmt(print(summary(process)))
  expect_true(any(grepl(
    "Births without a matched father: allowed", text, fixed = TRUE
  )))
  expect_false(any(grepl("Father matching", text, fixed = TRUE)))
  expect_false(any(grepl("Search weight", text, fixed = TRUE)))
})
