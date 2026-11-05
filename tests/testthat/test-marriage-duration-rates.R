make_duration_rate_history <- function() {
  persons <- data.frame(person_id = 1:4, sex = c("f", "m", "m", "m"),
    birth_time = -20, death_time = c(NA, NA, 4, NA),
    mother_id = NA_integer_, father_id = NA_integer_)
  unions <- data.frame(union_id = 1:3, wife_id = 1, husband_id = 2:4,
    start_time = c(1, 3, 5), end_time = c(2, 4, NA),
    end_reason = c("dissolution", "death_husband", NA))
  hammel_population(persons, unions, end_time = 6, groups = if ("group" %in% names(persons)) "group")
}

test_that("duration breaks select the previous ending clock and enable stratification", {
  pop <- make_duration_rate_history()
  result <- calculate_union_rates(pop, c(0, 6), c(0, 21, Inf),
    dur_breaks = c(0, .5, 1, 2, Inf))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2 * (2 + 2 * 4))
  expect_true(all(is.na(result$dur_start[result$previous_union == "none"])))
  expect_true(all(is.na(result$age_start[result$previous_union != "none"])))
  women <- subset(result, sex == "f")
  none <- subset(women, previous_union == "none")
  expect_equal(none$person_years, c(1, 0))
  expect_equal(none$unions, c(0, 1))
  for (status in c("separated", "widowed")) {
    rows <- subset(women, previous_union == status)
    expect_equal(rows$person_years, c(.5, .5, 0, 0))
    expect_equal(rows$unions, c(0, 0, 1, 0))
    expect_equal(rows$rate_union, c(0, 0, NA, NA))
  }
  expect_identical(calculate_union_rates(pop, c(0, 6), c(0, 21, Inf), dur_breaks = c(0, .5, 1, 2, Inf)), result)
  # Finite final duration bounds gain Inf.
  expect_identical(calculate_union_rates(pop, c(0, 6), c(0, 21, Inf),
    dur_breaks = c(0, .5, 1, 2)), result)
})

test_that("duration exposure uses history before the window and resets after each ending", {
  pop <- make_duration_rate_history()
  result <- calculate_union_rates(pop, c(2.5, 3, 4, 4.5, 5, 6), c(0, Inf),
    dur_breaks = c(0, .5, 1, Inf))
  dissolved <- subset(result, sex == "f" & previous_union == "separated")
  expect_equal(sum(dissolved$person_years), .5)
  expect_equal(dissolved$dur_start[dissolved$person_years > 0], .5)
  event <- subset(dissolved, unions > 0)
  expect_equal(event$time_start, 3)
  expect_equal(event$dur_start, 1)
  widowed <- subset(result, sex == "f" & previous_union == "widowed")
  expect_equal(sum(widowed$person_years), 1)
  expect_equal(widowed$dur_start[widowed$person_years > 0], c(0, .5))
  expect_equal(subset(widowed, unions > 0)$dur_start, 1)
  before <- calculate_union_rates(pop, c(2.5, 3), c(0, Inf), dur_breaks = c(0, Inf))
  expect_equal(sum(before$unions), 0) # Last observation endpoint excludes events.
})

test_that("group/sex/period duration cells preserve pooled counts and exposure", {
  pop <- make_duration_rate_history()
  pop$persons$group <- c(10, 20, 10, 20)
  pop$groups <- list(group = group())
  times <- c(0, 2, 4, 6)
  result <- calculate_union_rates(pop, times, c(0, Inf),
    dur_breaks = c(0, .5, 1, 2, Inf), by = "group")
  age <- calculate_union_rates(pop, times, c(0, Inf), by = c("group", "previous_union"))
  expect_type(result$group, "double")
  for (i in seq_len(nrow(age))) {
    rows <- result$group == age$group[i] & result$sex == age$sex[i] &
      result$time_start == age$time_start[i] & result$previous_union == age$previous_union[i]
    expect_equal(sum(result$person_years[rows]), age$person_years[i])
    expect_equal(sum(result$unions[rows]), age$unions[i])
  }
  # Age bounds apply only to 'none' rows; duration rows include all attained ages.
  limited <- calculate_union_rates(pop, c(0, 6), c(50, Inf),
    dur_breaks = c(0, Inf))
  expect_equal(sum(subset(limited, previous_union == "none")$person_years), 0)
  expect_equal(sum(subset(limited, previous_union == "separated")$person_years), 5)
  expect_equal(sum(subset(limited, previous_union == "widowed")$person_years), 1)
})

test_that("duration observation handles zero-duration episodes and empty outputs", {
  pop <- make_duration_rate_history()
  pop$unions$start_time <- c(1, 2, 2)
  pop$unions$end_time <- c(2, 2, NA)
  pop$unions$end_reason[2] <- "dissolution"
  counts <- calculate_union_rates(pop, c(0, 6), c(0, Inf), dur_breaks = c(0, .5, Inf))
  women <- subset(counts, sex == "f" & previous_union == "separated")
  expect_equal(women$unions, c(2, 0))
  expect_equal(women$person_years, c(0, 0))
  expect_true(all(is.na(women$rate_union)))
  pop$persons$group <- "A"
  pop$groups <- list(group = group())
  empty <- hammel_population(pop$persons[FALSE, ], end_time = 6, groups = if ("group" %in% names(pop$persons)) "group")
  grouped <- calculate_union_rates(empty, c(0, 6), c(0, Inf), dur_breaks = c(0, Inf), by = "group")
  expect_equal(nrow(grouped), 0)
  expect_type(grouped$dur_start, "double")
  expect_type(grouped$previous_union, "character")
  pooled <- calculate_union_rates(empty, c(0, 6), c(0, Inf), dur_breaks = c(0, Inf))
  expect_equal(nrow(pooled), 6)
  expect_true(all(pooled$person_years == 0))
})

test_that("invalid duration breaks fail and NULL preserves the age-only contract", {
  pop <- make_duration_rate_history()
  for (bad in list(NA, numeric(), c(1, 1), c(-1, 2), c(0, NA, 2), matrix(c(0, 2)), c(0, Inf, 3)))
    expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf), dur_breaks = bad), "dur_breaks")
  expect_identical(calculate_union_rates(pop, c(0, 6), c(0, Inf)),
    calculate_union_rates(pop, c(0, 6), c(0, Inf), dur_breaks = NULL))
  expect_false("dur_start" %in% names(calculate_union_rates(pop, c(0, 6), c(0, Inf))))
})

test_that("observed duration rates recover known female search hazards", {
  n <- 3000L
  model <- make_rate_model(2L * n)
  model$population$persons$sex[seq_len(n)] <- "f"
  model$population$unions <- tibble::tibble(union_id = seq_len(n),
    wife_id = seq_len(n), husband_id = n + seq_len(n), start_time = -5,
    end_time = 0, end_reason = "dissolution")
  rates <- data.frame(previous_union = c("none", rep("separated", 3), "widowed"),
    age_start = c(15, NA, NA, NA, NA), age_end = c(Inf, NA, NA, NA, NA),
    dur_start = c(NA, 0, 1, 2, 0), dur_end = c(NA, 1, 2, Inf, Inf),
    rate = c(0, .1, .3, .2, 0))
  model <- set_union_formation(model, "rate", rates = rates, kinship_generations = 0L)
  result <- simulate(model, 4, 771)
  observed <- calculate_union_rates(result, c(0, 4), c(0, Inf), dur_breaks = c(0, 1, 2, Inf))
  women <- subset(observed, sex == "f" & previous_union == "separated")
  target <- c(.1, .3, .2)
  expect_true(all(abs(women$rate_union - target) < 6 * sqrt(target / women$person_years)))
  expect_equal(sum(women$unions), nrow(result$unions) - n)
})
