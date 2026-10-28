union_count_history <- function() {
  data.frame(union_id = 1:3, start_time = c(0, 1, 2),
    end_time = c(3, 4, NA_real_),
    end_reason = c("dissolution", "death_husband", NA_character_))
}

union_count_population <- function() {
  persons <- data.frame(person_id = 1:6, sex = rep(c("f", "m"), 3),
    birth_time = -30, death_time = c(NA, NA, NA, 4, NA, NA),
    mother_id = NA_integer_, father_id = NA_integer_, group = c(20L, 10L, 10L, 20L, 20L, 10L))
  unions <- union_count_history()
  unions$wife_id <- c(1L, 3L, 5L)
  unions$husband_id <- c(2L, 4L, 6L)
  hammel_population(persons, unions, end_time = 5)
}

test_that("union counts split exact exposure by calendar time and duration", {
  counts <- calculate_union_counts(union_count_history(), c(0, 2, 5), c(0, 2, Inf))
  expect_s3_class(counts, "tbl_df")
  expect_named(counts, c("time_start", "time_end", "dur_start", "dur_end",
                        "union_years", "dissolutions", "rate_dissolution"))
  expect_equal(counts$time_start, c(0, 2, 0, 2))
  expect_equal(counts$dur_start, c(0, 0, 2, 2))
  expect_equal(counts$union_years, c(3, 3, 0, 3))
  expect_equal(counts$dissolutions, c(0, 0, 0, 1))
  expect_equal(counts$rate_dissolution, c(0, 0, NA, 1/3))
  pooled <- calculate_union_counts(union_count_history(), c(0, 5), c(0, Inf))
  expect_equal(pooled$union_years, 9)
  expect_equal(pooled$dissolutions, 1)
  expect_equal(pooled$rate_dissolution, 1/9)
  finer <- calculate_union_counts(union_count_history(), seq(0, 5, .25), seq(0, 5, .25))
  expect_equal(sum(finer$union_years), pooled$union_years)
  expect_equal(sum(finer$dissolutions), pooled$dissolutions)
  # Actual population and bare union records have identical pooled observations.
  expect_identical(calculate_union_counts(union_count_population(), c(0, 2, 5), c(0, 2, Inf)), counts)
})

test_that("events respect internal boundaries, censoring and zero-duration histories", {
  unions <- data.frame(union_id = 1:4, start_time = c(-1, 0, 1, 0),
    end_time = c(1, 2, 1, 1),
    end_reason = c("dissolution", "dissolution", "dissolution", "death_wife"))
  counts <- calculate_union_counts(unions, c(0, 1, 2), c(0, 1, 2, Inf))
  # Dissolution at the upper endpoint is excluded; zero-duration dissolution
  # belongs to the duration-zero cell at its event time.
  expect_equal(counts$dissolutions, c(0, 1, 0, 0, 0, 1))
  expect_equal(sum(counts$union_years), 4)
  late_long <- subset(counts, time_start == 1 & dur_start == 2)
  expect_equal(late_long$union_years, 0)
  expect_equal(late_long$dissolutions, 1)
  expect_true(is.na(late_long$rate_dissolution))
  # A missing end and a known future end both censor at the requested end.
  ongoing <- union_count_history()[3, ]
  future <- ongoing; future$end_time <- 10; future$end_reason <- "dissolution"
  expect_identical(calculate_union_counts(ongoing, c(0, 5), c(0, Inf)),
                   calculate_union_counts(future, c(0, 5), c(0, Inf)))
  # Beginning the window late uses attained duration, not zero duration.
  late <- calculate_union_counts(union_count_history(), c(2.5, 3.5), c(0, 1, 2, Inf))
  expect_equal(late$union_years, c(.5, 1, 1))
  expect_equal(late$dissolutions, c(0, 0, 1))
})

test_that("repeated and distinct overlapping union episodes each contribute once", {
  unions <- data.frame(union_id = 1:3, start_time = c(0, 1, 3),
    end_time = c(2, 2, 5), end_reason = "dissolution")
  counts <- calculate_union_counts(unions, c(0, 6), c(0, Inf))
  expect_equal(counts$union_years, 5)
  expect_equal(counts$dissolutions, 3)
})

test_that("grouped union counts follow women and preserve typed groups and totals", {
  population <- union_count_population()
  counts <- calculate_union_counts(population, c(0, 2, 5), c(0, 2, Inf), by_group = TRUE)
  expect_identical(unique(counts$group), c(20L, 10L))
  expect_type(counts$group, "integer")
  expect_equal(sum(counts$union_years[counts$group == 20L]), 6)
  expect_equal(sum(counts$union_years[counts$group == 10L]), 3)
  expect_equal(sum(counts$dissolutions[counts$group == 20L]), 1)
  expect_equal(sum(counts$dissolutions[counts$group == 10L]), 0)
  pooled <- calculate_union_counts(population, c(0, 2, 5), c(0, 2, Inf))
  expect_equal(unname(rowsum(counts$union_years, rep(1:4, 2))[, 1]), pooled$union_years)
  bare <- population$unions
  bare$group <- population$persons$group[match(bare$wife_id, population$persons$person_id)]
  expect_identical(calculate_union_counts(bare, c(0, 2, 5), c(0, 2, Inf), TRUE), counts)
  population$persons$group <- as.character(population$persons$group)
  expect_type(calculate_union_counts(population, c(0, 5), c(0, Inf), TRUE)$group, "character")
})

test_that("duration tails and empty grids follow Lexis conventions", {
  counts <- calculate_union_counts(union_count_history(), c(0, 5), c(0, 2))
  expect_equal(counts$dur_start, c(0, 2))
  expect_equal(counts$dur_end, c(2, Inf))
  empty <- union_count_history()[FALSE, ]
  counts <- calculate_union_counts(empty, c(0, 2, 5), c(0, 2))
  expect_s3_class(counts, "tbl_df")
  expect_equal(nrow(counts), 4)
  expect_equal(counts$union_years, rep(0, 4))
  expect_equal(counts$dissolutions, rep(0, 4))
  expect_true(all(is.na(counts$rate_dissolution)))
  empty$group <- integer()
  grouped <- calculate_union_counts(empty, c(0, 5), c(0, Inf), TRUE)
  expect_equal(nrow(grouped), 0)
  expect_type(grouped$group, "integer")
  expect_type(grouped$dissolutions, "double")
  # Population groups without union episodes still have explicit empty cells.
  population <- union_count_population()
  population$unions <- population$unions[FALSE, ]
  grouped <- calculate_union_counts(population, c(0, 5), c(0, Inf), TRUE)
  expect_equal(nrow(grouped), 2)
  expect_equal(grouped$union_years, c(0, 0))
})

test_that("union counts validate endpoints, breakpoints, histories and grouping", {
  population <- union_count_population()
  expect_error(calculate_union_counts(population, c(0, 6), c(0, Inf)), "endpoint")
  population["end_time"] <- list(NULL)
  expect_error(calculate_union_counts(population, c(0, 5), c(0, Inf)), "endpoint")
  unions <- union_count_history()
  for (breaks in list(c(0, 0), c(0, NA), c(0, Inf), 0, c(1, 0), matrix(c(0, 5))))
    expect_error(calculate_union_counts(unions, breaks, c(0, Inf)), "time_breaks")
  for (breaks in list(c(-1, 1), c(0, 0), c(0, NA), 0, c(1, 0), matrix(c(0, 5))))
    expect_error(calculate_union_counts(unions, c(0, 5), breaks), "dur_breaks")
  for (flag in list(NA, 1, "yes", c(TRUE, FALSE), matrix(TRUE)))
    expect_error(calculate_union_counts(unions, c(0, 5), c(0, Inf), flag), "by_group")
  expect_error(calculate_union_counts(unions, c(0, 5), c(0, Inf), TRUE), "group")
  expect_error(calculate_union_counts(rbind(unions, unions), c(0, 5), c(0, Inf)), "unique")
  expect_error(calculate_union_counts(unions[-1], c(0, 5), c(0, Inf)), "missing columns")
  bad <- unions; bad$end_time[1] <- -1
  expect_error(calculate_union_counts(bad, c(0, 5), c(0, Inf)), "precede")
  bad <- unions; bad$end_reason[1] <- NA_character_
  expect_error(calculate_union_counts(bad, c(0, 5), c(0, Inf)), "end reason")
  bad <- unions; bad$group <- c("A", NA, "B")
  expect_error(calculate_union_counts(bad, c(0, 5), c(0, Inf), TRUE), "group")
})
