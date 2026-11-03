marriage_history <- function() {
  persons <- data.frame(person_id = 1:4, sex = c("f", "m", "m", "m"),
    birth_time = -20, death_time = c(NA, NA, 4, NA),
    mother_id = NA_integer_, father_id = NA_integer_)
  unions <- data.frame(union_id = 1:3, wife_id = 1, husband_id = 2:4,
    start_time = c(1, 3, 5), end_time = c(2, 4, NA),
    end_reason = c("dissolution", "death_husband", NA))
  hammel_population(persons, unions, end_time = 6)
}

test_that("marriage counts assign exact exposure and starts by latest episode", {
  pop <- marriage_history()
  counts <- calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE)
  expect_s3_class(counts, "tbl_df")
  women <- subset(counts, sex == "f")
  expect_identical(women$previous_union, c("none", "separated", "widowed"))
  expect_equal(women$person_years, c(1, 1, 1))
  expect_equal(women$unions, c(1, 1, 1))
  expect_equal(women$rate_union, c(1, 1, 1))
  men <- subset(counts, sex == "m")
  expect_equal(men$person_years, c(9, 4, 0))
  expect_equal(men$unions, c(3, 0, 0))
  expect_equal(men$rate_union, c(1/3, 0, NA))
  # Pooling categories reproduces unpartnered exposure by sex.
  pooled <- collect_demographic_counts(pop, c(0, 6), c(0, Inf), by_union_type = TRUE)
  for (sex_value in c("m", "f"))
    expect_equal(sum(counts$person_years[counts$sex == sex_value]),
      subset(pooled, union_type == "unpartnered" & sex == sex_value)$person_years)
  expect_equal(sum(counts$unions), 2 * nrow(pop$unions))
  shuffled <- pop; shuffled$unions <- shuffled$unions[3:1, ]
  expect_identical(calculate_union_rates(shuffled, c(0, 6), c(0, Inf), by_previous_union = TRUE), counts)
})

test_that("marriage cells are half-open and use preformation categories", {
  pop <- marriage_history()
  counts <- calculate_union_rates(pop, 0:6, c(0, 21, 23, 25, Inf), by_previous_union = TRUE)
  starts <- subset(counts, sex == "f" & unions > 0)
  expect_equal(starts$time_start, c(1, 3, 5))
  expect_equal(starts$age_start, c(21, 23, 25))
  expect_equal(starts$previous_union, c("none", "separated", "widowed"))
  # Ending at the observation endpoint excludes the formation there.
  before <- calculate_union_rates(pop, c(2, 5), c(0, Inf), by_previous_union = TRUE)
  women <- subset(before, sex == "f")
  expect_equal(women$person_years, c(0, 1, 1))
  expect_equal(women$unions, c(0, 1, 0))
  # History before the window remains relevant to classification.
  after <- calculate_union_rates(pop, c(4, 6), c(0, Inf), by_previous_union = TRUE)
  expect_equal(subset(after, sex == "f")$unions, c(0, 0, 1))
})

test_that("touching and zero-duration unions change category without extra exposure", {
  pop <- marriage_history()
  pop$unions$start_time <- c(1, 2, 2)
  pop$unions$end_time <- c(2, 2, NA)
  pop$unions$end_reason[2] <- "dissolution"
  counts <- calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE)
  women <- subset(counts, sex == "f")
  expect_equal(women$person_years, c(1, 0, 0))
  expect_equal(women$unions, c(1, 2, 0))
  expect_true(is.na(women$rate_union[2]))
  pop$unions$start_time[2] <- 1.5
  expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE), "nonoverlapping")
})

test_that("male widowhood and repeated changes follow the latest ending", {
  pop <- marriage_history()
  pop$persons$sex <- c("m", "f", "f", "f")
  pop$unions$husband_id <- 1
  pop$unions$wife_id <- 2:4
  pop$unions$end_reason[2] <- "death_wife"
  counts <- calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE)
  men <- subset(counts, sex == "m")
  expect_equal(men$unions, c(1, 1, 1))
  expect_equal(men$person_years, c(1, 1, 1))
})

test_that("grouped marriage counts pool correctly and include empty categories", {
  pop <- marriage_history()
  pop$persons$group <- c(10, 20, 20, 10)
  counts <- calculate_union_rates(pop, c(0, 6), c(0, Inf), by_group = TRUE, by_previous_union = TRUE)
  expect_type(counts$group, "double")
  expect_equal(nrow(counts), 12)
  pooled <- aggregate(counts[c("person_years", "unions")],
    counts[c("sex", "previous_union")], sum)
  plain <- calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE)
  rows <- match(paste(plain$sex, plain$previous_union), paste(pooled$sex, pooled$previous_union))
  expect_equal(pooled$person_years[rows], plain$person_years)
  expect_equal(pooled$unions[rows], plain$unions)
  empty <- hammel_population(pop$persons[FALSE, ], end_time = 6)
  grouped <- calculate_union_rates(empty, c(0, 6), c(0, Inf), TRUE, by_previous_union = TRUE)
  expect_s3_class(grouped, "tbl_df")
  expect_equal(nrow(grouped), 0)
  expect_type(grouped$previous_union, "character")
  expect_type(grouped$rate_union, "double")
  zeros <- calculate_union_rates(empty, c(0, 6), c(0, Inf), by_previous_union = TRUE)
  expect_equal(zeros$person_years, rep(0, 6))
  expect_true(all(is.na(zeros$rate_union)))
})

test_that("marriage observation rejects invalid windows and ambiguous histories", {
  pop <- marriage_history()
  expect_error(calculate_union_rates(pop$persons, c(0, 6), c(0, Inf), by_previous_union = TRUE), "union histories")
  expect_error(calculate_union_rates(pop, c(0, 7), c(0, Inf), by_previous_union = TRUE), "endpoint")
  expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf), TRUE, by_previous_union = TRUE), "group")
  expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf), NA, by_previous_union = TRUE), "by_group")
  expect_error(calculate_union_rates(pop, c(0, 0), c(0, Inf), by_previous_union = TRUE), "time_breaks")
  pop$persons$death_time[2] <- 1.5
  expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf), by_previous_union = TRUE), "after a partner's death")
})

test_that("previous union is optional and pooled rates use summed exposure", {
  pop <- marriage_history()
  for (by_group in c(FALSE, TRUE)) {
    pop$persons$group <- c("A", "B", "A", "B")
    pooled <- calculate_union_rates(pop, c(0, 3, 6), c(0, 23, Inf), by_group)
    strata <- calculate_union_rates(pop, c(0, 3, 6), c(0, 23, Inf),
      by_group, by_previous_union = TRUE)
    expect_false("previous_union" %in% names(pooled))
    keys <- c("time_start", "time_end", "sex", "age_start", "age_end", if (by_group) "group")
    for (i in seq_len(nrow(pooled))) {
      rows <- Reduce(`&`, lapply(keys, function(k) strata[[k]] == pooled[[k]][i]))
      expect_equal(pooled$person_years[i], sum(strata$person_years[rows]))
      expect_equal(pooled$unions[i], sum(strata$unions[rows]))
      expected <- if (pooled$person_years[i] > 0)
        pooled$unions[i] / pooled$person_years[i] else NA_real_
      expect_equal(pooled$rate_union[i], expected)
    }
  }
  empty <- hammel_population(pop$persons[FALSE, ], end_time = 6)
  expect_equal(nrow(calculate_union_rates(empty, c(0, 6), c(0, Inf), TRUE)), 0)
  expect_equal(nrow(calculate_union_rates(empty, c(0, 6), c(0, Inf))), 2)
  for (bad in list(NULL, NA, 1, "yes", c(TRUE, FALSE), matrix(TRUE)))
    expect_error(calculate_union_rates(pop, c(0, 6), c(0, Inf),
      by_previous_union = bad), "by_previous_union")
})
