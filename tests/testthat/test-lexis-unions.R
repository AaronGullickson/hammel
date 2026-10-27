make_lexis_union_history <- function() {
  persons <- data.frame(person_id = 1:8, sex = c("f", "m", "m", rep("f", 5)),
    birth_time = c(-20, -30, -40, 0.25, 0.5, 1.5, 2.5, 3),
    death_time = c(NA, 1.5, rep(NA_real_, 6)),
    mother_id = c(rep(NA_integer_, 3), rep(1L, 5)),
    father_id = NA_integer_)
  unions <- data.frame(union_id = 1:2, wife_id = 1L, husband_id = 2:3,
    start_time = c(0.5, 2.5), end_time = c(1.5, NA),
    end_reason = c("death_husband", NA))
  hammel_population(persons, unions, end_time = 3)
}

test_that("long union counts match exact widowhood and remarriage examples", {
  population <- make_lexis_union_history()
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, 20, 30, 40, Inf), TRUE)
  expect_identical(counts$union_status, rep(c("overall", "partnered", "unpartnered"), 8))
  women <- subset(counts, sex == "f" & age_start == 20)
  expect_equal(women$person_years, c(3, 1.5, 1.5))
  expect_equal(women$unions, c(NA, NA, 2))
  expect_equal(women$rate_marriage, c(NA, NA, 2/1.5))
  expect_equal(women$births, c(4, 2, 2))
  expect_equal(women$rate_fertility, c(4/3, 2/1.5, 2/1.5))
  men <- subset(counts, sex == "m" & age_start >= 30)
  expect_equal(men$person_years, c(1.5, 1, 0.5, 3, 0.5, 2.5))
  expect_equal(men$unions, c(NA, NA, 1, NA, NA, 1))
  expect_equal(men$rate_marriage, c(NA, NA, 2, NA, NA, 0.4))
  expect_equal(men$deaths, c(1, 1, 0, 0, 0, 0))
  expect_equal(men$rate_mortality, c(1/1.5, 1, 0, 0, 0, 0))
  expect_true(all(is.na(counts$births[counts$sex == "m"])))
  expect_true(all(is.na(counts$rate_fertility[counts$sex == "m"])))
  total <- subset(counts, union_status == "overall")
  partnered <- subset(counts, union_status == "partnered")
  unpartnered <- subset(counts, union_status == "unpartnered")
  for (column in c("person_years", "births", "deaths"))
    expect_equal(partnered[[column]] + unpartnered[[column]], total[[column]])
  expect_equal(sum(unpartnered$unions), 4)
  plain <- calculate_lexis_counts(population, c(0, 3), c(0, 20, 30, 40, Inf))
  expect_identical(total[names(plain)], plain)
  expect_identical(calculate_lexis_counts(population, c(0, 3),
    c(0, 20, 30, 40, Inf), FALSE), plain)
})

test_that("union exposure and events respect age and period boundaries", {
  population <- make_lexis_union_history()
  counts <- calculate_lexis_counts(population, c(0, 0.5, 1.5, 2.5, 3),
    c(0, 20, 20.5, 21.5, 22.5, 30, 40, Inf), TRUE)
  partnered <- subset(counts, union_status == "partnered")
  unpartnered <- subset(counts, union_status == "unpartnered")
  expect_equal(sum(partnered$person_years), 3)
  expect_equal(sum(unpartnered$unions), 4)
  starts <- subset(unpartnered, sex == "f" & unions > 0)
  expect_equal(starts$time_start, c(0.5, 2.5))
  expect_equal(starts$age_start, c(20.5, 22.5))
  expect_equal(subset(partnered, sex == "f" & births > 0)$time_start, c(0.5, 2.5))
  expect_equal(subset(unpartnered, sex == "f" & births > 0)$time_start, c(0, 1.5))
  death <- subset(partnered, deaths > 0)
  expect_equal(death$time_start, 1.5)
  expect_equal(death$deaths, 1)
  expect_true(is.na(death$rate_mortality)) # No partnered exposure after the boundary death.
  expect_equal(attr(counts, "unassigned_births")$births, rep(0, 4))
  window <- calculate_lexis_counts(population, c(1, 2.5), c(0, Inf), TRUE)
  expect_equal(sum(window$unions, na.rm = TRUE), 0)
  expect_equal(sum(subset(window, union_status == "partnered")$person_years), 1)
})

test_that("overlaps, touching episodes and zero durations do not duplicate time", {
  population <- make_lexis_union_history()
  population$unions <- rbind(population$unions,
    transform(population$unions[1, ], union_id = 3, start_time = 1, end_time = 1.25,
      end_reason = "dissolution"),
    transform(population$unions[2, ], union_id = 4, start_time = 2, end_time = 2.5,
      end_reason = "dissolution"),
    transform(population$unions[2, ], union_id = 5, start_time = 1.75, end_time = 1.75,
      end_reason = "dissolution"))
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  expect_equal(subset(counts, union_status == "partnered")$person_years, c(2, 2))
  expect_equal(subset(counts, union_status == "unpartnered")$unions, c(5, 5))
})

test_that("optional summaries require histories and handle empty union tables", {
  population <- make_lexis_union_history()
  expect_error(calculate_lexis_counts(population$persons, c(0, 3), c(0, Inf), TRUE),
    "union histories")
  for (bad in list(NULL, NA, 1, "yes", c(TRUE, FALSE), matrix(TRUE)))
    expect_error(calculate_lexis_counts(population, c(0, 3), c(0, Inf), bad),
      "include_unions")
  population$unions <- population$unions[FALSE, ]
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  total <- subset(counts, union_status == "overall")
  partnered <- subset(counts, union_status == "partnered")
  unpartnered <- subset(counts, union_status == "unpartnered")
  expect_equal(partnered$person_years, c(0, 0))
  expect_true(all(is.na(partnered$rate_mortality)))
  expect_true(all(is.na(partnered$rate_fertility)))
  expect_equal(unpartnered$person_years, total$person_years)
  expect_equal(unpartnered$births, total$births)
  expect_equal(unpartnered$deaths, total$deaths)
  expect_equal(unpartnered$rate_marriage, c(0, 0))
  population <- hammel_population(population$persons[FALSE, ], end_time = 3)
  empty <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 6)
  expect_equal(empty$person_years, rep(0, 6))
  expect_equal(empty$births, c(rep(NA_real_, 3), rep(0, 3)))
  expect_equal(empty$unions, c(NA, NA, 0, NA, NA, 0))
  expect_true(all(is.na(empty$rate_marriage)))
  expect_true(all(is.na(empty$rate_mortality)))
  expect_true(all(is.na(empty$rate_fertility)))
})

test_that("death status uses the left limit at union boundaries", {
  population <- make_lexis_union_history()
  # Widow dies after her first union, before the second; remove that second union.
  population$unions <- population$unions[1, ]
  population$persons$death_time[1] <- 2
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  expect_equal(subset(counts, sex == "f")$deaths, c(1, 0, 1))
  expect_equal(subset(counts, sex == "m")$deaths, c(1, 1, 0))
  # A zero-duration union at death creates no preceding partnered status.
  population$unions$start_time <- 1.5
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  expect_equal(subset(counts, sex == "m")$deaths, c(1, 0, 1))
})

test_that("existing TFR and life-table grouping accepts long Lexis results", {
  counts <- calculate_lexis_counts(make_lexis_union_history(), c(0, 3),
    c(0, 20, 30, 40, Inf), TRUE)
  tfr <- calculate_tfr(counts)
  expect_identical(tfr$union_status, c("overall", "partnered", "unpartnered"))
  # Unobserved finite age bands make each status-specific total incomplete.
  expect_identical(tfr$tfr, rep(NA_real_, 3))
  expect_s3_class(tfr, "tbl_df")
  expect_identical(attr(tfr, "unassigned_births"), attr(counts, "unassigned_births"))
  life <- calculate_life_table(counts)
  expect_s3_class(life, "tbl_df")
  expect_equal(nrow(life), nrow(counts))
  for (status in c("overall", "partnered", "unpartnered")) {
    plain <- calculate_life_table(subset(counts, union_status == status,
      select = -union_status))
    observed <- subset(life, union_status == status, select = -union_status)
    expect_equal(observed, plain)
  }
})
