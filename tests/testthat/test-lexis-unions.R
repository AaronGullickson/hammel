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

test_that("union counts match exact widowhood and remarriage examples", {
  population <- make_lexis_union_history()
  counts <- calculate_lexis_counts(population, c(0, 3), c(0, 20, 30, 40, Inf),
    include_unions = TRUE)
  women <- subset(counts, sex == "f" & age_start == 20)
  expect_equal(women$person_years, 3)
  expect_equal(women$person_years_partnered, 1.5)
  expect_equal(women$person_years_unpartnered, 1.5)
  expect_equal(women$unions, 2)
  expect_equal(women$births, 4)
  expect_equal(women$births_partnered, 2)
  expect_equal(women$births_unpartnered, 2)
  men <- subset(counts, sex == "m" & age_start >= 30)
  expect_equal(men$person_years_partnered, c(1, 0.5))
  expect_equal(men$person_years_unpartnered, c(0.5, 2.5))
  expect_equal(men$unions, c(1, 1))
  expect_true(all(is.na(counts$births_partnered[counts$sex == "m"])))
  expect_true(all(is.na(counts$births_unpartnered[counts$sex == "m"])))
  expect_equal(counts$person_years_partnered + counts$person_years_unpartnered,
    counts$person_years)
  expect_equal(counts$births_partnered + counts$births_unpartnered, counts$births)
  expect_equal(sum(counts$unions), 4)
  plain <- calculate_lexis_counts(population, c(0, 3), c(0, 20, 30, 40, Inf))
  expect_identical(counts[names(plain)], plain)
})

test_that("union exposure and events respect age and period boundaries", {
  population <- make_lexis_union_history()
  counts <- calculate_lexis_counts(population, c(0, 0.5, 1.5, 2.5, 3),
    c(0, 20, 20.5, 21.5, 22.5, 30, 40, Inf), include_unions = TRUE)
  expect_equal(sum(counts$person_years_partnered), 3)
  expect_equal(sum(counts$unions), 4)
  starts <- subset(counts, sex == "f" & unions > 0)
  expect_equal(starts$time_start, c(0.5, 2.5))
  expect_equal(starts$age_start, c(20.5, 22.5))
  births <- subset(counts, sex == "f" & births > 0)
  expect_equal(births$births_partnered, c(0, 1, 0, 1))
  expect_equal(births$births_unpartnered, c(1, 0, 1, 0))
  expect_equal(attr(counts, "unassigned_births")$births, rep(0, 4))
  window <- calculate_lexis_counts(population, c(1, 2.5), c(0, Inf), TRUE)
  expect_equal(sum(window$unions), 0)
  expect_equal(sum(window$person_years_partnered), 1)
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
  expect_equal(counts$person_years_partnered, c(2, 2))
  expect_equal(counts$unions, c(5, 5))
  expect_equal(counts$person_years_partnered + counts$person_years_unpartnered,
    counts$person_years)
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
  expect_equal(counts$person_years_partnered, c(0, 0))
  expect_equal(counts$person_years_unpartnered, counts$person_years)
  expect_equal(counts$births_unpartnered, counts$births)
  expect_equal(counts$unions, c(0, 0))
  population <- hammel_population(population$persons[FALSE, ], end_time = 3)
  empty <- calculate_lexis_counts(population, c(0, 3), c(0, Inf), TRUE)
  expect_s3_class(empty, "tbl_df")
  expect_equal(empty$person_years_partnered, c(0, 0))
  expect_equal(empty$births_partnered, c(NA, 0))
})
