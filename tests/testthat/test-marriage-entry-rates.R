marriage_entry_history <- function() {
  persons <- tibble::tibble(person_id = 1:2, sex = c("f", "m"),
    birth_time = -20, death_time = NA_real_, mother_id = NA_integer_,
    father_id = NA_integer_, group = c("A", "B"))
  unions <- tibble::tibble(union_id = 1:2, wife_id = 1L, husband_id = 2L,
    start_time = c(1, 3), end_time = c(3, NA),
    end_reason = c("conversion", NA), union_type = c("cohab", "marriage"),
    previous_cohab_id = c(NA_integer_, 1L))
  hammel_population(persons, unions, end_time = 5)
}

test_that("marriage entry includes conversion and nonmarried exposure", {
  pop <- marriage_entry_history()
  marriage <- calculate_marriage_rates(pop, c(0,5), c(0,Inf))
  union <- calculate_union_rates(pop, c(0,5), c(0,Inf))
  expect_s3_class(marriage, "tbl_df")
  expect_equal(marriage$person_years, c(3,3))
  expect_equal(marriage$marriages, c(1,1))
  expect_equal(marriage$rate_marriage, c(1/3,1/3))
  expect_equal(union$person_years, c(1,1))
  expect_equal(union$unions, c(1,1))
  grouped <- calculate_marriage_rates(pop, c(0,5), c(0,Inf), by = "group")
  expect_equal(sum(grouped$marriages), 2)
  expect_equal(sum(grouped$person_years), 6)
  expect_equal(grouped$group, rep(c("A", "B"), each = 2))
  expect_true(all(is.na(grouped$rate_marriage[grouped$person_years == 0])))
})

test_that("marriage conversion boundaries and age cells are half open", {
  pop <- marriage_entry_history()
  x <- calculate_marriage_rates(pop, 0:5, c(0,21,23,Inf))
  events <- subset(x, marriages > 0)
  expect_equal(events$time_start, c(3,3))
  expect_equal(events$age_start, c(23,23))
  expect_equal(sum(x$person_years), 6)
  expect_equal(sum(calculate_marriage_rates(pop,c(0,3),c(0,Inf))$marriages), 0)
  after <- calculate_marriage_rates(pop,c(3,5),c(0,Inf))
  expect_equal(after$marriages,c(1,1))
  expect_equal(after$person_years,c(0,0))
  expect_true(all(is.na(after$rate_marriage)))
  cohab <- pop
  cohab$unions <- cohab$unions[1, ]
  cohab$unions$end_time <- NA_real_
  cohab$unions$end_reason <- NA_character_
  x <- calculate_marriage_rates(cohab,c(0,5),c(0,Inf))
  expect_equal(x$person_years,c(5,5))
  expect_equal(x$marriages,c(0,0))
})

test_that("direct and legacy marriages retain existing exposure accounting", {
  pop <- marriage_entry_history()
  pop$unions <- pop$unions[2, ]
  pop$unions$previous_cohab_id <- NA_integer_
  direct <- calculate_marriage_rates(pop,c(0,5),c(0,Inf))
  expect_equal(direct$person_years,c(3,3))
  expect_equal(direct$marriages,c(1,1))
  pop$unions$union_type <- NULL
  pop$unions$previous_cohab_id <- NULL
  expect_identical(calculate_marriage_rates(pop,c(0,5),c(0,Inf)),direct)
  pop$unions <- pop$unions[0, ]
  pop$persons$death_time <- c(2,NA)
  x <- calculate_marriage_rates(pop,c(0,5),c(0,Inf))
  expect_equal(x$person_years,c(5,2))
  expect_equal(x$marriages,c(0,0))
  expect_error(calculate_marriage_rates(pop,c(0,6),c(0,Inf)),"endpoint")
  expect_error(calculate_marriage_rates(pop,c(0,5),c(0,Inf), by = NA), "`by` must")
  pop$persons <- pop$persons[0, ]
  pop$ancestry <- pop$ancestry[0, ]
  expect_s3_class(calculate_marriage_rates(pop,c(0,5),c(0,Inf)),"tbl_df")
})


test_that("cohabitation endings preserve nonmarried exposure", {
  for (reason in c("dissolution", "death_husband")) {
    pop <- marriage_entry_history()
    pop$unions <- pop$unions[1, ]
    pop$unions$end_reason <- reason
    if (reason == "death_husband") pop$persons$death_time[2] <- 3
    counts <- calculate_marriage_rates(pop, c(0,5), c(0,Inf))
    expect_equal(counts$person_years,
      if (reason == "death_husband") c(3,5) else c(5,5))
    expect_equal(counts$marriages,c(0,0))
  }
})

test_that("marriage previous categories include cohabiting exposure and conversion", {
  pop <- marriage_entry_history()
  # A prior legal marriage ended before the observation window.
  prior <- pop$unions[2, ]
  prior$union_id <- 3L
  prior$start_time <- -3
  prior$end_time <- -1
  prior$end_reason <- "dissolution"
  prior$previous_cohab_id <- NA_integer_
  pop$unions <- rbind(prior, pop$unions)
  for (grouped in c(FALSE, TRUE)) {
    pooled <- calculate_marriage_rates(pop,c(0,5),c(0,Inf), by = if (isTRUE(grouped)) "group")
    split <- calculate_marriage_rates(pop,c(0,5),c(0,Inf), by = c(if (isTRUE(grouped)) "group", "previous_union"))
    expect_equal(sum(split$person_years),sum(pooled$person_years))
    expect_equal(sum(split$marriages),sum(pooled$marriages))
    expect_equal(sum(split$person_years[split$previous_union == "separated"]),6)
    expect_equal(sum(split$marriages[split$previous_union == "separated"]),2)
    expect_equal(sum(split$person_years[split$previous_union != "separated"]),0)
  }
  expect_error(calculate_marriage_rates(pop,c(0,5),c(0,Inf), by = NA), "`by` must")
})
