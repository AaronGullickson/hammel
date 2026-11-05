conversion_rate_history <- function() {
  persons <- tibble::tibble(person_id=1:2,sex=c("f","m"),birth_time=c(-20,-25),
    death_time=NA_real_,mother_id=NA_integer_,father_id=NA_integer_,group=c(20L,10L))
  unions <- tibble::tibble(union_id=1:3,wife_id=1L,husband_id=2L,
    start_time=c(0,2,3),end_time=c(2,3,6),
    end_reason=c("conversion","dissolution","dissolution"),
    union_type=c("cohab","marriage","cohab"),previous_cohab_id=c(NA_integer_,1L,NA_integer_))
  hammel_population(persons,unions,end_time=8)
}

test_that("conversion rates use cohabitation exposure and duration", {
  pop <- conversion_rate_history()
  x <- calculate_cohab_convert_rate(pop,c(0,8),dur_breaks=c(0,1,2))
  expect_s3_class(x,"tbl_df")
  expect_named(x,c("time_start","time_end","dur_start","dur_end","union_years","conversions","rate_cohab_convert"))
  expect_equal(x$dur_end,c(1,2,Inf))
  expect_equal(x$union_years,c(2,2,1))
  expect_equal(x$conversions,c(0,0,1))
  expect_equal(x$rate_cohab_convert,c(0,0,1))
  pooled <- calculate_cohab_convert_rate(pop,c(0,8),dur_breaks=c(0,Inf))
  expect_equal(pooled$union_years,5)
  expect_equal(pooled$conversions,1)
  expect_equal(pooled$rate_cohab_convert,.2)
  grouped <- calculate_cohab_convert_rate(pop,c(0,8),dur_breaks=c(0,Inf),by = "group_w")
  expect_identical(grouped$group_w,c(20L,10L))
  expect_equal(grouped$union_years,c(5,0))
  expect_equal(grouped$conversions,c(1,0))
  expect_equal(grouped$rate_cohab_convert,c(.2,NA))
})

test_that("conversion age rates follow woman age and conserve totals", {
  pop <- conversion_rate_history()
  x <- calculate_cohab_convert_rate(pop,c(0,8),age_breaks=c(20,22,25,Inf))
  expect_equal(x$union_years,c(2,2,1))
  expect_equal(x$conversions,c(0,1,0))
  expect_equal(x$rate_cohab_convert,c(0,.5,0))
  expect_false(any(c("dur_start","dur_end") %in% names(x)))
  expect_equal(sum(x$union_years),5)
  for (grouped in c(FALSE,TRUE)) {
    age <- calculate_cohab_convert_rate(pop,c(0,2,4,8),age_breaks=c(0,22,24,Inf), by = if (isTRUE(grouped)) "group_w")
    dur <- calculate_cohab_convert_rate(pop,c(0,2,4,8),dur_breaks=c(0,1,2,Inf), by = if (isTRUE(grouped)) "group_w")
    expect_equal(sum(age$union_years),sum(dur$union_years))
    expect_equal(sum(age$conversions),sum(dur$conversions))
  }
})

test_that("conversion observation censors and assigns boundary events", {
  pop <- conversion_rate_history()
  before <- calculate_cohab_convert_rate(pop,c(0,2),dur_breaks=c(0,Inf))
  expect_equal(before$union_years,2)
  expect_equal(before$conversions,0)
  at <- calculate_cohab_convert_rate(pop,c(2,3),dur_breaks=c(0,2,Inf))
  expect_equal(at$union_years,c(0,0))
  expect_equal(at$conversions,c(0,1))
  expect_true(all(is.na(at$rate_cohab_convert)))
  pop$unions$end_time[3] <- NA_real_; pop$unions$end_reason[3] <- NA_character_
  expect_equal(calculate_cohab_convert_rate(pop,c(0,8),dur_breaks=c(0,Inf))$union_years,7)
  pop$unions$end_time[3] <- 5; pop$unions$end_reason[3] <- "death_husband"
  pop$persons$death_time[2] <- 5
  expect_equal(calculate_cohab_convert_rate(pop,c(0,8),dur_breaks=c(0,Inf))$union_years,4)
})

test_that("conversion observation validates axes, histories, flags and empty data", {
  pop <- conversion_rate_history()
  expect_error(calculate_cohab_convert_rate(pop,c(0,8)),"exactly one")
  expect_error(calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,Inf), age_breaks = c(0,Inf)),"exactly one")
  expect_error(calculate_cohab_convert_rate(pop,c(0,9), dur_breaks = c(0,Inf)),"endpoint")
  expect_error(calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,1,1)),"dur_breaks")
  expect_error(calculate_cohab_convert_rate(pop,c(0,8),age_breaks=c(-1,1)),"age_breaks")
  expect_error(calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,Inf), by = NA),"`by` must")
  expect_error(calculate_cohab_convert_rate(pop$unions,c(0,8), dur_breaks = c(0,Inf)),"hammel_population")
  pop$unions <- pop$unions[0, ]
  x <- calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,Inf))
  expect_equal(x$union_years,0)
  expect_equal(x$conversions,0)
  expect_true(is.na(x$rate_cohab_convert))
  pop$unions$union_type <- NULL
  expect_error(calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,Inf)),"typed")
  pop$unions$union_type <- character()
  pop$persons <- pop$persons[0, ]; pop$ancestry <- pop$ancestry[0, ]
  x <- calculate_cohab_convert_rate(pop,c(0,8), dur_breaks = c(0,Inf), by = "group_w")
  expect_s3_class(x,"tbl_df")
  expect_equal(nrow(x),0)
  expect_type(x$group_w,"integer")
})
