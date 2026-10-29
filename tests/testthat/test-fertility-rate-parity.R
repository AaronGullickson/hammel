parity_observation_history <- function() {
  persons <- tibble::tibble(person_id = c(10, 20, 30, 40, 50),
    sex = c("f", "m", "f", "m", "f"),
    birth_time = c(-20, -20, -1, 1, 2),
    death_time = c(NA, NA, -.5, NA, NA),
    mother_id = c(NA, NA, 10, 10, 10), father_id = NA_real_,
    group = c("A", "B", "B", "B", "A"), parity = 99L)
  hammel_population(persons, end_time = 3)
}

test_that("separate rates contain only their process summaries", {
  pop <- parity_observation_history()
  m <- calculate_mortality_rates(pop, c(0,3), c(0,Inf))
  f <- calculate_fertility_rates(pop, c(0,3), c(0,Inf))
  expect_false(any(c("births", "rate_fertility") %in% names(m)))
  expect_false(any(c("deaths", "rate_mortality") %in% names(f)))
  expect_true(all(f$sex == "f"))
  expect_null(attr(m, "unassigned_births"))
  expect_s3_class(attr(f, "unassigned_births"), "tbl_df")
})

test_that("parity exposure uses pre-window births and pre-birth event parity", {
  pop <- parity_observation_history()
  x <- calculate_fertility_rates(pop, c(0,3), c(15,Inf), parity_breaks = 0:3)
  expect_equal(x$person_years, c(0,1,1,1))
  expect_equal(x$births, c(0,1,1,0))
  expect_equal(x$rate_fertility, c(NA,1,1,0))
  expect_equal(x$parity_end, c(1,2,3,Inf))
  pooled <- calculate_fertility_rates(pop, c(0,3), c(15,Inf))
  expect_equal(sum(x$person_years), pooled$person_years)
  expect_equal(sum(x$births), pooled$births)
  broad <- calculate_fertility_rates(pop, c(0,3), c(15,Inf), parity_breaks = c(0,2,Inf))
  expect_equal(broad$person_years, c(1,2))
  expect_equal(broad$births, c(1,1))
})

test_that("parity pools across group, time, age and union type exactly", {
  pop <- parity_observation_history()
  pop$unions <- tibble::tibble(union_id=1:2, wife_id=10, husband_id=20,
    start_time=c(0,1), end_time=c(1,2),
    end_reason=c("conversion","dissolution"), union_type=c("cohab","marriage"),
    previous_cohab_id=c(NA,1))
  # Both births are birth-first ties: cohab at conversion, marriage at dissolution.
  x <- calculate_fertility_rates(pop,c(0,1,2,3),c(0,20,22,Inf),
    by_group=TRUE,by_union_type=TRUE,parity_breaks=0:3)
  pooled <- calculate_fertility_rates(pop,c(0,1,2,3),c(0,20,22,Inf),
    by_group=TRUE,by_union_type=TRUE)
  cols <- c("time_start","time_end","age_start","age_end","group","union_type")
  observed <- dplyr::summarise(dplyr::group_by(x,dplyr::across(dplyr::all_of(cols))),
    person_years=sum(person_years),births=sum(births),.groups="drop")
  expected <- dplyr::select(pooled,dplyr::all_of(c(cols,"person_years","births")))
  observed <- dplyr::arrange(observed,!!!rlang::syms(cols))
  expected <- dplyr::arrange(expected,!!!rlang::syms(cols))
  attr(expected,"unassigned_births") <- NULL
  expect_equal(observed,expected)
  expect_equal(sum(x$births[x$union_type=="cohab"]),1)
  expect_equal(sum(x$births[x$union_type=="marriage"]),1)
})

test_that("parity boundaries and empty histories are validated", {
  pop <- parity_observation_history()
  for (bad in list(c(1,2),c(0,.5,Inf),c(0,0,Inf),c(0,NA,Inf),c(-1,0,Inf)))
    expect_error(calculate_fertility_rates(pop,c(0,3),c(0,Inf),parity_breaks=bad),"parity_breaks")
  empty <- hammel_population(pop$persons[FALSE,],end_time=3)
  x <- calculate_fertility_rates(empty,c(0,3),c(0,Inf),by_group=TRUE,parity_breaks=c(0,Inf))
  expect_s3_class(x,"tbl_df")
  expect_equal(nrow(x),0L)
  expect_type(x$rate_fertility,"double")
})

test_that("simultaneous births progress parity and final endpoints exclude events", {
  pop <- parity_observation_history()
  pop$persons$birth_time[4:5] <- 1
  x <- calculate_fertility_rates(pop,c(0,1,3),c(15,Inf),parity_breaks=0:3)
  # Two births at time 1 count in the later period at parities one and two;
  # parity two has events but no exposure, so its rate is missing.
  expect_equal(x$births,c(0,0,0,1,0,1,0,0))
  expect_equal(x$person_years,c(0,0,1,0,0,0,0,2))
  expect_true(is.na(x$rate_fertility[6]))
  y <- calculate_fertility_rates(pop,c(0,1),c(15,Inf),parity_breaks=0:3)
  expect_equal(sum(y$births),0)
  expect_equal(sum(y$person_years),1)
})
