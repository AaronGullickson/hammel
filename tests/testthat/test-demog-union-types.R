typed_demog_history <- function() {
  persons <- tibble::tibble(person_id=1:5,sex=c("f","m","f","f","f"),
    birth_time=c(-20,-20,1,3,4),death_time=rep(NA_real_,5),
    mother_id=c(NA_integer_,NA_integer_,1L,1L,1L),father_id=NA_integer_,
    group=c("A","B","B","B","B"))
  unions <- tibble::tibble(union_id=1:2,wife_id=1L,husband_id=2L,
    start_time=c(1,3),end_time=c(3,4),end_reason=c("conversion","dissolution"),
    union_type=c("cohab","marriage"),previous_cohab_id=c(NA_integer_,1L))
  hammel_population(persons,unions,end_time=5)
}

test_that("typed demographic exposure and births partition totals", {
  pop <- typed_demog_history()
  for (grouped in c(FALSE,TRUE)) {
    x <- collect_demographic_counts(pop,c(0,5),c(0,20,Inf),
      by_union_type=TRUE,by_group=grouped)
    expect_identical(unique(x$union_type),c("unpartnered","cohab","marriage"))
    total <- collect_demographic_counts(pop,c(0,5),c(0,20,Inf),by_group=grouped)
    parts <- lapply(unique(x$union_type),function(type) subset(x,union_type==type))
    for (column in c("person_years","births","deaths"))
      expect_equal(Reduce(`+`,lapply(parts,`[[`,column)),total[[column]])
    mother <- subset(x,sex == "f" & age_start == 20)
    if (grouped) mother <- subset(mother,group == "A")
    expect_equal(mother$person_years,c(2,2,1))
    expect_equal(mother$births,c(0,2,1))
    expect_equal(mother$rate_fertility,c(0,1,1))
  }
  x <- collect_demographic_counts(pop,c(0,1,3,4,5),c(0,20,23,24,Inf),TRUE)
  tied <- subset(x,sex == "f" & births > 0)
  expect_identical(tied$union_type,c("cohab","cohab","marriage"))
  expect_equal(tied$time_start,c(1,3,4))
  expect_equal(tied$births,c(1,1,1))
  expect_equal(subset(tied,time_start==3)$person_years,0)
  expect_equal(subset(tied,time_start==4)$person_years,0)
})

test_that("typed death strata use the left limit", {
  pop <- typed_demog_history()
  pop$persons <- pop$persons[1:2, ]
  pop$ancestry <- pop$ancestry[pop$ancestry$person_id %in% 1:2, ]
  pop$persons$death_time[2] <- 4
  pop$unions$end_reason[2] <- "death_husband"
  x <- collect_demographic_counts(pop,c(0,4,5),c(0,Inf),TRUE)
  death <- subset(x,deaths>0)
  expect_equal(death$union_type,"marriage")
  expect_equal(death$time_start,4)
  expect_equal(death$person_years,0)
  expect_equal(sum(x$deaths),1)
})

test_that("type discovery handles single types and empty typed histories", {
  pop <- typed_demog_history()
  pop$unions <- pop$unions[1, ]
  pop$unions$end_reason <- "dissolution"
  x <- collect_demographic_counts(pop,c(0,5),c(0,Inf),TRUE)
  expect_identical(unique(x$union_type),c("unpartnered","cohab"))
  pop$unions <- pop$unions[0, ]
  x <- collect_demographic_counts(pop,c(0,5),c(0,Inf),TRUE)
  expect_identical(unique(x$union_type),"unpartnered")
  expect_equal(x$person_years,collect_demographic_counts(pop,c(0,5),c(0,Inf))$person_years)
  pop$persons <- pop$persons[0, ]; pop$ancestry <- pop$ancestry[0, ]
  x <- collect_demographic_counts(pop,c(0,5),c(0,Inf),TRUE,TRUE)
  expect_equal(nrow(x),0)
  expect_type(x$union_type,"character")
})

test_that("typed overlapping histories are rejected", {
  pop <- typed_demog_history()
  pop$unions$start_time[2] <- 2
  pop$unions$end_reason[1] <- "dissolution"
  pop$unions$previous_cohab_id[2] <- NA_integer_
  expect_error(collect_demographic_counts(pop,c(0,5),c(0,Inf),TRUE),"nonoverlapping")
})
