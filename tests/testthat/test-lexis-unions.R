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


test_that("union type strata conserve demographic counts without overall rows", {
  pop <- make_lexis_union_history()
  x <- demographic_counts(pop,c(0,3),c(0,20,30,40,Inf), by = "union_type")
  # union_type varies slowest: rows follow time, age and sex within each type.
  expect_identical(x$union_type,rep(c("unpartnered","partnered"),each = 8))
  pooled <- demographic_counts(pop,c(0,3),c(0,20,30,40,Inf))
  unpartnered <- subset(x,union_type == "unpartnered")
  partnered <- subset(x,union_type == "partnered")
  for (column in c("person_years","births","deaths"))
    expect_equal(unpartnered[[column]] + partnered[[column]],pooled[[column]])
  women <- subset(x,sex == "f" & age_start == 20)
  expect_equal(women$person_years,c(1.5,1.5))
  expect_equal(women$births,c(2,2))
  expect_equal(women$rate_fertility,c(2/1.5,2/1.5))
  expect_true(all(is.na(x$births[x$sex == "m"])))
  expect_identical(attr(x,"unassigned_births"),attr(pooled,"unassigned_births"))
})

test_that("untyped boundaries and overlap merging retain event semantics", {
  pop <- make_lexis_union_history()
  counts <- demographic_counts(pop,c(0,.5,1.5,2.5,3),
    c(0,20,20.5,21.5,22.5,30,40,Inf), by = "union_type")
  partnered <- subset(counts,union_type == "partnered")
  unpartnered <- subset(counts,union_type == "unpartnered")
  expect_equal(sum(partnered$person_years),3)
  expect_equal(subset(partnered,sex == "f" & births > 0)$time_start,c(.5,2.5))
  expect_equal(subset(unpartnered,sex == "f" & births > 0)$time_start,c(0,1.5))
  expect_equal(subset(partnered,deaths > 0)$time_start,1.5)
  pop$unions$end_reason[1] <- "dissolution"
  counts <- demographic_counts(pop,c(1.5,2.5),c(0,Inf), by = "union_type")
  x <- subset(counts,sex == "f" & union_type == "partnered")
  expect_equal(x$births,1)
  expect_equal(x$person_years,0)
  pop$unions <- rbind(pop$unions,transform(pop$unions[1, ],union_id=3,husband_id=3,start_time=1))
  duplicate <- demographic_counts(pop,c(1.5,2.5),c(0,Inf), by = "union_type")
  expect_identical(duplicate$births,counts$births)
})

test_that("union type summaries validate flags and empty histories", {
  pop <- make_lexis_union_history()
  expect_error(demographic_counts(pop$persons,c(0,3),c(0,Inf), by = "union_type"),"union histories")
  for (bad in list(NA,1,"yes",TRUE,c("union_type","union_type"),matrix("union_type")))
    expect_error(demographic_counts(pop,c(0,3),c(0,Inf), by = bad),"by")
  pop$unions <- pop$unions[0, ]
  x <- demographic_counts(pop,c(0,3),c(0,Inf), by = "union_type")
  expect_equal(subset(x,union_type == "partnered")$person_years,c(0,0))
  expect_true(all(is.na(subset(x,union_type == "partnered")$rate_fertility)))
  pop <- hammel_population(pop$persons[0, ],end_time=3)
  x <- demographic_counts(pop,c(0,3),c(0,Inf), by = "union_type")
  expect_s3_class(x,"tbl_df")
  expect_equal(nrow(x),4)
  expect_equal(x$person_years,rep(0,4))
  expect_equal(x$births,c(NA,0,NA,0))
})

test_that("derived measures group by observed union type", {
  x <- demographic_counts(make_lexis_union_history(),c(0,3),c(0,20,30,40,Inf), by = "union_type")
  tfr <- calculate_tfr(x)
  expect_identical(tfr$union_type,c("unpartnered","partnered"))
  expect_identical(tfr$tfr,rep(NA_real_,2))
  life <- calculate_life_table(x)
  expect_equal(nrow(life),nrow(x))
  for (type in unique(x$union_type)) {
    plain <- calculate_life_table(subset(x,union_type == type,select=-union_type))
    expect_equal(subset(life,union_type == type,select=-union_type),plain)
  }
})
