unpartnered_birth_history <- function() {
  times <- c(0,1,2,2.5,3,4,4.5,5,6,7,7.5,8)
  persons <- tibble::tibble(person_id=1:15,sex=c("f","m","m",rep("f",12)),
    birth_time=c(-30,-30,-30,times),death_time=c(NA,4,rep(NA_real_,13)),
    mother_id=c(rep(NA_integer_,3),rep(1L,12)),
    father_id=c(rep(NA_integer_,3),ifelse(times %in% c(1,2,3),2L,3L)),
    group=rep(c(20L,10L,30L),5),custom=letters[1:15])
  unions <- tibble::tibble(union_id=1:4,wife_id=1L,husband_id=c(2L,2L,3L,3L),
    start_time=c(1,3,5,6),end_time=c(2,4,6,7),
    end_reason=c("dissolution","death_husband","conversion","dissolution"),
    union_type=c("marriage","cohab","cohab","marriage"),
    previous_cohab_id=c(NA_integer_,NA_integer_,NA_integer_,3L))
  hammel_population(persons,unions,end_time=9, groups = if ("group" %in% names(persons)) "group")
}

test_that("unpartnered births exclude all partnership types with exact birth ties", {
  pop <- unpartnered_birth_history()
  before <- pop
  x <- extract_unpartnered_births(pop)
  expect_s3_class(x,"tbl_df")
  expect_equal(x$birth_time,c(0,2.5,4,4.5,7.5,8))
  expect_equal(x$mother_id,rep(1L,6))
  expect_equal(x$father_id,rep(3L,6)) # Repeated pair retained per birth.
  expected <- augment_persons(pop)
  expect_identical(x,expected[match(x$person_id,expected$person_id), ])
  expect_type(x$group_mother,"integer")
  expect_identical(pop,before)
  pop["end_time"] <- list(NULL)
  expect_error(extract_unpartnered_births(pop), "end_time")
})

test_that("missing parents, overlapping histories and untyped histories are handled", {
  pop <- unpartnered_birth_history()
  pop$persons$father_id[4] <- NA_integer_
  pop$persons$mother_id[7] <- NA_integer_
  expect_equal(extract_unpartnered_births(pop)$birth_time,c(4,4.5,7.5,8))
  pop <- unpartnered_birth_history()
  duplicate <- pop$unions[1, ]
  duplicate$union_id <- 5L
  pop$unions <- rbind(pop$unions,duplicate)
  expect_equal(nrow(extract_unpartnered_births(pop)),6)
  pop$unions$union_type <- NULL
  pop$unions$previous_cohab_id <- NULL
  # Conversion needs a typed link; use a real ending for the untyped fixture.
  pop$unions$end_reason[3] <- "dissolution"
  expect_equal(extract_unpartnered_births(pop)$birth_time,c(0,2.5,4,4.5,7.5,8))
  pop$persons$group <- NULL
  pop$groups <- NULL
  x <- extract_unpartnered_births(pop)
  expect_false(any(c("group_mother","group_father") %in% names(x)))
})

test_that("birth extractor preserves empty schema and validates histories", {
  pop <- unpartnered_birth_history()
  pop$unions <- pop$unions[0, ]
  expect_equal(nrow(extract_unpartnered_births(pop)),12)
  pop$persons$father_id <- NA_integer_
  x <- extract_unpartnered_births(pop)
  expect_equal(nrow(x),0)
  expect_named(x,names(augment_persons(pop)))
  pop$persons <- pop$persons[0, ]
  pop$ancestry <- pop$ancestry[0, ]
  expect_s3_class(extract_unpartnered_births(pop),"tbl_df")
  expect_type(extract_unpartnered_births(pop)$person_id,"integer")
  expect_error(extract_unpartnered_births(list()),"hammel_population")
  pop <- unpartnered_birth_history()
  pop$persons$father_id[4] <- 999L
  expect_error(extract_unpartnered_births(pop),"father_id")
})
