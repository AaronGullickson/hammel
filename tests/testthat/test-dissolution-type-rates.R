test_that("dissolution type strata partition episode exposure and events", {
  unions <- tibble::tibble(union_id = 1:4, start_time = c(0,2,0,1),
    end_time = c(2,4,3,NA), end_reason = c("conversion","dissolution","death_wife",NA),
    union_type = c("cohab","marriage","cohab","marriage"),
    group_w = c("A","B","A","A"), group_h = c("B","A","A","A"))
  for (grouped in c(FALSE,TRUE)) {
    pooled <- calculate_dissolution_rates(unions,c(0,5),c(0,Inf), by = if (isTRUE(grouped)) c("group_w", "group_h"))
    split <- calculate_dissolution_rates(unions,c(0,5),c(0,Inf), by = c(if (isTRUE(grouped)) c("group_w", "group_h"), "union_type"))
    expect_s3_class(split,"tbl_df")
    expect_equal(sum(split$union_years),sum(pooled$union_years))
    expect_equal(sum(split$dissolutions),sum(pooled$dissolutions))
    expect_equal(sum(split$union_years[split$union_type == "cohab"]),5)
    expect_equal(sum(split$union_years[split$union_type == "marriage"]),6)
    expect_equal(sum(split$dissolutions[split$union_type == "cohab"]),0)
    expect_equal(sum(split$dissolutions[split$union_type == "marriage"]),1)
  }
  x <- calculate_dissolution_rates(unions,c(0,5),c(0,1,Inf),by = "union_type")
  expect_equal(x$union_years,c(2,3,2,4))
  expect_equal(x$dissolutions,c(0,0,0,1))
  expect_equal(x$rate_dissolution,c(0,0,0,.25))
  expect_error(calculate_dissolution_rates(unions,c(0,5),c(0,Inf), by = NA),"`by` must")
  bare <- unions; bare$union_type <- NULL
  expect_error(calculate_dissolution_rates(bare,c(0,5),c(0,Inf),by = "union_type"),"union_type")
  for (bad in list(c(NA,"marriage","cohab","marriage"),rep("partnered",4))) {
    bare$union_type <- bad
    expect_error(calculate_dissolution_rates(bare,c(0,5),c(0,Inf),by = "union_type"),"cohab or marriage")
  }
  empty <- calculate_dissolution_rates(unions[0, ],c(0,5),c(0,Inf),by = "union_type")
  expect_identical(empty$union_type,c("cohab","marriage"))
  expect_equal(empty$union_years,c(0,0))
  expect_true(all(is.na(empty$rate_dissolution)))
  empty_group <- calculate_dissolution_rates(unions[0, ],c(0,5),c(0,Inf),
    by = c("group_w", "group_h", "union_type"))
  expect_equal(nrow(empty_group),0)
  expect_type(empty_group$union_type,"character")
})

test_that("typed population dissolution exposure resets on conversion", {
  pop <- hammel_population(
    tibble::tibble(person_id=1:2,sex=c("f","m"),birth_time=-20,
      death_time=NA_real_,mother_id=NA_integer_,father_id=NA_integer_,group=c("A","B")),
    tibble::tibble(union_id=1:2,wife_id=1L,husband_id=2L,start_time=c(1,3),
      end_time=c(3,NA),end_reason=c("conversion",NA),union_type=c("cohab","marriage"),
      previous_cohab_id=c(NA_integer_,1L)),end_time=5, groups = "group")
  x <- calculate_dissolution_rates(pop,c(0,5),c(0,1,Inf),by = c("group_w", "group_h", "union_type"))
  expect_equal(sum(x$union_years[x$union_type == "cohab"]),2)
  expect_equal(sum(x$union_years[x$union_type == "marriage"]),2)
  expect_equal(sum(x$dissolutions),0)
  expect_equal(sum(x$union_years[x$dur_start == 0]),2)
})
