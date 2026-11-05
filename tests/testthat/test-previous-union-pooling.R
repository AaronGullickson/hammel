previous_pool_history <- function() {
  persons <- tibble::tibble(person_id=1:6,sex=c("f","f",rep("m",4)),
    birth_time=-20,death_time=c(NA,NA,NA,-1,NA,NA),mother_id=NA_integer_,
    father_id=NA_integer_,group=rep(c(10L,20L),3))
  unions <- tibble::tibble(union_id=1:5,wife_id=c(1,2,1,1,2),husband_id=c(3,4,5,5,6),
    start_time=c(-5,-5,0,1,2),end_time=c(-1,-1,1,NA,NA),
    end_reason=c("dissolution","death_husband","conversion",NA,NA),
    union_type=c("marriage","marriage","cohab","marriage","marriage"),
    previous_cohab_id=c(NA_integer_,NA_integer_,NA_integer_,3L,NA_integer_))
  hammel_population(persons,unions,end_time=3)
}

test_that("pooled previous-union summaries sum exposure and events before exact intervals", {
  pop <- previous_pool_history()
  for (marriage in c(FALSE,TRUE)) for (grouped in c(FALSE,TRUE)) {
    fun <- if(marriage) calculate_marriage_rates else calculate_union_rates
    for (duration in if(marriage) list(NULL) else list(NULL,c(0,1,Inf))) {
      for (typed in if(marriage) FALSE else c(FALSE,TRUE)) {
        args <- list(population=pop,time_breaks=0:3,age_breaks=c(0,21,Inf),
          by = c("previous_union", if (grouped) "group", if (typed) "union_type"),conf.level=.95)
        if(!marriage) args <- c(args,list(dur_breaks=duration))
        detailed <- do.call(fun,args)
        pooled <- do.call(fun,c(args,list(pool_dissolved=TRUE)))
        expect_setequal(unique(detailed$previous_union),c("none","separated","widowed"))
        expect_setequal(unique(pooled$previous_union),c("none","dissolved"))
        events <- if(marriage) "marriages" else "unions"
        rate <- if(marriage) "rate_marriage" else "rate_union"
        expected_keys <- detailed[setdiff(names(detailed),c("person_years",events,rate,"conf_low","conf_high"))]
        expected_keys$previous_union[expected_keys$previous_union!="none"] <- "dissolved"
        for(i in seq_len(nrow(pooled))) {
          rows <- rep(TRUE,nrow(detailed))
          for(key in names(expected_keys)) {
            value <- pooled[[key]][i]
            rows <- rows & if(is.na(value)) is.na(expected_keys[[key]]) else
              !is.na(expected_keys[[key]]) & expected_keys[[key]]==value
          }
          count <- sum(detailed[[events]][rows])
          exposure <- sum(detailed$person_years[rows])
          expect_equal(pooled[[events]][i],count)
          expect_equal(pooled$person_years[i],exposure)
          expect_equal(pooled[[rate]][i],if(exposure>0) count/exposure else NA_real_)
          ci <- if(exposure>0) as.numeric(stats::poisson.test(count,T=exposure)$conf.int) else c(NA_real_,NA_real_)
          expect_equal(c(pooled$conf_low[i],pooled$conf_high[i]),ci)
        }
        if(grouped) expect_type(pooled$group,"integer")
      }
    }
  }
})

test_that("pooling requires previous strata and validates its option", {
  pop <- previous_pool_history()
  for(fun in list(calculate_union_rates,calculate_marriage_rates)) {
    for(bad in list(NA,1,"yes",c(TRUE,FALSE)))
      expect_error(fun(pop,c(0,3),c(0,Inf),pool_dissolved=bad),"pool_dissolved")
    expect_error(fun(pop,c(0,3),c(0,Inf),pool_dissolved=TRUE),"requires.*previous_union")
    x <- fun(pop,c(0,3),c(0,Inf),by = "previous_union",pool_dissolved=TRUE)
    expect_setequal(x$previous_union,c("none","dissolved"))
    expect_false("conf_low" %in% names(x))
    empty <- hammel_population(pop$persons[FALSE,],end_time=3)
    x <- fun(empty,c(0,3),c(0,Inf),by = "previous_union",pool_dissolved=TRUE,conf.level=.95)
    expect_true(all(x$person_years==0 & is.na(x$conf_low) & is.na(x$conf_high)))
  }
})

 test_that("duration stratification permits pooling without a redundant previous flag", {
  x <- calculate_union_rates(previous_pool_history(),c(0,3),c(0,Inf),
    dur_breaks=c(0,1,Inf),pool_dissolved=TRUE)
  expect_setequal(x$previous_union,c("none","dissolved"))
 })
