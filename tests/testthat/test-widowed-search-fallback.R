test_that("pooled previous-union inputs equal detailed schedules for both ending causes", {
  for (duration in c(FALSE,TRUE)) for (reason in c("dissolution","death_husband")) {
    rates <- data.frame(previous_union=c("none","dissolved"),age_start=15,age_end=Inf,rate=c(0,100))
    if(duration) {
      rates$age_start[2] <- rates$age_end[2] <- NA_real_
      rates$dur_start <- c(NA,0)
      rates$dur_end <- c(NA,Inf)
    }
    detailed <- rbind(rates[1,],transform(rates[2,],previous_union="separated"),
      transform(rates[2,],previous_union="widowed"))
    m <- make_rate_model(4)
    m$population$unions <- tibble::tibble(union_id=1,wife_id=1,husband_id=2,
      start_time=-5,end_time=-1,end_reason=reason)
    if(reason=="death_husband") m$population$persons$death_time[2] <- -1
    m <- set_union_formation(m,"rate",rates=rates,kinship_generations=0)
    actual <- simulate(m,.5,seed=42)
    expected <- simulate(set_union_formation(m,"rate",rates=detailed,kinship_generations=0),.5,seed=42)
    expect_equal(actual,expected)
    expect_equal(nrow(actual$unions),2L)
    expect_equal(m$processes$union_formation$rates,tibble::as_tibble(rates))
  }
})

test_that("search rates and cohabitation probabilities accept the same complete alternatives", {
  pooled <- data.frame(previous_union=c("none","dissolved"),age_start=15,age_end=Inf,rate=c(0,100))
  detailed <- rbind(pooled[1,],transform(pooled[2,],previous_union="separated"),
    transform(pooled[2,],previous_union="widowed"))
  for (x in list(pooled,detailed)) {
    p <- x
    p$prob <- .5
    p$rate <- NULL
    m <- set_union_formation(make_rate_model(4),"rate",rates=x,cohab_prob=p,
      union_convert_rate=data.frame(dur_start=0,dur_end=Inf,rate=0))
    expect_s3_class(m,"hammel_model")
    expect_silent(simulate(m,0,seed=1))
  }
  expect_error(union_formation("rate",rates=detailed[-3,]),"complete|all three")
  expect_error(union_formation("rate",rates=rbind(pooled,detailed[3,])),"mix")
  expect_error(union_formation("rate",rates=rbind(pooled,detailed[2:3,])),"mix")
  calendar <- rbind(transform(pooled,group="A",time_start=0,time_end=1),
    transform(detailed,group="A",time_start=1,time_end=2),
    transform(detailed,group="B",time_start=0,time_end=2))
  expect_silent(union_formation("rate",rates=calendar))
  expect_equal(sort(unique(expand_previous_union_rates(calendar)$previous_union)),
    sort(c("none","separated","widowed")))

  # Execute pooled and detailed inputs across differing calendar/group bands.
  p <- calendar
  p$prob <- ifelse(p$previous_union=="none",0,1)
  p$rate <- NULL
  duplicate_pooled <- function(x) {
    pooled_rows <- x$previous_union=="dissolved"
    rbind(x[!pooled_rows,],
      transform(x[pooled_rows,],previous_union="separated"),
      transform(x[pooled_rows,],previous_union="widowed"))
  }
  m <- make_rate_model(4)
  m$population$persons$group <- c("A","A","B","B")
  m$population$groups <- list(group = group())
  m$population$unions <- tibble::tibble(union_id=1,wife_id=1,husband_id=2,
    start_time=-5,end_time=-1,end_reason="death_husband")
  m$population$persons$death_time[2] <- -1
  attach <- function(rates,prob) set_union_formation(m,"rate",rates=rates,
    cohab_prob=prob,union_convert_rate=data.frame(dur_start=0,dur_end=Inf,rate=0),
    kinship_generations=0)
  actual <- simulate(attach(calendar,p),2,seed=42)
  expected <- simulate(attach(duplicate_pooled(calendar),duplicate_pooled(p)),2,seed=42)
  expect_equal(actual,expected)
  expect_equal(actual$unions$union_type[2],"cohab")

  # Equal effective categories cannot switch clocks when representation changes.
  inconsistent <- calendar
  inconsistent$dur_start <- inconsistent$dur_end <- NA_real_
  rows <- inconsistent$previous_union=="dissolved"
  inconsistent$age_start[rows] <- inconsistent$age_end[rows] <- NA_real_
  inconsistent$dur_start[rows] <- 0
  inconsistent$dur_end[rows] <- Inf
  expect_error(union_formation("rate",rates=inconsistent),"same age or duration clock")
})
