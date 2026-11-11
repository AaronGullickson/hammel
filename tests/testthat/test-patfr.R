patfr_schedule <- function(rates=c(.2,.4,.6), lower=20, upper=21) {
  tibble::tibble(age_start=lower,age_end=upper,
    parity_start=seq_along(rates)-1L,
    parity_end=c(seq_len(length(rates)-1L),Inf),rate=rates)
}

test_that("PATFR reproduces hand-calculated HFD midpoint births", {
  x <- calculate_patfr(patfr_schedule())
  # Unit cohort, first-year q1=2/11, L0=10/11, L1=1/11,
  # so first births=2/11, second births=2/55, third+ births=0.
  expect_equal(x$patfr,12/55,tolerance=1e-14)
  expect_equal(x$age_start,20)
  expect_equal(x$age_end,21)
  expect_s3_class(x,"tbl_df")
  expect_identical(calculate_tfr(patfr_schedule()),x)
  expect_named(x,c("age_start","age_end","patfr"))
  expect_equal(calculate_patfr(patfr_schedule(c(0,1,1)))$patfr,0)
})

test_that("PATFR obeys first-birth probability and repeated-birth benchmarks", {
  x <- patfr_schedule(c(.2,0),lower=15,upper=45)
  q <- .2/(1+.2/2)
  expect_equal(calculate_patfr(x)$patfr,1-(1-q)^30,tolerance=1e-13)
  # Identical rates at all parities give constant total birth intensity.
  for (rate in c(0,.1,.7)) {
    x <- patfr_schedule(rep(rate,5),lower=15,upper=45)
    expect_equal(calculate_patfr(x)$patfr,30*rate,tolerance=1e-12)
  }
  # Pooled parity zero+ is a recurrent event process, without a parity cap.
  x <- patfr_schedule(.4,lower=15,upper=45)
  expect_equal(calculate_patfr(x)$patfr,12)
  x <- patfr_schedule(c(.8,10))
  expect_gt(calculate_patfr(x)$patfr,2)
})

test_that("age and parity aggregation preserve constant-schedule results", {
  broad <- patfr_schedule(rep(.2,5),lower=15,upper=45)
  annual <- do.call(rbind,lapply(15:44,function(age) transform(broad,age_start=age,age_end=age+1)))
  expect_equal(calculate_patfr(annual),calculate_patfr(broad))
  expect_equal(calculate_patfr(annual[nrow(annual):1, ]),calculate_patfr(broad))
  pooled <- tibble::tibble(age_start=15,age_end=45,
    parity_start=c(0,3),parity_end=c(3,Inf),rate=.2)
  expect_equal(calculate_patfr(pooled)$patfr,6)
  expect_equal(calculate_patfr(transform(pooled,age_start=15.25,age_end=45.75))$patfr,6.1)
})

test_that("PATFR separates all TFR strata and preserves metadata and group types", {
  base <- patfr_schedule(rep(.1,3),lower=15,upper=20)
  grid <- expand.grid(time_start=c(2000,2005),group=c(20L,10L),
    union_type=c("cohab","marriage"),stringsAsFactors=FALSE)
  rows <- lapply(seq_len(nrow(grid)),function(i) {
    x <- base
    x$rate <- i*.01
    x$time_start <- grid$time_start[i]
    x$time_end <- grid$time_start[i]+5
    x$group <- grid$group[i]
    x$union_type <- grid$union_type[i]
    x
  })
  rates <- do.call(rbind,rows)
  rates$sex <- "f"
  male <- rates
  male$sex <- "m"
  male$rate <- 1
  rates <- rbind(rates,male)
  unknown <- tibble::tibble(time_start=2000,time_end=2010,births=2)
  attr(rates,"unassigned_births") <- unknown
  x <- calculate_patfr(rates)
  expect_equal(x$patfr,seq_len(nrow(grid))*.05)
  expect_identical(x$group,grid$group)
  expect_identical(x$union_type,grid$union_type)
  expect_identical(x$time_start,grid$time_start)
  expect_identical(attr(x,"unassigned_births"),unknown)
  expect_identical(calculate_tfr(rates),x)
  names(rates)[names(rates)=="union_type"] <- "union_status"
  expect_identical(calculate_tfr(rates)$union_status,grid$union_type)
})

test_that("missing PATFR rates invalidate only their included stratum", {
  rates <- patfr_schedule()
  other <- rates
  rates$group <- "A"
  other$group <- "B"
  rates$rate[2] <- NA_real_
  x <- calculate_patfr(rbind(rates,other))
  expect_true(is.na(x$patfr[1]))
  expect_equal(x$patfr[2],12/55)
  rates$rate <- NA_real_
  expect_true(is.na(calculate_patfr(rates)$patfr))
})

test_that("PATFR omits open ages and warns on known nonzero omitted rates", {
  rates <- patfr_schedule()
  tail <- transform(rates,age_start=21,age_end=Inf,rate=0)
  expect_warning(x <- calculate_patfr(rbind(rates,tail)),NA)
  expect_equal(x$patfr,12/55)
  tail$rate[1] <- .01
  expect_warning(y <- calculate_patfr(rbind(rates,tail)),"PATFR.*underestimated")
  expect_equal(y,x)
  expect_warning(calculate_tfr(rbind(rates,tail)),"PATFR.*underestimated")
  tail$rate <- NA_real_
  expect_warning(z <- calculate_patfr(rbind(rates,tail)),NA)
  expect_equal(z,x)
  tail$rate <- .01
  expect_warning(
    expect_error(calculate_patfr(tail), "finite end"),
    "PATFR.*underestimated"
  )
})

test_that("PATFR validates age/parity coverage and invalid values", {
  rates <- patfr_schedule()
  expect_error(calculate_patfr(rates[-1, ]),"start at zero")
  expect_error(calculate_patfr(rates[-3, ]),"open final")
  expect_error(calculate_patfr(rbind(rates,rates[1, ])),"exactly once")
  bad <- rates
  bad$parity_start[2] <- 1.5
  expect_error(calculate_patfr(bad),"integer")
  bad <- rates
  bad$parity_start[2] <- 2
  bad$parity_end[2] <- 3
  expect_error(calculate_patfr(bad),"exactly once")
  bad <- rates
  bad$age_end[2] <- 20.5
  expect_error(calculate_patfr(bad),"exactly once")
  bad <- rates
  bad$rate[2] <- -1
  expect_error(calculate_patfr(bad),"nonnegative")
  bad <- rates
  bad$rate[2] <- Inf
  expect_error(calculate_patfr(bad),"finite")
  bad <- rates
  bad$rate[1] <- 3
  expect_error(calculate_patfr(bad),"probabilities")
  bad <- rates
  bad$parity_end <- NULL
  expect_error(calculate_tfr(bad),"parity_end")
})

test_that("PATFR accepts observation rate names and rejects empty tables", {
  rates <- patfr_schedule()
  names(rates)[names(rates)=="rate"] <- "rate_fertility"
  expect_equal(calculate_patfr(rates)$patfr,12/55)
  rates$group <- 20L
  expect_error(calculate_patfr(rates[0, ]), "at least one age band")
  expect_error(calculate_tfr(rates[0, ]), "at least one age band")
})

test_that("fertility processes dispatch to PATFR and mortality is rejected", {
  process <- fertility(data.frame(age_start=15,age_end=20,rate=.1))
  process$rates <- patfr_schedule()
  expect_identical(calculate_tfr(process),calculate_patfr(process$rates))
  expect_identical(calculate_patfr(process),calculate_patfr(process$rates))
  death <- mortality(data.frame(sex=c("m","f"),age_start=0,age_end=Inf,rate=.1))
  expect_error(calculate_patfr(death),"fertility process")
})
