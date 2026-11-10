test_that("fractional fertility epochs retain the run-start accounting grid", {
  # Women aged 20 with a child born at the run start and a 10-year interval
  # stay ineligible; rates change at fractional times. Unrealized births are
  # counted on the one-year grid from the run start.
  n <- 2000
  for(origin in c(2000,2000.25)) {
    m <- make_rate_model(3,search_rate=0) |>
      set_population(spacing_population(n, age = 20, sons = 0, end_time = origin))
    rates <- data.frame(age_start=15,age_end=50,rate=.1)
    epochs <- do.call(rbind,lapply(1:4,function(i) {
      x <- rates
      x$rate <- c(.1,.2,.4,.3)[i]
      x$time_start <- origin+c(0,.5,1.25,2.2)[i]
      x$time_end <- origin+c(.5,1.25,2.2,Inf)[i]
      rbind(transform(x,union_type="partnered"),transform(x,union_type="unpartnered"))
    }))
    m <- set_fertility(m,epochs,birth_interval=10)
    out <- run_simulation(m,3.4,seed=42)
    u <- out$unrealized_births
    expect_true(all(u$period_end<=origin+ceiling(u$period_end-origin)))
    # Records do not bridge any run-anchored annual boundary.
    expect_true(all(floor(u$period_start-origin+1e-9)==ceiling(u$period_end-origin-1e-9)-1))
    observed <- calculate_fertility_rates(out,origin+c(0,1,2,3,3.4),c(0,Inf),unrealized=TRUE)
    expect_equal(count_run_events(out, "births"),0)
    expect_equal(observed$births,c(n,0,0,0)) # the supplied children were born at the run start
    for (k in 1:4) expect_count(observed$unrealized_births[k], c(.15,.35,.32,.12)[k] * n)
    expect_false(anyNA(observed$rate_fertility_adjusted))
  }
})

test_that("execution population coverage checks ignore unreachable historical and future slices", {
  m <- make_rate_model(4,search_rate=0)
  m$population$persons$group <- c("A","A","B","B")
  m$population$groups <- list(group = group())
  m$population$end_time <- 2000
  m$population$persons$birth_time <- m$population$persons$birth_time+2000
  extend <- function(x) {
    rbind(transform(x,group="A",time_start=1900,time_end=1950),
      transform(x,group="A",time_start=1950,time_end=2010),
      transform(x,group="B",time_start=1950,time_end=2010),
      transform(x,group="B",time_start=2010,time_end=2100))
  }
  mortality <- extend(m$processes$mortality$rates)
  mortality$age_start[mortality$time_start==1900] <- 1
  mortality$age_end[mortality$time_start==2010] <- 50
  m <- set_mortality(m,mortality)
  m <- set_fertility(m,extend(m$processes$fertility$rates))
  m <- set_union_formation(m,"rate",rates=extend(m$processes$union_formation$rates),
    cohab_prob=transform(extend(data.frame(age_start=0,age_end=Inf,rate=0)),prob=0)[,
      c("age_start","age_end","group","time_start","time_end","prob")],
    union_convert_rate=transform(extend(data.frame(dur_start=0,dur_end=Inf,rate=0)),group_w=group,group=NULL))
  pairs <- directional_pairs(data.frame(group_pair=c("A:A","A:B","B:B"),dur_start=0,dur_end=Inf,rate=0))
  dissolution <- rbind(transform(pairs[1,],time_start=1900,time_end=1950),
    transform(pairs,time_start=1950,time_end=2010),
    transform(pairs[3,],time_start=2010,time_end=2100))
  m <- set_union_dissolution(m,dissolution)
  expect_s3_class(run_simulation(m,10,seed=42),"hammel_population")
  expect_s3_class(run_simulation(m,0,seed=42),"hammel_population")
  expect_error(run_simulation(m,10.1,seed=42),"Mortality schedules")
  full_age <- m
  fixed <- m$processes$mortality$rates
  fixed$age_end[fixed$time_start==2010] <- Inf
  full_age <- set_mortality(full_age,fixed)
  expect_error(run_simulation(full_age,10.1,seed=42),"missing schedules")
  old <- m
  old$population$end_time <- 1940
  old$population$persons$birth_time <- old$population$persons$birth_time-60
  expect_error(run_simulation(old,1,seed=42),"Mortality schedules")
  # Unreachable malformed rows still fail construction.
  invalid <- m$processes$mortality$rates
  invalid$age_end[1] <- -1
  expect_error(set_mortality(m,invalid))
})
