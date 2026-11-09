parity_rates <- function(values=c(.2,.3,.4),lower=15,upper=50) {
  tibble::tibble(age_start=lower,age_end=upper,parity_start=seq_along(values)-1L,
    parity_end=c(seq_len(length(values)-1L),Inf),rate=values)
}

parity_model <- function(values=c(.2,.3,.4),interval=.75,n=100) {
  quiet_unreachable(make_model(n,sex="f",age=20) |>
    set_fertility(parity_rates(values),birth_interval=interval))
}

test_that("parity schedules require complete partitions in every stratum", {
  rates <- parity_rates()
  expect_s3_class(quiet_unreachable(fertility(rates)),"hammel_process")
  expect_error(fertility(rates[-1, ]),"cover")
  expect_error(fertility(rates[-2, ]),"cover")
  expect_error(fertility(rates[-3, ]),"cover")
  for (column in c("parity_start","parity_end")) {
    bad <- rates
    bad[[column]] <- NULL
    expect_error(fertility(bad),column)
  }
  for (bad_value in c(-1,.5,Inf,NA_real_)) {
    bad <- rates
    bad$parity_start[1] <- bad_value
    expect_error(fertility(bad))
  }
  bad <- rates
  bad$parity_end[1] <- 2
  expect_error(fertility(bad),"cover")
  expect_error(fertility(rbind(rates,rates[1, ])),"Age intervals")
  grouped <- rbind(transform(rates,group="A"),transform(rates,group="B"))
  expect_s3_class(quiet_unreachable(fertility(grouped)),"hammel_process")
  expect_error(fertility(grouped[-5, ]),"cover")
  typed <- rbind(transform(rates,union_type="partnered"),transform(rates,union_type="unpartnered"))
  expect_s3_class(quiet_unreachable(fertility(typed)),"hammel_process")
  expect_error(fertility(typed[-5, ]),"cover")
  expect_s3_class(fertility(parity_rates(.2)),"hammel_process")
})

test_that("birth history initializes parity without custom column overrides", {
  model <- parity_model(c(0,20,0),interval=0,n=2)
  people <- model$population$persons
  # Woman 1 has one recorded child, already deceased. Founder 2 has none.
  child <- people[1, ]
  child$person_id <- 3L
  child$birth_time <- -1
  child$death_time <- -.5
  child$mother_id <- 1L
  people <- rbind(people,child)
  people$parity <- 99L
  model <- set_population(model,hammel_population(people,end_time=0, groups = if ("group" %in% names(people)) "group"))
  result <- run_simulation(model,1,42)
  born <- subset(result$persons,birth_time>=0)
  expect_equal(nrow(born),1)
  expect_equal(born$mother_id,1)
  expect_equal(result$persons$parity[1:3],rep(99L,3)) # Custom field is not authoritative.
  expect_equal(run_simulation(set_population(model,result),1,
    43)$simulation$births,0)
})

test_that("successful births increment parity and retain actual birth spacing", {
  model <- parity_model(c(20,20,0),interval=.75,n=50)
  # Large positive prior rates are infeasible for compensation, so use zero
  # spacing for this high-hazard birth-count limit; test spacing separately.
  model <- set_fertility(model,parity_rates(c(20,20,0)),birth_interval=0)
  result <- run_simulation(model,1,42)
  counts <- table(result$persons$mother_id)
  expect_true(all(counts<=2))
  expect_equal(sum(counts),100)
  spaced <- parity_model(c(.6,.6,.6),interval=.75,n=300)
  x <- run_simulation(spaced,8,100)
  times <- split(x$persons$birth_time[!is.na(x$persons$mother_id)],
    x$persons$mother_id[!is.na(x$persons$mother_id)])
  differences <- unlist(lapply(times,function(t) diff(sort(t))))
  expect_true(length(differences)>0)
  expect_true(all(differences>=.75-1e-12))
})

test_that("cancelled births and initial nulliparity do not advance parity", {
  model <- parity_model(c(.8,0),n=20)
  model <- set_union_formation(model,"on_demand",cancel_birth_if_no_match=TRUE)
  x <- run_simulation(model,3,42)
  expect_equal(x$simulation$births,0)
  model <- set_population(parity_model(c(.8,0),n=20),x)
  y <- run_simulation(model,3,43)
  expect_gt(y$simulation$births,0)
  expect_true(all(table(y$persons$mother_id)<=1))
})

test_that("absent parity and zero-spacing equal-parity rates give the same fertility", {
  # Actual-parity cells (D101) split the same total intensity into more clocks,
  # so histories differ by seed but agree in distribution.
  n <- 2000
  plain <- make_model(n,sex="f",birth_rate=.2,birth_interval=0)
  strat <- set_fertility(plain,parity_rates(rep(.2,3)),birth_interval=0)
  founders <- function(x) sum(x$persons$mother_id %in% seq_len(n))
  a <- founders(run_simulation(plain,4,42))
  b <- founders(run_simulation(strat,
    4,43))
  expect_lt(abs(a - b), 4 * sqrt(2 * .2 * 4 * n))
  expect_false("parity" %in% names(run_simulation(strat,1,42)$persons))
})

test_that("parity composes with grouped typed fertility, conversion and stratified rates", {
  set.seed(33)
  pop <- generate_uniform_population(500,50,0,groups = list(group = c(A=1,B=1)))
  base <- parity_rates(c(.3,.2,.1))
  rates <- rbind(transform(base,union_type="unpartnered"),
    transform(base,union_type="cohab",rate=rate*.8),
    transform(base,union_type="marriage",rate=rate*1.1))
  b_factor <- c(unpartnered=1.1,cohab=.9,marriage=1.2)
  stratified <- do.call(rbind,lapply(names(b_factor),function(type)
    stratify_rates(rates[rates$union_type==type,],list(group=c(A=1,B=b_factor[[type]])))))
  model <- hammel_model() |> set_population(pop) |>
    set_mortality(data.frame(sex=c("m","f"),age_start=0,age_end=Inf,rate=.01)) |>
    (\(m) quiet_unreachable(set_fertility(m,stratified)))() |>
    set_union_formation("rate",rates=data.frame(age_start=15,age_end=Inf,rate=.4),
      cohab_prob=data.frame(age_start=15,age_end=Inf,prob=.5),
      union_convert_rate=data.frame(dur_start=0,dur_end=Inf,rate=.2),
      kinship_generations=0L,diagnostics=TRUE) |>
    set_union_dissolution(data.frame(dur_start=0,dur_end=Inf,rate=.1))
  explicit <- rbind(transform(rates,group="A"),transform(rates,group="B"))
  factor <- c(unpartnered=1.1,cohab=.9,marriage=1.2)[explicit$union_type]
  explicit$rate <- explicit$rate * ifelse(explicit$group=="B",factor,1)
  equivalent <- quiet_unreachable(set_fertility(model,explicit))
  x <- run_simulation(model,8,44)
  expect_identical(run_simulation(equivalent,8,44),x)
  expect_gt(x$simulation$births,0)
  expect_gt(x$simulation$conversions,0)
  expect_false(anyNA(x$unions$wife_id))
})

test_that("grouped zero-spacing parity schedules match shared rates in distribution", {
  set.seed(123)
  pop <- generate_uniform_population(400,50,0,groups = list(group = c(A=1,B=1)))
  base <- data.frame(age_start=15,age_end=50,union_type=c("unpartnered","partnered"),rate=.2)
  model <- hammel_model() |> set_population(pop) |>
    set_mortality(data.frame(sex=c("m","f"),age_start=0,age_end=Inf,rate=.01)) |>
    set_fertility(base,birth_interval=0) |>
    set_union_formation("rate",rates=data.frame(age_start=15,age_end=Inf,rate=.3),kinship_generations=0L)
  rates <- do.call(rbind,lapply(0:2,function(k) transform(base,parity_start=k,parity_end=if(k==2) Inf else k+1)))
  # Different cell structure, same intensity: compare birth counts statistically.
  a <- run_simulation(set_fertility(model,rates,birth_interval=0),8,
    44)$simulation$births
  b <- run_simulation(model,8,45)$simulation$births
  expect_lt(abs(a - b), 4 * sqrt(a + b))
})

test_that("history parity is recovered even when the preceding run was unstratified", {
  base <- make_model(100,sex="f",age=20,birth_rate=.3,birth_interval=0)
  burned <- run_simulation(base,5,42)
  observed <- table(burned$persons$mother_id)
  model <- set_population(parity_model(c(0,0,0),interval=.75,n=0),burned)
  expect_equal(run_simulation(model,1,43)$simulation$births,0)
  # Only women with exactly one recorded birth may give a further birth.
  model <- set_fertility(model,parity_rates(c(0,20,0)),birth_interval=0)
  x <- run_simulation(model,1,44)
  mothers <- x$persons$mother_id[x$persons$birth_time>=5]
  mothers <- mothers[!is.na(mothers)]
  expect_gt(length(mothers),0)
  expect_true(all(as.character(mothers) %in% names(observed)[observed==1]))
})

test_that("edited parity schedules are revalidated and broad bands run", {
  model <- parity_model()
  model$processes$fertility$rates <- model$processes$fertility$rates[-1, ]
  expect_error(run_simulation(model,0),"cover")
  # Very broad parity bands need no per-parity storage until women reach them.
  rates <- data.frame(age_start=15,age_end=Inf,parity_start=c(0,1000000),
    parity_end=c(1000000,Inf),rate=c(.1,.2))
  broad <- quiet_unreachable(set_fertility(make_model(20,sex="f",age=20),rates))
  expect_gt(run_simulation(broad,5,3)$simulation$births,0)
})
