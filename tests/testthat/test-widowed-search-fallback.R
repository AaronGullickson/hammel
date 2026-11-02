test_that("omitted widowed rates match explicit copies for age and duration", {
  for (duration in c(FALSE, TRUE)) {
    rates <- data.frame(previous_union=c("none","dissolved"),
      age_start=c(15,15), age_end=Inf, rate=c(0,100))
    if (duration) {
      rates$age_start[2] <- rates$age_end[2] <- NA_real_
      rates$dur_start <- c(NA,0)
      rates$dur_end <- c(NA,Inf)
    }
    explicit <- rbind(rates, transform(rates[2,], previous_union="widowed"))
    m <- make_rate_model(4)
    m$population$unions <- tibble::tibble(union_id=1, wife_id=1, husband_id=2,
      start_time=-5, end_time=-1, end_reason="death_husband")
    m$population$persons$death_time[2] <- -1
    m <- set_union_formation(m,"rate",rates=rates,kinship_generations=0)
    actual <- simulate(m,.5,seed=42)
    expected <- simulate(set_union_formation(m,"rate",rates=explicit,
      kinship_generations=0),.5,seed=42)
    expect_equal(actual,expected)
    expect_equal(nrow(actual$unions),2L)
    expect_equal(actual$unions$end_reason[1],"death_husband")
    expect_equal(m$processes$union_formation$rates,tibble::as_tibble(rates))
    # Edited stored tables are revalidated and receive the same execution fallback.
    m$processes$union_formation$rates <- tibble::as_tibble(rates)
    expect_equal(simulate(m,.5,seed=42),actual)
  }
})

test_that("fallback is local to groups and active calendar slices", {
  base <- data.frame(previous_union=c("none","dissolved"),
    age_start=15,age_end=Inf,rate=c(0,100))
  a <- rbind(transform(base,group="A",time_start=0,time_end=.4),
    transform(base,group="A",time_start=.4,time_end=1))
  # Explicit zero in A's first slice suppresses fallback only in that slice.
  a <- rbind(a,transform(a[2,],previous_union="widowed",rate=0))
  b <- rbind(transform(base,group="B",time_start=0,time_end=1),
    transform(base[2,],group="B",previous_union="widowed",rate=0,
      time_start=0,time_end=1))
  rates <- rbind(a,b)
  explicit <- rbind(rates,transform(a[4,],previous_union="widowed"))
  for (group in c("A","B")) {
    m <- make_rate_model(4)
    m$population$persons$group <- rep(c("A","B"),each=2)
    m$population$persons$group[1] <- group
    m$population$unions <- tibble::tibble(union_id=1,wife_id=1,husband_id=2,
      start_time=-5,end_time=-1,end_reason="death_husband")
    m$population$persons$death_time[2] <- -1
    m <- set_union_formation(m,"rate",rates=rates,kinship_generations=0)
    actual <- simulate(m,1,seed=42)
    expected <- simulate(set_union_formation(m,"rate",rates=explicit,
      kinship_generations=0),1,seed=42)
    expect_equal(actual,expected)
    expect_equal(nrow(actual$unions),if(group=="A") 2L else 1L)
    if(group=="A") expect_gte(actual$unions$start_time[2],.4)
  }
})

test_that("explicit widowed gaps and effective clock consistency are preserved", {
  rates <- data.frame(previous_union=c("none","dissolved","widowed"),
    age_start=c(15,15,30),age_end=Inf,rate=c(0,100,100))
  m <- make_rate_model(4)
  m$population$unions <- tibble::tibble(union_id=1,wife_id=1,husband_id=2,
    start_time=-5,end_time=-1,end_reason="death_husband")
  m$population$persons$death_time[2] <- -1
  m <- set_union_formation(m,"rate",rates=rates,kinship_generations=0)
  expect_equal(nrow(simulate(m,1,seed=42)$unions),1L)
  bad <- rbind(transform(rates,group="A"),transform(rates[-3,],group="B"))
  bad$dur_start <- bad$dur_end <- NA_real_
  idx <- bad$group=="B" & bad$previous_union=="dissolved"
  bad$age_start[idx] <- bad$age_end[idx] <- NA_real_
  bad$dur_start[idx] <- 0
  bad$dur_end[idx] <- Inf
  expect_error(union_formation("rate",rates=bad),"same age or duration clock")
  calendar <- rates
  calendar$dur_start <- c(NA,0,NA)
  calendar$dur_end <- c(NA,Inf,NA)
  calendar$age_start[2] <- calendar$age_end[2] <- NA_real_
  calendar <- rbind(transform(calendar,time_start=0,time_end=1),
    transform(calendar[-3,],time_start=1,time_end=2))
  # Explicit age-based widowhood first, then inherited duration-based widowhood.
  expect_error(union_formation("rate",rates=calendar),"across calendar periods")
  p <- transform(rates[-3,],prob=.5)
  p$rate <- NULL
  expect_error(set_union_formation(m,"rate",rates=rates,cohab_prob=p),"all three")
})
