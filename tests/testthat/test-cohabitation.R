cohab_fixture <- function(n = 2L, prob = 1, conversion = data.frame(dur_start = 0, dur_end = Inf, rate = 0)) {
  model <- make_rate_model(n, search_rate = 20)
  set_union_formation(model, "rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 20),
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = prob),
    union_convert_rate = conversion, kinship_generations = 0L)
}

initial_cohabs <- function(n = 1L, conversion = .3, breakup = .2) {
  model <- cohab_fixture(2L * n,
    conversion = data.frame(dur_start = 0, dur_end = Inf, rate = conversion))
  model$population$persons$sex <- rep(c("f", "m"), n)
  model$population$unions <- tibble::tibble(union_id = seq_len(n),
    wife_id = seq.int(1L, 2L*n, 2L), husband_id = seq.int(2L, 2L*n, 2L),
    start_time = 0, end_time = NA_real_, end_reason = NA_character_, union_type = "cohab")
  model$processes$union_formation$rates$rate <- 0
  set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf,
    rate = c(breakup, 0), union_type = c("cohab", "marriage")))
}

test_that("cohabitation validates probability, conversion, type coverage and rate mode", {
  model <- cohab_fixture()
  expect_s3_class(model$processes$union_formation$cohab_prob, "tbl_df")
  expect_error(set_union_formation(model, "on_demand", cohab_prob = data.frame(age_start=0,age_end=Inf,prob=.5)), "rate-based")
  for (prob in c(-1, 1.1, NA, Inf)) {
    bad <- model; bad$processes$union_formation$cohab_prob$prob <- prob
    expect_error(simulate(bad, 0), "prob")
  }
  bad <- model; bad$processes$union_formation$cohab_prob$rate <- .1
  expect_error(simulate(bad, 0), "not.*rate")
  expect_error(set_union_formation(model, "rate", rates = model$processes$union_formation$rates,
    cohab_matching = age_gap_matching_model(0,5)), "require.*cohab_prob")
  bad <- model; bad$processes$union_formation$union_convert_rate$previous_union <- "none"
  expect_error(simulate(bad, 0), "Conversion rates cannot be stratified by `previous_union`")
  expect_error(fertility(data.frame(age_start=15,age_end=50,rate=.1,union_status="partnered")), "replaced")
  expect_error(fertility(data.frame(age_start=15,age_end=50,rate=.1,union_type=c("unpartnered","cohab"))), "union_type")
  expect_error(union_dissolution(data.frame(dur_start=0,dur_end=Inf,rate=.1,union_type="cohab")), "both cohab and marriage")
  expect_warning(simulate(cohab_fixture(conversion=NULL), .1, 42), "no conversions")
})

test_that("probability extremes and gaps assign entry types without simultaneous partners", {
  for (p in c(0, 1)) {
    model <- cohab_fixture(100, p)
    model$population$persons$sex[1:50] <- "f"
    result <- simulate(model, 1, 42)
    expect_equal(nrow(result$unions), 50)
    expect_true(all(result$unions$union_type == if (p) "cohab" else "marriage"))
    expect_equal(anyDuplicated(result$unions$husband_id), 0L)
  }
  model <- cohab_fixture()
  model$processes$union_formation$cohab_prob$age_start <- 50
  expect_equal(simulate(model, .5, 42)$unions$union_type, "marriage")
  expect_identical(simulate(model, .5, 42), simulate(model, .5, 42))
  model <- cohab_fixture(1)
  result <- simulate(model, .5, 42)
  expect_equal(nrow(result$unions), 0)
})

test_that("conversion preserves partners, resets dissolution duration and links episodes", {
  model <- initial_cohabs(conversion = 1e5, breakup = 0)
  model$processes$union_formation$union_convert_rate$dur_start <- 2
  early <- simulate(model, 1, 42)
  expect_equal(nrow(early$unions), 1)
  late <- simulate(set_population(model, early), 2, 43)
  expect_equal(nrow(late$unions), 2)
  expect_equal(late$unions$end_reason[1], "conversion")
  expect_equal(late$unions$start_time[2], late$unions$end_time[1])
  expect_equal(late$unions$previous_cohab_id[2], late$unions$union_id[1])
  expect_equal(late$unions$wife_id, c(1,1))
  expect_equal(late$unions$husband_id, c(2,2))
  expect_equal(late$simulation$conversions, 1)
  # Current partnership is continuous and conversion is not a second entry.
  rates <- calculate_union_rates(late, c(0,3), c(0,Inf))
  expect_equal(sum(rates$unions), 2)
  counts <- demographic_counts(late, c(0,3), c(0,Inf), by = "union_type")
  expect_equal(sum(subset(counts, union_type != "unpartnered")$person_years), 6)
  schedule <- data.frame(union_type=c("cohab","marriage"),dur_start=c(0,1),dur_end=Inf,rate=c(0,1e5))
  model <- set_union_dissolution(model, schedule)
  result <- simulate(model, 3.5, 43)
  expect_gte(result$unions$end_time[2] - result$unions$start_time[2], 1)
})

test_that("conversion and cohab breakup follow competing annual hazards", {
  n <- 4000L
  result <- simulate(initial_cohabs(n), 2, 415)
  expected <- c(conversion=.3/.5*(1-exp(-1)), dissolution=.2/.5*(1-exp(-1)))
  reasons <- result$unions$end_reason[seq_len(n)]
  observed <- c(sum(reasons=="conversion",na.rm=TRUE),sum(reasons=="dissolution",na.rm=TRUE))
  expect_true(all(abs(observed-n*expected)<6*sqrt(n*expected*(1-expected))))
  expect_equal(result$simulation$conversions, observed[1])
  expect_equal(result$simulation$dissolutions, observed[2])
})

test_that("fertility follows type and conversion retains actual-birth spacing", {
  model <- initial_cohabs(conversion = 1e5, breakup = 0)
  model$processes$union_formation$union_convert_rate$dur_start <- 1
  model <- set_fertility(model, data.frame(age_start=15,age_end=50,
    union_type=c("unpartnered","cohab","marriage"),rate=c(0,20,0)),birth_interval=0)
  result <- simulate(model, 2, 517)
  kids <- subset(result$persons,birth_time>=0)
  expect_gt(nrow(kids), 0)
  expect_true(all(kids$father_id==2))
  expect_true(all(kids$birth_time<=result$unions$start_time[2]))
  expect_equal(nrow(result$unions), 2)
  # A shared partnered schedule remains valid with either type.
  model <- set_fertility(model,data.frame(age_start=15,age_end=50,
    union_type=c("unpartnered","partnered"),rate=c(0,.2)))
  expect_s3_class(simulate(model,1,42),"hammel_population")
})

test_that("cohabitation endings restore legal status by default and support the toggle", {
  for (counts in c(FALSE, TRUE)) {
    model <- initial_cohabs(conversion=0,breakup=1e5)
    model$processes$union_formation$cohabitation_counts_as_previous_union <- counts
    result <- simulate(model,.1,42)
    info <- previous_union_info_cpp(result$persons,result$unions,counts)
    expect_equal(info$previous_union,c(if(counts)1L else 0L,if(counts)1L else 0L))
    model <- initial_cohabs(conversion=0,breakup=0)
    model$processes$union_formation$cohabitation_counts_as_previous_union <- counts
    model <- set_mortality(model,data.frame(sex=c("f","m"),age_start=0,age_end=Inf,rate=c(1e5,0)))
    result <- simulate(model,.1,42)
    info <- previous_union_info_cpp(result$persons,result$unions,counts)
    expect_equal(info$previous_union[2],if(counts)2L else 0L)
    # The surviving man is available for a new woman immediately.
    model <- cohab_fixture(3)
    model$population$persons$sex <- c("f","m","f")
    model$population$unions <- initial_cohabs()$population$unions
    model$population$persons$birth_time[1] <- -80
    model <- set_mortality(model,data.frame(sex=c("f","f","m"),age_start=c(0,80,0),age_end=c(80,Inf,Inf),rate=c(0,1e5,0)))
    result <- simulate(model,1,42)
    expect_true(any(result$unions$wife_id==3 & result$unions$husband_id==2))
  }
})

test_that("type-specific group strata built by binding stratified subsets", {
  model <- cohab_fixture(40)
  model$population$persons$sex[1:20] <- "f"
  model$population$persons$group <- rep(c("A","B"),20)
  model$population$groups <- list(group = group())
  typed <- data.frame(age_start=15,age_end=50,union_type=c("unpartnered","cohab","marriage"),rate=c(.03,.15,.1))
  factors <- list(unpartnered=c(A=1,B=.8),cohab=c(A=1.2,B=.9),marriage=c(A=1,B=1.1))
  stratified <- do.call(rbind, lapply(names(factors), function(type)
    stratify_rates(typed[typed$union_type==type,],list(group=factors[[type]]))))
  implicit <- set_fertility(model,stratified,birth_interval=.1)
  explicit <- stratified
  explicit$rate <- typed$rate[match(explicit$union_type,typed$union_type)] *
    mapply(function(type,group) factors[[type]][[group]], explicit$union_type, explicit$group)
  reference <- set_fertility(model,as.data.frame(explicit),birth_interval=.1)
  expect_identical(simulate(implicit,2,617),simulate(reference,2,617))
  expect_error(stratify_rates(typed,list(union_type=c(cohab=1))),"already has")
})

test_that("cohabitation and marriage use distinct matching and diagnostic weights", {
  model <- cohab_fixture(600,.5)
  model$population$persons$sex[1:300] <- "f"
  model$population$persons$group <- rep(c("A","B"),300)
  model$population$groups <- list(group = group())
  model$processes$union_formation$matching <- matching_model(~1,numeric(),group_exogamy = NA)
  model$processes$union_formation$cohab_matching <- matching_model(~1,numeric(),group_exogamy=c("A:B"=10))
  model$processes$union_formation$diagnostics <- TRUE
  result <- simulate(model,.3,713)
  u <- get_detailed_unions(result)
  marriage <- subset(u,union_type=="marriage")
  cohab <- subset(u,union_type=="cohab")
  expect_gt(nrow(marriage),50)
  expect_gt(nrow(cohab),50)
  expect_true(all(marriage$group_w==marriage$group_h))
  expect_gt(mean(cohab$group_w!=cohab$group_h),.9)
  expect_true(all(is.finite(cohab$.log_weight_A)))
  expect_equal(anyDuplicated(u$husband_id),0L)
})

test_that("cohab probabilities select group and previous category and require coverage", {
  model <- cohab_fixture(200)
  model$population$persons$sex[1:100] <- "f"
  model$population$persons$group <- rep(c("A","B"),100)
  model$population$groups <- list(group = group())
  p <- expand.grid(group=c("A","B"),previous_union=c("none","separated","widowed"),stringsAsFactors=FALSE)
  p$age_start <- 15; p$age_end <- Inf
  p$prob <- as.numeric(p$group=="B" & p$previous_union=="none")
  model$processes$union_formation$cohab_prob <- p
  result <- simulate(model,.5,813)
  u <- get_detailed_unions(result)
  expect_true(all(u$union_type[u$group_w=="B"]=="cohab"))
  expect_true(all(u$union_type[u$group_w=="A"]=="marriage"))
  bad <- model; bad$processes$union_formation$cohab_prob <- p[-1,]
  expect_error(simulate(bad,0),"all three")
})

test_that("conversion age clocks and pending birth ties retain the conversion", {
  model <- initial_cohabs(conversion=0,breakup=0)
  model$processes$union_formation$union_convert_rate <- data.frame(age_start=21,age_end=Inf,rate=1e5)
  result <- simulate(model,2,413)
  expect_gte(result$unions$start_time[2],1)
  # Process clocks draw thresholds in refresh order: fertility, then (when its
  # rate is positive) dissolution, then conversion; zero-rate processes draw
  # none. With one unit per process each due time is its threshold over the
  # rate. Ties are placed at .5, away from her integer-age fertility cells.
  seed <- 71L
  set.seed(seed); draws <- rexp(3)
  model <- initial_cohabs(conversion=draws[2]/.5,breakup=0)
  model <- set_fertility(model,data.frame(age_start=15,age_end=50,
    union_type=c("unpartnered","cohab","marriage"),rate=c(0,draws[1]/.5,0)),birth_interval=0)
  tied <- simulate(model,1.5,seed)
  kids <- subset(tied$persons,birth_time>=0)
  expect_equal(nrow(kids),1)
  expect_equal(kids$birth_time,.5,tolerance=0)
  expect_equal(tied$unions$end_time[1],.5,tolerance=0)
  expect_equal(tied$unions$end_reason[1],"conversion")
  # Dissolution at the same time wins, preventing conversion.
  model$processes$union_dissolution$rates$rate[1] <- draws[2]/.5
  model$processes$union_formation$union_convert_rate$rate <- draws[3]/.5
  ended <- simulate(model,1.5,seed)
  expect_equal(ended$unions$end_reason[1],"dissolution")
  expect_equal(ended$simulation$conversions,0)
})

test_that("earlier legal history governs subsequent searches after cohabitant death", {
  for (counts in c(FALSE,TRUE)) {
    model <- cohab_fixture(3)
    model$population$persons$birth_time[3] <- -80
    model$population$unions <- tibble::tibble(union_id=1:2,wife_id=1,husband_id=2:3,
      start_time=c(-10,-2),end_time=c(-5,NA),end_reason=c("dissolution",NA),union_type=c("marriage","cohab"))
    model$processes$union_formation$rates$rate <- 0
    model$processes$union_formation$cohabitation_counts_as_previous_union <- counts
    model <- set_mortality(model,data.frame(sex=c("f","m","m"),age_start=c(0,0,80),
      age_end=c(Inf,80,Inf),rate=c(0,0,1e5)))
    result <- simulate(model,.1,42)
    info <- previous_union_info_cpp(result$persons,result$unions,counts)
    expect_equal(info$previous_union[1],if(counts)2L else 1L)
    expect_equal(info$legal_previous[1],1L)
    model <- set_population(model,result)
    model$processes$union_formation$rates <- data.frame(age_start=15,age_end=Inf,
      previous_union=c("none","separated","widowed"),rate=c(0,1e5,0))
    continued <- simulate(model,.1,43)
    expect_equal(nrow(continued$unions),if(counts)2L else 3L)
  }
})

test_that("typed dissolution retains symmetric pair hazards", {
  model <- initial_cohabs(3,conversion=0,breakup=0)
  model$population$persons$group <- c("A","B","A","A","B","B")
  model$population$groups <- list(group = group())
  rates <- expand.grid(group_w=c("A","B"),group_h=c("A","B"),union_type=c("cohab","marriage"),stringsAsFactors=FALSE)
  rates$dur_start <- 0; rates$dur_end <- Inf
  rates$rate <- ifelse(rates$union_type=="cohab" & rates$group_w!=rates$group_h,1e5,0)
  model <- set_union_dissolution(model,rates)
  result <- simulate(model,.1,42)
  expect_equal(result$simulation$dissolutions,1)
  expect_equal(result$unions$end_reason[1],"dissolution")
  model <- set_union_dissolution(model,transform(rates,rate=ifelse(group_w!=group_h,0,rate)))
  expect_equal(simulate(model,.1,42)$simulation$dissolutions,0)
})

test_that("union entry types share exposure and exclude conversion", {
  history <- function() hammel_population(
    data.frame(person_id = 1:4, sex = c("f", "m", "m", "m"),
      birth_time = -20, death_time = c(NA, NA, 4, NA),
      mother_id = NA_integer_, father_id = NA_integer_, group = c("A", "B", "A", "B")),
    data.frame(union_id = 1:3, wife_id = 1, husband_id = 2:4,
      start_time = c(1, 3, 5), end_time = c(2, 4, NA),
      end_reason = c("dissolution", "death_husband", NA)), end_time = 6, groups = "group")
  pop <- history()
  pop$unions$union_type <- c("cohab", "marriage", "cohab")
  for (previous in c(FALSE, TRUE)) {
    for (duration in list(NULL, c(0, 1, Inf))) {
      pooled <- calculate_union_rates(pop, c(0, 6), c(0, Inf), dur_breaks = duration, by = if (isTRUE(previous)) "previous_union")
      typed <- calculate_union_rates(pop, c(0, 6), c(0, Inf), dur_breaks = duration, by = c(if (isTRUE(previous)) "previous_union", "union_type"))
      cohab <- subset(typed, union_type == "cohab")
      marriage <- subset(typed, union_type == "marriage")
      expect_equal(cohab$person_years, pooled$person_years)
      expect_equal(marriage$person_years, pooled$person_years)
      expect_equal(cohab$unions + marriage$unions, pooled$unions)
      expect_equal(cohab$rate_union + marriage$rate_union, pooled$rate_union)
    }
  }
  grouped <- calculate_union_rates(pop, c(0,6), c(0,Inf), by = c("group", "previous_union", "union_type"))
  expect_equal(sum(grouped$unions), 6)
  expect_equal(sum(grouped$person_years[grouped$union_type == "cohab"]),
    sum(grouped$person_years[grouped$union_type == "marriage"]))
  expect_error(calculate_union_rates(history(), c(0,6), c(0,Inf),
    by = "union_type"), "typed")
  expect_error(calculate_union_rates(pop, c(0,6), c(0,Inf), by = NA), "`by` must")
  converted <- simulate(initial_cohabs(conversion = 100, breakup = 0), 1, 42)
  typed <- calculate_union_rates(converted, c(0,1), c(0,Inf), by = "union_type")
  expect_equal(sum(typed$unions), 2) # Existing cohab starts at the window boundary.
  expect_equal(sum(typed$unions[typed$union_type == "marriage"]), 0)
})
