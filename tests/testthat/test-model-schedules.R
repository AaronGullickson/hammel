model_calendar <- function(first, second = first, boundary = .5, start = 0, end = Inf) {
  tibble::tibble(time_start = c(start, boundary), time_end = c(boundary, end),
    model = list(first, second))
}
model_choice_fixture <- function(role = "matching", specification, static_general = NULL,
                                 birth = FALSE, on_demand = FALSE) {
  m <- make_rate_model(3, search_rate = if (birth) 0 else 1)
  m$population$persons$birth_time <- c(-20, -20, -40)
  if (birth) {
    r <- data.frame(age_start = 15, age_end = 50, rate = 20)
    if (!on_demand) r <- rbind(transform(r, union_type = "unpartnered"),
      transform(r, union_type = "partnered"))
    m <- set_fertility(m, r, birth_interval = 2)
  }
  args <- list(model = m, initiation = if (on_demand) "on_demand" else "rate",
    matching = if (is.null(static_general)) create_age_gap_matching_model(0, .1) else static_general,
    cancel_birth_if_no_match = FALSE, kinship_generations = 0)
  if (!on_demand) args$rates <- m$processes$union_formation$rates
  args[role] <- list(specification)
  if (role == "cohab_matching") {
    args$cohab_prob <- data.frame(age_start = 0, age_end = Inf, prob = 1)
    args$union_convert_rate <- data.frame(dur_start = 0, dur_end = Inf, rate = 0)
  }
  do.call(set_union_formation, args)
}
model_inheritance_fixture <- function(inheritance, extra = NULL) {
  m <- make_model(2, sex = "f", birth_rate = 100, birth_interval = 0)
  m$population$persons$sex <- c("f", "m")
  m$population$persons$group <- c("A", "B")
  m$population$groups <- list(group = group(ancestry = TRUE))
  m$population$unions <- tibble::tibble(union_id = 1, wife_id = 1, husband_id = 2,
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  set_inheritance(set_fertility(m, m$processes$fertility$rates, birth_interval = 0, ffab = 1),
                  inheritance, extra = extra)
}

test_that("model arguments accept direct objects and inspectable model-list tibbles", {
  matching <- create_age_gap_matching_model(0, 5)
  inheritance <- create_group_inheritance("mother")
  for (role in c("matching", "father_matching", "cohab_matching")) {
    p <- model_calendar(matching)
    attr(p, "source") <- list(example = "synthetic model calendar")
    m <- model_choice_fixture(role, as.data.frame(p))
    expect_s3_class(m$processes$union_formation[[role]], "tbl_df")
    expect_identical(m$processes$union_formation[[role]]$model[[1]], matching)
    expect_s3_class(model_choice_fixture(role, matching), "hammel_model")
  }
  p <- model_calendar(inheritance)
  m <- model_inheritance_fixture(as.data.frame(p))
  expect_s3_class(m$groups$group$inheritance, "tbl_df")
  expect_identical(m$groups$group$inheritance$model[[1]], inheritance)
  expect_s3_class(model_inheritance_fixture("mother"), "hammel_model")
})

test_that("all scheduled models and interval bounds validate before execution", {
  mm <- create_age_gap_matching_model(0, 5)
  s <- model_calendar(mm)
  bad_entries <- list(NULL, "mother", create_group_inheritance("mother"), s, list())
  for (entry in bad_entries) {
    bad <- s; bad$model[2] <- list(entry)
    expect_error(model_choice_fixture("matching", bad), "model")
  }
  bad_intervals <- list(s[setdiff(names(s), "model")], s[setdiff(names(s), "time_end")],
    transform(s, time_start = NA_real_), transform(s, time_start = Inf),
    transform(s, time_end = time_start), transform(s, time_end = NA_real_),
    transform(s, time_start = as.character(time_start)), s[FALSE, ], rbind(s, s[1, ]))
  for (bad in bad_intervals) expect_error(model_choice_fixture("matching", bad))
  bad <- s; bad$model[[2]]$coefficients[1] <- NA_real_
  expect_error(model_choice_fixture("matching", bad), "matching\\$model\\[\\[2\\]\\]")
  bad <- model_calendar(create_group_inheritance("mother")); bad$model[2] <- list(mm)
  expect_error(model_inheritance_fixture(bad), "inheritance\\$model")
  bad$model[2] <- list("father")
  expect_error(model_inheritance_fixture(bad), "complete model")
  edited <- model_choice_fixture("matching", s)
  edited$processes$union_formation$matching$model[[2]]$coefficients[1] <- NA_real_
  expect_error(simulate(edited, 0), "matching\\$model")
})

test_that("coverage applies to all four arguments and zero-year starts", {
  mm <- create_age_gap_matching_model(0, 5)
  for (role in c("matching", "father_matching", "cohab_matching", "inheritance")) {
    source <- if (role == "inheritance") create_group_inheritance("mother") else mm
    s <- model_calendar(source, boundary = .5, start = .1, end = 1)
    m <- if (role == "inheritance") model_inheritance_fixture(s) else model_choice_fixture(role, s)
    expect_error(simulate(m, 0), paste0(role, ".*coverage gap"))
    m$population$end_time <- .1
    expect_s3_class(simulate(m, .9, seed = 3), "hammel_population")
    expect_error(simulate(m, 1, seed = 3), paste0(role, ".*coverage gap"))
    s$time_start[1] <- 0
    s$time_start[2] <- .7
    m <- if (role == "inheritance") model_inheritance_fixture(s) else model_choice_fixture(role, s)
    expect_error(simulate(m, 1, seed = 3), "coverage gap")
    expect_s3_class(simulate(m, .4, seed = 3), "hammel_population")
  }
})

test_that("general and cohabitation matching use the new model at a tied pending search", {
  first <- create_age_gap_matching_model(0, .1)
  second <- create_age_gap_matching_model(20, .1)
  # The first draw is the search clock's (zero death and fertility rates draw
  # none); one searching woman at rate 1 makes it her search time.
  seed <- 783; set.seed(seed); boundary <- rexp(1)
  for (role in c("matching", "cohab_matching")) {
    m <- model_choice_fixture(role, model_calendar(first, second, boundary))
    out <- simulate(m, boundary + .01, seed = seed)
    expect_equal(nrow(out$unions), 1)
    expect_equal(out$unions$start_time, boundary, tolerance = 0)
    expect_equal(out$unions$husband_id, 3)
    expect_equal(nrow(simulate(m, boundary, seed = seed)$unions), 0)
  }
})

test_that("father and on-demand matching select preferences at birth execution", {
  first <- create_age_gap_matching_model(0, .1)
  second <- create_age_gap_matching_model(20, .1)
  seed <- 183; set.seed(seed); boundary <- rexp(1) / 20
  for (on_demand in c(FALSE, TRUE)) {
    role <- if (on_demand) "matching" else "father_matching"
    m <- model_choice_fixture(role, model_calendar(first, second, boundary),
      birth = TRUE, on_demand = on_demand)
    out <- simulate(m, boundary + .01, seed = seed)
    born <- out$persons[out$persons$birth_time >= 0, ]
    expect_equal(nrow(born), 1)
    expect_equal(born$birth_time, boundary, tolerance = 0)
    expect_equal(born$father_id, 3)
  }
})

test_that("NULL father and cohabitation matching follow general calendars dynamically", {
  first <- create_age_gap_matching_model(0, .1)
  second <- create_age_gap_matching_model(20, .1)
  for (role in c("father_matching", "cohab_matching")) {
    birth <- role == "father_matching"
    # The first draw is the fertility clock's (rate 20) or the search clock's (rate 1).
    seed <- 183; set.seed(seed); boundary <- if (birth) rexp(1) / 20 else rexp(1)
    general <- model_calendar(first, second, boundary)
    fallback <- model_choice_fixture(role, NULL, general, birth = birth)
    fixed <- model_choice_fixture(role, first, general, birth = birth)
    a <- simulate(fallback, boundary + .01, seed = seed)
    b <- simulate(fixed, boundary + .01, seed = seed)
    if (birth) {
      expect_equal(tail(a$persons$father_id, 1), 3)
      expect_equal(tail(b$persons$father_id, 1), 2)
    } else {
      expect_equal(a$unions$husband_id, 3)
      expect_equal(b$unions$husband_id, 2)
    }
  }
})

test_that("inheritance switches newborn identity without relabeling older people or ancestry", {
  inheritance <- model_calendar(create_group_inheritance("mother"),
    create_group_inheritance("father"), boundary = .25)
  m <- model_inheritance_fixture(inheritance)
  out <- simulate(m, .5, seed = 981)
  kids <- out$persons[out$persons$birth_time >= 0, ]
  expect_gt(sum(kids$birth_time < .25), 5)
  expect_gt(sum(kids$birth_time >= .25), 5)
  expect_true(all(kids$group == ifelse(kids$birth_time < .25, "A", "B")))
  expect_equal(out$persons$group[1:2], c("A", "B"))
  expect_equal(out$ancestry$group$A[-c(1,2)], rep(.5, nrow(kids)))
  expect_equal(out$ancestry$group$B[-c(1,2)], rep(.5, nrow(kids)))
  expect_equal(out$unions[, names(m$population$unions)], m$population$unions)
})

test_that("inheritance applies the new model at an exact tied birth", {
  seed <- 981; set.seed(seed); boundary <- rexp(1) / 100
  s <- model_calendar(create_group_inheritance("mother"), create_group_inheritance("father"), boundary)
  out <- simulate(model_inheritance_fixture(s), boundary + 1e-7, seed = seed)
  expect_equal(out$simulation$births, 1)
  expect_equal(tail(out$persons$birth_time, 1), boundary, tolerance = 0)
  expect_equal(tail(out$persons$group, 1), "B")
})

test_that("future custom identities use a common registry and require complete model coverage", {
  late <- create_group_inheritance("custom", fun = function(...) "C")
  s <- model_calendar(create_group_inheritance("mother"), late, .25)
  m <- model_inheritance_fixture(s, extra = "C")
  out <- simulate(m, .5, seed = 32)
  kids <- out$persons[out$persons$birth_time >= 0, ]
  expect_true(all(kids$group == ifelse(kids$birth_time < .25, "A", "C")))
  expect_equal(out$ancestry$group$A[-c(1,2)], rep(.5, nrow(kids)))
  bad <- m
  rates <- bad$processes$mortality$rates
  bad <- set_mortality(bad, rbind(transform(rates, group = "A"), transform(rates, group = "B")))
  expect_error(simulate(bad, .1), "missing schedules.*C")
  bad <- m
  bad$processes$union_formation$matching <- create_matching_model(~ female_group,
    c(female_groupB = 0), group_levels = c("A", "B"))
  expect_error(simulate(bad, .1), "group_levels")
})

test_that("inheritance boundaries do not call custom functions without a birth", {
  fun <- create_group_inheritance("custom", fun = function(...) stop("must not execute"))
  m <- model_inheritance_fixture(model_calendar(create_group_inheritance("mother"), fun, .1))
  m$processes$fertility$rates$rate <- 0
  expect_s3_class(simulate(m, .5, seed = 1), "hammel_population")
})

test_that("weighted calibration switches from hard within-group matching and back", {
  early <- create_age_gap_matching_model(0, 5, prevent_exogamy = TRUE)
  late <- create_matching_model(~ 1, numeric(), delta = c(A = log(3)))
  for (reverse in c(FALSE, TRUE)) {
    before <- if (reverse) late else early
    after <- if (reverse) early else late
    m <- make_rate_model(400, search_rate = .4)
    m$population$persons$sex <- rep(c("f", "m"), each = 200)
    m$population$persons$group <- rep(c("A", "B"), 200)
    m$population$groups <- list(group = group())
    m$processes$union_formation$diagnostics <- TRUE
    m$processes$union_formation$matching <- model_calendar(before, after, 1e-6)
    reference <- m; reference$processes$union_formation$matching <- after
    expected <- simulate(reference, .5, seed = 921)
    expect_true(all(expected$unions$start_time > 1e-6))
    expect_gt(nrow(expected$unions), 10)
    actual <- simulate(m, .5, seed = 921)
    expect_equal(actual$unions, expected$unions, tolerance = 1e-9)
  }
})

test_that("matching formulas can switch between native and R formula paths", {
  first <- create_age_gap_matching_model(0, .1)
  second <- create_matching_model(~ abs(age_difference), c("abs(age_difference)" = 100))
  # The first draw is the search clock's (zero death and fertility rates draw
  # none); one searching woman at rate 1 makes it her search time.
  seed <- 783; set.seed(seed); boundary <- rexp(1)
  out <- simulate(model_choice_fixture("matching", model_calendar(first, second, boundary)),
    boundary + .01, seed = seed)
  expect_equal(out$unions$husband_id, 3)
})

test_that("scheduled models preserve RNG isolation and handle zero and empty runs", {
  mm <- create_age_gap_matching_model(0, 5)
  m <- model_choice_fixture("matching", model_calendar(mm))
  set.seed(72); old <- .Random.seed
  out <- simulate(m, 1, seed = 83)
  expect_identical(.Random.seed, old)
  expect_identical(simulate(m, 1, seed = 83), out)
  expect_equal(nrow(simulate(m, 0, seed = 83)$unions), 0)
  empty <- make_model(0)
  empty$processes$union_formation$matching <- model_calendar(mm)
  expect_equal(nrow(simulate(empty, 1, seed = 83)$persons), 0)
  expect_equal(nrow(simulate(make_model(2), .5, seed = 1)$persons), 2)
})

test_that("inheritance schedules update hierarchy ranks and passing parameters", {
  early <- create_group_inheritance("hierarchy", order = c("A", "B"), passing = 0)
  for (late in list(create_group_inheritance("hierarchy", order = c("B", "A")),
      create_group_inheritance("hierarchy", order = c("A", "B"), passing = 1))) {
    m <- model_inheritance_fixture(model_calendar(early, late, .25))
    out <- simulate(m, .5, seed = 981)
    kids <- out$persons[out$persons$birth_time >= 0, ]
    expect_true(all(kids$group == ifelse(kids$birth_time < .25, "A", "B")))
  }
})

test_that("different custom inheritance intervals share identity indices and callback RNG", {
  for (labels in list(c("A", "B", "C", "D"), c(10L, 20L, 30L, 40L))) {
    first <- create_group_inheritance("custom",
      fun = local({value <- labels[3]; function(...) value}))
    second <- create_group_inheritance("custom",
      fun = local({value <- labels[4]; function(...) value}))
    m <- model_inheritance_fixture(model_calendar(first, second, .25))
    m$population$persons$group <- labels[1:2]
    m <- set_groups(m, group = group(categories = labels, inheritance = model_calendar(first, second, .25)))
    out <- simulate(m, .5, seed = 981)
    kids <- out$persons[out$persons$birth_time >= 0, ]
    expect_equal(kids$group, ifelse(kids$birth_time < .25, labels[3], labels[4]))
    expect_identical(typeof(kids$group), typeof(labels))
  }
  late <- create_group_inheritance("custom",
    fun = function(...) sample(c("A", "C"), 1))
  m <- model_inheritance_fixture(model_calendar(create_group_inheritance("mother"), late, .1), extra = "C")
  set.seed(12); before <- .Random.seed
  a <- simulate(m, .5, seed = 51)
  expect_identical(.Random.seed, before)
  expect_identical(simulate(m, .5, seed = 51), a)
  kids <- a$persons[a$persons$birth_time >= .1, ]
  expect_setequal(unique(kids$group), c("A", "C"))
})

test_that("constant schedules preserve direct-model behavior and existing unions", {
  mm <- create_age_gap_matching_model(0, 5)
  direct <- model_choice_fixture("matching", mm)
  scheduled <- model_choice_fixture("matching", model_calendar(mm))
  expect_identical(simulate(scheduled, 1, seed = 3), simulate(direct, 1, seed = 3))
  m <- model_inheritance_fixture(create_group_inheritance("mother"))
  m$processes$fertility$rates$rate <- 0
  m$processes$union_formation$matching <- model_calendar(mm,
    create_age_gap_matching_model(0, 5, prevent_exogamy = TRUE), .1)
  out <- simulate(m, .5, seed = 3)
  expect_equal(out$unions[, names(m$population$unions)], m$population$unions)
})

test_that("restarts select the active inheritance model and preserve earlier assignments", {
  m <- model_inheritance_fixture(model_calendar(create_group_inheritance("mother"),
    create_group_inheritance("father"), .25))
  early <- simulate(m, .25, seed = 981)
  out <- simulate(set_population(m, early), .25, seed = 982)
  old <- nrow(early$persons)
  expect_equal(out$persons$group[seq_len(old)], early$persons$group)
  expect_gt(nrow(out$persons), old)
  expect_true(all(out$persons$group[-seq_len(old)] == "B"))
})

test_that("model choices update before shared rate and probability changes", {
  early <- create_age_gap_matching_model(0, .1)
  late <- create_age_gap_matching_model(20, .1)
  m <- model_choice_fixture("cohab_matching", model_calendar(early, late, .3))
  m$processes$union_formation$rates <- data.frame(age_start = 0, age_end = Inf,
    rate = c(0, 1e5), time_start = c(0, .3), time_end = c(.3, Inf))
  m$processes$union_formation$cohab_prob <- data.frame(age_start = 0, age_end = Inf,
    prob = c(0, 1), time_start = c(0, .3), time_end = c(.3, Inf))
  out <- simulate(m, .31, seed = 23)
  expect_equal(out$unions$union_type, "cohab")
  expect_equal(out$unions$husband_id, 3)
  expect_gte(out$unions$start_time, .3)
})

test_that("inheritance-only boundaries preserve matching execution without births", {
  m <- make_rate_model(100,search_rate=.4)
  m$population$persons$sex <- rep(c("f","m"),each=50)
  m$population$persons$group <- rep(c("A","B"),50)
  m$population$groups <- list(group = group())
  m <- set_inheritance(set_fertility(m,m$processes$fertility$rates), create_group_inheritance("mother"))
  scheduled <- set_inheritance(set_fertility(m,m$processes$fertility$rates), model_calendar(create_group_inheritance("mother"),
      create_group_inheritance("father"),boundary=.25))
  actual <- simulate(scheduled,2,seed=812)
  expected <- simulate(m,2,seed=812)
  expect_equal(actual$persons,expected$persons)
  expect_equal(actual$unions,expected$unions)
  expect_equal(actual$simulation$births,0)
})
