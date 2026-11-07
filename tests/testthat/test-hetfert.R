het_model <- function(n = 100L, specification = het_fert(.5, .4), rate = .2,
                      interval = .75, ffab = 0) {
  p <- tibble::tibble(person_id = seq_len(n), sex = rep("f", n),
    birth_time = rep(-20.1, n), death_time = rep(NA_real_, n),
    mother_id = rep(NA_real_, n), father_id = rep(NA_real_, n))
  hammel_model() |> set_population(hammel_population(p, end_time = 0, groups = if ("group" %in% names(p)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = rate),
      hetfert = specification, birth_interval = interval, ffab = ffab) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
}

test_that("het_fert validates parameters and edited process specifications", {
  expect_s3_class(het_fert(.5, .4), "hammel_het_fert")
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "1"))
    expect_error(het_fert(bad, .4), "cv")
  for (bad in list(-.1, 1.1, NA_real_, Inf, c(0, 1), "1"))
    expect_error(het_fert(.5, bad), "alpha")
  expect_error(het_model(specification = list(cv = .5, alpha = .4)), "het_fert")
  m <- het_model()
  m$processes$fertility$hetfert$alpha <- 2
  expect_error(simulate(m, 1, 1), "alpha")
  set.seed(1); before <- .Random.seed
  het_fert(.5, .5)
  expect_identical(.Random.seed, before)
})

test_that("founder draws have the specified moments, obey seed policy and include zero runs", {
  m <- het_model(20000L)
  set.seed(15); before <- .Random.seed
  p <- simulate(m, 0, 42)
  expect_identical(.Random.seed, before)
  expect_identical(p, simulate(m, 0, 42))
  expect_false("hetfert" %in% names(m$population$persons))
  expect_equal(mean(p$persons$hetfert), 1, tolerance = .015)
  expect_equal(sd(p$persons$hetfert), .5, tolerance = .025)
  expect_true(all(is.finite(p$persons$hetfert) & p$persons$hetfert > 0))
  expect_identical(p$simulation$hetfert, het_fert(.5, .4))
  expect_s3_class(p$persons, "tbl_df")
  set.seed(15); simulate(m, 0)
  expect_false(identical(.Random.seed, before))
  empty <- simulate(het_model(0), 0, 42)
  expect_type(empty$persons$hetfert, "double")
  expect_equal(nrow(empty$persons), 0)
})

test_that("supplied female weights are authoritative and male weights unused", {
  m <- het_model(3)
  m$population$persons$sex[3] <- "m"
  m$population$persons$hetfert <- c(.25, 4, NA)
  p <- simulate(m, 0, 42)
  expect_identical(p$persons$hetfert, c(.25, 4, NA_real_))
  for (bad in c(NA, 0, -1, Inf)) {
    edited <- m
    edited$population$persons$hetfert[1] <- bad
    expect_error(simulate(edited, 0, 42), "strictly positive")
  }
  m$population$persons$hetfert <- c("1", "2", NA)
  expect_error(simulate(m, 0, 42), "numeric column")
  m <- het_model(3, specification = NULL)
  m$population$persons$hetfert <- c(NA, 0, -1)
  expect_false("hetfert" %in% names(simulate(m, 0, 42)$persons))
})

test_that("zero CV retains uniform histories, and NULL removes custom weights", {
  m <- het_model(150, specification = NULL, interval = 0)
  plain <- simulate(m, 3, 19)
  weighted <- m
  weighted$processes$fertility$hetfert <- het_fert(0, .5)
  p <- simulate(weighted, 3, 19)
  expect_true(all(p$persons$hetfert[p$persons$sex == "f"] == 1))
  expect_true(all(is.na(p$persons$hetfert[p$persons$sex == "m"])))
  p$persons$hetfert <- NULL
  p$simulation$hetfert <- NULL
  expect_identical(p, plain)
  m$population$persons$hetfert <- rep(100, 150)
  p <- simulate(m, 3, 19)
  p$persons$hetfert <- NULL
  expect_identical(p, plain)
})

test_that("sampling uses relative weights while keeping fixed cell intensity", {
  m <- het_model(2, het_fert(.5, .4), rate = 10, interval = 0)
  m$population$persons$hetfert <- c(1, 9)
  p <- simulate(m, 20, 51)
  children <- p$persons$mother_id[!is.na(p$persons$mother_id)]
  n <- length(children)
  expect_lt(abs(n - 400), 4 * sqrt(400))
  expect_lt(abs(sum(children == 2) - .9 * n), 4 * sqrt(n * .9 * .1))
  expect_true(all(is.na(p$persons$hetfert[p$persons$sex == "m"])))
})

test_that("daughters use actual maternal weights with exact endpoint inheritance", {
  m <- het_model(150, het_fert(.7, 1), rate = .8, interval = 0, ffab = 1)
  m$population$persons$hetfert <- rep(c(.25, 4), 75)
  p <- simulate(m, 2, 42)
  children <- !is.na(p$persons$mother_id)
  mothers <- match(p$persons$mother_id[children], p$persons$person_id)
  expect_gt(sum(children), 100)
  expect_identical(p$persons$hetfert[children], p$persons$hetfert[mothers])
  m$processes$fertility$hetfert <- het_fert(0, .5)
  p <- simulate(m, 1, 42)
  children <- !is.na(p$persons$mother_id)
  mothers <- match(p$persons$mother_id[children], p$persons$person_id)
  expect_equal(p$persons$hetfert[children], sqrt(p$persons$hetfert[mothers]))
})

test_that("newborn innovation has the specified normal-scale distribution", {
  for (alpha in c(0, .6)) {
    m <- het_model(6000, het_fert(.5, alpha), rate = .6, interval = 0, ffab = 1)
    m$population$persons$hetfert <- rep(c(.25, 4), 3000)
    p <- simulate(m, 1, 62)
    children <- !is.na(p$persons$mother_id)
    mothers <- match(p$persons$mother_id[children], p$persons$person_id)
    variance <- log1p(.5^2); mu <- -variance / 2
    innovation <- (log(p$persons$hetfert[children]) - mu -
      alpha * (log(p$persons$hetfert[mothers]) - mu)) / sqrt(variance * (1 - alpha^2))
    expect_lt(abs(mean(innovation)), 4 / sqrt(length(innovation)))
    expect_lt(abs(var(innovation) - 1), 4 * sqrt(2 / (length(innovation) - 1)))
    # Reproducing mothers are selected: checking innovation does not claim a
    # stationary distribution among daughters at intermediate alpha.
  }
})

test_that("weights survive restarts and spacing removes/reintroduces weighted members", {
  old <- options(hammel.check_cells = TRUE); on.exit(options(old), add = TRUE)
  m <- het_model(100, het_fert(1, .8), rate = .8, interval = .75)
  first <- simulate(m, 3, 71)
  again <- simulate(set_population(m, first), 2, 72)
  expect_identical(head(again$persons$hetfert, nrow(first$persons)), first$persons$hetfert)
  times <- split(again$persons$birth_time[!is.na(again$persons$mother_id)],
    again$persons$mother_id[!is.na(again$persons$mother_id)])
  expect_true(all(unlist(lapply(times, function(t) diff(sort(t)))) >= .75 - 1e-12))
  expect_identical(first, simulate(m, 3, 71))
  # Parent sex histories, rather than survival, control retention of values.
  first$persons$death_time[1] <- first$end_time
  expect_identical(simulate(set_population(m, first), 0, 1)$persons$hetfert, first$persons$hetfert)
})

test_that("weighted hazards integrate with grouped partner/father matching and union types", {
  old <- options(hammel.check_cells = TRUE); on.exit(options(old), add = TRUE)
  set.seed(81)
  pop <- create_baseline_population(1000, max_age = 45)
  pop$persons$group <- rep(c("A", "B"), 500)
  pop$groups <- list(group = group())
  rates <- expand.grid(union_type = c("unpartnered", "cohab", "marriage"),
    parity_start = 0:1, stringsAsFactors = FALSE)
  rates$parity_end <- ifelse(rates$parity_start == 0, 1, Inf)
  rates$age_start <- 15; rates$age_end <- 45
  rates$rate <- .15
  model <- hammel_model() |> set_population(pop) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .01)) |>
    set_fertility(stratify_rates(rates, list(group = c(A = 1, B = 1.4))), hetfert = het_fert(1, .5)) |>
    set_union_formation("rate", matching = create_age_gap_matching_model(0, 5,
      group_exogamy = c("A:B" = -1)),
      rates = data.frame(age_start = 15, age_end = Inf, rate = .3),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = .5),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .3)) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .2))
  p <- simulate(model, 8, 82)
  expect_gt(p$simulation$births, 0)
  expect_gt(p$simulation$conversions, 0)
  expect_gt(p$simulation$dissolutions, 0)
  expect_true(all(is.finite(p$persons$hetfert[p$persons$sex == "f"])))
  # The fertility-driven on-demand channel uses the same weighted hazard sums.
  overall <- rates[rates$union_type == "unpartnered", ]
  overall$union_type <- NULL
  model <- model |> set_fertility(overall, hetfert = het_fert(.7, .4)) |>
    set_union_formation("on_demand", matching = create_age_gap_matching_model(0, 5,
      group_exogamy = c("A:B" = -1)))
  expect_gt(simulate(model, 5, 83)$simulation$births, 0)
})

test_that("removing a dominant eligible weight preserves tiny remaining weights", {
  old <- options(hammel.check_cells = TRUE); on.exit(options(old), add = TRUE)
  m <- het_model(2, het_fert(.5, .5), rate = 10, interval = .75)
  m$population$persons$hetfert <- c(1e100, 1e-100)
  p <- simulate(m, 3, 91)
  # The small-weight mother must become selectable while the other is in spacing.
  expect_gt(sum(p$persons$mother_id %in% 1), 1)
  expect_gt(sum(p$persons$mother_id %in% 2), 1)
  for (id in 1:2) {
    times <- p$persons$birth_time[p$persons$mother_id %in% id]
    expect_true(all(diff(sort(times)) >= .75 - 1e-12))
  }
})

test_that("failed generated-weight initialization restores a supplied seed", {
  m <- het_model(100, het_fert(1e308, .5))
  set.seed(10); before <- .Random.seed
  expect_error(simulate(m, 0, 42), "strictly positive")
  expect_identical(.Random.seed, before)
})

test_that("disabled runs remove weights and re-enabling redraws all women", {
  enabled <- het_model(40, het_fert(.5, .4), rate = .8, interval = 0, ffab = 1)
  first <- simulate(enabled, 1, 101)
  disabled <- set_population(enabled, first) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .8),
      birth_interval = 0, ffab = 1)
  zero <- simulate(disabled, 0, 102)
  expect_false("hetfert" %in% names(zero$persons))
  expect_null(zero$simulation$hetfert)
  second <- simulate(disabled, 1, 102)
  expect_gt(nrow(second$persons), nrow(first$persons))
  expect_false("hetfert" %in% names(second$persons))
  expect_true("hetfert" %in% names(first$persons))
  expect_identical(disabled$population$persons$hetfert, first$persons$hetfert)
  restarted <- simulate(set_population(enabled, second), 0, 103)
  expect_true(all(is.finite(restarted$persons$hetfert)))
  expect_true(all(restarted$persons$hetfert > 0))
  expect_false(identical(head(restarted$persons$hetfert, nrow(first$persons)),
    first$persons$hetfert))
  empty <- het_model(0, specification = NULL)
  empty$population$persons$hetfert <- numeric()
  expect_false("hetfert" %in% names(simulate(empty, 0, 1)$persons))
})

test_that("weight validation handles missing sex and logical NA placeholders clearly", {
  p <- het_model(2)$population$persons
  p$hetfert <- c(1, 2); p$sex[1] <- NA_character_
  expect_error(validate_fertility_weights(p), "nonmissing.*sex")
  p$sex[] <- "m"; p$hetfert <- c(NA, NA)
  expect_no_error(validate_fertility_weights(p))
  m <- set_population(het_model(2), hammel_population(p, end_time = 0, groups = if ("group" %in% names(p)) "group"))
  out <- simulate(m, 0, 1)
  expect_identical(out$persons$hetfert, c(NA_real_, NA_real_))
  p$sex[1] <- "f"
  expect_error(validate_fertility_weights(p), "Remove the column")
  p <- p[FALSE, ]
  expect_no_error(validate_fertility_weights(p))
})

test_that("daughters' weights follow the lognormal spread for cv above one", {
  # alpha = 0: every daughter draws afresh in the engine, whose variance uses
  # the large-cv form log1p(cv^2) = 2 log(cv) + log1p(1 / cv^2).
  m <- het_model(2000, het_fert(2, 0), rate = .5, interval = 0, ffab = 1)
  p <- simulate(m, 2, seed = 91)$persons
  w <- p$hetfert[!is.na(p$mother_id)]
  expect_gt(length(w), 1000)
  expect_true(all(is.finite(w) & w > 0))
  target <- sqrt(log1p(2^2))
  expect_lt(abs(sd(log(w)) - target), 4 * target / sqrt(2 * length(w)))
  expect_lt(abs(mean(log(w)) + target^2 / 2), 4 * target / sqrt(length(w)))
})
