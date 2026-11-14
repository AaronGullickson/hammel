strata_mortality_rates <- function() {
  data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
}

strata_model <- function() {
  model <- make_rate_model(6, search_rate = 0.4, birth_interval = 0.75)
  model$population$persons$sex <- rep(c("f", "m"), 3)
  model$population$persons$group <- c("A", "B", "C", "A", "B", "C")
  model$population$groups <- list(group = group())
  model$processes$fertility$rates$rate <- c(0.2, 0.1)
  model
}

# Build the same table by hand, row block by category.
explicit_strata <- function(rates, factors) {
  out <- do.call(rbind, lapply(names(factors), function(g) {
    transform(as.data.frame(rates), rate = rate * factors[[g]], group = g)
  }))
  rownames(out) <- NULL
  out
}

test_that("stratify_rates copies rows per category and multiplies rates", {
  rates <- strata_mortality_rates()
  attr(rates, "source") <- list(example = "synthetic")
  out <- stratify_rates(rates, list(group = c(A = 1.2, B = 1)))
  expect_s3_class(out, "tbl_df")
  expect_identical(out$group, c("A", "A", "B", "B"))
  expect_equal(out$rate, c(0.012, 0.012, 0.01, 0.01))
  expect_identical(out$sex, c("m", "f", "m", "f"))
  expect_identical(attr(out, "source"), attr(rates, "source"))
  # Crossing multiplies factors, with the first column varying slowest.
  crossed <- stratify_rates(rates[1, ], list(group = c(A = 2, B = 1), union_type = c(u = 1, p = 3)))
  expect_identical(crossed$group, c("A", "A", "B", "B"))
  expect_identical(crossed$union_type, c("u", "p", "u", "p"))
  expect_equal(crossed$rate, 0.01 * c(2, 6, 1, 3))
  # Already-stratified tables can be stratified further; calendar bands are kept.
  timed <- transform(rates, time_start = 0, time_end = Inf)
  twice <- stratify_rates(stratify_rates(timed, list(group = c(A = 1))), list(union_type = c(u = 1)))
  expect_identical(twice$union_type, c("u", "u"))
  expect_identical(twice$time_start, c(0, 0))
  empty <- stratify_rates(rates[FALSE, ], list(group = c(A = 1)))
  expect_equal(nrow(empty), 0)
  expect_true("group" %in% names(empty))
})

test_that("stratify_rates validates its inputs", {
  rates <- strata_mortality_rates()
  for (bad in list(NULL, list(), list(c(A = 1)), list(group = c(A = 1), group = c(B = 1)),
                   list(group = numeric()), list(group = c(1, 2)), list(group = c(A = NA_real_)),
                   list(group = c(A = Inf)), list(group = c(A = -1)), list(group = c(A = 1, A = 2)),
                   list(group = setNames(1, "")), list(group = c(A = "1")), list(group = list(A = 1)),
                   c(group = 1), data.frame(group = 1)))
    expect_error(stratify_rates(rates, bad), "by")
  expect_error(stratify_rates(rates, list(rate = c(A = 1))), "band, rate or calendar")
  expect_error(stratify_rates(rates, list(age_start = c(A = 1))), "band, rate or calendar")
  expect_error(stratify_rates(rates, list(`1x` = c(A = 1))), "valid column names")
  # Any other name may be a group dimension, checked at simulation.
  expect_named(stratify_rates(rates, list(education = c(low = 1))), c(names(rates), "education"))
  expect_error(stratify_rates(rates, list(sex = c(f = 1))), "already has")
  expect_error(stratify_rates(rates$rate, list(group = c(A = 1))), "data frame")
  expect_error(stratify_rates(data.frame(age_start = 0), list(group = c(A = 1))), "rate")
  big <- transform(rates, rate = 1e308)
  expect_error(stratify_rates(big, list(group = c(A = 10))), "finite")
  expect_no_error(stratify_rates(rates, list(group = c(A = 0, B = 2L))))
})

test_that("setters reject stratified columns their process does not support", {
  base <- data.frame(age_start = 15, age_end = Inf, rate = 0.1)
  expect_error(mortality(stratify_rates(strata_mortality_rates(), list(union_type = c(u = 1)))),
               "Mortality rates cannot be stratified by `union_type`")
  expect_error(fertility(stratify_rates(base, list(sex = c(f = 1)))), "Fertility")
  expect_error(union_formation("rate", rates = stratify_rates(base, list(sex = c(f = 1)))),
               "Union-search")
  expect_error(set_union_formation(hammel_model(), "rate", rates = base,
    cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 0.5),
    union_convert_rate = stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = 0.1),
                                        list(previous_union = c(none = 1)))), "Conversion")
  expect_error(union_dissolution(stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = 0.1),
                                                list(sex = c(f = 1)))), "Union-dissolution")
  expect_error(group_transition(stratify_rates(data.frame(origin = "A", destination = "B",
    age_start = 0, age_end = Inf, rate = 0.1), list(previous_union = c(none = 1)))), "Group-transition")
  # The check uses the table's columns, however they were added.
  expect_error(mortality(transform(strata_mortality_rates(), union_type = "cohab")),
               "Mortality rates cannot be stratified by `union_type`")
  expect_error(union_dissolution(data.frame(sex = "f", dur_start = 0, dur_end = Inf, rate = 0.1)),
               "Union-dissolution")
  # Supported columns pass.
  expect_s3_class(mortality(stratify_rates(data.frame(age_start = 0, age_end = Inf, rate = 0.01),
    list(group = c(A = 1), sex = c(m = 1, f = 0.9)))), "hammel_mortality")
  expect_s3_class(fertility(stratify_rates(base, list(group = c(A = 1), union_type = c(partnered = 1, unpartnered = 0.5)))),
                  "hammel_fertility")
  expect_s3_class(union_formation("rate", rates = stratify_rates(base,
    list(previous_union = c(none = 1, separated = 2, widowed = 1)))), "hammel_union_formation")
  expect_s3_class(group_transition(stratify_rates(data.frame(age_start = 0, age_end = Inf, rate = 0.1),
    list(origin = c(A = 1), destination = c(B = 1), sex = c(f = 1, m = 2)))), "hammel_group_transition")
})

test_that("stratified tables equal explicit schedules for each process and together", {
  factors <- c(A = 1.1, B = 0.9, C = 1)
  for (selected in list("mortality", "fertility", "union_formation",
                        c("mortality", "fertility", "union_formation"))) {
    stratified <- explicit <- strata_model()
    for (name in selected) {
      rates <- stratified$processes[[name]]$rates
      stratified$processes[[name]]$rates <- stratify_rates(rates, list(group = factors))
      explicit$processes[[name]]$rates <- explicit_strata(rates, factors)
    }
    expect_identical(run_simulation(stratified, 20, seed = 654), run_simulation(
      explicit, 20, seed = 654))
  }
})

test_that("zero factors suppress only the selected risk", {
  model <- make_model(2, death_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$population$groups <- list(group = group())
  model$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
                                                    list(group = c(A = 0, B = 1)))
  result <- run_simulation(model, 1, seed = 654)
  expect_true(is.na(result$persons$death_time[1]))
  expect_false(is.na(result$persons$death_time[2]))
  model <- make_rate_model(2, search_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$population$groups <- list(group = group())
  search <- model$processes$union_formation$rates
  model$processes$union_formation$rates <- stratify_rates(search, list(group = c(A = 0, B = 1)))
  expect_equal(nrow(run_simulation(model, 1, seed = 654)$unions), 0)
  model$processes$union_formation$rates <- stratify_rates(search, list(group = c(A = 1, B = 0)))
  expect_equal(nrow(run_simulation(model, 1, seed = 654)$unions), 1)
})

test_that("every group must be listed and groups must be declared", {
  model <- strata_model()
  model$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
                                                    list(group = c(A = 1.2, B = 1)))
  expect_error(run_simulation(model, 0), "missing schedules.*C")
  model$population$persons$group <- NULL
  model$population$groups <- NULL
  expect_error(run_simulation(model, 0), "not a declared group dimension")
})

test_that("stratified tables support numeric group codes and restarts", {
  model <- strata_model()
  model$population$persons$group <- match(model$population$persons$group, c("A", "B", "C"))
  model$population$groups <- list(group = group())
  factors <- c("1" = 1.1, "2" = 0.9, "3" = 1)
  explicit <- model
  for (name in c("mortality", "fertility", "union_formation")) {
    rates <- model$processes[[name]]$rates
    model$processes[[name]]$rates <- stratify_rates(rates, list(group = factors))
    explicit$processes[[name]]$rates <- explicit_strata(rates, factors)
  }
  result <- run_simulation(model, 5, seed = 654)
  expect_identical(result, run_simulation(explicit, 5, seed = 654))
  expect_identical(run_simulation(set_population(model, result), 5, seed = 123),
                   run_simulation(set_population(explicit, result), 5, seed =
                     123))
})

test_that("rate summary bands collapse only for identical stratum sets", {
  rates <- tibble::tibble(
    race = c("A", "A", "B", "B"),
    age_start = c(15, 20, 20, 15),
    age_end = c(20, 25, 25, 20),
    parity_start = c(0, 1, 1, 0),
    parity_end = c(1, Inf, Inf, 1),
    rate = c(0, 0.2, 0.4, 0.1)
  )
  out <- summarize_rate_strata(rates)
  expect_identical(out$stratified_by, "race")
  expect_identical(out$strata, tibble::tibble(race = c("A", "B")))
  expect_false(out$bands_by_stratum)
  expect_named(out$bands, c(
    "age_start", "age_end", "parity_start", "parity_end"
  ))
  expect_equal(nrow(out$bands), 2L)
  expect_identical(out, summarize_rate_strata(rbind(rates, rates[1, ])))

  # Marginal age and parity sets match, but their joint coverage differs.
  rates$parity_start[3:4] <- c(0, 1)
  rates$parity_end[3:4] <- c(1, Inf)
  out <- summarize_rate_strata(rates)
  expect_true(out$bands_by_stratum)
  expect_identical(out$bands$race, rates$race)
  expect_equal(nrow(out$bands), 4L)
})

test_that("rate summaries retain mixed clocks and calendar strata", {
  rates <- tibble::tibble(
    previous_union = c("none", "dissolved"),
    age_start = c(18, NA_real_),
    age_end = c(Inf, NA_real_),
    dur_start = c(NA_real_, 0),
    dur_end = c(NA_real_, 10),
    time_start = c(2000, 2010),
    time_end = c(2010, Inf),
    rate = c(0.1, 0)
  )
  out <- summarize_rate_strata(rates)
  expect_true(out$bands_by_stratum)
  expect_equal(out$bands, rates[setdiff(names(rates), "rate")])
  expect_s3_class(out$strata, "tbl_df")
  expect_null(summarize_rate_strata(NULL))

  pooled <- summarize_rate_strata(rates[2:8])
  expect_length(pooled$stratified_by, 0L)
  expect_false(pooled$bands_by_stratum)
  expect_s3_class(pooled$bands, "tbl_df")

  empty <- summarize_rate_strata(rates[FALSE, ])
  expect_equal(nrow(empty$bands), 0L)
  expect_false(empty$bands_by_stratum)
})

test_that("rate metadata printing separates shared and differing bands", {
  rates <- tibble::tibble(
    sex = c("m", "f"),
    age_start = c(0, 0),
    age_end = c(Inf, Inf),
    rate = c(0.1, 0.2)
  )
  metadata <- summarize_rate_strata(rates)
  text <- capture.output(value <- print_rate_strata(metadata))
  expect_identical(value, metadata)
  expect_true(any(grepl("sex: m, f", text, fixed = TRUE)))
  expect_equal(sum(grepl("[0, Inf)", text, fixed = TRUE)), 1L)
  expect_false(any(grepl("Bands by stratum", text, fixed = TRUE)))

  rates$age_start[2] <- 5
  text <- capture.output(print_rate_strata(summarize_rate_strata(rates)))
  expect_true(any(grepl("sex = m", text, fixed = TRUE)))
  expect_true(any(grepl("sex = f", text, fixed = TRUE)))
  expect_true(any(grepl("[5, Inf)", text, fixed = TRUE)))
  expect_length(capture.output(print_rate_strata(NULL)), 0L)
})

test_that("printed marginal bands retain parity and calendar intervals", {
  rates <- tibble::tibble(
    age_start = c(15, 20),
    age_end = c(20, 25),
    parity_start = c(0, 1),
    parity_end = c(1, Inf),
    time_start = c(2000, 2010),
    time_end = c(2010, Inf),
    rate = 0.1
  )
  text <- capture.output(print_rate_strata(summarize_rate_strata(rates)))
  expect_true(any(grepl("None", text, fixed = TRUE)))
  expect_true(any(grepl("Age bands (years)", text, fixed = TRUE)))
  expect_true(any(grepl("Calendar bands (years)", text, fixed = TRUE)))
  expect_true(any(grepl("[1, Inf)", text, fixed = TRUE)))
  expect_true(any(grepl("[2010, Inf)", text, fixed = TRUE)))
})

test_that("large crossed schedules print compact marginal band sequences", {
  grid <- expand.grid(age = 15:49, year = 1970:2020)
  rates <- tibble::tibble(
    age_start = grid$age,
    age_end = grid$age + 1,
    time_start = grid$year,
    time_end = grid$year + 1,
    rate = 0.1
  )
  metadata <- summarize_rate_strata(rates)
  text <- capture.output(value <- print_rate_strata(metadata))
  expect_identical(value, metadata)
  expect_equal(nrow(metadata$bands), 35L * 51L)
  expect_lt(length(text), 15L)
  expect_true(any(grepl("35 bands; width 1", text, fixed = TRUE)))
  expect_true(any(grepl("51 bands; width 1", text, fixed = TRUE)))
  expect_true(any(grepl("[49, 50)", text, fixed = TRUE)))
})

test_that("irregular and open-ended band lists are bounded and explicit", {
  start <- c(0:8, 12)
  end <- c(1:9, Inf)
  text <- format_rate_bands(start, end)
  expect_match(text, "10 bands; showing 5", fixed = TRUE)
  expect_match(text, "[12, Inf)", fixed = TRUE)
  expect_false(grepl("width", text, fixed = TRUE))
  expect_identical(format_rate_bands(0, Inf), "[0, Inf)")

  # Equal widths with gaps must not be described as a regular sequence.
  text <- format_rate_bands(seq(0, 18, 2), seq(1, 19, 2))
  expect_match(text, "showing 5", fixed = TRUE)
})

test_that("calendar periods sharing band structures print together", {
  rates <- tibble::tibble(
    time_start = rep(c(2000, 2001, 2002, 2003), each = 2),
    time_end = rep(c(2001, 2002, 2003, 2004), each = 2),
    age_start = c(15, 20, 20, 15, 15, 25, 15, 20),
    age_end = c(20, 25, 25, 20, 25, 35, 20, 25),
    rate = seq_len(8) / 100
  )
  metadata <- summarize_rate_strata(rates)
  text <- capture.output(value <- print_rate_strata(metadata))
  expect_identical(value, metadata)
  expect_true(any(grepl("Bands by calendar period", text, fixed = TRUE)))
  expect_equal(sum(grepl("Calendar bands (years)", text, fixed = TRUE)), 2L)
  expect_equal(sum(grepl("Age bands (years)", text, fixed = TRUE)), 2L)
  expect_true(any(grepl(
    "[2000, 2001), [2001, 2002), [2003, 2004)", text, fixed = TRUE
  )))
  expect_true(any(grepl("[15, 25), [25, 35)", text, fixed = TRUE)))
})

test_that("calendar structure grouping preserves joint parity coverage", {
  rates <- tibble::tibble(
    time_start = rep(c(2000, 2010), each = 2),
    time_end = rep(c(2010, Inf), each = 2),
    age_start = rep(c(15, 20), 2),
    age_end = rep(c(20, 25), 2),
    parity_start = c(0, 1, 1, 0),
    parity_end = c(1, Inf, Inf, 1),
    rate = 0.1
  )
  text <- capture.output(print_rate_strata(summarize_rate_strata(rates)))
  expect_true(any(grepl("Bands by calendar period", text, fixed = TRUE)))
  expect_equal(sum(grepl("Parity bands", text, fixed = TRUE)), 2L)
})
