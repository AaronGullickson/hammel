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
  out <- do.call(rbind, lapply(names(factors), function(g)
    transform(as.data.frame(rates), rate = rate * factors[[g]], group = g)))
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
  expect_error(group_transition(stratify_rates(data.frame(group_origin = "A", group_dest = "B",
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
    list(group_origin = c(A = 1), group_dest = c(B = 1), sex = c(f = 1, m = 2)))), "hammel_group_transition")
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
    expect_identical(simulate(stratified, 20, seed = 654), simulate(explicit, 20, seed = 654))
  }
})

test_that("zero factors suppress only the selected risk", {
  model <- make_model(2, death_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$population$groups <- list(group = group())
  model$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
                                                    list(group = c(A = 0, B = 1)))
  result <- simulate(model, 1, seed = 654)
  expect_true(is.na(result$persons$death_time[1]))
  expect_false(is.na(result$persons$death_time[2]))
  model <- make_rate_model(2, search_rate = 1e8)
  model$population$persons$group <- c("A", "B")
  model$population$groups <- list(group = group())
  search <- model$processes$union_formation$rates
  model$processes$union_formation$rates <- stratify_rates(search, list(group = c(A = 0, B = 1)))
  expect_equal(nrow(simulate(model, 1, seed = 654)$unions), 0)
  model$processes$union_formation$rates <- stratify_rates(search, list(group = c(A = 1, B = 0)))
  expect_equal(nrow(simulate(model, 1, seed = 654)$unions), 1)
})

test_that("every group must be listed and groups must be declared", {
  model <- strata_model()
  model$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
                                                    list(group = c(A = 1.2, B = 1)))
  expect_error(simulate(model, 0), "missing schedules.*C")
  model$population$persons$group <- NULL
  model$population$groups <- NULL
  expect_error(simulate(model, 0), "not a declared group dimension")
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
  result <- simulate(model, 5, seed = 654)
  expect_identical(result, simulate(explicit, 5, seed = 654))
  expect_identical(simulate(set_population(model, result), 5, seed = 123),
                   simulate(set_population(explicit, result), 5, seed = 123))
})
