status_rates <- function(partnered = 0, unpartnered = 0) {
  data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(partnered, unpartnered))
}

test_that("conditional fertility validates labels, completeness and per-status bands", {
  rates <- status_rates(0.2, 0.05)
  expect_s3_class(fertility(rates)$rates, "tbl_df")
  expect_identical(as.data.frame(fertility(rates)$rates), rates)
  for (labels in list(c("partnered", "total"), c("partnered", NA),
                      c("partnered", "partnered"), c("married", "unmarried"),
                      factor(c("partnered", "unpartnered")))) {
    bad <- rates
    bad$union_type <- labels
    expect_error(fertility(bad), "union_type")
  }
  expect_error(fertility(rates[1, ]), "union_type")
  expect_error(fertility(rbind(rates, rates[1, ])), "overlap")
  different <- rates
  different$age_start <- c(20, 15)
  different$age_end <- c(45, 40)
  expect_s3_class(fertility(different), "hammel_fertility")
  model <- make_model()
  edited <- fertility(rates)
  edited$rates$union_type[2] <- "total"
  expect_error(set_process(model, edited), "union_type")
  model$processes$fertility <- edited
  expect_error(simulate(model, 0), "union_type")
})

test_that("on-demand rejects status schedules, rate overall warns exactly once", {
  model <- set_fertility(make_model(), status_rates())
  expect_error(simulate(model, 0), "overall fertility schedule")
  expect_error(simulate(model, 1), "overall fertility schedule")
  model <- make_rate_model(1, search_rate = 0)
  model <- set_fertility(model, data.frame(age_start = 15, age_end = 50, rate = 20),
                         birth_interval = 0)
  warnings <- character()
  result <- withCallingHandlers(simulate(model, 1, seed = 52), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 1)
  expect_match(warnings, "overall fertility schedule.*union_type")
  expect_gt(result$simulation$births, 0)
  expect_true(all(is.na(result$persons$father_id)))
  expect_equal(nrow(result$unions), 0)
  expect_warning(simulate(model, 0), "overall fertility schedule")
})

test_that("unpartnered births select fathers without triggering unions", {
  model <- set_fertility(make_rate_model(search_rate = 0), status_rates(0, 20), 0)
  expect_no_warning(result <- simulate(model, 1, seed = 52))
  expect_gt(result$simulation$births, 0)
  expect_equal(nrow(result$unions), 0)
  expect_true(all(result$persons$father_id[result$persons$birth_time >= 0] == 2))
  expect_identical(result, simulate(model, 1, seed = 52))
})

test_that("formation and widowhood switch fertility schedules without bypassing spacing", {
  # An initial recorded birth at time zero sets a spacing boundary at 0.5.
  model <- make_rate_model(3, search_rate = 100)
  p <- model$population$persons
  p$birth_time[3] <- 0
  p$mother_id[3] <- 1
  p$sex[3] <- "m"
  model <- set_population(model, hammel_population(p, end_time = 0, groups = if ("group" %in% names(p)) "group"))
  # A high partnered rate makes births near-certain once spacing ends, so the
  # check does not depend on the seed and several births must respect spacing.
  model <- set_fertility(model, status_rates(5, 0), 0.5)
  result <- simulate(model, 2, seed = 52)
  expect_lt(result$unions$start_time[1], 0.5)
  children <- result$persons[result$persons$person_id > 3, ]
  expect_gt(nrow(children), 0)
  expect_true(all(children$birth_time >= 0.5))
  expect_true(all(diff(c(0, children$birth_time)) >= 0.5))
  expect_true(all(children$father_id == 2))

  # Partnered fertility is zero. Widowhood enables births, with no new search.
  model <- make_rate_model(search_rate = 0)
  model$population$persons$birth_time[2] <- -30
  model$population$unions <- tibble::tibble(union_id = 1, wife_id = 1,
    husband_id = 2, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  model <- set_fertility(model, status_rates(0, 20), 0)
  model <- set_mortality(model, data.frame(sex = c("f", "m", "m"),
    age_start = c(0, 0, 30.3), age_end = c(Inf, 30.3, Inf), rate = c(0, 0, 1e8)))
  result <- simulate(model, 1, seed = 52)
  children <- result$persons[result$persons$birth_time >= 0, ]
  expect_gt(nrow(children), 0)
  expect_true(all(children$birth_time > result$persons$death_time[2]))
  expect_true(all(is.na(children$father_id)))
  continued <- simulate(set_population(model, result), 0.2, seed = 53)
  expect_identical(continued$unions, result$unions)
  expect_gt(continued$simulation$births, 0)
})

test_that("fertility compensation and different age bands remain separate by status", {
  model <- make_rate_model(search_rate = 0)
  model$population$persons$birth_time[1] <- -20
  # Both schedules would be infeasible if mistakenly integrated together.
  model <- set_fertility(model, status_rates(0.8, 0.8), 0.75)
  expect_no_error(simulate(model, 1, seed = 52))
  rates <- status_rates(20, 20)
  rates$age_start <- c(30, 15)
  rates$age_end <- c(40, 19)
  model <- set_fertility(model, rates, 0)
  expect_equal(simulate(model, 1, seed = 52)$simulation$births, 0)
})

test_that("constant status-specific fertility gives the corresponding birth counts", {
  model <- make_rate_model(6000, search_rate = 0)
  model$population$persons$sex[1:3000] <- "f"
  model$population$unions <- tibble::tibble(union_id = 1:1500,
    wife_id = 1:1500, husband_id = 3001:4500, start_time = -1,
    end_time = NA_real_, end_reason = NA_character_)
  model <- set_fertility(model, status_rates(0.3, 0.1), 0)
  result <- simulate(model, 1, seed = 137)
  children <- result$persons[result$persons$birth_time >= 0, ]
  partnered <- children$mother_id <= 1500
  expect_lt(abs(sum(partnered) - 450), 6 * sqrt(450))
  expect_lt(abs(sum(!partnered) - 150), 6 * sqrt(150))
  expect_true(all(children$father_id[partnered] == children$mother_id[partnered] + 3000))
  expect_true(all(children$father_id[!partnered] %in% 4501:6000))
})
