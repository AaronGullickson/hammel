# Matching calibration takes search intensities from the search cells and
# cohabitation probabilities from the probability tracker. With
# hammel.check_cells, every search solve compares the cell-derived marriage and
# cohabitation totals with an explicit scan of searching women, recomputing each
# woman's probability from her schedule.

test_that("cell search intensities match explicit scans through every state change", {
  set.seed(31)
  pop <- create_baseline_population(3000, max_age = 60, groups = list(group = c(A = .5, B = .5)))
  mortality <- data.frame(sex = rep(c("m", "f"), 2), age_start = 0, age_end = Inf,
    rate = c(.01, .01, .03, .03), time_start = rep(c(0, 2), each = 2), time_end = rep(c(2, Inf), each = 2))
  # Never-partnered women search by age; separated women by time since their
  # union ended (with a band edge at one year); widowed women by age.
  search <- data.frame(previous_union = c("none", "none", "none", "separated", "separated", "widowed"),
    age_start = c(0, 15, 25, NA, NA, 0), age_end = c(15, 25, Inf, NA, NA, Inf),
    dur_start = c(NA, NA, NA, 0, 1, NA), dur_end = c(NA, NA, NA, 1, Inf, NA),
    rate = c(0, .3, .1, .05, .4, .05))
  model <- hammel_model() |>
    set_population(pop) |>
    set_mortality(mortality) |>
    set_fertility(data.frame(union_type = c("partnered", "unpartnered"), age_start = 15,
      age_end = 45, rate = c(.15, .03)), birth_interval = .75) |>
    set_union_formation("rate", rates = search,
      matching = matching_model(~ age_difference, c(age_difference = 0),
        group_exogamy = c("A:B" = -1)),
      cohab_prob = data.frame(age_start = c(0, 25), age_end = c(25, Inf), prob = c(.6, .3)),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .1),
      kinship_generations = 0) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .1))
  withr::local_options(hammel.check_cells = TRUE)
  out <- run_simulation(model, 6, seed = 32)
  starts <- out$unions$start_time >= 0
  expect_gt(sum(starts), 500)
  expect_gt(sum(out$unions$union_type[starts] == "cohab"), 100)
  expect_gt(sum(out$unions$end_reason %in% "dissolution"), 50)
  # Remarriage after dissolution exercises the duration clock.
  wives <- out$unions$wife_id[starts]
  expect_gt(sum(duplicated(wives)), 10)
  # The same model with search weights uses the same calibration path.
  weighted <- set_union_formation(model, "rate", rates = search,
    matching = matching_model(~ age_difference, c(age_difference = 0),
      group_exogamy = c("A:B" = -1)),
    cohab_prob = data.frame(age_start = c(0, 25), age_end = c(25, Inf), prob = c(.6, .3)),
    union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .1),
    kinship_generations = 0,
    # The baseline declares its categories in the order given (A, B).
    weights = weight_model(~ group, c(groupB = log(2))))
  expect_s3_class(run_simulation(weighted, 3, seed = 33), "hammel_population")
})
