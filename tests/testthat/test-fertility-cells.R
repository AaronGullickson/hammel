# Cell-accounted fertility under birth spacing (D101).

cell_model <- function(n = 2000, rate = .2, spacing = .75, age = 20,
                       rates = data.frame(age_start = 15, age_end = 50, rate = rate)) {
  people <- data.frame(person_id = seq_len(n), sex = "f", birth_time = -age,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
  hammel_model() |>
    set_population(hammel_population(people, end_time = 0, groups = if ("group" %in% names(people)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(rates, birth_interval = spacing) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
}

# Births to the founders and their woman-years, ignoring descendants.
founder_rate <- function(pop, n, years) {
  births <- sum(pop$persons$mother_id %in% seq_len(n))
  births / (n * years)
}

test_that("older specifications carrying spacing_method are accepted", {
  model <- cell_model(10)
  old <- model
  old$processes$fertility$spacing_method <- "cell"
  expect_null(set_process(old, old$processes$fertility)$processes$fertility$spacing_method)
  old$processes$fertility$spacing_method <- "midpoint"
  expect_warning(attached <- set_process(old, old$processes$fertility), "no longer available")
  expect_null(attached$processes$fertility$spacing_method)
  expect_error(fertility(data.frame(age_start = 15, age_end = 50, rate = .1),
                         spacing_method = "cell"), "unused argument")
})

test_that("simulations report unrealized births as a typed tibble", {
  result <- simulate(cell_model(20), 0, seed = 1)
  out <- result$simulation$unrealized_births
  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 0L)
  expect_identical(names(out), c("period_start", "period_end", "age_start", "age_end",
                                 "unrealized_births"))
})

test_that("cells without an eligible woman count unrealized births", {
  # Women whose previous birth (a son) was at -0.1 with a 2-year interval are
  # ineligible until 1.9: their cells have members but no eligible woman, and
  # each would-be birth is counted as unrealized.
  n <- 2000
  model <- hammel_model() |>
    set_population(spacing_population(n)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .5),
                  birth_interval = 2) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  result <- simulate(model, 1.5, seed = 4)
  out <- result$simulation$unrealized_births
  expect_equal(nrow(result$persons), 2 * n)
  expect_equal(out$period_start, c(0, 1))
  expect_equal(out$period_end, c(1, 1.5))
  expect_true(all(out$unrealized_births == round(out$unrealized_births)))
  expect_count(out$unrealized_births[1], .5 * n)
  expect_count(out$unrealized_births[2], .25 * n)
  # Cells are 1-year age sub-bands: they are 30 at the start and 31 after a year.
  expect_equal(out$age_start, c(30, 31))
  expect_equal(out$age_end, c(31, 32))
})

test_that("cells split ages at whole years and pooled parity bands by actual parity", {
  # Women aged 20.5 at actual parity 4 in a pooled [2,Inf) band, last birth
  # at -0.1 with a 2-year interval: members and none eligible.
  n <- 1000
  rates <- data.frame(age_start = 15, age_end = Inf, parity_start = c(0, 2),
                      parity_end = c(2, Inf), rate = c(.1, .4))
  model <- hammel_model() |>
    set_population(spacing_population(n, 20.5, c(-4, -3, -2, -0.1))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    (\(m) quiet_unreachable(set_fertility(m, rates, birth_interval = 2)))() |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  expect_warning(fertility(rates, birth_interval = 2),
                 "parity 2 at ages \\[15, 19\\) \\(earliest eligible age 19\\)")
  out <- simulate(model, 1, seed = 7)$simulation$unrealized_births
  expect_equal(out$parity_start, c(4, 4))
  expect_equal(out$parity_end, c(5, 5))
  expect_equal(out$age_start, c(20, 21))
  expect_equal(out$age_end, c(21, 22))
  expect_count(out$unrealized_births[1], .2 * n)
  expect_count(out$unrealized_births[2], .2 * n)
})

test_that("zero spacing generates births at the input rate", {
  n <- 3000; years <- 2
  pop <- simulate(cell_model(n, rate = .1, spacing = 0), years, seed = 21)
  births <- sum(pop$persons$mother_id %in% seq_len(n))
  expected <- .1 * n * years
  expect_lt(abs(births - expected), 4 * sqrt(expected))
  expect_identical(nrow(pop$simulation$unrealized_births), 0L)
})

test_that("unstratified cell accounting recovers the input rate from a no-history start", {
  # Nobody is initially in spacing; cell accounting uses the actual eligible
  # share from the first instant, so the first year is not inflated.
  n <- 10000
  first_year <- function(pop) sum(pop$persons$mother_id %in% seq_len(n)) / n
  expect_lt(abs(first_year(simulate(cell_model(n, rate = .3), 1, seed = 31)) - .3),
            4 * sqrt(.3 / n))
  n <- 3000
  longer <- simulate(cell_model(n, rate = .3), 3, seed = 32)
  expect_lt(abs(founder_rate(longer, n, 3) - .3), 4 * sqrt(.3 / (3 * n)))
})

test_that("parity cells recover the pooled upper parity rate in an initialized cohort", {
  n <- 3000
  rates <- data.frame(age_start = 15, age_end = Inf, parity_start = 0:2,
    parity_end = c(1, 2, Inf), rate = c(.2, .3, .4))
  pop <- simulate(quiet_unreachable(cell_model(n, rates = rates)), 10, seed = 41)
  observed <- calculate_fertility_rates(pop, age_breaks = c(15, Inf), time_breaks = c(0, 10),
                                        parity_breaks = c(0, 1, 2))
  top <- observed[observed$parity_start == 2, ]
  expect_lt(abs(top$rate_fertility - .4), 4 * sqrt(.4 / top$person_years))
  out <- pop$simulation$unrealized_births
  expect_true(all(c("parity_start", "parity_end") %in% names(out)))
  # Start-up: nobody can be at parity one and eligible before .75 years.
  expect_gt(sum(out$unrealized_births[out$parity_start == 1 & out$period_start == 0]), 0)
})

test_that("cell mode runs with groups, union types, cohabitation, matching and restarts", {
  # Internal option: rescan cell counts at every event and calibration totals
  # at every fertility-driven search.
  old <- options(hammel.check_cells = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(5)
  base <- create_baseline_population(1500, max_age = 60)
  base$persons$group <- sample(c("A", "B"), 1500, TRUE)
  base$groups <- list(group = group())
  fert <- data.frame(age_start = c(15, 25, 35), age_end = c(25, 35, 50), rate = c(.1, .12, .04))
  sf <- rbind(transform(fert, union_type = "unpartnered", rate = rate * .3),
              transform(fert, union_type = "cohab", rate = rate * .8),
              transform(fert, union_type = "marriage"))
  mort <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .01)
  model <- hammel_model() |> set_population(base) |> set_mortality(mort) |>
    set_fertility(stratify_rates(sf, list(group = c(A = 1, B = 1.4)))) |>
    set_union_formation("rate",
      create_matching_model(~ age_difference, c(age_difference = -.1), group_exogamy = c("A:B" = -1)),
      rates = data.frame(age_start = 15, age_end = Inf, rate = .15),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = .4),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = .2)) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = c(.1, .03),
                                     union_type = c("cohab", "marriage")))
  first <- simulate(model, 5, seed = 51)
  expect_identical(first, simulate(model, 5, seed = 51))
  expect_gt(first$simulation$births, 0)
  expect_true(all(c("group", "union_type") %in% names(first$simulation$unrealized_births)))
  restarted <- simulate(set_population(model, first), 3, seed = 52)
  expect_equal(restarted$simulation$start_time, 5)
  expect_true(all(restarted$simulation$unrealized_births$period_start >= 5))
  # Grouped on-demand mode exercises cell-derived calibration totals.
  od <- hammel_model() |> set_population(base) |> set_mortality(mort) |>
    set_fertility(fert) |>
    set_union_formation("on_demand",
      create_matching_model(~ age_difference, c(age_difference = -.1), group_exogamy = c("A:B" = -1)))
  expect_gt(simulate(od, 5, seed = 53)$simulation$births, 0)
})

test_that("cancelled attempts keep eligibility and do not create births", {
  pop <- simulate(set_union_formation(cell_model(50, rate = .5), "on_demand",
                                      cancel_birth_if_no_match = TRUE), 2, seed = 61)
  expect_identical(pop$simulation$births, 0)
  expect_identical(nrow(pop$simulation$unrealized_births), 0L)
})

test_that("unreachable parity/age cells are reported from the schedule and interval", {
  # Constant rates from 15: parity k is first eligible at 15 + k * .75.
  rates <- data.frame(age_start = 15, age_end = 40, parity_start = 0:5,
                      parity_end = c(1:5, Inf), rate = .3)
  cells <- hammel:::unreachable_parity_cells(rates, .75)
  expect_equal(unique(cells$parity), 2:5)
  expect_equal(cells$earliest_eligible[cells$parity == 5][1], 18.75)
  expect_equal(max(cells$age_end[cells$parity == 5]), 18)
  expect_warning(fertility(rates), "parity 5 at ages \\[15, 18\\) \\(earliest eligible age 18.75\\)")
  # Partly reachable cells are not reported: parity 1 is eligible from 15.75.
  expect_false(1 %in% cells$parity)
  # Zero or no spacing, and schedules without parity, raise no warning.
  expect_no_warning(fertility(rates, birth_interval = 0))
  expect_no_warning(fertility(data.frame(age_start = 15, age_end = 40, rate = .3)))
  # A later start of positive rates delays eligibility at higher parities.
  # Parity 2 is then reachable (eligible from 16.5); parity 3 is eligible from
  # 20.75, so [20,21) is partly reachable; parity 4 (from 21.5) is not.
  later <- transform(rates, age_start = ifelse(parity_start >= 2, 20, 15))
  late <- hammel:::unreachable_parity_cells(later, .75)
  expect_equal(unique(late$parity), 4:5)
  expect_true(all(late$age_start >= 20))
  # A zero schedule at one parity makes higher parities unreachable at any age.
  gap <- transform(rates, rate = ifelse(parity_start == 1, 0, .3))
  out <- hammel:::unreachable_parity_cells(gap, .75)
  expect_true(all(is.infinite(out$earliest_eligible)))
  expect_warning(fertility(gap), "unreachable at any age")
  # Groups are checked separately and named.
  grouped <- rbind(transform(rates, group = "A"), transform(later, group = "B"))
  expect_warning(fertility(grouped), "group A, parity 2")
  by_group <- hammel:::unreachable_parity_cells(grouped, .75)
  expect_equal(unique(by_group$parity[by_group$group == "A"]), 2:5)
  expect_equal(unique(by_group$parity[by_group$group == "B"]), 4:5)
})

test_that("upper-open bands stay one cell only when interval x rate is very small", {
  # Women aged 30 whose last birth was at -0.1 with a 2-year interval.
  open_case <- function(rate, n) {
    model <- hammel_model() |>
      set_population(spacing_population(n)) |>
      set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
      set_fertility(data.frame(age_start = 15, age_end = Inf, rate = rate),
                    birth_interval = 2) |>
      set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
    simulate(model, 1.5, seed = 8)$simulation$unrealized_births
  }
  small <- open_case(.0004, 50000) # 2 * .0004 = .0008 < .001: one open cell
  expect_equal(small$age_start, c(15, 15))
  expect_equal(small$age_end, c(Inf, Inf))
  expect_count(small$unrealized_births[1], .0004 * 50000)
  expect_count(small$unrealized_births[2], .0002 * 50000)
  large <- open_case(.5, 2000)    # split into yearly cells
  expect_equal(large$age_start, c(30, 31))
  expect_equal(large$age_end, c(31, 32))
  expect_count(large$unrealized_births[1], .5 * 2000)
  expect_count(large$unrealized_births[2], .25 * 2000)
})

test_that("released eligibility does not consume blocked time or force a birth", {
  # Review reproduction: blocked 0.65 years at rate .5 gives .325 expected
  # unrealized births per woman before release, and release must not trigger
  # an immediate birth.
  n <- 2000
  model <- hammel_model() |>
    set_population(spacing_population(n)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .5),
                  birth_interval = .75) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  result <- simulate(model, .9, seed = 2)
  expect_count(result$simulation$unrealized_births$unrealized_births, .5 * .65 * n)
  # After release at .65 the women give births at the input rate.
  released <- sum(result$persons$birth_time >= .65 & result$persons$mother_id %in% seq_len(n))
  expect_count(released, .5 * .25 * n)
  # One woman, many intervals (sons only, so no descendants reproduce): none
  # may sit exactly at the minimum spacing.
  one <- data.frame(person_id = 1, sex = "f", birth_time = -15, death_time = NA_real_,
                    mother_id = NA_real_, father_id = NA_real_)
  long <- hammel_model() |>
    set_population(hammel_population(one, end_time = 0, groups = if ("group" %in% names(one)) "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 1000, rate = 2),
                  birth_interval = .75, ffab = 0) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(long, 300, seed = 3)
  times <- sort(pop$persons$birth_time[pop$persons$mother_id %in% 1])
  gaps <- diff(times)
  expect_gt(length(gaps), 200)
  expect_identical(sum(abs(gaps - .75) < 1e-12), 0L)
  # Excess over the interval is exponential with the input-equivalent rate.
  expect_lt(abs(mean(gaps - .75) - 1 / 2), 4 * (1 / 2) / sqrt(length(gaps)))
})

test_that("the unreachable-cell check terminates for infinite unreachable ranges", {
  rates <- data.frame(age_start = 15, age_end = Inf, parity_start = 0:1,
                      parity_end = c(1, Inf), rate = c(0, .1))
  cells <- hammel:::unreachable_parity_cells(rates, .75)
  expect_equal(nrow(cells), 1)
  expect_equal(cells$age_end, Inf)
  expect_warning(fertility(rates), "parity 1 at ages \\[15, Inf\\) \\(unreachable at any age\\)")
})

test_that("calculate_fertility_rates attributes unrealized births exactly or returns NA", {
  # Women aged 30 in spacing until 1.9: unrealized .5 per woman in [0,1),
  # .25 in [1,1.5).
  n <- 2000
  model <- hammel_model() |>
    set_population(spacing_population(n)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .5),
                  birth_interval = 2) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 1.5, seed = 4)
  r <- calculate_fertility_rates(pop, c(0, 1, 1.5), c(15, 50), unrealized = TRUE)
  r <- r[r$age_start == 15, ]
  expect_count(r$unrealized_births[1], .5 * n)
  expect_count(r$unrealized_births[2], .25 * n)
  expect_true(all(abs(r$rate_fertility_adjusted - .5) < 4 * sqrt(.5 / r$person_years)))
  expect_equal(r$rate_fertility, c(0, 0))
  # Without the argument the schema is unchanged.
  expect_false("unrealized_births" %in% names(calculate_fertility_rates(pop, c(0, 1.5), c(15, 50))))
  # A break that cuts through a recorded one-year period makes the rows unknown.
  cut <- calculate_fertility_rates(pop, c(0, .5, 1.5), c(15, 50), unrealized = TRUE)
  expect_true(all(is.na(cut$unrealized_births[cut$age_start == 15])))
  # Ages outside the cells get zero; a break inside a cell's year gives NA.
  ages <- calculate_fertility_rates(pop, c(0, 1.5), c(15, 30.5, 50), unrealized = TRUE)
  expect_true(all(is.na(ages$unrealized_births[ages$age_start %in% c(15, 30.5)])))
  expect_equal(ages$unrealized_births[ages$age_start == 50], 0)
  # Restarts keep only the latest run: earlier periods are unknown.
  later <- simulate(set_population(model, pop), 1, seed = 5)
  both <- calculate_fertility_rates(later, c(0, 1.5, 2.5), c(15, 50), unrealized = TRUE)
  expect_true(is.na(both$unrealized_births[both$time_start == 0 & both$age_start == 15]))
  expect_false(is.na(both$unrealized_births[both$time_start == 1.5 & both$age_start == 15]))
  # Overall schedules cannot be split by union type.
  typed <- calculate_fertility_rates(pop, c(0, 1.5), c(15, 50), by = "union_type",
                                     unrealized = TRUE)
  expect_true(all(is.na(typed$unrealized_births[typed$age_start == 15])))
  # A population that was never simulated has no unrealized record.
  expect_error(calculate_fertility_rates(model$population, c(-1, 0), c(15, 50),
                                         unrealized = TRUE), "unrealized-births record")
  expect_error(calculate_fertility_rates(pop, c(0, 1.5), c(15, 50), unrealized = NA), "TRUE or FALSE")
})

test_that("unrealized births pool into observed parity bands and restore cohort rates", {
  # Actual parity 4 inside a pooled [2,Inf) band: .2 per woman in each
  # one-year age cell.
  n <- 1000
  rates <- data.frame(age_start = 15, age_end = Inf, parity_start = c(0, 2),
                      parity_end = c(2, Inf), rate = c(.1, .4))
  model <- hammel_model() |>
    set_population(spacing_population(n, 20.5, c(-4, -3, -2, -0.1))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    (\(m) quiet_unreachable(set_fertility(m, rates, birth_interval = 2)))() |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 1, seed = 7)
  r <- calculate_fertility_rates(pop, c(0, 1), c(15, Inf), parity_breaks = c(0, 2),
                                 unrealized = TRUE)
  expect_count(r$unrealized_births[r$parity_start == 2], .4 * n)
  expect_equal(r$unrealized_births[r$parity_start == 0], 0)
  pooled <- calculate_fertility_rates(pop, c(0, 1), c(15, Inf), unrealized = TRUE)
  expect_count(pooled$unrealized_births[pooled$age_start == 15], .4 * n)
  expect_lt(abs(pooled$rate_fertility_adjusted[pooled$age_start == 15] - .4), 4 * sqrt(.4 / n))
  # Initialized cohort: adding back start-up unrealized births restores parity 1.
  n <- 3000
  cohort <- data.frame(age_start = 15, age_end = Inf, parity_start = 0:2,
                       parity_end = c(1, 2, Inf), rate = c(.2, .3, .4))
  out <- simulate(quiet_unreachable(cell_model(n, rates = cohort)), 10, seed = 41)
  adj <- calculate_fertility_rates(out, c(0, 10), c(15, Inf), parity_breaks = c(0, 1, 2),
                                   unrealized = TRUE)
  adj <- adj[adj$time_start == 0, ]
  expect_true(all(abs(adj$rate_fertility_adjusted - c(.2, .3, .4)) <
                  4 * sqrt(c(.2, .3, .4) / adj$person_years)))
  expect_lt(adj$rate_fertility[2], adj$rate_fertility_adjusted[2])
})

test_that("unrealized births follow group and union-type strata and refuse parity splits", {
  # Group-A women, unpartnered, in spacing until 1.9; their sons are group B.
  n <- 2000
  rates <- data.frame(age_start = 15, age_end = 50, rate = c(.2, .5),
                      union_type = c("partnered", "unpartnered"))
  model <- hammel_model() |>
    set_population(spacing_population(n, groups = c("A", "B"))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(rates, birth_interval = 2) |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
                        cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 1, seed = 6)
  r <- calculate_fertility_rates(pop, c(0, 1), c(15, 50), by = c("group", "union_type"), unrealized = TRUE)
  r <- r[r$age_start == 15, ]
  cell <- r$group == "A" & r$union_type == "unpartnered"
  expect_count(r$unrealized_births[cell], .5 * n)
  expect_true(all(r$unrealized_births[!cell] == 0))
  expect_lt(abs(r$rate_fertility_adjusted[cell] - .5), 4 * sqrt(.5 / n))
  # This schedule has no parity columns, so parity rows cannot be attributed.
  p <- calculate_fertility_rates(pop, c(0, 1), c(15, 50), parity_breaks = c(0, 1),
                                 unrealized = TRUE)
  expect_true(all(is.na(p$unrealized_births[p$age_start == 15])))
  expect_equal(p$unrealized_births[p$age_start == 50], c(0, 0))
})

test_that("a cut period makes only the overlapping parity rows unknown", {
  # Actual parity 4 in a pooled [2,Inf) band; a break at .5 cuts the record.
  rates <- data.frame(age_start = 15, age_end = Inf, parity_start = c(0, 2),
                      parity_end = c(2, Inf), rate = c(.1, .4))
  model <- hammel_model() |>
    set_population(spacing_population(200, 20.5, c(-4, -3, -2, -0.1))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    (\(m) quiet_unreachable(set_fertility(m, rates, birth_interval = 2)))() |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 1, seed = 7)
  r <- calculate_fertility_rates(pop, c(0, .5, 1), c(15, Inf), parity_breaks = c(0, 2),
                                 unrealized = TRUE)
  r <- r[r$age_start == 15, ]
  expect_true(all(is.na(r$unrealized_births[r$parity_start == 2])))
  expect_equal(r$unrealized_births[r$parity_start == 0], c(0, 0))
})
