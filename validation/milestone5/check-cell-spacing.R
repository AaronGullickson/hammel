# Offline engine validation of cell-accounted fertility (D101) against midpoint
# compensation (D028/D099). Synthetic schedules only; no rates are downloaded.
# This reports recovery and bias; it sets no scientific acceptance threshold.
# Run from the repository root:
#   Rscript validation/milestone5/check-cell-spacing.R [cores]
# Writes cell-spacing-replications.csv and cell-spacing-summary.csv here.
#
# Scenarios (10 seeds per method, seeds 9701-9710):
#   cohort     initialized cohort, 3,000 women aged 20 at parity 0, ten years,
#              no mortality; inputs .2/.3/.4 (2+), spacing .75.
#   young      generational, 300 founders aged 20, same parity inputs on ages
#              15-40, no mortality before 60; observed in years [80,100) at
#              ages [18,19) and [20,35). Pooled (2+) and actual (2,3,4+) parity.
#   restart    unstratified, 5,000-person baseline run 40 years under a TFR ~2
#              schedule, then restarted with a TFR ~7.7 schedule and spacing
#              1.25 (the milestone 1 Niger stress case); annual years 35-46.
#   union      rate-mode unions with partnered/unpartnered schedules, search,
#              dissolution and widowhood; years [40,60) by union type.
#   groups     on-demand, two groups with multiplier B = 1.4, unmatched births
#              kept; years [40,60) by group.
# Recovery ratio = observed births / sum over cells of input rate x observed
# person-years, so 1 means the inputs are recovered for that composition.
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[[1]]) else 4L
pkgload::load_all(quiet = TRUE)
seeds <- 9701:9710
methods <- c("midpoint", "cell")
out_dir <- "validation/milestone5"

mort_zero <- data.frame(sex = rep(c("m", "f"), each = 2), age_start = c(0, 60),
                        age_end = c(60, Inf), rate = c(0, 1.26))
mort_life <- data.frame(sex = rep(c("m", "f"), each = 4),
  age_start = c(0, 1, 50, 70), age_end = c(1, 50, 70, Inf),
  rate = c(.005, .002, .02, .08))
bands <- seq(15, 45, 5)
low <- data.frame(age_start = bands, age_end = bands + 5,
                  rate = c(.02, .08, .11, .1, .06, .02, .005))
niger <- transform(low, rate = c(.2, .35, .33, .28, .22, .12, .04))
parity_rates <- function(end) data.frame(age_start = 15, age_end = end,
  parity_start = 0:2, parity_end = c(1, 2, Inf), rate = c(.2, .3, .4))
founders <- function(n, age = 20) hammel_population(data.frame(person_id = seq_len(n),
  sex = "f", birth_time = -age, death_time = NA_real_, mother_id = NA_integer_,
  father_id = NA_integer_), end_time = 0)

# Births and person-years by age band (and optional strata) with input rates.
recovery <- function(pop, times, input, by = NULL, multiplier = NULL) {
  r <- calculate_fertility_rates(pop, time_breaks = times, age_breaks = c(bands, 50),
    by_union_type = identical(by, "union_type"), by_group = identical(by, "group"))
  key <- if (is.null(by)) rep("all", nrow(r)) else as.character(r[[by]])
  r$input <- input$rate[match(r$age_start, input$age_start)]
  if (identical(by, "union_type")) {
    r$input <- input$rate[match(paste(r$age_start, r$union_type),
                                paste(input$age_start, input$union_type))]
  }
  if (!is.null(multiplier)) r$input <- r$input * multiplier[key]
  r$expected <- r$input * r$person_years
  agg <- aggregate(cbind(births, expected) ~ time_start + key,
                   data = transform(r, key = key), FUN = sum)
  data.frame(time_start = agg$time_start, stratum = agg$key,
             value = agg$births / agg$expected, births = agg$births)
}

run_cohort <- function(method, seed) {
  model <- hammel_model() |> set_population(founders(3000)) |>
    set_mortality(transform(mort_zero, rate = 0)) |>
    set_fertility(parity_rates(Inf), spacing_method = method) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 10, seed)
  r <- calculate_fertility_rates(pop, time_breaks = c(0, 10), age_breaks = c(15, Inf),
                                 parity_breaks = c(0, 1, 2))
  r <- r[r$time_start == 0, ]
  u <- pop$simulation$unrealized_births
  lost <- if (is.null(u)) rep(NA, 3) else
    vapply(c(0, 1, 2), function(k) sum(u$unrealized_births[u$parity_start == k]), 0)
  data.frame(scenario = "cohort", method, seed, measure = "rate",
    window = "[0,10)", ages = "[15,Inf)", stratum = paste0("parity ", c("0", "1", "2+")),
    target = c(.2, .3, .4), value = r$rate_fertility, births = r$births,
    exposure = r$person_years, unrealized = lost)
}

run_young <- function(method, seed) {
  model <- hammel_model() |> set_population(founders(300)) |>
    set_mortality(mort_zero) |>
    set_fertility(parity_rates(40), spacing_method = method) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 100, seed)
  one <- function(ages, breaks, labels, targets) {
    r <- calculate_fertility_rates(pop, time_breaks = c(80, 100), age_breaks = ages,
                                   parity_breaks = breaks)
    # A finite final age break gains an open band; keep the requested ages.
    r <- r[r$age_start == ages[1] & r$time_start == 80, ]
    data.frame(scenario = "young", method, seed, measure = "rate", window = "[80,100)",
      ages = sprintf("[%g,%g)", ages[1], ages[2]), stratum = labels, target = targets,
      value = r$rate_fertility, births = r$births, exposure = r$person_years,
      unrealized = NA)
  }
  rbind(one(c(18, 19), c(0, 1, 2), c("parity 0", "parity 1", "parity 2+"), c(.2, .3, .4)),
        one(c(18, 19), c(0, 1, 2, 3, 4), paste("actual parity", c(0:3, "4+")), c(.2, .3, .4, .4, .4)),
        one(c(20, 35), c(0, 1, 2), c("parity 0", "parity 1", "parity 2+"), c(.2, .3, .4)))
}

baseline <- function(n, seed, groups = FALSE) {
  set.seed(seed)
  pop <- create_baseline_population(n, max_age = 70)
  if (groups) pop$persons$group <- sample(c("A", "B"), n, TRUE)
  pop
}

run_restart <- function(method, seed) {
  model <- hammel_model() |> set_population(baseline(5000, seed)) |>
    set_mortality(mort_life) |> set_fertility(low, spacing_method = method) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  first <- simulate(model, 40, seed)
  second <- simulate(set_population(model, first) |>
    set_fertility(niger, birth_interval = 1.25, spacing_method = method), 6, seed + 1000)
  a <- recovery(first, 35:40, low)
  b <- recovery(second, 40:46, niger)
  ab <- rbind(a, b)
  data.frame(scenario = "restart", method, seed, measure = "recovery ratio",
    window = sprintf("[%g,%g)", ab$time_start, ab$time_start + 1), ages = "[15,50)",
    stratum = ifelse(ab$time_start < 40, "TFR ~2 schedule", "TFR ~7.7, spacing 1.25"),
    target = 1, value = ab$value, births = ab$births, exposure = NA, unrealized = NA)
}

run_union <- function(method, seed) {
  status <- rbind(transform(low, union_type = "partnered", rate = rate * 1.6),
                  transform(low, union_type = "unpartnered", rate = rate * .3))
  model <- hammel_model() |> set_population(baseline(6000, seed)) |>
    set_mortality(mort_life) |> set_fertility(status, spacing_method = method) |>
    set_union_formation("rate", rates = data.frame(age_start = 15, age_end = Inf, rate = .2),
                        cancel_birth_if_no_match = FALSE) |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .06))
  pop <- simulate(model, 60, seed)
  r <- recovery(pop, c(40, 60), status, by = "union_type")
  data.frame(scenario = "union", method, seed, measure = "recovery ratio",
    window = "[40,60)", ages = "[15,50)", stratum = r$stratum, target = 1,
    value = r$value, births = r$births, exposure = NA, unrealized = NA)
}

run_groups <- function(method, seed) {
  multiplier <- c(A = 1, B = 1.4)
  model <- hammel_model() |> set_population(baseline(6000, seed, groups = TRUE)) |>
    set_mortality(mort_life) |>
    set_fertility(low, group_prop_hazard = multiplier, spacing_method = method) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(model, 60, seed)
  r <- recovery(pop, c(40, 60), low, by = "group", multiplier = multiplier)
  data.frame(scenario = "groups", method, seed, measure = "recovery ratio",
    window = "[40,60)", ages = "[15,50)", stratum = paste("group", r$stratum), target = 1,
    value = r$value, births = r$births, exposure = NA, unrealized = NA)
}

jobs <- expand.grid(scenario = c("cohort", "young", "restart", "union", "groups"),
                    method = methods, seed = seeds, stringsAsFactors = FALSE)
runners <- list(cohort = run_cohort, young = run_young, restart = run_restart,
                union = run_union, groups = run_groups)
results <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  job <- jobs[j, ]
  suppressWarnings(runners[[job$scenario]](job$method, job$seed))
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) stop(paste(unique(unlist(results[failed])), collapse = "\n"))
replications <- do.call(rbind, results)
write.csv(replications, file.path(out_dir, "cell-spacing-replications.csv"), row.names = FALSE)

key <- c("scenario", "measure", "window", "ages", "stratum", "target", "method")
groups <- split(replications, replications[key], drop = TRUE)
summary <- do.call(rbind, lapply(groups, function(d) {
  ok <- is.finite(d$value)
  cbind(d[1, key], n_valid = sum(ok), mean = mean(d$value[ok]),
        mcse = sd(d$value[ok]) / sqrt(sum(ok)),
        mean_unrealized = mean(d$unrealized))
}))
summary <- summary[order(summary$scenario, summary$window, summary$ages, summary$stratum,
                         summary$method), ]
write.csv(summary, file.path(out_dir, "cell-spacing-summary.csv"), row.names = FALSE)
print(summary[c("scenario", "window", "ages", "stratum", "target", "method", "n_valid",
                "mean", "mcse", "mean_unrealized")], row.names = FALSE, digits = 4)
