# General demographic fixtures explicitly preserve unmatched births. Tests of
# cancellation and constructors separately exercise the production default TRUE.
make_model <- function(n = 1000, sex = "m", age = 20, death_rate = 0,
           birth_rate = 0, birth_interval = 0.75) {
 persons <- data.frame(person_id = seq_len(n), sex = rep(sex, n),
  birth_time = rep(-age, n), death_time = rep(NA_real_, n),
  mother_id = rep(NA_integer_, n), father_id = rep(NA_integer_, n))
 hammel_model() |>
  set_population(hammel_population(persons, end_time = 0)) |>
  set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
              age_end = Inf, rate = death_rate)) |>
  set_fertility(data.frame(age_start = 15, age_end = 50, rate = birth_rate),
         birth_interval = birth_interval) |>
  set_union_formation("on_demand",
    create_matching_model(~ age_difference, c(age_difference = 0)),
    cancel_birth_if_no_match = FALSE)
}

make_rate_model <- function(n = 2, search_rate = 1, birth_rate = 0, birth_interval = 0) {
  model <- make_model(n, birth_rate = birth_rate, birth_interval = birth_interval)
  if (n) model$population$persons$sex[1] <- "f"
  rates <- model$processes$fertility$rates
  model <- set_fertility(model, rbind(
    transform(rates, union_type = "partnered"),
    transform(rates, union_type = "unpartnered", rate = 0)), birth_interval)
  set_union_formation(model, "rate",
    create_matching_model(~ age_difference, c(age_difference = 0)),
    rates = data.frame(age_start = 15, age_end = Inf, rate = search_rate),
    cancel_birth_if_no_match = FALSE)
}

# Simple synthetic parity schedules often give positive rates in young cells
# that no eligible woman can reach; silence only that advisory warning.
quiet_unreachable <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("no woman can be eligible", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

# Models that differ only in how equal rates are stratified create different
# cells, so their histories differ draw by draw while agreeing in distribution.
# Compare union and dissolution counts over replications.
expect_equivalent_unions <- function(a, b, years, seeds = 1:20) {
  counts <- function(model) vapply(seeds, function(seed) {
    unions <- simulate(model, years, seed)$unions
    c(nrow(unions), sum(unions$end_reason %in% "dissolution"))
  }, numeric(2))
  x <- counts(a)
  y <- counts(b)
  for (k in 1:2) {
    se <- sqrt(var(x[k, ]) / length(seeds) + var(y[k, ]) / length(seeds))
    expect_lt(abs(mean(x[k, ]) - mean(y[k, ])), 4 * se + 1e-9)
  }
}
