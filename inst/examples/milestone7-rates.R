# Calendar-varying annual rates with fractional years (D118).
# All numbers below are illustrative, not empirical rates or defaults.
library(hammel)

set.seed(42)
population <- create_baseline_population(50000, max_age = 60,
                                         end_time = 2015.5)

mortality_rates <- data.frame(
  sex = rep(c("m", "f"), 2),
  age_start = 0, age_end = Inf,
  time_start = rep(c(2015.5, 2017.25), each = 2),
  time_end = rep(c(2017.25, Inf), each = 2),
  rate = c(.012, .010, .009, .008)
)
fertility_rates <- data.frame(
  age_start = 20, age_end = 40,
  time_start = c(2015.5, 2018.1),
  time_end = c(2018.1, Inf),
  rate = c(.10, .15)
)

model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates, birth_interval = .75) |>
  set_union_formation("on_demand")

# Starts at 2015.5, finishes at 2019.8; events at the endpoint are excluded.
result <- simulate(model, years = 4.3, seed = 123)
result$simulation
calculate_fertility_rates(result, time_breaks = c(2015.5, 2018.1, 2019.8),
                          age_breaks = c(20, 40))
calculate_mortality_rates(result, time_breaks = c(2015.5, 2018.1, 2019.8),
                          age_breaks = c(0, Inf))
