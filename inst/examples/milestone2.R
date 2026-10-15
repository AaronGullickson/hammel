# Updated for the required union process in milestone 2.
# For the original milestone 1 behavior, use the accepted report's pinned commit.
# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)

set.seed(42)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates) |>
  set_union_formation("on_demand",
    create_matching_model(~ age_difference, c(age_difference = 0)))

timing <- system.time(result <- simulate(model, years = 100, seed = 42))
print(result$simulation)
cat("Living persons:", sum(is.na(result$persons$death_time)), "\n")
cat("Simulation elapsed seconds:", timing[["elapsed"]], "\n")
