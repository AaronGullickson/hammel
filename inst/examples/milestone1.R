# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
library(hammel)

set.seed(42)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality_hmd("USA", year = 2000) |>
  set_fertility_hfd("USA", year = 2000)

timing <- system.time(result <- simulate(model, years = 100, seed = 42))
print(result$simulation)
cat("Living persons:", sum(is.na(result$persons$death_time)), "\n")
cat("Simulation elapsed seconds:", timing[["elapsed"]], "\n")

library(dplyr)
library(ggplot2)
p <- result$persons |> 
  mutate(
    living = is.na(death_time),
    age = if_else(
      living, 
      result$end_time - birth_time,
      death_time - birth_time
    )
  ) |>
  as_tibble() |>
  arrange(age)

p |>
  filter(living) |>
  ggplot(aes(x = age)) +
  geom_histogram(breaks = seq(from = 0, to = 100, by = 5), col = "grey20", fill = "grey80")

# Results can initialize a subsequent model with the same or new processes:
continued_model <- hammel_model() |>
  set_population(result) |>
  set_mortality(model$processes$mortality$rates) |>
  set_fertility(model$processes$fertility$rates)
