# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)
# Save once with saveRDS(); use readRDS() here for subsequent sessions.

set.seed(42)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)

timing <- system.time(result <- simulate(model, years = 100, seed = 42))
print(result$simulation)
cat("Living persons:", sum(is.na(result$persons$death_time)), "\n")
cat("Simulation elapsed seconds:", timing[["elapsed"]], "\n")

# get lexis diagram counts in 5 year time intervals and one year age intervals
lexis_counts <- result |> calculate_lexis_counts(seq(2000, 2100, 10), seq(0, 100, 1))

ggplot(lexis_counts, aes(x = age_start, y = asmr))+
  geom_line(aes(group = factor(time_start), color = factor(time_start)))+
  geom_line(data = model$processes$mortality$rates, aes(y = rate), linewidth = 2, alpha = 0.5)+
  facet_wrap(~sex)+
  scale_y_log10()+
  theme_bw()

model$processes$fertility |> calculate_tfr()
lexis_counts |> calculate_tfr()

model$processes$mortality |> 
  calculate_life_table() |>
  filter(age_start == 0) |>
  select(age_start, sex, ex)

lexis_counts |> 
  calculate_life_table() |>
  filter(age_start == 0) |>
  select(time_start, time_end, age_start, sex, ex)

result |> calculate_population_size(seq(from = 2000, to = 2100, by = 10)) |>
  ggplot(aes(x = time, y = population))+
  geom_point()+
  geom_line()

# Results can initialize a subsequent model with the same or new processes:
result2 <- model |>
  set_population(result) |>
  simulate(100, seed = 13)

result2 |> calculate_population_size(seq(from = 2000, to = 2200, by = 10)) |>
  ggplot(aes(x = time, y = population))+
  geom_point()+
  geom_line()

