# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)
library(purrr)

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)
# Save once with saveRDS(); use readRDS() here for subsequent sessions.

set.seed(42)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)

results <- map(1:100, function(x) {simulate(model, years = 100, seed = x)})

lexii <- map(results, calculate_lexis_counts, seq(2000, 2100, 10), seq(0, 100, 1))

tfr <- map(lexii, calculate_tfr) |>
  bind_rows()

life_exp <- map(lexii, function(lexis) {
  lexis |> 
    calculate_life_table() |>
    filter(age_start == 0) |>
    select(time_start, time_end, age_start, sex, ex)
}) |> 
  bind_rows()

tfr_summary <- tfr |>
  group_by(time_start) |>
  summarize(tfr_median = median(tfr),
            tfr_mean = mean(tfr))

tfr_real <- fertility_rates |> calculate_tfr()

ggplot(tfr, aes(x = time_start, y = tfr))+
  geom_point(alpha = 0.2)+
  geom_line(data = tfr_summary, aes(y = tfr_median), color = "navy")+
  geom_line(data = tfr_summary, aes(y = tfr_mean), color = "tomato")+
  geom_hline(yintercept = tfr_real$tfr, linetype = 2)
