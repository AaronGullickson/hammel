# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)
library(purrr)

# Setup -----------------------------------------------------------------------

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)
# Save once with saveRDS(); use readRDS() here for subsequent sessions.

set.seed(42)
population <- create_baseline_population(100000, max_age = 100, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)

# Run simulation --------------------------------------------------------------

results <- map(1:100, function(x) {simulate(model, years = 100, seed = x)})

lexii <- map(
  results, 
  calculate_lexis_counts, c(2000:2009, seq(2010, 2100, 10)), seq(0, 100, 1)
)

# Assess fertility rates ------------------------------------------------------

tfr <- map(lexii, calculate_tfr) |>
  bind_rows()

tfr_summary <- tfr |>
  mutate(time_mid = (time_start + time_end) / 2) |>
  group_by(time_mid) |>
  summarize(
    median = median(tfr),
    mean = mean(tfr),
    p05 = quantile(tfr, 0.05),
    p25 = quantile(tfr, 0.25),
    p75 = quantile(tfr, 0.75),
    p95 = quantile(tfr, 0.95),
    .groups = "drop"
  )

tfr_real <- fertility_rates |> calculate_tfr()

ggplot(tfr_summary, aes(x = time_mid))+
  geom_ribbon(aes(ymin = p05, ymax = p95), color = NA, fill = "grey80")+
  geom_ribbon(aes(ymin = p25, ymax = p75), color = NA, fill = "grey")+
  geom_line(aes(y = median), color = "tomato")+
  geom_line(aes(y = mean), color = "navy")+
  geom_hline(yintercept = tfr_real$tfr, linetype = 2)+
  theme_bw()

asfr <- map(lexii, function(lexis) {
  lexis |>
    filter(sex == "f", age_start > 12 & age_end < 56) |>
    mutate(time_mid = (time_start + time_end) / 2) |>
    select(time_mid, age_start, age_end, asfr)
}) |>
  bind_rows()

asfr_summary <- asfr |>
  group_by(time_mid, age_start) |>
  summarize(
    median = median(asfr),
    mean = mean(asfr),
    p05 = quantile(asfr, 0.05),
    p25 = quantile(asfr, 0.25),
    p75 = quantile(asfr, 0.75),
    p95 = quantile(asfr, 0.95),
    .groups = "drop"
  )

ggplot(asfr_summary, aes(x = age_start))+
  geom_line(aes(y = mean, group = time_mid))+
  geom_line(
    data = fertility_rates, 
    aes(y = rate), 
    linewidth = 2,
    color = "red", 
    alpha = 0.4
  )+
  theme_bw()


# Assess mortality rates ------------------------------------------------------

life_exp <- map(lexii, function(lexis) {
  lexis |> 
    calculate_life_table() |>
    filter(age_start == 0) |>
    select(time_start, time_end, age_start, sex, ex)
}) |> 
  bind_rows()

life_exp_summary <- life_exp |>
  mutate(time_mid = (time_start + time_end) / 2) |>
  group_by(time_mid, sex) |>
  summarize(
    median = median(ex, na.rm = TRUE),
    mean = mean(ex, na.rm = TRUE),
    p05 = quantile(ex, 0.05, na.rm = TRUE),
    p25 = quantile(ex, 0.25, na.rm = TRUE),
    p75 = quantile(ex, 0.75, na.rm = TRUE),
    p95 = quantile(ex, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

ex_real <- mortality_rates |> 
  calculate_life_table() |>
  filter(age_start == 0) |>
  select(sex, ex)

# this comparison won't be exact because I truncate the open ended interval to 100 not 110
# because sparseness can lead to NAs in those cells
ggplot(life_exp_summary, aes(x = time_mid))+
  geom_ribbon(aes(ymin = p05, ymax = p95), color = NA, fill = "grey80")+
  geom_ribbon(aes(ymin = p25, ymax = p75), color = NA, fill = "grey")+
  geom_line(aes(y = median), color = "tomato")+
  geom_line(aes(y = mean), color = "navy")+
  geom_hline(data = ex_real, aes(yintercept = ex), linetype = 2)+
  facet_wrap(~sex)+
  theme_bw()

# Lets redo lexis curve for ASMR comparison out to 110 and larger time chunks
lexii <- map(
  results, 
  calculate_lexis_counts, seq(2000, 2100, 20), seq(0, 110, 1)
)

asmr <- map(lexii, function(lexis) {
  lexis |>
    mutate(time_mid = (time_start + time_end) / 2) |>
    select(time_mid, age_start, age_end, sex, asmr)
}) |>
  bind_rows()

asmr_summary <- asmr |>
  group_by(time_mid, age_start, sex) |>
  summarize(
    median = median(asmr, na.rm = TRUE),
    mean = mean(asmr),
    p05 = quantile(asmr, 0.05, na.rm = TRUE),
    p25 = quantile(asmr, 0.25, na.rm = TRUE),
    p75 = quantile(asmr, 0.75, na.rm = TRUE),
    p95 = quantile(asmr, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(asmr_summary, aes(x = age_start))+
  geom_line(aes(y = mean, group = time_mid))+
  geom_line(
    data = mortality_rates, 
    aes(y = rate), 
    linewidth = 2,
    color = "red", 
    alpha = 0.4
  )+
  facet_wrap(~sex)+
  scale_y_log10()+
  theme_bw()
