library(here)
#devtools::load_all(".")
devtools::install()
library(hammel)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(gt)

# Get rates ------------------------------------------------------------

mortality_rates_usa2000 <- get_mortality_hmd("USA", year = 2000)
fertility_rates_overall_usa <- read_csv(
  here("validation", "data", "usa_2000_fertility_overall.csv")
)
fertility_rates_by_status_usa <- read_csv(
  here("validation", "data", "usa_2000_fertility_by_union_status.csv")
) |>
  rename(union_type = union_status)
mar_rates_female_usa <- read_csv(
  here("validation", "data", "usa_2001_marriage_rates_female.csv")
)

# coarse dissolution rates
# derived from NSFG https://www.cdc.gov/nchs/data/nhsr/nhsr049.pdf
dissolution_rates <- tibble(
  dur_start = c(0, 5, 10, 15),
  dur_end   = c(5, 10, 15, 20),
  rate = -log(
    c(0.80, 0.68, 0.60, 0.52) /
    c(1.00, 0.80, 0.68, 0.60)
  ) / 5
)

# Dissoluation rates -----------------------------------------------------

# start with simple one group case

pop_base <- create_baseline_population(50000, 50, 2000)

model <- hammel_model() |>
  set_population(pop_base) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_female_usa) |>
  set_union_dissolution(dissolution_rates)

system.time(result <- simulate(model, 100, 42))

table(result$unions$end_reason, exclude = NULL)

union_counts <- result |>
  calculate_dissolution_rates(
    time_breaks = seq(2025, 2100, 75),
    dur_breaks = c(0, 5, 10, 15, 20, Inf)
  )

union_counts |>
  filter(dur_start < 20) |>
  ggplot(aes(x = dur_start + 2.5, y = rate_dissolution))+
  geom_point()+
  geom_line(data = dissolution_rates, aes(y = rate), linetype = 2)+
  theme_bw()

# now try separate groups
pop_base_grp2 <- create_baseline_population(50000, 50, 2000, group_weights = c("A" = 1, "B" = 1))

model <- hammel_model() |>
  set_population(pop_base_grp2) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_female_usa) |>
  set_union_dissolution(dissolution_rates, group_prop_hazard = c("A:B" = 1.1))

system.time(result <- simulate(model, 100, 42))

table(result$unions$end_reason, exclude = NULL)

union_counts <- result |>
  calculate_dissolution_rates(
    time_breaks = seq(2025, 2100, 75),
    dur_breaks = c(0, 5, 10, 15, 20, Inf),
    by_group = TRUE
  )

dissolution_rates_grp2 <- bind_rows(
  dissolution_rates |> mutate(group_pair = "A:A"),
  dissolution_rates |> mutate(group_pair = "B:B"),
  dissolution_rates |> mutate(group_pair = "A:B", rate = rate * 1.1)
)

union_counts |>
  filter(dur_start < 20) |>
  ggplot(aes(x = dur_start + 2.5, y = rate_dissolution, 
             group = group_pair, color = group_pair))+
  geom_point()+
  geom_line(data = dissolution_rates_grp2, aes(y = rate), linetype = 2)+
  theme_bw()

model <- hammel_model() |>
  set_population(pop_base_grp2) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_female_usa) |>
  set_union_dissolution(dissolution_rates_grp2)

system.time(result <- simulate(model, 100, 42))

table(result$unions$end_reason, exclude = NULL)

union_counts <- result |>
  calculate_dissolution_rates(
    time_breaks = seq(2025, 2100, 75),
    dur_breaks = c(0, 5, 10, 15, 20, Inf),
    by_group = TRUE
  )

union_counts |>
  filter(dur_start < 20) |>
  ggplot(aes(x = dur_start + 2.5, y = rate_dissolution, 
             group = group_pair, color = group_pair))+
  geom_point()+
  geom_line(data = dissolution_rates_grp2, aes(y = rate), linetype = 2)+
  theme_bw()

# Remarriage rates -----------------------------------------------------------

# vary remarriage rates but all are age specific

mar_rates_by_prev <- bind_rows(
  mar_rates_female_usa |> mutate(previous_union = "none"),
  mar_rates_female_usa |> mutate(previous_union = "dissolved", rate = rate * 1.1),
  mar_rates_female_usa |> mutate(previous_union = "widowed", rate = rate * 0.9)
)

model <- hammel_model() |>
  set_population(pop_base) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_by_prev) |>
  set_union_dissolution(dissolution_rates)

system.time(result <- simulate(model, 100, 42))

table(result$unions$end_reason, exclude = NULL)

mar_counts <- calculate_marriage_rates(
  result,
  time_breaks = seq(2025, 2100, 75),
  age_breaks = seq(15, 70, 5), 
  by_previous_union = TRUE
)

mar_counts |>
  filter(sex == "f") |>
  ggplot(aes(x = age_start + 2.5, y = rate_marriage, group = previous_union, color = previous_union))+
  geom_point()+
  geom_line(data = mar_rates_by_prev, aes(y = rate), linetype = 2)+
  theme_bw()

# use duration specific rates - for testing make up some plausible values

mar_rates_dissolved <- tribble(
  ~dur_start, ~dur_end, ~rate,
  0,          1,        0.005,
  1,          5,        0.040,
  5,          15,       0.070,
  15,         Inf,      0.010
)

mar_rates_widowed <- tribble(
  ~dur_start, ~dur_end, ~rate,
  0,          1,        0.005,
  1,          5,        0.060,
  5,          15,       0.030,
  15,         Inf,      0.010
)

mar_rates_by_dur <- bind_rows(
  mar_rates_female_usa |> mutate(previous_union = "none"),
  mar_rates_dissolved |> mutate(previous_union = "dissolved"),
  mar_rates_widowed |> mutate(previous_union = "widowed")
) |>
  # don't have to resort them, but for human eyes its nice
  select(previous_union, starts_with("age"), starts_with("dur"), rate)

model <- hammel_model() |>
  set_population(pop_base) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_by_dur) |>
  set_union_dissolution(dissolution_rates)

system.time(result <- simulate(model, 100, 42))

table(result$unions$end_reason, exclude = NULL)

mar_counts <- calculate_marriage_rates(
  result,
  time_breaks = seq(2025, 2100, 75),
  age_breaks = seq(15, 70, 5),
  dur_breaks = c(0, 1, 5, 15)
)

mar_counts |>
  filter(sex == "f" & previous_union != "none") |>
  ggplot(aes(x = dur_start + 2.5, y = rate_marriage, group = previous_union, color = previous_union))+
  geom_point()+
  geom_line(data = mar_rates_by_dur |> filter(previous_union != "none"), 
            aes(y = rate), linetype = 2)+
  theme_bw()

mar_counts |>
  filter(sex == "f" & previous_union == "none") |>
  ggplot(aes(x = age_start + 2.5, y = rate_marriage))+
  geom_point()+
  geom_line(data = mar_rates_by_dur |> filter(previous_union == "none"), 
            aes(y = rate), linetype = 2)+
  theme_bw()


# Check Cohabitation -----------------------------------------------------------------------------------

# what all do we need here?
# - probablity schedule by age (and possibly previous_union and group) for cohabitations
# - fertility rates with cohabs and marriages separated rather than a single partnered
# - dissolution rates that are stratified by union_type
# - conversion rates from cohabitation to marriage, duration or age dependent

# make up some plausible probabilities for testing
prob_cohab <- tribble(
  ~age_start, ~age_end, ~prob,
  15,         20,       0.4,
  20,         25,       0.25,
  25,         30,       0.15,
  30,         35,       0.1,
  35,         40,       0.05,
  40,         Inf,      0.01
)

fertility_rates_by_union_type <- bind_rows(
  # unpartnered is the same
  fertility_rates_by_status_usa |> filter(union_type == "unpartnered"),
  # married uses the partnered version
  fertility_rates_by_status_usa |> filter(union_type == "partnered") |> mutate(union_type = "marriage"),
  # cohab is half the married rate
  fertility_rates_by_status_usa |> filter(union_type == "partnered") |> mutate(union_type = "cohab", rate = 0.5 * rate)
)

dissolution_rates_by_union_type <- bind_rows(
  # double the rate of dissolution for cohabitations
  dissolution_rates |> mutate(union_type = "cohab", rate = 2 * rate),
  dissolution_rates |> mutate(union_type = "marriage")
)

# lets do duration dependent conversion rates from cohab to marriage
convert_rate <- tribble(
  ~dur_start, ~dur_end, ~rate,
  0,          5,        0.15,
  5,          10,       0.30,
  10,         Inf,      0.05
)


model <- hammel_model() |>
  set_population(pop_base) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_union_type) |>
  set_union_formation(
    "rate", 
    rates = mar_rates_by_prev, 
    cohab_prob = prob_cohab, 
    union_convert_rate = convert_rate
  ) |>
  set_union_dissolution(dissolution_rates_by_union_type)

system.time(result <- simulate(model, 100, 42))

table(result$unions$union_type, result$unions$end_reason)

mar_counts <- calculate_marriage_rates(
  result,
  time_breaks = seq(2025, 2100, 75),
  age_breaks = seq(15, 70, 5),
  by_previous_union = TRUE
)

union_counts <- calculate_union_rates(
  result,
  time_breaks = seq(2025, 2100, 75),
  age_breaks = seq(15, 70, 5),
  by_union_type = TRUE, 
  by_previous_union = TRUE
)

union_counts |>
  filter(sex == "f") |>
  ggplot(aes(x = age_start + 2.5, y = rate_union, color = union_type))+
  geom_point()+
  geom_line()+
  facet_wrap(~previous_union)+
  theme_bw()

# input mar rates should match overall union rates
calculate_union_rates(
  result,
  time_breaks = seq(2025, 2100, 75),
  age_breaks = seq(15, 70, 5),
  by_previous_union = TRUE 
) |>
  filter(sex == "f") |>
  ggplot(aes(x = age_start + 2.5, y = rate_union, group = previous_union, color = previous_union))+
  geom_point()+
  geom_line(data = mar_rates_by_prev, aes(y = rate), linetype = 2)


dissolution_counts <- result |>
  calculate_dissolution_rates(
    time_breaks = seq(2025, 2100, 75),
    dur_breaks = c(0, 5, 10, 15, 20, Inf),
    by_union_type = TRUE
  )

dissolution_counts |>
  filter(dur_start < 20) |>
  ggplot(aes(x = dur_start, y = rate_dissolution, group = union_type, color = union_type))+
  geom_point()+
  geom_line(data = dissolution_rates_by_union_type, aes(y = rate), linetype = 2)+
  theme_bw()

calculate_demog_rates(
  result, 
  time_breaks = seq(2025, 2100, 75), 
  age_breaks = seq(15, 40, 5), 
  by_union_type = TRUE
) |>
  filter(sex == "f") |>
  ggplot(aes(x = age_start, y = rate_fertility, group = union_type, color = union_type))+
  geom_point()+
  geom_line(data = fertility_rates_by_union_type, aes(y = rate), linetype = 2)+
  theme_bw()

devtools::load_all(debug = FALSE, recompile = TRUE)
