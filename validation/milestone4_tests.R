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
)
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

# vary remarriage rates

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
