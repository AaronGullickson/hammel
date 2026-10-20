# load libraries --------------------------------------------------

devtools::load_all()
library(ggplot2)
library(dplyr)
library(here)
library(readr)
library(tidyr)

# Get rates --------------------------------------------------------

mortality_rates_usa2000 <- get_mortality_hmd("USA", year = 2000)
mortality_rates_usa1950 <- get_mortality_hmd("USA", year = 1950)

mortality_rates_grouped <- bind_rows(
  mortality_rates_usa2000 |> mutate(group = "A"),
  mortality_rates_usa1950 |> mutate(group = "B")
)

fertility_rates_overall_usa <- read_csv(
  here("validation", "data", "usa_2000_fertility_overall.csv")
)
fertility_rates_by_status_usa <- read_csv(
  here("validation", "data", "usa_2000_fertility_by_union_status.csv")
)
fertility_rates_overall_england <- read_csv(
  here("validation", "data", "england_wales_2015_overall_fertility.csv")
)
fertility_rates_by_status_england <- read_csv(
  here("validation", "data", "england_wales_2015_fertility_by_union_status.csv")
)

mar_rates_female_usa <- read_csv(
  here("validation", "data", "usa_2001_marriage_rates_female.csv")
)
mar_rates_england <- read_csv(
  here("validation", "data", "england_wales_2015_female_union_rates.csv")
)

calculate_life_table(mortality_rates_grouped) |>
  filter(age_start == 0) |>
  select(sex, group, ex)

calculate_tfr(fertility_rates_overall_usa)

# Test simulations -----------------------------------------------------------

pop_baseline <- create_baseline_population(20000, max_age = 70, end_time = 2000, 
                                           group_weights = c("A" = 2, "B" = 1))

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation("on_demand") |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100))

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

lexis_counts <- calculate_lexis_counts(result, seq(2000, 2100, 10), 0:100, by_group = TRUE)

lexis_counts |> calculate_tfr()

lexis_counts |> calculate_life_table() |>
  filter(age_start == 0, !is.na(ex)) |>
  select(time_start, time_end, sex, group, ex) |>
  print(n = 40)

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_grouped) |>
  set_union_formation("on_demand") |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100))

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

lexis_counts <- calculate_lexis_counts(result, seq(2000, 2100, 10), 0:100, by_group = TRUE)

lexis_counts |> calculate_tfr()

lexis_counts |> calculate_life_table() |>
  filter(age_start == 0, !is.na(ex)) |>
  select(time_start, time_end, sex, group, ex) |>
  print(n = 40)

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa, group_prop_hazard = c("A" = 1.1, "B" = 0.9)) |>
  set_mortality(mortality_rates_usa2000, group_prop_hazard = c("B" = 1.8)) |>
  set_union_formation("on_demand") |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100))

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

lexis_counts <- calculate_lexis_counts(result, seq(2000, 2100, 10), 0:85, by_group = TRUE)

lexis_counts |> calculate_tfr()

lexis_counts |> calculate_life_table() |>
  filter(age_start == 0, !is.na(ex)) |>
  select(time_start, time_end, sex, group, ex) |>
  print(n = 40)
