# Updated for the required union process in milestone 2.
# For the original milestone 1 behavior, use the accepted report's pinned commit.
# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)
library(here)
library(readr)

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)
# Use data on England and Wales that provides marital and non-marital rates

marriage_rates <- read_csv(
  here("validation", "data", "england_wales_2015_female_union_rates.csv")
)
fertility_rates_status <- read_csv(
  here("validation", "data", "england_wales_2015_fertility_by_union_status.csv")
)
fertility_rates <- read_csv(
  here("validation", "data", "england_wales_2015_overall_fertility.csv")
)

fertility_rates |> calculate_tfr()
fertility_rates_status |> calculate_tfr()

set.seed(42)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates_status) |>
  set_union_formation(
    "rate",
    create_matching_model(
      ~ age_difference + I(age_difference^2),
      c(age_difference = 0.072, "I(age_difference^2)" = -0.014)
    ),
    rates = marriage_rates,
    kinship_generations = 2
  )

timing <- system.time(result <- simulate(model, years = 100, seed = 42))
print(result$simulation)
cat("Living persons:", sum(is.na(result$persons$death_time)), "\n")
cat("Simulation elapsed seconds:", timing[["elapsed"]], "\n")

husbands <- result$persons |>
  filter(sex == "m") |>
  rename(
    husband_id = person_id,
    birth_time_h = birth_time,
    death_time_h = death_time
  ) |>
  select(husband_id, birth_time_h, death_time_h)

wives <- result$persons |>
  filter(sex == "f") |>
  rename(
    wife_id = person_id,
    birth_time_w = birth_time,
    death_time_w = death_time
  ) |>
  select(wife_id, birth_time_w, death_time_w)

unions <- result$unions |>
  left_join(husbands) |>
  left_join(wives) |>
  mutate(
    age_mar_h = start_time - birth_time_h,
    age_mar_w = start_time - birth_time_w,
    age_diff = age_mar_h - age_mar_w,
    decade = floor(start_time / 10) * 10
  )

temp <- unions |>
  group_by(decade) |>
  summarize(age_diff = mean(age_diff))

ggplot(unions, aes(x = factor(decade), y = age_diff))+
  geom_boxplot()+
  geom_jitter(alpha = 0.02)+
  geom_hline(yintercept = 0.072 / (2 * 0.014), color = "red", linetype = 2)+
  theme_bw()


# how many remain unmarried?
result$persons |>
  mutate(
    ever_married = person_id %in% c(result$unions$husband_id, result$unions$wife_id),
    birth_decade = floor(birth_time / 10) * 10
  ) |>
  group_by(birth_decade, sex) |>
    summarize(prop_ever_married = mean(ever_married)) |>
  ggplot(aes(x = birth_decade, y = prop_ever_married, group = sex, color = sex))+
  geom_line()+
  geom_point()+
  theme_bw()


lexis <- calculate_lexis_counts(result, seq(2000, 2100, 5), 0:100)

lexis |> calculate_tfr()

