# Updated for the required union process in milestone 2.
# For the original milestone 1 behavior, use the accepted report's pinned commit.
# Requires HMDHFDplus and HMD/HFD credentials in the R environment.
#library(hammel)
devtools::load_all()
library(ggplot2)
library(dplyr)
library(here)
library(readr)
library(tidyr)

mortality_rates <- get_mortality_hmd("USA", year = 2000)
fertility_rates <- get_fertility_hfd("USA", year = 2000)
# Use data on England and Wales that provides marital and non-marital rates

marriage_rates <- read_csv(
  here("validation", "data", "england_wales_2015_female_union_rates.csv")
)
fertility_rates_status <- read_csv(
  here("validation", "data", "england_wales_2015_fertility_by_union_status.csv")
) |>
  rename(union_type = union_status)
fertility_rates <- read_csv(
  here("validation", "data", "england_wales_2015_overall_fertility.csv")
)

fertility_rates |> calculate_tfr()
fertility_rates_status |> calculate_tfr()

set.seed(39)
population <- create_baseline_population(50000, max_age = 70, end_time = 2000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates_status) |>
  #set_fertility(fertility_rates)
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

#result <- population
#timing <- system.time(for(i in 1:100) {
#  model <- model |> set_population(result)
#  result <- simulate(model, years = 1, seed = 27)
#})

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


lexis <- calculate_demog_rates(result, seq(2000, 2100, 5), seq(15, 65, 5), by_union_type = TRUE)
marriages <- calculate_union_rates(result, seq(2000, 2100, 5), seq(15, 65, 5)) |>
  group_by(time_start, time_end, age_start, age_end, sex) |>
  summarize(rate_union = if (sum(person_years) > 0)
    sum(unions) / sum(person_years) else NA_real_, .groups = "drop")


lexis |> calculate_tfr()


asfr <- lexis |>
  filter(sex == "f") |>
  select(time_start, time_end, age_start, age_end, union_type,
         rate = rate_fertility)


asfr |>
  filter(age_start >= 15 & age_start <=45, time_start > 2020) |>
  group_by(union_type, age_start) |>
  summarize(rate = mean(rate), .groups = "drop") |>
  ggplot(aes(x = age_start, y = rate, group = union_type, color = union_type))+
  geom_point()+
  geom_line(linetype = 2)+
  geom_line(data = fertility_rates_status)+
  theme_bw()

lexis |>
  filter(sex == "f", union_type == "total", age_start <= 45) |>
  group_by(age_start) |>
  summarize(rate = mean(rate_fertility, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = age_start, y = rate))+
  geom_point()+
  geom_line(linetype = 2)+
  geom_line(data = fertility_rates)+
  theme_bw()

marriages |>
  filter(sex == "f") |>
  group_by(age_start) |>
  summarize(rate = mean(rate_union), .groups = "drop") |>
  ggplot(aes(x = age_start, y = rate))+
  geom_point()+
  geom_line(linetype = 2)+
  geom_line(data = marriage_rates)+
  theme_bw()
  
marriages |>
  group_by(age_start, sex) |>
  summarize(rate = mean(rate_union), .groups = "drop") |>
  ggplot(aes(x = age_start, y = rate, group = sex, color = sex))+
  geom_point()+
  geom_line()+
  theme_bw()

tibble(
  age_diff = seq(from = -20, to = 20, by = 0.1),
  logit = 0.025 * age_diff - 0.005 * age_diff^2
) |>
  ggplot(aes(x = age_diff, y = logit))+
  geom_line()+
  geom_vline(xintercept = 0.072 / (2 * 0.014), linetype = 2)+
  theme_bw()

# if h is halving time of the utility and d is the preferred age difference where utility is maximized, then:
# b_quadratic = - log(2)/h^2
# b_linear = -2 * d * b_quadratic
