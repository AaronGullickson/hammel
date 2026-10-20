# load libraries --------------------------------------------------

devtools::load_all()
library(ggplot2)
library(dplyr)
library(here)
library(readr)
library(tidyr)

# Set inputs ------------------------------------------------------------------

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

pop_baseline <- create_baseline_population(20000, max_age = 70, end_time = 2000, 
                                           group_weights = c("A" = 2, "B" = 1))

# Test group specific rates ---------------------------------------------------

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

# Test exogamy parameters -----------------------------------------------------

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

result |> get_detailed_unions() |>
  mutate(time_span = floor(start_time / 100) * 100) |>
  count(time_span, group_w, group_h) |>
  pivot_wider(id_cols = time_span, names_from = c(group_w, group_h), values_from = n) |>
  mutate(log_odds = log(A_B) + log(B_A) - log(A_A) - log(B_B))

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A:B" = -4, "B:A" = -2)
    )
  ) |>
  simulate(years = 100)

result |> get_detailed_unions() |>
  mutate(time_span = floor(start_time / 100) * 100) |>
  count(time_span, group_w, group_h) |>
  pivot_wider(id_cols = time_span, names_from = c(group_w, group_h), values_from = n) |>
  mutate(log_odds = log(A_B) + log(B_A) - log(A_A) - log(B_B))

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5, prevent_exogamy = TRUE,
    )
  ) |>
  simulate(years = 100)

result |> get_detailed_unions() |>
  mutate(time_span = floor(start_time / 100) * 100) |>
  count(time_span, group_w, group_h) |>
  pivot_wider(id_cols = time_span, names_from = c(group_w, group_h), values_from = n)

result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5, prevent_exogamy = c("A:B"),
    )
  ) |>
  simulate(years = 100)

result |> get_detailed_unions() |>
  mutate(time_span = floor(start_time / 100) * 100) |>
  count(time_span, group_w, group_h) |>
  pivot_wider(id_cols = time_span, names_from = c(group_w, group_h), values_from = n)

# Test inheritance rules

# defaults to mother
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(fertility_rates_overall_usa) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# fathers should be similar
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("father")
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# randome parent
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("random_parent")
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# same-sex
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("same_sex")
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# opposite-sex
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("opposite_sex")
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# ancestry
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("ancestry")
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# one drop rule with no passing
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("hierarchy", order = c("B", "A"))
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# one drop rule with 5% passing
result <- hammel_model() |>
  set_population(pop_baseline) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = create_group_inheritance("hierarchy", order = c("B", "A"), passing = 0.05)
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -4)
    )
  ) |>
  simulate(years = 100)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

# custom function
inheritance <- create_group_inheritance(
  "custom",
  groups = c("A", "B", "AB"),
  fun = function(child_ancestry, mother_group, father_group, ...) {
    # Same identity or missing father: follow mother
    if (is.na(father_group) || mother_group == father_group) {
      return(mother_group)
    }

    # Entirely unknown ancestry: follow mother
    known <- 1 - child_ancestry[["unknown"]]
    if (known <= 0) return(mother_group)
    
    # if we are here then assign as a function of group A 
    # ancestry
    prop_A <- child_ancestry[["group_A"]] / known

    logits <- c(
      A  = -5.138 + 11.88 * prop_A - 15.3 * prop_A^2 + 9.553 * prop_A^3,
      B  = 0.949 - 9.39 * prop_A + 11.935 * prop_A^2 - 8.483 * prop_A^3,
      AB = 0
    )
    weights <- exp(logits - max(logits))

    sample(names(weights), size = 1, prob = weights)
  }
)

pop_even <- create_baseline_population(100000, 70, 2000, group_weights = c("A" = 1, "B" = 1))

result <- hammel_model() |>
  set_population(pop_even) |>
  set_fertility(
    fertility_rates_overall_usa,
    inheritance = inheritance 
  ) |>
  set_mortality(mortality_rates_usa2000) |>
  set_union_formation(
    "on_demand", 
    matching = create_age_gap_matching_model(
      0, 5,
      group_exogamy = c("A::B" = -2)
    )
  ) |>
  simulate(years = 200)

calculate_population_size(result, c(2000, 2050, 2100), by_group = TRUE) |>
  group_by(time) |>
  mutate(prop = population / sum(population)) |>
  ungroup()

result$persons |>
  filter(is.na(death_time)) |>
  left_join(result$ancestry) |>
  group_by(group) |>
  summarize(prop_A = mean(group_A), prop_B = mean(group_B), .groups = "drop")

result$persons |>
  left_join(result$ancestry) |>
  left_join(
    result$persons |> 
      select(person_id, group) |> 
      rename(father_id = person_id, group_father = group)
  ) |>
  left_join(
    result$persons |> 
      select(person_id, group) |> 
      rename(mother_id = person_id, group_mother = group)
  ) |>
  filter(!is.na(group_father) & !is.na(group_mother) &
    group_father != group_mother) |>
  mutate(prop_a_group = cut(group_A, seq(0, 1, 0.1))) |>
  group_by(prop_a_group, group) |>
  summarize(n = n()) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  ggplot(aes(x = prop_a_group, y = prop, group = group, color = group))+
  geom_line()+
  geom_point(aes(size = n))+
  theme_bw()
