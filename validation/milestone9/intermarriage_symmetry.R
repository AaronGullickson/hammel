# Sex symmetry of race intermarriage with a second group dimension
# (milestone 9, stage 6). Race carries a matching pair effect; region has
# none and is either independent of race or strongly associated with it.
# Compares A-wife/B-husband with B-wife/A-husband unions formed after year
# 10, overall and among couples whose partners share a region. Synthetic
# rates; no data are downloaded. Takes a few minutes.
library(hammel)
library(dplyr)

mortality <- data.frame(sex = rep(c("m", "f"), each = 3), age_start = rep(c(0, 1, 50), 2),
                        age_end = rep(c(1, 50, Inf), 2), rate = rep(c(.006, .0015, .045), 2))

run <- function(setting, seed) {
  set.seed(seed)
  groups <- switch(setting,
    one = data.frame(race = c("A", "B"), weight = c(.85, .15)),
    independent = data.frame(race = c("A", "B", "A", "B"), region = c("N", "N", "S", "S"),
                             weight = c(.425, .075, .425, .075)),
    associated = data.frame(race = c("A", "B", "A", "B"), region = c("N", "N", "S", "S"),
                            weight = c(.75, .02, .10, .13)))
  model <- hammel_model() |>
    set_population(create_baseline_population(30000, groups = groups)) |>
    set_mortality(mortality) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = .065)) |>
    set_union_formation("on_demand", kinship_generations = 1,
      matching = create_age_gap_matching_model(2, 5, group_exogamy = list(race = c("A:B" = -1))))
  if (setting != "one") model <- set_groups(model, race = group(), region = group(inheritance = "mother"))
  unions <- get_detailed_unions(simulate(model, 60, seed = seed)) |> filter(start_time >= 10)
  shared <- if (setting != "one") mutate(filter(unions, region_w == region_h), region = region_w)
  bind_rows(mutate(unions, region = "all"), shared) |>
    group_by(region) |>
    summarise(AB = sum(race_w == "A" & race_h == "B"), BA = sum(race_w == "B" & race_h == "A"),
              unions = n(), .groups = "drop") |>
    mutate(setting = setting)
}

replicates <- c(one = 4, independent = 4, associated = 10)
results <- bind_rows(lapply(names(replicates), function(s)
  bind_rows(lapply(seq_len(replicates[[s]]), function(k) run(s, k)))))
results |>
  group_by(setting, region) |>
  summarise(AB = sum(AB), BA = sum(BA), unions = sum(unions), .groups = "drop") |>
  mutate(ratio = AB / BA, se_log = sqrt(1 / AB + 1 / BA),
         low = ratio * exp(-1.96 * se_log), high = ratio * exp(1.96 * se_log)) |>
  print(n = Inf)
