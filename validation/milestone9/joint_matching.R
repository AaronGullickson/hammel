# Joint matching across two dimensions (milestone 9, step 2).
#
# 1. Sex symmetry within combinations. Race {A 85%, B 15%} and region {N, S}
#    both carry symmetric pair effects (A:B = -1, N:S = -1). For each pair of
#    different combinations (such as A|N and B|S), unions with the wife in one
#    and the husband in the other should be about as common as the reverse.
#    Region is either independent of race or associated with it.
# 2. Calibration cost. Run time as the number of matching combinations grows,
#    compared with one dimension having the same number of categories.
#
# Synthetic rates; no data are downloaded. Takes several minutes.
library(hammel)
library(dplyr)

mortality <- data.frame(sex = rep(c("m", "f"), each = 3), age_start = rep(c(0, 1, 50), 2),
                        age_end = rep(c(1, 50, Inf), 2), rate = rep(c(.006, .0015, .045), 2))

base_model <- function(population, matching) {
  hammel_model() |>
    set_population(population) |>
    set_mortality(mortality) |>
    set_fertility(data.frame(age_start = 15, age_end = 45, rate = .065)) |>
    set_union_formation("on_demand", kinship_generations = 1, matching = matching)
}

# --- 1. Symmetry -------------------------------------------------------------

symmetry_run <- function(setting, seed) {
  set.seed(seed)
  groups <- data.frame(race = c("A", "B", "A", "B"), region = c("N", "N", "S", "S"),
    weight = switch(setting, independent = c(.425, .075, .425, .075),
                             associated = c(.75, .02, .10, .13)))
  population <- create_baseline_population(30000, groups = groups)
  matching <- create_age_gap_matching_model(2, 5,
    group_exogamy = list(race = c("A:B" = -1), region = c("N:S" = -1)))
  model <- set_groups(base_model(population, matching), race = group(), region = group(inheritance = "mother"))
  get_detailed_unions(simulate(model, 60, seed = seed)) |>
    filter(start_time >= 10) |>
    transmute(w = paste(race_w, region_w, sep = "|"), h = paste(race_h, region_h, sep = "|"))
}

symmetry <- bind_rows(lapply(c("independent", "associated"), function(setting) {
  unions <- bind_rows(lapply(1:6, function(seed) symmetry_run(setting, seed)))
  counts <- count(unions, w, h)
  pairs <- counts |>
    filter(w < h) |>
    left_join(rename(counts, h = w, w = h, n_reverse = n), by = c("w", "h")) |>
    mutate(n_reverse = coalesce(n_reverse, 0L), setting = setting)
  pairs
})) |>
  mutate(ratio = n / n_reverse, se_log = sqrt(1 / n + 1 / n_reverse),
         low = ratio * exp(-1.96 * se_log), high = ratio * exp(1.96 * se_log))
print(symmetry, n = Inf)

# --- 1b. Symmetry with separate markets (D138) --------------------------------
# Race effect with region prohibited: each region is its own market. Compare
# A-wife/B-husband with B-wife/A-husband within each region and pooled.

markets_run <- function(seed) {
  set.seed(seed)
  groups <- data.frame(race = c("A", "B", "A", "B"), region = c("N", "N", "S", "S"),
                       weight = c(.75, .02, .10, .13))
  population <- create_baseline_population(30000, groups = groups)
  matching <- create_age_gap_matching_model(2, 5,
    group_exogamy = list(race = c("A:B" = -1), region = NA))
  model <- set_groups(base_model(population, matching), race = group(), region = group(inheritance = "mother"))
  unions <- get_detailed_unions(simulate(model, 60, seed = seed)) |> filter(start_time >= 10)
  bind_rows(mutate(unions, region = "all"), mutate(unions, region = region_w)) |>
    group_by(region) |>
    summarise(AB = sum(race_w == "A" & race_h == "B"), BA = sum(race_w == "B" & race_h == "A"),
              cross_region = sum(region_w != region_h), .groups = "drop")
}
markets <- bind_rows(lapply(1:6, markets_run)) |>
  group_by(region) |>
  summarise(AB = sum(AB), BA = sum(BA), cross_region = sum(cross_region), .groups = "drop") |>
  mutate(ratio = AB / BA, se_log = sqrt(1 / AB + 1 / BA),
         low = ratio * exp(-1.96 * se_log), high = ratio * exp(1.96 * se_log))
print(markets)

# --- 2. Calibration cost -----------------------------------------------------

timing_run <- function(race_k, region_k, joint) {
  set.seed(1)
  races <- LETTERS[seq_len(race_k)]
  regions <- paste0("R", seq_len(region_k))
  if (joint) {
    groups <- expand.grid(race = races, region = regions, stringsAsFactors = FALSE)
    groups$weight <- 1
    population <- create_baseline_population(20000, groups = groups)
    pairs <- function(labels) setNames(rep(-1, length(labels) - 1), paste0(labels[1], ":", labels[-1]))
    exogamy <- list(race = pairs(races))
    if (region_k > 1) exogamy$region <- pairs(regions)
    model <- set_groups(base_model(population, create_age_gap_matching_model(2, 5, group_exogamy = exogamy)),
                        race = group(), region = group())
  } else {
    k <- race_k * region_k
    labels <- paste0("G", seq_len(k))
    population <- create_baseline_population(20000, groups = list(group = setNames(rep(1, k), labels)))
    model <- base_model(population, create_age_gap_matching_model(2, 5,
      group_exogamy = setNames(rep(-1, k - 1), paste0("G1:", labels[-1]))))
  }
  elapsed <- system.time(simulate(model, 20, seed = 2))[["elapsed"]]
  data.frame(combinations = race_k * region_k, joint = joint, seconds = elapsed)
}

cost <- bind_rows(lapply(list(c(2, 1), c(2, 2), c(2, 4), c(4, 4), c(4, 8)), function(k)
  bind_rows(timing_run(k[1], k[2], TRUE), timing_run(k[1], k[2], FALSE))))
print(cost)

# Prohibitions: 4 races x 4 regions (16 combinations).
prohibition_run <- function(exogamy, label) {
  set.seed(1)
  groups <- expand.grid(race = LETTERS[1:4], region = paste0("R", 1:4), stringsAsFactors = FALSE)
  groups$weight <- 1
  model <- set_groups(base_model(create_baseline_population(20000, groups = groups),
    create_age_gap_matching_model(2, 5, group_exogamy = exogamy)), race = group(), region = group())
  data.frame(setting = label, seconds = system.time(suppressWarnings(simulate(model, 20, seed = 2)))[["elapsed"]])
}
race <- c("A:B" = -1, "A:C" = -1, "A:D" = -1)
print(bind_rows(
  prohibition_run(list(race = race, region = c("R1:R2" = -1, "R1:R3" = -1, "R1:R4" = -1)), "effects in both"),
  prohibition_run(list(race = race, region = NA), "race effects, region prohibited"),
  prohibition_run(list(race = NA, region = NA), "both prohibited")))
