library(here)
devtools::load_all(".")
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(gt)


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

pop_base <- create_baseline_population(50000, 50, 2000)

model <- hammel_model() |>
  set_population(pop_base) |>
  set_mortality(mortality_rates_usa2000) |>
  set_fertility(fertility_rates_by_status_usa) |>
  set_union_formation("rate", rates = mar_rates_female_usa) |>
  set_union_dissolution(dissolution_rates)

result <- simulate(model, 100, 42)

table(result$unions$end_reason, exclude = NULL)

result |>
  calculate_union_counts(
    time_breaks = seq(2025, 2100, 75),
    dur_breaks = c(0, 5, 10, 15, 20, Inf)
  )
