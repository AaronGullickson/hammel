# Create dataset of crude immigration rates by European region over time in 
# the US

# Load packages ---------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(readr)

# Load raw data ---------------------------------------------------------------

dir_source <- here("data-raw", "usa_historical_migration")

immig_counts <- read_csv(here(dir_source, "us_immigration_c89_c119.csv"))

pop_counts <- read_csv(here(dir_source, "us_resident_population_a7.csv"))

# Calculate rates -------------------------------------------------------------


usa_historical_migration <- immig_counts |>
  # remove missing poland cases with zeroes
  mutate(c96 = if_else(is.na(c96), 0, c96)) |>
  mutate(
    # c92 is Ireland
    irish = c92,
    # c95 is Germany
    german = c95,
    # NW Europe + Germany: c91 GB, c93 Scandanavia, c94 Other NW Europe,
    # c95 Germany
    nw_europe = c91 + c93 + c94,
    # southern Europe: Italy c100, other Southern c101
    south_europe = c100 + c101,
    # eastern Europe: c96 Poland, c97 Other Central, c 98 USSR and Baltic,
    # c99 other Eastern
    east_europe = c96 + c97 + c98 + c99,
    # all others in one group, c89 is total
    other = c89 - irish -german - nw_europe - south_europe - east_europe
  ) |>
  select(year, irish, german, ends_with("_europe"), other) |>
  pivot_longer(cols = -year,
               names_to = "stream", values_to = "count") |>
  left_join(pop_counts, by = "year") |>
  mutate(
    rate = count / resident_population,
    time_start = year,
    time_end = year + 1
  ) |>
  select(stream, time_start, time_end, rate)

# ggplot2::ggplot(usa_historical_migration, 
#                 ggplot2::aes(x = time_start + 0.5, y = rate, color = stream))+
#   ggplot2::geom_line()

usethis::use_data(usa_historical_migration, overwrite = TRUE)
