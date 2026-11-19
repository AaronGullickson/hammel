library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

immig_counts <- read_csv(
  here("validation", "data", "usa_historical_migration",
    "us_immigration_c89_c119.csv"
  )
)

immig_counts <- immig_counts |>
  # remove missing poland cases with zeroes
  mutate(c96 = if_else(is.na(c96), 0, c96)) |>
  mutate(
    # c92 is Ireland
    irish = c92,
    # NW Europe + Germany: c91 GB, c93 Scandanavia, c94 Other NW Europe, c95 Germany
    nw_europe = c91 + c93 + c94 + c95,
    # southern Europe: Italy c100, other Southern c101
    south_europe = c100 + c101,
    # eastern Europe: c96 Poland, c97 Other Central, c 98 USSR and Baltic, c99 other Eastern
    east_europe = c96 + c97 + c98 + c99
  ) |>
  select(year, irish, ends_with("_europe"))

immig_counts |>
  pivot_longer(cols = c("irish", ends_with("_europe"))) |>
  ggplot(aes(x = year, y = value, color = name))+
  geom_line()

 
pop_counts <- read_csv(
  here("validation", "data", "usa_historical_migration",
    "us_resident_population_a7.csv"
  )
)

immig_rates <- immig_counts |>
  pivot_longer(cols = c("irish", ends_with("_europe")),
               names_to = "stream", values_to = "count") |>
  left_join(pop_counts, by = "year") |>
  mutate(
    rate = count / resident_population,
    time_start = year,
    time_end = year + 1
  ) |>
  select(stream, time_start, time_end, rate)

ggplot(immig_rates, aes(x = time_start + 0.5, y = rate, color = stream))+
  geom_line()

pal_stream <- c(
  "Irish" = "green",
  "NW Europe" = "dodgerblue",
  "Southern Europe" = "yellow",
  "Eastern Europe" = "firebrick"
)

immig_rates |>
  mutate(
    stream = factor(
      stream,
      levels = c("irish", "nw_europe", "south_europe", "east_europe"),
      labels = c("Irish", "NW Europe", "Southern Europe", "Eastern Europe")
    )
) |>
  ggplot(aes(x = time_start + 0.5, y = rate, fill = stream))+
  geom_area()+
  scale_fill_manual(values = pal_stream)

saveRDS(
  immig_rates,
  file = here("validation", "data", "usa_historical_migration", "immig_rates.rds")
)
