# load libraries --------------------------------------------------

devtools::load_all()
library(ggplot2)
library(dplyr)
library(here)
library(readr)
library(tidyr)

# Get rates --------------------------------------------------------

mortality_rates_usa <- get_mortality_hmd("USA", year = 2000)

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
