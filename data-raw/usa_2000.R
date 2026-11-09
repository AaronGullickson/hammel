# Build the bundled United States example rates, around 2000, from the
# prepared CSVs in validation/data/usa_2000_2001. Run from the package root
# after the preparation scripts there. Requires usethis. Sources, methods and
# caveats are in that folder's README files and in R/data.R.

root <- "validation/data/usa_2000_2001"
read_rates <- function(file) {
  tibble::as_tibble(read.csv(file.path(root, file)))
}

usa_2000_mortality <- read_rates("usa_1999_2001_mortality.csv")
usa_2000_fertility <- read_rates("usa_2000_fertility_overall.csv")
usa_2000_fertility_union <- read_rates("usa_2000_fertility_by_union_status.csv")
# Current fertility input names this column `union_type`.
names(usa_2000_fertility_union)[
  names(usa_2000_fertility_union) == "union_status"
] <- "union_type"
usa_2000_marriage <- read_rates("usa_2001_marriage_rates_female.csv")
usa_2000_migration <- read_rates("usa_2000_regional_migration.csv")

stopifnot(
  nrow(usa_2000_mortality) == 44L,
  nrow(usa_2000_fertility) == 6L,
  nrow(usa_2000_fertility_union) == 12L,
  nrow(usa_2000_marriage) == 12L,
  nrow(usa_2000_migration) == 132L
)

usethis::use_data(
  usa_2000_mortality,
  usa_2000_fertility,
  usa_2000_fertility_union,
  usa_2000_marriage,
  usa_2000_migration,
  overwrite = TRUE
)
