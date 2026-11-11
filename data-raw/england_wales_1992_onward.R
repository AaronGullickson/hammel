# Bundle the local inputs used by the second demographic-transition article.
# Run offline from the package root. Original sources, attribution and methods
# remain in validation/data/england_wales_1992_onward and R/data.R.

root <- "validation/data/england_wales_1992_onward"
checksums <- c(
  fertility_by_union_type_1992.csv = "5a1b0d77f79b5fd6c5a46d6fd5ee78c6",
  female_marriage_rates_1992_2019.csv = "84ea067121a85422e87ac2f5c588960c",
  female_premarital_cohabitation_1994_2019.csv =
    "4406b5ababfa8b549fd8a75bc37994b8",
  great_britain_cohabitation_cohort_outcomes.csv =
    "4a635adb573ed0ca3f7f3accb58c0229",
  marital_divorce_death_cumulative_2005.csv =
    "f42f7504c293d15c005646a2d0a8caec",
  marital_divorce_death_cumulative_2010.csv =
    "cad9a54f8b8c5aa5958a1a31d63ce918"
)
stopifnot(identical(
  unname(tools::md5sum(file.path(root, names(checksums)))),
  unname(checksums)
))
read_source <- function(file) {
  input <- read.csv(
    file.path(root, file),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  tibble::as_tibble(input)
}

# Keep the vignette's finite denominator-age convention, not an open tail.
fertility <- read_source("fertility_by_union_type_1992.csv")
names(fertility)[names(fertility) == "denominator_age_start"] <- "age_start"
names(fertility)[names(fertility) == "denominator_age_end"] <- "age_end"
fertility <- fertility[c("union_type", "age_start", "age_end", "rate")]

# These remain observed marriage rates, not tuned union-search rates.
marriage <- read_source("female_marriage_rates_1992_2019.csv")
marriage <- marriage[marriage$marriage_order != "all_marriages", ]
stopifnot(all(marriage$marriage_order %in% c("first_marriage", "remarriage")))
marriage$time_start <- marriage$year
marriage$time_end <- marriage$year + 1L
marriage$previous_union <- ifelse(
  marriage$marriage_order == "first_marriage", "none", "dissolved"
)
marriage <- marriage[c(
  "time_start", "time_end", "age_start", "age_end", "previous_union", "rate"
)]

england_wales_1992_onward <- list(
  fertility = fertility,
  marriage = marriage,
  premarital_cohabitation = read_source(
    "female_premarital_cohabitation_1994_2019.csv"
  ),
  cohabitation_outcomes = read_source(
    "great_britain_cohabitation_cohort_outcomes.csv"
  ),
  marital_divorce_death_2005 = read_source(
    "marital_divorce_death_cumulative_2005.csv"
  ),
  marital_divorce_death_2010 = read_source(
    "marital_divorce_death_cumulative_2010.csv"
  )
)
stopifnot(
  identical(vapply(england_wales_1992_onward, nrow, integer(1)), c(
    fertility = 18L,
    marriage = 672L,
    premarital_cohabitation = 373L,
    cohabitation_outcomes = 3L,
    marital_divorce_death_2005 = 52L,
    marital_divorce_death_2010 = 60L
  )),
  sum(is.na(marriage$rate)) == 48L,
  all(is.finite(fertility$rate) & fertility$rate >= 0),
  all(marriage$rate >= 0, na.rm = TRUE)
)
save(
  england_wales_1992_onward,
  file = "data/england_wales_1992_onward.rda",
  compress = "xz",
  version = 2
)
