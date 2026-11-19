# Rebuild the bundled U.S. mortality, fertility, marriage and migration rates.
# Run from the package root with Poppler installed.

library(here)
library(readr)
library(readxl)
library(tibble)
library(usethis)

root <- here("data-raw", "usa_2000")

# MORTALITY ------

sources <- c(
  m = "usa_life_table_1999_2001_males.xls",
  f = "usa_life_table_1999_2001_females.xls"
)
titles <- c(
  m = "Table 2. Life table for males: United States, 1999-2001",
  f = "Table 3. Life table for females: United States, 1999-2001"
)
starts <- c(0, 1, seq(5, 100, 5))
ends <- c(starts[-1], Inf)

abridge <- function(sex) {
  x <- read_excel(
    here(root, sources[[sex]]),
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  stopifnot(
    identical(x[[1]][1], titles[[sex]]),
    identical(unname(unlist(x[3, 2:6])), c("nqx", "lx", "ndx", "nLx", "Tx"))
  )
  # Single-year rows follow the "Years" heading; earlier rows are in days.
  rows <- seq(which(x[[1]] == "Years") + 1L, length.out = 110L)
  ages <- strsplit(x[[1]][rows], "-", fixed = TRUE)
  age <- as.numeric(vapply(ages, `[`, character(1), 1L))
  stopifnot(identical(age, as.numeric(0:109)), identical(
    as.numeric(vapply(ages, `[`, character(1), 2L)), as.numeric(1:110)
  ))
  lx <- as.numeric(x[[3]][rows])
  dx <- as.numeric(x[[4]][rows])
  big_l <- as.numeric(x[[5]][rows])
  big_t <- as.numeric(x[[6]][rows])
  stopifnot(
    all(is.finite(c(lx, dx, big_l, big_t))),
    lx[1] == 1e5,
    isTRUE(all.equal(head(lx, -1) - dx[-110], lx[-1], tolerance = 1e-6))
  )
  # Closed bands: deaths over person-years lived in the band. Open band:
  # survivors at 100 over their remaining person-years (l100 / T100).
  rate <- vapply(seq_along(starts), function(i) {
    if (is.infinite(ends[i])) {
      return(lx[age == starts[i]] / big_t[age == starts[i]])
    }
    band <- age >= starts[i] & age < ends[i]
    sum(dx[band]) / sum(big_l[band])
  }, numeric(1))
  data.frame(sex = sex, age_start = starts, age_end = ends, rate = rate)
}

mortality <- rbind(abridge("f"), abridge("m"))
stopifnot(
  nrow(mortality) == 44L,
  all(is.finite(mortality$rate)),
  all(mortality$rate > 0)
)

# FERTILITY ------

pdf <- here(root, "usa_fertility_revised_nchs_2003.pdf")
txt <- tempfile(fileext = ".txt")
status <- system2("pdftotext", c("-layout", shQuote(pdf), shQuote(txt)))
stopifnot(status == 0L)
lines <- readLines(txt, warn = FALSE)
unlink(txt)
extract <- function(table, marital_status) {
  start <- grep(paste0("^Table ", table, "\\. Birth rates"), lines)[1]
  section <- lines[start:(start + 30L)]
  row <- grep("^2000 ", section, value = TRUE)[1]
  # The year is followed by dotted leaders, then nine decimal rate columns.
  numbers <- regmatches(row, gregexpr("[0-9]+\\.[0-9]+", row))[[1]]
  values <- as.numeric(numbers)
  stopifnot(length(values) == 9L, all(is.finite(values)))
  data.frame(
    marital_status = marital_status,
    age_start = seq(15, 40, 5),
    age_end = seq(20, 45, 5),
    rate_per_1000 = values[c(2, 5:9)],
    rate = values[c(2, 5:9)] / 1000
  )
}
fertility <- rbind(extract(8, "married"), extract(7, "unmarried"))
stopifnot(nrow(fertility) == 12L,
  isTRUE(all.equal(fertility$rate_per_1000,
    c(287.8, 216.1, 161.6, 115.8, 48.1, 9.6, 39, 72.2, 58.5, 39.3, 19.7, 5))))
union_fertility <- data.frame(
  union_type = ifelse(
    fertility$marital_status == "married", "partnered", "unpartnered"
  ),
  fertility[c("age_start", "age_end", "rate")]
)

# Table 1, first (all-races) panel. Its final retained rate is strictly
# ages 40-44, unlike the 40+ numerator used in Tables 7 and 8.
start <- grep("^Table 1\\. Number of births", lines)[1]
section <- lines[start:(start + 25L)]
stopifnot(any(grepl("All races", section, fixed = TRUE)))
row <- grep("^2000 ", section, value = TRUE)[1]
values <- as.numeric(regmatches(row, gregexpr("[0-9]+\\.[0-9]+", row))[[1]])
# Decimal columns: crude, general, total fertility, 10-14, 15-19,
# 15-17, 18-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49.
stopifnot(length(values) == 13L, all(is.finite(values)),
  isTRUE(all.equal(values[c(5, 8:12)], c(47.7, 109.7, 113.5, 91.2, 39.7, 8.0))))
overall <- data.frame(age_start = seq(15, 40, 5), age_end = seq(20, 45, 5),
  rate = values[c(5, 8:12)] / 1000)
age_columns <- c("age_start", "age_end")
partnered <- union_fertility$union_type == "partnered"
stopifnot(identical(
  overall[age_columns],
  union_fertility[partnered, age_columns]
))

# MARRIAGE ------

x <- read_excel(
  here(root, "usa_marriage_transitions_2001_census.xls"),
  col_names = FALSE,
  col_types = "text",
  .name_repair = "minimal"
)
stopifnot(
  grepl("Unmarried at beginning", x[[2]][8], fixed = TRUE),
  grepl("year prior", x[[1]][39], fixed = TRUE)
)
# Use female annual transitions per 1,000 directly as rates per person-year.
rows <- 25:36
marriage <- data.frame(
  age_start = seq(15, 70, 5),
  age_end = c(seq(20, 70, 5), Inf),
  rate = as.numeric(x[[3]][rows]) / 1000
)
stopifnot(
  all(is.finite(marriage$rate)),
  all(marriage$rate >= 0 & marriage$rate < 1),
  all(head(marriage$age_end, -1) == tail(marriage$age_start, -1))
)

# MIGRATION ------

# Ballpark age-specific migration rates between US census regions, around
# 2000. Run from the package root. See README.md.
#
# Inputs, transcribed from U.S. Census Bureau, "Geographical Mobility: March
# 1999 to March 2000" (P20-538, Schachter 2001):
# - Figure 3: region-to-region migrants, thousands. The figure's values were
#   assigned to flows using Table D's regional in- and out-migrant totals,
#   which they reproduce exactly (checked below).
# - Table B: percent moving to a different state, by age (proxy for the age
#   pattern of moves between regions).
# Region populations are Census 2000 counts (thousands).

regions <- c("NE", "MW", "S", "W")
flows <- matrix(
  c(
    NA, 83, 405, 127,
    72, NA, 364, 204,
    183, 416, NA, 432,
    108, 223, 489, NA
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(regions, regions)
)
# Table D out-migrant and in-migrant totals.
stopifnot(
  rowSums(flows, na.rm = TRUE) == c(615, 640, 1031, 820),
  colSums(flows, na.rm = TRUE) == c(363, 722, 1258, 763)
)
population <- c(NE = 53594, MW = 64393, S = 100237, W = 63198)

# Table B, percent moving to a different state, by age; total 3.1. The first
# band (ages 1-4 in the source) is applied from age 0.
interstate <- data.frame(
  age_start = c(0, 5, 10, 20, 25, 30, 35, 45, 55, 65, 85),
  age_end = c(5, 10, 20, 25, 30, 35, 45, 55, 65, 85, Inf),
  percent = c(4.0, 3.3, 3.1, 6.2, 6.0, 4.0, 3.0, 1.9, 1.9, 1.1, 1.0)
)
interstate$multiplier <- interstate$percent / 3.1

pairs <- expand.grid(
  origin = regions,
  destination = regions,
  stringsAsFactors = FALSE
)
pairs <- pairs[pairs$origin != pairs$destination, ]
pair_flows <- flows[cbind(pairs$origin, pairs$destination)]
pairs$all_ages <- pair_flows / population[pairs$origin]
migration <- merge(
  pairs, interstate[c("age_start", "age_end", "multiplier")], by = NULL
)
migration$rate <- signif(migration$all_ages * migration$multiplier, 4)
rows <- order(migration$origin, migration$destination, migration$age_start)
columns <- c("origin", "destination", "age_start", "age_end", "rate")
migration <- migration[rows, columns]

# BUNDLED DATA ------

# Preserve the precision and column types of the original CSV-based build.
bundle_table <- function(table) {
  csv <- capture.output(write.csv(table, row.names = FALSE))
  columns <- read_csv(
    I(paste(csv, collapse = "\n")),
    col_types = cols(.default = col_character())
  )
  as_tibble(lapply(columns, type.convert, as.is = TRUE))
}
usa_2000_mortality <- bundle_table(mortality)
usa_2000_fertility <- bundle_table(overall)
usa_2000_fertility_union <- bundle_table(union_fertility)
usa_2000_marriage <- bundle_table(marriage)
usa_2000_migration <- bundle_table(migration)

stopifnot(
  nrow(usa_2000_mortality) == 44L,
  nrow(usa_2000_fertility) == 6L,
  nrow(usa_2000_fertility_union) == 12L,
  nrow(usa_2000_marriage) == 12L,
  nrow(usa_2000_migration) == 132L
)

use_data(
  usa_2000_mortality,
  usa_2000_fertility,
  usa_2000_fertility_union,
  usa_2000_marriage,
  usa_2000_migration,
  overwrite = TRUE
)
