# Rebuild the bundled England and Wales observations from original sources.
# Run from the package root with readxl and tibble installed.

root <- "data-raw/england_wales_1992_onward"
read_sheet <- function(file, sheet) {
  data <- readxl::read_excel(
    file.path(root, file),
    sheet = sheet,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  as.data.frame(data, check.names = FALSE)
}
number <- function(x) suppressWarnings(as.numeric(x))
age_bounds <- function(label) {
  nums <- number(regmatches(label, gregexpr("[0-9]+", label))[[1]])
  if (length(nums) == 2L) {
    return(c(nums[1], nums[2] + 1))
  }
  if (length(nums) == 1L && grepl("over|\\+", label)) {
    return(c(nums, Inf))
  }
  stop("Unrecognised age label: ", label)
}

# Published Table B1 uses finite population denominator age bands.
fertility_1992 <- matrix(
  c(
    286.7, 192, 22,
    196.3, 127.9, 31,
    169.5, 107.9, 32.7,
    103.3, 101.2, 26.3,
    35.5, 62.3, 13.3,
    5.7, 16.6, 2.9
  ),
  ncol = 3,
  byrow = TRUE
)
fertility <- tibble::tibble(
  union_type = rep(c("marriage", "cohab", "unpartnered"), each = 6),
  age_start = rep(c(16L, 20L, 25L, 30L, 35L, 40L), 3),
  age_end = rep(c(20L, 25L, 30L, 35L, 40L, 45L), 3),
  rate = as.vector(fertility_1992) / 1000
)

# Annual female population denominators: old demographic estimates through
# 2001; survey-based revised estimates 2002 onward. Keep the method break.
pop <- list()
for (year in 1992:2001) {
  d <- read_sheet("ons_marital_population_1971_2001.xls", as.character(year))
  rows <- which(grepl("^[0-9]+$|^[0-9]+ and over$", d[[1]]))
  for (i in rows) {
    a <- number(sub(" and over", "", d[i, 1]))
    for (status in c("Single", "Widowed", "Divorced")) {
      col <- c(Single = 9, Widowed = 11, Divorced = 12)[status]
      pop[[length(pop) + 1L]] <- data.frame(
        year = year,
        marital_status = status,
        age_start = a,
        age_end = if (grepl("over", d[i, 1])) Inf else a + 1,
        population = number(d[i, col]) * 1000
      )
    }
  }
}
d <- read_sheet("ons_marital_population_2002_2019.xlsx", "1 E&W Marital Status")
start <- which(d[[1]] == "FEMALES")
hdr <- start + 1L
year_cols <- which(number(unlist(d[hdr, ])) %in% 2002:2019)
for (i in seq.int(hdr + 1L, nrow(d))) {
  status <- d[i, 1]
  label <- d[i, 2]
  if (is.na(status) || !status %in% c("Single", "Widowed", "Divorced") ||
        is.na(label) || !grepl("^[0-9]+[-+]|^[0-9]+$", label)) {
    next
  }
  bounds <- age_bounds(label)
  for (j in year_cols) {
    pop[[length(pop) + 1L]] <- data.frame(
      year = number(d[hdr, j]),
      marital_status = status,
      age_start = bounds[1],
      age_end = bounds[2],
      population = number(d[i, j])
    )
  }
}
pop <- do.call(rbind, pop)
stopifnot(setequal(unique(pop$year), 1992:2019), !anyDuplicated(pop[1:4]))

# Counts and denominators use exactly the same age bands. No interpolation,
# population smoothing, imputation of suppressed cells, or split of open ages.
d <- read_sheet("ons_marriages_2019.xlsx", "3a")
header <- which(d[[1]] == "Year")
labels <- unlist(d[header, ])
cols <- which(
  grepl("^Women", labels) &
    grepl("Under 20|[0-9]+ to [0-9]+|80 and over", labels)
)
marriage_statuses <- list(
  first_marriage = "Single",
  remarriage = c("Divorced", "Widowed"),
  all_marriages = c("Single", "Divorced", "Widowed")
)
out <- list()
for (year in 1992:2019) {
  for (j in cols) {
    label <- sub("Women[\r\n]+", "", d[header, j])
    bounds <- if (label == "Under 20") c(16, 20) else age_bounds(label)
    for (kind in names(marriage_statuses)) {
      statuses <- marriage_statuses[[kind]]
      rows <- which(d[[1]] == as.character(year) & d[[2]] %in% statuses)
      count <- sum(number(unlist(d[rows, j])))
      matches <- pop$year == year & pop$marital_status %in% statuses &
        pop$age_start >= bounds[1] & pop$age_end <= bounds[2]
      p <- pop[matches, , drop = FALSE]
      stopifnot(nrow(p) > 0L)
      denom <- sum(p$population) # NA stays NA; a suppressed cell is not zero.
      out[[length(out) + 1L]] <- data.frame(
        year = year,
        age_start = bounds[1],
        age_end = bounds[2],
        marriage_order = kind,
        marriages = count,
        rate = if (is.finite(denom) && denom > 0) count / denom else NA_real_
      )
    }
  }
}
out <- do.call(rbind, out)
stopifnot(nrow(out) == 28L * 12L * 3L, all(out$rate >= 0, na.rm = TRUE))
# Counts by band must reconcile with published status-specific age-known totals.
for (year in 1992:2019) {
  for (kind in names(marriage_statuses)) {
    statuses <- marriage_statuses[[kind]]
    rows <- which(d[[1]] == as.character(year) & d[[2]] %in% statuses)
    unstated <- number(d[rows, 40])
    # ONS note 10: age-not-stated cases were distributed over ages in 1961-93.
    if (year <= 1993) {
      unstated[] <- 0
    }
    expected <- sum(number(d[rows, 22])) - sum(unstated) -
      sum(number(d[rows, 23]))
    counts <- out$marriages[out$year == year & out$marriage_order == kind]
    stopifnot(sum(counts) == expected)
  }
}

# These remain observed marriage rates, not tuned union-search rates.
marriage <- tibble::as_tibble(out)
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

# ONS Table 11a measures shared addresses at marriage.
d <- read_sheet("ons_marriages_2019.xlsx", "11a")
rows <- which(number(d[[1]]) %in% 1994:2019)
premarital_cohabitation <- tibble::tibble(
  year = as.integer(number(d[rows, 1])),
  age_groups = d[rows, 2],
  prop_cohabit = number(d[rows, 6]) / 100
)

# Table 3 gives competing-risk outcomes for women under 45 at union start.
cohabitation_outcomes <- tibble::tibble(
  cohort_start = c(1990L, 1995L, 2000L),
  cohort_end = c(1995L, 2000L, 2005L),
  marriage_5 = c(.463, .399, .395),
  separation_5 = c(.232, .305, .354),
  continuing_5 = c(.305, .296, .252),
  marriage_10 = c(.561, .494, NA),
  separation_10 = c(.320, .381, NA),
  continuing_10 = c(.120, .124, NA)
)
stopifnot(all(abs(rowSums(cohabitation_outcomes[3:5]) - 1) < .002))

life <- read.csv(file.path(root, "marital_life_table_2005.csv"))
survivor_diff <- head(life$surviving_marriages, -1) -
  head(life$divorce_events, -1) - head(life$death_events, -1) -
  tail(life$surviving_marriages, -1)
stopifnot(
  nrow(life) == 52L,
  identical(life$dur_start, c(0:50, 55L)),
  all(life$dur_end > life$dur_start),
  all(life$prob_divorce >= 0 & life$prob_death >= 0),
  all(life$prob_divorce + life$prob_death < 1),
  all(abs(survivor_diff) <= 1)
)
# Columns g/h in the paper apply to the END of each interval. In particular
# source duration 4 gives incidence by anniversary 5, not anniversary 4.
# Integer life-table event counts retain more precision than the 0.1%-rounded
# published cumulative columns (which hide small nonzero late divorces).
initial_marriages <- life$surviving_marriages[1]
cumulative <- tibble::tibble(
  year = 2005L,
  anniversary = life$dur_end,
  cumulative_divorce = cumsum(life$divorce_events) / initial_marriages,
  cumulative_death = cumsum(life$death_events) / initial_marriages
)
divorce_diff <- cumulative$cumulative_divorce -
  life$published_cumulative_divorce_end
death_diff <- cumulative$cumulative_death - life$published_cumulative_death_end
stopifnot(all(abs(divorce_diff) <= .00051), all(abs(death_diff) <= .00051))

divorce_2010 <- file.path(root, "ons_divorce_death_cumulative_2010.csv")
d <- read.csv(divorce_2010, skip = 6)
cumulative_2010 <- tibble::tibble(
  year = 2010L,
  anniversary = as.integer(number(d[[1]])),
  cumulative_divorce = number(d[[2]]) / 100,
  cumulative_death = number(d[[3]]) / 100
)

england_wales_1992_onward <- list(
  fertility = fertility,
  marriage = marriage,
  premarital_cohabitation = premarital_cohabitation,
  cohabitation_outcomes = cohabitation_outcomes,
  marital_divorce_death_2005 = cumulative,
  marital_divorce_death_2010 = cumulative_2010
)
# Preserve the precision and column types of the original CSV-based build.
england_wales_1992_onward <- lapply(england_wales_1992_onward, function(table) {
  csv <- capture.output(write.csv(table, row.names = FALSE, na = ""))
  tibble::as_tibble(read.csv(text = csv, fileEncoding = "UTF-8"))
})
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
