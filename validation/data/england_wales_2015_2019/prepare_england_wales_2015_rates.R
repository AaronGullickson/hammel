# Run from the package root. Requires readxl; no network access.
root <- "validation/data/england_wales_2015_2019"
path <- file.path(root, "england_wales_births_by_parents_2015_ons.xls")
x <- readxl::read_excel(path, sheet = "Table 6", col_names = FALSE,
                        col_types = "text", .name_repair = "minimal")
stopifnot(grepl("per 1,000", x[[11]][12], fixed = TRUE),
          grepl("single, widowed and divorced", x[[11]][56], fixed = TRUE))
labels <- unlist(x[9, 12:18], use.names = FALSE)
stopifnot(identical(labels, c("Under 20", "20 to 24", "25 to 29", "30 to 34",
                             "35 to 39", "40 to 44", "45 and over")))
sections <- list(total = 14:32, partnered = 36:54, unpartnered = 58:76)
long <- do.call(rbind, lapply(names(sections), function(status) {
  rows <- sections[[status]]
  rows <- rows[grepl("^[12][0-9]{3}( |$)", x[[1]][rows]) %in% TRUE]
  do.call(rbind, lapply(rows, function(i) {
    # Trailing 1 in year labels is a footnote, not part of the year.
    year <- as.integer(substr(x[[1]][i], 1, 4))
    raw <- as.numeric(unlist(x[i, 12:18], use.names = FALSE))
    data.frame(year = year, union_status = status,
      source_age_group = labels, age_start = seq(15, 45, 5),
      age_end = seq(20, 50, 5),
      births = as.numeric(unlist(x[i, 3:9], use.names = FALSE)),
      rate_per_1000 = round(raw, 1), rate = round(raw, 1) / 1000)
  }))
}))
rownames(long) <- NULL
long <- long[order(long$year, long$union_status, long$age_start), ]
stopifnot(nrow(long) == 17L * 3L * 7L, all(is.finite(long$rate)),
          all(long$rate >= 0), !anyDuplicated(long[c("year", "union_status", "age_start")]))
# Published married and unmarried birth counts sum to the overall counts.
counts <- split(long$births, long$union_status)
stopifnot(all(counts$partnered + counts$unpartnered == counts$total))
write.csv(long, file.path(root, "england_wales_fertility_ons_table6.csv"), row.names = FALSE)
conditional <- subset(long, year == 2015 & union_status != "total",
                       select = c(union_status, age_start, age_end, rate))
overall <- subset(long, year == 2015 & union_status == "total",
                   select = c(age_start, age_end, rate))
stopifnot(nrow(conditional) == 14L, nrow(overall) == 7L,
          conditional$rate[conditional$union_status == "partnered" &
                            conditional$age_start == 20] == 0.2027)
write.csv(conditional, file.path(root, "england_wales_2015_fertility_by_union_status.csv"),
          row.names = FALSE)
write.csv(overall, file.path(root, "england_wales_2015_overall_fertility.csv"), row.names = FALSE)
# The archived marriage chart includes 2015. Preserve its original age groups.
m <- read.csv(file.path(root, "england_wales_marriage_rates_ons_figure2b.csv"),
              skip = 6, check.names = FALSE)
values <- as.numeric(unlist(m[m$Year == 2015, -1], use.names = FALSE))
stopifnot(length(values) == 11L, all(is.finite(values)), all(values >= 0))
marriage <- data.frame(age_start = c(16, seq(20, 65, 5)),
                       age_end = c(seq(20, 65, 5), Inf), rate = values / 1000)
write.csv(marriage, file.path(root, "england_wales_2015_female_union_rates.csv"), row.names = FALSE)
