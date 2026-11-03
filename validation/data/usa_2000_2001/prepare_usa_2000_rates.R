# Run from the package root. Requires readxl and the pdftotext command.
# Rebuild offline from the archived official sources; no interpolation.
root <- "validation/data/usa_2000_2001"
pdf <- file.path(root, "usa_fertility_revised_nchs_2003.pdf")
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
  data.frame(year = 2000L, marital_status = marital_status,
    source_table = table, source_age_group = c("15-19", "20-24", "25-29",
      "30-34", "35-39", "40-44"), age_start = seq(15, 40, 5),
    age_end = seq(20, 45, 5), rate_per_1000 = values[c(2, 5:9)],
    rate = values[c(2, 5:9)] / 1000)
}
fertility <- rbind(extract(8, "married"), extract(7, "unmarried"))
stopifnot(nrow(fertility) == 12L,
  isTRUE(all.equal(fertility$rate_per_1000,
    c(287.8, 216.1, 161.6, 115.8, 48.1, 9.6, 39, 72.2, 58.5, 39.3, 19.7, 5))))
write.csv(fertility, file.path(root, "usa_2000_fertility_by_marital_status.csv"),
          row.names = FALSE)
model <- data.frame(union_status = ifelse(fertility$marital_status == "married",
  "partnered", "unpartnered"), fertility[c("age_start", "age_end", "rate")])
write.csv(model, file.path(root, "usa_2000_fertility_by_union_status.csv"),
          row.names = FALSE)

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
stopifnot(identical(overall[c("age_start", "age_end")],
  model[model$union_status == "partnered", c("age_start", "age_end")]))
write.csv(overall, file.path(root, "usa_2000_fertility_overall.csv"), row.names = FALSE)

x <- readxl::read_excel(file.path(root, "usa_marriage_transitions_2001_census.xls"),
  col_names = FALSE, col_types = "text", .name_repair = "minimal")
stopifnot(grepl("Unmarried at beginning", x[[2]][8], fixed = TRUE),
          grepl("year prior", x[[1]][39], fixed = TRUE))
marriage <- do.call(rbind, lapply(c("m", "f"), function(sex) {
  rows <- if (sex == "m") 11:22 else 25:36
  data.frame(year = 2001L, sex = sex,
    source_age_group = sub("^\\.", "", x[[1]][rows]),
    age_start = seq(15, 70, 5), age_end = c(seq(20, 70, 5), Inf),
    unmarried_at_start_thousands = as.numeric(x[[2]][rows]),
    transitions_per_1000 = as.numeric(x[[3]][rows]),
    margin_of_error_90_per_1000 = as.numeric(x[[4]][rows]),
    annual_transition_proportion = as.numeric(x[[3]][rows]) / 1000)
}))
stopifnot(nrow(marriage) == 24L,
  all(is.finite(marriage$transitions_per_1000)),
  all(marriage$annual_transition_proportion >= 0 &
      marriage$annual_transition_proportion < 1),
  !anyDuplicated(marriage[c("sex", "age_start")]))
write.csv(marriage, file.path(root, "usa_2001_marriage_transitions.csv"),
          row.names = FALSE)

# User-selected approximation: use annual transitions / 1000 directly as
# per-person-year search rates, without a probability-to-hazard transformation.
# Female rates configure the current engine; male rates are comparison data.
for (sex in c("f", "m")) {
  selected <- marriage[marriage$sex == sex, ]
  rates <- data.frame(age_start = selected$age_start, age_end = selected$age_end,
    rate = selected$transitions_per_1000 / 1000)
  stopifnot(nrow(rates) == 12L, all(is.finite(rates$rate)), all(rates$rate >= 0),
    all(head(rates$age_end, -1) == tail(rates$age_start, -1)),
    identical(rates$rate, selected$annual_transition_proportion))
  label <- if (sex == "f") "female" else "male"
  write.csv(rates, file.path(root, paste0("usa_2001_marriage_rates_", label, ".csv")),
    row.names = FALSE)
}
