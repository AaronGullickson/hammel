# Run from the repository root. Offline preparation; no interpolation.
# Numerical transcriptions and limitations: usa_1980_2019_marital_rates.README.md.
# NCHS: women/1980 column, Table 4 p.8; wife/1980 column, Table 3 p.19.
# NCFMR: printed age-specific labels in the indicated female series.
root <- "validation/data"
rows <- function(measure, year, starts, ends, values, source, denominator) {
  stopifnot(length(starts) == length(ends), length(starts) == length(values),
    all(is.finite(starts)), all(ends > starts), all(is.finite(values)),
    all(values >= 0), all(head(ends, -1) == tail(starts, -1)), is.infinite(tail(ends, 1)))
  labels <- ifelse(is.finite(ends), paste0(starts, "-", ends - 1), paste0(starts, "+"))
  tibble::tibble(measure = measure, year = year, sex = "f", source_age_group = labels,
    age_start = starts, age_end = ends, denominator = denominator,
    rate_per_1000 = values, rate = values / 1000, source_id = source)
}
endpoints <- list(
  rows("first_marriage", 1980L, c(15,18,20,25,30,35,45,65), c(18,20,25,30,35,45,65,Inf),
    c(19.8,87.3,119.8,101.6,56.3,27.1,7.0,0.9), "NCHS_MVSR_32_11_Table4", "never_married_women"),
  rows("remarriage", 1980L, c(15,25,30,35,45,65), c(25,30,35,45,65,Inf),
    c(231.0,178.6,129.0,78.5,19.7,2.3), "NCHS_MVSR_32_11_Table4", "divorced_or_widowed_women"),
  rows("divorce", 1980L, c(15,20,25,30,35,45,55,65), c(20,25,30,35,45,55,65,Inf),
    c(45.9,52.2,43.3,33.2,23.2,10.1,3.8,1.6), "NCHS_Series21_43_Table3", "married_women"),
  rows("first_marriage", 2017L, c(15,20,25,30,35,40,45,50), c(20,25,30,35,40,45,50,Inf),
    c(7.8,42.9,85.7,83.0,55.0,37.0,29.6,12.9), "NCFMR_FP19_05_Figure1", "never_married_women"),
  rows("remarriage", 2018L, c(20,25,35,45,55,65), c(25,35,45,55,65,Inf),
    c(120.5,92.5,60.6,35.1,13.2,2.8), "NCFMR_FP20_13_Figure3", "divorced_or_widowed_women"),
  rows("divorce", 2019L, c(15,25,35,45,55,65), c(25,35,45,55,65,Inf),
    c(26.3,21.3,21.4,17.5,11.3,5.8), "NCFMR_FP21_16_Figure2_Women", "married_women")
)
all_rates <- do.call(rbind, endpoints)
stopifnot(inherits(all_rates, "tbl_df"), nrow(all_rates) == 42L,
  !anyDuplicated(all_rates[c("measure", "year", "sex", "age_start")]),
  isTRUE(all.equal(all_rates$rate * 1000, all_rates$rate_per_1000)),
  all_rates$rate[all_rates$measure == "first_marriage" & all_rates$year == 1980 &
    all_rates$age_start == 20] == .1198,
  all_rates$rate[all_rates$measure == "divorce" & all_rates$year == 2019 &
    all_rates$age_start == 65] == .0058,
  min(endpoints[[5]]$age_start) == 20)
for (x in endpoints) {
  name <- sprintf("usa_%s_%s_rates_female.csv", x$year[1], x$measure[1])
  write.csv(x[c("age_start", "age_end", "rate")], file.path(root, name), row.names = FALSE)
}
write.csv(all_rates, file.path(root, "usa_1980_2019_marital_rates.csv"), row.names = FALSE)
message("Prepared 42 published female rate cells across six endpoint schedules; no interpolation.")
