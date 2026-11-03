# Rebuild from the archived ONS Figure 2b download, without network access.
# Run from the package root: Rscript validation/data/england_wales_2015_2019/prepare_england_wales_2019_union_rates.R
source_path <- "validation/data/england_wales_2015_2019/england_wales_marriage_rates_ons_figure2b.csv"
source <- read.csv(source_path, skip = 6, check.names = FALSE)
rates <- as.numeric(unlist(source[source$Year == 2019, -1], use.names = FALSE))
stopifnot(length(rates) == 11L, all(is.finite(rates)), all(rates >= 0))
# Finite intervals [16,20), [20,25), ..., [60,65); open [65,Inf).
midpoints <- c(18, seq(22.5, 62.5, 5))
ages <- 16:64
single <- approx(midpoints, rates[1:10], xout = ages + 0.5, rule = 2)$y / 1000
output <- data.frame(age_start = c(ages, 65), age_end = c(ages + 1, Inf),
                     rate = c(single, rates[11] / 1000))
stopifnot(nrow(output) == 50L, all(is.finite(output$rate)),
          all(output$rate >= 0), all(output$age_end > output$age_start),
          all(head(output$age_end, -1) == tail(output$age_start, -1)),
          abs(output$rate[output$age_start == 27] - 0.0414) < 1e-12,
          abs(output$rate[output$age_start == 32] - 0.0532) < 1e-12)
write.csv(output, "validation/data/england_wales_2015_2019/england_wales_2019_female_union_rates.csv",
          row.names = FALSE)
