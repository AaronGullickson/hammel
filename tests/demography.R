library(hammel)

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
persons <- data.frame(person_id = 1:3, sex = c("f", "m", "f"),
  birth_time = c(-0.5, 0.25, 0.5), death_time = c(0.75, NA, NA),
  mother_id = c(NA_integer_, NA_integer_, 1L), father_id = NA_integer_)
population <- hammel_population(persons, end_time = 1)
exposure <- calculate_lexis_counts(population, c(0, 1), c(0, 1, Inf))
stopifnot(isTRUE(all.equal(exposure$person_years[exposure$sex == "f"], c(1, 0.25))),
          isTRUE(all.equal(exposure$person_years[exposure$sex == "m"], c(0.75, 0))),
          sum(exposure$person_years) == 2)
mx <- calculate_lexis_counts(population, c(0, 1), c(0, 1, Inf))
stopifnot(mx$deaths[mx$sex == "f" & mx$age_start == 1] == 1,
          mx$asmr[mx$sex == "f" & mx$age_start == 1] == 4,
          is.na(mx$asmr[mx$sex == "m" & mx$age_start == 1]))
asfr <- calculate_lexis_counts(population, c(0, 1), c(0, 1, 2))
stopifnot(identical(asfr$births[asfr$sex == "f"], c(0, 1, 0)),
          all(is.na(asfr$births[asfr$sex == "m"])),
          all(is.na(asfr$asfr[asfr$sex == "m"])),
          asfr$asfr[asfr$sex == "f"][2] == 4,
          attr(asfr, "unassigned_births")$births == 1)

# Fractional periods and age widths split exposure without calendar-year loops.
split <- calculate_lexis_counts(population, c(0, 0.5, 1), c(0, 0.5, 1, Inf))
stopifnot(isTRUE(all.equal(sum(split$person_years), sum(exposure$person_years))))

# Count events on lower boundaries, excluding the final time/age boundary.
p <- persons[1, ]; p$birth_time <- 0; p$death_time <- 1
boundary <- calculate_lexis_counts(p, c(0, 1, 2), c(0, 1, Inf))
stopifnot(sum(boundary$deaths) == 1,
          boundary$deaths[boundary$sex == "f" & boundary$time_start == 1 &
                            boundary$age_start == 1] == 1)
excluded <- calculate_lexis_counts(p, c(0, 1), c(0, 1, Inf))
stopifnot(sum(excluded$deaths) == 0, sum(excluded$person_years) == 1)
expect_error(calculate_lexis_counts(population, c(0, 2), c(0, Inf)))

# Independent histories with repeated IDs must be analyzed separately.
stacked <- rbind(persons, persons)
expect_error(calculate_lexis_counts(stacked, c(0, 1), c(0, Inf)))

# Unequal age widths, missing bands, and period grouping.
rates <- data.frame(age_start = c(15, 20, 30), age_end = c(20, 30, 50),
                    asfr = c(0.1, 0.2, 0.05))
stopifnot(calculate_tfr(rates)$tfr == 3.5)
bad <- rates; bad$asfr[2] <- NA_real_
stopifnot(calculate_tfr(bad)$tfr == 1.5)
expect_error(calculate_tfr(rates[c(1, 3), ]))
bad <- rates; bad$age_end[3] <- Inf
stopifnot(calculate_tfr(bad)$tfr == 2.5)
rates$time_start <- 0; rates$time_end <- 1
later <- rates; later$time_start <- 1; later$time_end <- 2; later$asfr <- later$asfr * 2
stopifnot(identical(calculate_tfr(rbind(rates, later))$tfr, c(3.5, 7)))

# Exact piecewise constant hazard life tables: constant hazard gives e=1/m.
rates <- data.frame(age_start = c(0, 1, 5), age_end = c(1, 5, Inf), mx = 0.1)
lt <- calculate_life_table(rates)
stopifnot(all(abs(lt$ex - 10) < 1e-12),
          abs(lt$lx[3] - exp(-0.5)) < 1e-12,
          abs(sum(lt$dx) - 1) < 1e-12,
          abs(sum(lt$Lx) - 10) < 1e-12)
scaled <- calculate_life_table(rates, radix = 100000)
stopifnot(isTRUE(all.equal(scaled$ex, lt$ex)),
          isTRUE(all.equal(scaled$lx, lt$lx * 100000)))
tiny <- rates; tiny$mx <- 1e-12
lt <- calculate_life_table(tiny)
stopifnot(all(abs(lt$ex/1e12 - 1) < 1e-12), abs(lt$ax[1] - 0.5) < 1e-10)
zero <- rates; zero$mx <- 0
lt <- calculate_life_table(zero)
stopifnot(all(is.infinite(lt$ex)), all(lt$qx == 0), all(lt$lx == 1),
          identical(lt$Lx[1:2], c(1, 4)))
missing <- rates; missing$mx[2] <- NA_real_
stopifnot(all(is.na(calculate_life_table(missing)$ex)))
closed <- rates; closed$age_end[3] <- 100
expect_error(calculate_life_table(closed))
expect_error(calculate_life_table(rates[-1, ]))
expect_error(calculate_life_table(rates, radix = 0))
female <- rates; female$sex <- "f"
male <- rates; male$sex <- "m"; male$mx <- 0.2
tables <- calculate_life_table(rbind(female, male))
stopifnot(all(abs(tables$ex[tables$sex == "f"] - 10) < 1e-12),
          all(abs(tables$ex[tables$sex == "m"] - 5) < 1e-12))
high <- data.frame(age_start = c(0, 1), age_end = c(1, Inf), mx = 40)
stopifnot(calculate_life_table(high)$lx[2] > 0,
          all(abs(calculate_life_table(high)$ex - 1/40) < 1e-12))

# An empty population gives zero exposure, zero counts, and missing rates.
empty <- hammel_population(persons[FALSE, ], end_time = 1)
empty_rates <- calculate_lexis_counts(empty, c(0, 1), c(0, 1, Inf))
stopifnot(all(empty_rates$person_years == 0), all(empty_rates$deaths == 0),
          all(is.na(empty_rates$asmr)))

# The shared summary feeds both calculations directly.
summary <- calculate_lexis_counts(population, c(0, 1), c(0, 1, 2, Inf))
stopifnot(calculate_tfr(summary)$tfr == 4,
          nrow(calculate_life_table(summary)) == nrow(summary))
stopifnot(identical(names(calculate_life_table(summary)),
  c("time_start", "time_end", "sex", "age_start", "age_end", "mx",
    "n", "ax", "qx", "lx", "lx_end", "dx", "Lx", "Tx", "ex")))
stopifnot(identical(calculate_tfr(summary),
                    calculate_tfr(summary[is.finite(summary$age_end), ])))

# A finite final breakpoint starts an automatically appended open age band.
stopifnot(identical(calculate_lexis_counts(population, c(0, 1), c(0, 1)),
                    calculate_lexis_counts(population, c(0, 1), c(0, 1, Inf))))
