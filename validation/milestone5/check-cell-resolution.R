# Offline diagnostic: how cell resolution affects within-band recovery under
# cell-accounted fertility (D101). The generational "young" scenario of
# check-cell-spacing.R (300 founders aged 20, constant parity inputs .2/.3/.4
# on ages 15-40, no mortality before 60, seeds 9701-9710), observed in years
# [80,100). The same rates are supplied as one 25-year band or as 1- or
# 0.25-year bands, with pooled [2,Inf) or separate actual-parity bands. Cells
# are split to at most one year and by actual parity (revised decision 5), so
# the 25-year, 1-year and per-parity inputs should give identical results;
# 0.25-year inputs give 0.25-year cells. Run from the repository root:
#   Rscript validation/milestone5/check-cell-resolution.R
# Prints means and MCSEs; writes cell-resolution-summary.csv.
pkgload::load_all(quiet = TRUE)
mort <- data.frame(sex = rep(c("m", "f"), each = 2), age_start = c(0, 60), age_end = c(60, Inf), rate = c(0, 1.26))
founders <- hammel_population(data.frame(person_id = 1:300, sex = "f", birth_time = -20,
  death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_), end_time = 0)
make <- function(width, split_parity) {
  ages <- seq(15, 40 - width, width)
  pb <- if (split_parity) list(s = 0:4, e = c(1:4, Inf), r = c(.2, .3, .4, .4, .4)) else
    list(s = 0:2, e = c(1, 2, Inf), r = c(.2, .3, .4))
  do.call(rbind, lapply(seq_along(pb$s), function(k) data.frame(age_start = ages,
    age_end = ages + width, parity_start = pb$s[k], parity_end = pb$e[k], rate = pb$r[k])))
}
run <- function(width, split, seed) {
  m <- hammel_model() |> set_population(founders) |> set_mortality(mort) |>
    set_fertility(make(width, split)) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  pop <- simulate(m, 100, seed)
  f <- function(ages, br) { r <- calculate_fertility_rates(pop, time_breaks = c(80, 100), age_breaks = ages, parity_breaks = br)
    r[r$age_start == ages[1] & r$time_start == 80, "rate_fertility"][[1]] }
  c(f(c(18, 19), c(0, 1, 2)), f(c(18, 19), c(0, 1, 2, 3, 4)), f(c(20, 35), c(0, 1, 2)),
    sum(pop$simulation$unrealized_births$unrealized_births) / pop$simulation$births)
}
lab <- c("18:p0","18:p1","18:p2+", "18:a0","18:a1","18:a2","18:a3","18:a4+", "20-35:p0","20-35:p1","20-35:p2+","unrealized share")
rows <- list()
for (cfg in list(c(25, 0), c(1, 0), c(1, 1), c(.25, 1))) {
  res <- parallel::mclapply(9701:9710, function(s) run(cfg[1], cfg[2] == 1, s), mc.cores = 5)
  M <- do.call(rbind, res)
  cat(sprintf("\nwidth %g, split parity %s\n", cfg[1], cfg[2] == 1))
  print(rbind(mean = round(colMeans(M), 4), mcse = round(apply(M, 2, sd) / sqrt(10), 4)) |> `colnames<-`(lab))
  rows[[length(rows) + 1]] <- data.frame(cell_width = cfg[1], actual_parity_cells = cfg[2] == 1,
    measure = lab, target = c(.2, .3, .4, .2, .3, .4, .4, .4, .2, .3, .4, NA),
    mean = colMeans(M), mcse = apply(M, 2, sd) / sqrt(10), n = nrow(M))
}
write.csv(do.call(rbind, rows), "validation/milestone5/cell-resolution-summary.csv", row.names = FALSE)
