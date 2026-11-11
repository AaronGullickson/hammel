# Benchmark count_kin() by population size, kinship depth and `by`.
# Run from the package root against an optimized install, e.g.
#   R CMD INSTALL --preclean --library=<lib> .
#   Rscript validation/kinship/benchmark-count-kin.R <lib>
# A devtools::load_all() build compiles C++ without optimization and runs
# about twice as slowly.

lib <- commandArgs(trailingOnly = TRUE)[1]
library(hammel, lib.loc = lib)

seconds <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - start
}

times <- seq(0, 100, by = 10)
rows <- list()
for (n in c(5000, 20000, 50000, 100000)) {
  model <- hammel_model() |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("rate", rates = usa_2000_marriage) |>
    set_population(generate_uniform_population(n, end_time = 0))
  population <- suppressWarnings(run_simulation(model, 100, seed = 1))
  for (generations in 2:3) {
    for (by in list(NULL, "sex")) {
      elapsed <- median(replicate(3, seconds(count_kin(
        population, times,
        generations = generations, by = by
      ))))
      rows[[length(rows) + 1L]] <- data.frame(
        start = n, persons = nrow(population$persons),
        generations = generations,
        by = if (is.null(by)) "none" else by, seconds = elapsed
      )
    }
  }
}
print(do.call(rbind, rows), row.names = FALSE)
