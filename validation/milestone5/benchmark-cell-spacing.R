# Offline timing comparison of midpoint compensation and cell-accounted
# fertility (D101) on identical synthetic inputs; one fresh process per run.
# Install an optimized build first (pkgload::load_all() compiles with -O0 and
# library assertions, which distorts the comparison), then run from the
# repository root:
#   R CMD INSTALL --no-multiarch --library=LIB .
#   Rscript validation/milestone5/benchmark-cell-spacing.R LIB
# Writes cell-spacing-benchmark.csv. Describes this machine only.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 4) {
  # Worker: library, workload, method, seed -> elapsed seconds and births.
  suppressPackageStartupMessages(library(hammel, lib.loc = args[1]))
  workload <- args[2]; method <- args[3]; seed <- as.integer(args[4])
  set.seed(seed)
  pop <- create_baseline_population(50000, max_age = 60)
  mort <- data.frame(sex = rep(c("m", "f"), each = 4), age_start = c(0, 1, 50, 70),
    age_end = c(1, 50, 70, Inf), rate = c(.005, .002, .02, .08))
  # HFD-shaped: single years 12-54 and an upper-open 55+ band.
  ages <- 12:54
  base <- c(.11 * exp(-((ages - 28) / 7)^2), .0001)
  starts <- c(ages, 55); ends <- c(ages + 1, Inf)
  fert <- if (workload == "single-year ages") {
    data.frame(age_start = starts, age_end = ends, rate = base)
  } else {
    do.call(rbind, lapply(0:4, function(k) data.frame(age_start = starts, age_end = ends,
      parity_start = k, parity_end = if (k == 4) Inf else k + 1,
      rate = base * c(.9, 1.3, .7, .6, .6)[k + 1])))
  }
  births <- suppressWarnings(fertility(fert, spacing_method = method))
  model <- hammel_model() |> set_population(pop) |> set_mortality(mort) |>
    set_process(births) |>
    set_union_formation("on_demand", cancel_birth_if_no_match = FALSE)
  simulate(model, 1, seed) # warm native code
  gc()
  t <- suppressWarnings(system.time(result <- simulate(model, 100, seed))[["elapsed"]])
  cat(t, result$simulation$births, "\n")
  quit(save = "no")
}
lib <- if (length(args)) args[1] else stop("Supply the library holding an optimized hammel build.")
rows <- list()
for (workload in c("single-year ages", "single-year ages x 5 parities"))
  for (seed in c(42, 71, 105)) {
    order <- if (seed == 71) c("cell", "midpoint") else c("midpoint", "cell")
    for (method in order) {
      out <- system2("Rscript", c("validation/milestone5/benchmark-cell-spacing.R",
        shQuote(lib), shQuote(workload), method, seed), stdout = TRUE)
      v <- as.numeric(strsplit(trimws(tail(out, 1)), " ")[[1]])
      rows[[length(rows) + 1]] <- data.frame(workload, seed, method, seconds = v[1], births = v[2])
    }
  }
res <- do.call(rbind, rows)
write.csv(res, "validation/milestone5/cell-spacing-benchmark.csv", row.names = FALSE)
print(aggregate(seconds ~ workload + method, res, median))
