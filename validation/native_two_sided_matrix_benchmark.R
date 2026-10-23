# Compare native bulk calculations with explicit vectorized R matrix math.
# Run after an optimized build, or build cleanly here. No rate downloads.
# TZ=America/Los_Angeles Rscript validation/native_two_sided_matrix_benchmark.R
pkgbuild::clean_dll()
pkgbuild::compile_dll(debug = FALSE, force = TRUE, quiet = TRUE)
devtools::load_all(compile = FALSE, quiet = TRUE)
output_dir <- Sys.getenv("HAMMEL_MATRIX_BENCH_DIR", "/tmp/hammel-native-matrix-benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Both implementations use identical quadratic utilities, group effects, and
# additive directional probabilities. Shift each column to avoid underflow.
draw_r_matrix <- function(fb, mb, fg, mg, draws) {
  vapply(seq_len(draws), function(k) {
    gap <- outer(fb, mb, "-")
    utility <- (-log(2) / 25) * gap^2 - 2.5 * outer(fg, mg, "!=")
    columns <- exp(sweep(utility, 2, apply(utility, 2, max), "-"))
    b <- columns[1, ] / colSums(columns)
    a <- exp(utility[1, ] - max(utility[1, ])); a <- a / sum(a)
    weights <- a + b
    which(runif(1) * sum(weights) < cumsum(weights))[1]
  }, integer(1))
}
prepared <- prepare_matching_model(create_age_gap_matching_model(0, 5,
  group_exogamy = c("A::B" = -2.5)), c("A", "B"), TRUE)
rows <- list()
for (size in c(50L, 250L, 1000L)) {
  set.seed(721 + size)
  fb <- -runif(size, 15, 65); fb[1] <- -25
  mb <- -runif(size, 15, 65)
  fg <- as.integer(runif(size) < .25); mg <- as.integer(runif(size) < .25)
  draws <- if (size <= 250) 100L else 20L
  for (trial in 1:3) {
    answers <- list()
    for (method in if (trial %% 2) c("native_bulk", "R_matrix") else c("R_matrix", "native_bulk")) {
      invisible(gc()); set.seed(730 + trial)
      elapsed <- system.time(answer <- if (method == "native_bulk")
        inspect_two_sided_matching_cpp(fb, mb, fg, mg, prepared, draws)$choices else
        draw_r_matrix(fb, mb, fg, mg, draws))[["elapsed"]]
      answers[[method]] <- answer
      rows[[length(rows) + 1L]] <- data.frame(size, trial, method, draws, elapsed,
        milliseconds_per_market = 1000 * elapsed / draws)
    }
    stopifnot(identical(answers$native_bulk, answers$R_matrix))
  }
}
result <- do.call(rbind, rows)
write.csv(result, file.path(output_dir, "matrix_comparison.csv"), row.names = FALSE)
saveRDS(list(time = Sys.time(), session = sessionInfo(),
  hashes = tools::md5sum(c("src/simulate.cpp", "validation/native_two_sided_matrix_benchmark.R"))),
  file.path(output_dir, "metadata.rds"))
print(aggregate(milliseconds_per_market ~ size + method, result, median))
