# Synthetic performance comparison, separate from milestone 3 scientific
# assessment. Run from the package root with:
# TZ=America/Los_Angeles Rscript validation/two_sided_matching_benchmark.R
# Optional output directory: HAMMEL_BENCH_DIR=/tmp/hammel-two-sided-benchmark
# This uses local synthetic rates only and never renders or edits the report.

pkgbuild::compile_dll(debug = FALSE, force = TRUE, quiet = TRUE)
devtools::load_all(quiet = TRUE)
output_dir <- Sys.getenv("HAMMEL_BENCH_DIR", "/tmp/hammel-two-sided-benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matching <- create_age_gap_matching_model(0, 5, group_exogamy = c("A::B" = -2.5))
prepared <- prepare_matching_model(matching, c("A", "B"), TRUE)
rows <- list()
for (size in c(10L, 50L, 100L, 300L, 1000L)) {
  set.seed(721 + size)
  fb <- -runif(size, 15, 65); fb[1] <- -25
  mb <- -runif(size, 15, 65)
  fg <- as.integer(runif(size) < 0.25); mg <- as.integer(runif(size) < 0.25)
  draws <- if (size <= 100) 200L else if (size <= 300) 50L else 20L
  for (trial in 1:3) {
    answers <- list()
    for (mode in if (trial %% 2) c("on_demand", "eager") else c("eager", "on_demand")) {
      invisible(gc())
      set.seed(730 + trial)
      elapsed <- system.time(answer <- inspect_sequential_matching_cpp(
        fb, mb, fg, mg, prepared, draws, precompute = mode == "eager"))[["elapsed"]]
      answers[[mode]] <- answer
      rows[[length(rows) + 1L]] <- data.frame(size, trial, mode, draws,
        elapsed, milliseconds_per_market = 1000 * elapsed / draws,
        mean_steps = mean(answer$steps), mean_choice_pairs = mean(answer$scored_pairs))
    }
    stopifnot(identical(answers$on_demand, answers$eager))
  }
  message("Finished fixed-market size ", size)
}
micro <- tibble::as_tibble(do.call(rbind, rows))
write.csv(micro, file.path(output_dir, "fixed_markets.csv"), row.names = FALSE)

# These workloads include sampling women, screening, event scheduling, births,
# and mortality. The population size and events can differ between methods;
# retain them alongside runtime rather than claiming identical work.
set.seed(741)
population <- create_baseline_population(2000, max_age = 80, end_time = 0,
  group_weights = c(A = 3, B = 1))
death_rates <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
birth_rates <- data.frame(age_start = 15, age_end = 45, rate = 0.08)
base <- hammel_model() |> set_population(population) |>
  set_mortality(death_rates) |> set_fertility(birth_rates)
rows <- list()
for (sample_size in c(50L, 200L, 1000L)) {
  for (trial in 1:3) {
    methods <- if (trial %% 2) c("female_choice", "two_sided") else c("two_sided", "female_choice")
    for (method in methods) {
      model <- set_union_formation(base, "on_demand", matching, sample_size = sample_size,
        matching_method = method)
      invisible(gc())
      elapsed <- system.time(result <- simulate(model, 20, seed = 750 + trial))[["elapsed"]]
      rows[[length(rows) + 1L]] <- data.frame(sample_size, trial, method, years = 20,
        initial_size = nrow(population$persons), elapsed,
        final_history_size = nrow(result$persons), births = result$simulation$births,
        deaths = result$simulation$deaths, unions = nrow(result$unions))
    }
  }
  message("Finished simulation sample size ", sample_size)
}
simulations <- tibble::as_tibble(do.call(rbind, rows))
write.csv(simulations, file.path(output_dir, "simulations.csv"), row.names = FALSE)

# A pointwise-equivalent R formula measures callback overhead without changing
# the preference curve. Smaller workloads keep the script usable interactively.
fallback <- create_matching_model(~ I(abs(age_difference)^2),
  c("I(abs(age_difference)^2)" = -log(2) / 25), group_exogamy = c("A::B" = -2.5))
rows <- list()
for (sample_size in c(50L, 200L)) {
  for (trial in 1:3) {
    for (scorer in if (trial %% 2) c("native", "R_formula") else c("R_formula", "native")) {
      model <- set_union_formation(base, "on_demand",
        if (scorer == "native") matching else fallback,
        sample_size = sample_size, matching_method = "two_sided")
      invisible(gc())
      elapsed <- system.time(result <- simulate(model, 5, seed = 760 + trial))[["elapsed"]]
      rows[[length(rows) + 1L]] <- data.frame(sample_size, trial, scorer, years = 5,
        elapsed, births = result$simulation$births, unions = nrow(result$unions))
    }
  }
  message("Finished custom-formula sample size ", sample_size)
}
callbacks <- tibble::as_tibble(do.call(rbind, rows))
write.csv(callbacks, file.path(output_dir, "callbacks.csv"), row.names = FALSE)

metadata <- list(time = Sys.time(), commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  code_hashes = tools::md5sum(c("src/simulate.cpp", "R/union_formation.R", "R/simulate.R")),
  build = "pkgbuild::compile_dll(debug = FALSE): installed R compiler flags",
  session = sessionInfo())
saveRDS(metadata, file.path(output_dir, "metadata.rds"))
print(aggregate(milliseconds_per_market ~ size + mode, micro, median))
print(aggregate(elapsed ~ sample_size + method, simulations, median))
print(aggregate(elapsed ~ sample_size + scorer, callbacks, median))
message("Benchmark results: ", output_dir)
