# Offline performance comparison for D082. Run from the repository root.
# Rscript validation/benchmark-calibration-intensity.R make LIB INPUT.rds
# Rscript validation/benchmark-calibration-intensity.R run LIB INPUT.rds OUTPUT.rds
# Install the baseline commit and working tree into separate optimized libraries.
# Rscript validation/benchmark-calibration-intensity.R compare OLDLIB INPUT.rds NEWLIB OUTPUT.rds
# Rscript validation/benchmark-calibration-intensity.R memory LIB INPUT.rds OUTPUT.rds
# All schedules below are synthetic workload inputs, not empirical estimates.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
if (args[1] == "compare") {
  stopifnot(length(args) == 5L)
  input <- readRDS(args[3])
  files <- character()
  for (i in seq_along(input$seeds)) {
    part <- input; part$seeds <- input$seeds[i]
    input_file <- tempfile(fileext = ".rds"); saveRDS(part, input_file)
    order <- if (i %% 2L) c("old", "new") else c("new", "old")
    for (version in order) {
      output <- tempfile(fileext = ".rds")
      library_path <- if (version == "old") args[2] else args[4]
      status <- system2(file.path(R.home("bin"), "Rscript"),
        shQuote(c("validation/benchmark-calibration-intensity.R", "run",
          library_path, input_file, output)))
      stopifnot(status == 0L)
      measured <- readRDS(output)
      measured$results$version <- version
      files <- c(files, output)
      saveRDS(measured, output)
    }
  }
  measured <- lapply(files, readRDS)
  saveRDS(list(results = do.call(rbind, lapply(measured, `[[`, "results")),
    sessions = lapply(measured, `[[`, "session"),
    input_md5 = unname(tools::md5sum(args[3]))), args[5])
  quit(status = 0L)
}
library(hammel, lib.loc = args[2])
if (args[1] == "memory") {
  stopifnot(length(args) == 4L, file.exists("/proc/self/status"))
  input <- readRDS(args[3]); gc()
  result <- simulate(input$models$grouped_previous, input$years, input$seeds[1])
  hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  saveRDS(list(peak_kib = as.numeric(gsub("[^0-9]", "", hwm)),
    births = result$simulation$births, unions = nrow(result$unions)), args[4])
  quit(status = 0L)
}
if (args[1] == "make") {
  set.seed(280126)
  population <- create_baseline_population(50000, 50, 2000,
    group_weights = c(A = 3, B = 1))
  ungrouped <- population$persons
  ungrouped$group <- NULL
  ungrouped <- hammel_population(ungrouped, end_time = 2000)
  ages <- 0:100
  mortality <- data.frame(sex = rep(c("f", "m"), each = length(ages)),
    age_start = ages, age_end = c(1:100, Inf),
    rate = rep(.00015 * exp(.085 * ages), 2))
  overall <- data.frame(age_start = seq(15, 40, 5), age_end = seq(20, 45, 5),
    rate = c(.03, .09, .12, .1, .04, .008))
  fertility <- rbind(transform(overall, union_type = "partnered", rate = rate * 1.5),
    transform(overall, union_type = "unpartnered", rate = rate * .5))
  search <- data.frame(age_start = seq(15, 65, 5), age_end = c(seq(20, 65, 5), Inf),
    rate = c(.02, .08, .14, .1, .05, .05, .04, .025, .01, .007, .002))
  previous <- rbind(transform(search, previous_union = "none"),
    transform(search, previous_union = "dissolved", rate = rate * .8),
    transform(search, previous_union = "widowed", rate = rate * .5))
  dissolution <- data.frame(dur_start = c(0, 5, 10, 15, 20),
    dur_end = c(5, 10, 15, 20, Inf), rate = c(.045, .033, .025, .029, .02))
  make_model <- function(population, rates = search, initiation = "rate") {
    model <- hammel_model() |>
      set_population(population) |>
      set_mortality(mortality) |>
      set_fertility(if (initiation == "rate") fertility else overall) |>
      set_union_formation(initiation, rates = if (initiation == "rate") rates else NULL,
        matching = if ("group" %in% names(population$persons))
          create_age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -2)) else
          create_age_gap_matching_model(0, 5)) |>
      set_union_dissolution(dissolution)
    model
  }
  models <- list(ungrouped_rate = make_model(ungrouped),
    grouped_rate = make_model(population), grouped_previous = make_model(population, previous),
    grouped_on_demand = make_model(population, initiation = "on_demand"))
  saveRDS(list(models = models, years = 100, seeds = c(42L, 71L, 105L)), args[3])
} else {
  stopifnot(args[1] == "run", length(args) == 4L)
  input <- readRDS(args[3])
  rows <- list()
  # Warm native code and allocators without including warm-up in reported times.
  invisible(simulate(input$models[[1]], 1, seed = 1))
  for (seed in input$seeds) for (case in names(input$models)) {
    gc()
    elapsed <- system.time(result <- simulate(input$models[[case]], input$years, seed))[["elapsed"]]
    unions <- get_detailed_unions(result)
    a_share <- if ("group" %in% names(result$persons)) {
      living <- is.na(result$persons$death_time)
      mean(result$persons$group[living] == "A")
    } else NA_real_
    pairing_ratio <- if ("group_w" %in% names(unions))
      sum(unions$group_w == "A" & unions$group_h == "B") /
        sum(unions$group_w == "B" & unions$group_h == "A") else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(case, seed, elapsed,
      births = result$simulation$births, deaths = result$simulation$deaths,
      unions = nrow(result$unions), dissolutions = result$simulation$dissolutions,
      a_share, pairing_ratio)
    message(case, " seed=", seed, " elapsed=", elapsed)
    rm(result, unions)
  }
  saveRDS(list(results = do.call(rbind, rows), session = sessionInfo(),
    library = normalizePath(args[2]), input_md5 = unname(tools::md5sum(args[3]))), args[4])
}
