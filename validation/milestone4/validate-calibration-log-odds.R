# D082 matching recovery. Offline synthetic inputs; run from repository root.
# Rscript validation/milestone4/validate-calibration-log-odds.R LIB OUTPUT_PREFIX [REPS] [CHOICES] [age|duration]
# Uses recorded engine weight offsets and independently reconstructed male pools.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
library(hammel, lib.loc = args[1])
stopifnot(requireNamespace("survival", quietly = TRUE))
suppressPackageStartupMessages(library(survival))
reps <- if (length(args) >= 3L) as.integer(args[3]) else 20L
choices <- if (length(args) >= 4L) as.integer(args[4]) else 1500L
clock <- if (length(args) >= 5L) match.arg(args[5], c("age", "duration")) else "age"

# Uniform sampling of alternatives excludes the chosen man, then includes him
# with certainty. Every possible chosen man has the same subset-selection
# probability conditional on a resulting sampled set, requiring no correction.
choice_data <- function(result, sample_size, alternatives = 50L) {
  p <- result$persons; u <- result$unions
  male <- which(p$sex == "m")
  ids <- p$person_id
  wives <- match(u$wife_id, ids)
  husbands <- match(u$husband_id, ids)
  selected <- sample.int(nrow(u), min(sample_size, nrow(u)))
  sets <- lapply(seq_along(selected), function(k) {
    row <- selected[k]; now <- u$start_time[row]
    # Immediately before this formation: exclude prior ongoing unions, but
    # include the chosen man whose new episode starts exactly at now.
    partnered <- husbands[u$start_time < now & (is.na(u$end_time) | u$end_time > now)]
    eligible <- male[p$birth_time[male] + 15 <= now &
      (is.na(p$death_time[male]) | p$death_time[male] > now) & !male %in% partnered]
    chosen <- husbands[row]
    stopifnot(chosen %in% eligible)
    pool <- eligible[eligible != chosen]
    sampled <- if (length(pool)) pool[sample.int(length(pool), min(alternatives, length(pool)))] else integer()
    candidates <- c(chosen, sampled)
    gf <- p$group[wives[row]]; gm <- p$group[candidates]
    data.frame(set = k, chosen = seq_along(candidates) == 1L,
      age_diff = p$birth_time[wives[row]] - p$birth_time[candidates],
      AB = as.numeric(gf == "A" & gm == "B"),
      BA = as.numeric(gf == "B" & gm == "A"),
      log_weight = ifelse(gm == "A", u$.log_weight_A[row], u$.log_weight_B[row]))
  })
  do.call(rbind, sets)
}

rows <- list()
for (replication in seq_len(reps)) {
  set.seed(282600L + replication)
  population <- create_baseline_population(20000, 50, 2000, group_weights = c(A = 3, B = 1))
  raw <- data.frame(age_start = c(15, 20, 25, 30, 35, 45, 60),
    age_end = c(20, 25, 30, 35, 45, 60, Inf), rate = c(.02, .08, .14, .1, .06, .03, .005))
  search <- rbind(transform(raw, previous_union = "none"),
    transform(raw, previous_union = "dissolved", rate = rate * 1.4),
    transform(raw, previous_union = "widowed", rate = rate * .7))
  if (clock == "duration") {
    none <- transform(raw, previous_union = "none", dur_start = NA_real_, dur_end = NA_real_)
    duration <- data.frame(age_start = NA_real_, age_end = NA_real_,
      dur_start = c(0, 1, 3, 8), dur_end = c(1, 3, 8, Inf))
    search <- rbind(none,
      transform(duration, previous_union = "dissolved", rate = c(0, .15, .1, .05)),
      transform(duration, previous_union = "widowed", rate = c(0, .07, .04, .02)))
  }
  for (case in c("group_only", "age_dependent")) {
    matching <- if (case == "group_only")
      create_matching_model(~ 1, numeric(), group_exogamy = c("A:B" = -3)) else
      create_age_gap_matching_model(2, 5, group_exogamy = c("A:B" = -3))
    model <- hammel_model() |>
      set_population(population) |>
      set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .01)) |>
      set_fertility(data.frame(age_start = 15, age_end = 50,
        union_type = c("partnered", "unpartnered"), rate = c(.12, .04))) |>
      set_union_formation("rate", matching, rates = search, diagnostics = TRUE,
        kinship_generations = 0L, group_prop_hazard = c(B = 1.2)) |>
      set_union_dissolution(data.frame(dur_start = c(0, 5, 15),
        dur_end = c(5, 15, Inf), rate = c(.06, .03, .015)))
    result <- simulate(model, 60, seed = 283600L + replication)
    set.seed(284600L + replication)
    data <- choice_data(result, choices)
    # Reconstruct eligibility without kin screening because it was explicitly
    # disabled above, rather than silently including ineligible alternatives.
    fit <- if (case == "group_only")
      survival::clogit(chosen ~ AB + BA + offset(log_weight) + strata(set),
        data = data, method = "efron") else
      survival::clogit(chosen ~ age_diff + I(age_diff^2) + AB + BA +
        offset(log_weight) + strata(set), data = data, method = "efron")
    b <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    log_odds <- unname(b["AB"] + b["BA"])
    log_odds_se <- sqrt(sum(covariance[c("AB", "BA"), c("AB", "BA")]))
    detailed <- get_detailed_unions(result)
    table <- table(factor(detailed$group_w, c("A", "B")),
                   factor(detailed$group_h, c("A", "B")))
    pooled_log_odds <- log(table[1, 2]) + log(table[2, 1]) - log(table[1, 1]) - log(table[2, 2])
    rows[[length(rows) + 1L]] <- data.frame(replication, case, clock,
      choices = length(unique(data$set)), unions = nrow(result$unions),
      AB = unname(b["AB"]), BA = unname(b["BA"]), log_odds, log_odds_se,
      linear = if (case == "age_dependent") unname(b["age_diff"]) else NA_real_,
      quadratic = if (case == "age_dependent") unname(b["I(age_diff^2)"]) else NA_real_,
      pooled_log_odds)
    message(case, " rep=", replication, " recovered log OR=", round(log_odds, 4))
  }
}
results <- do.call(rbind, rows)
write.csv(results, paste0(args[2], "-replications.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(results, results$case), function(x) {
  do.call(rbind, lapply(c("AB", "BA", "log_odds", "linear", "quadratic", "pooled_log_odds"), function(parameter) {
    values <- x[[parameter]]
    data.frame(case = x$case[1], parameter, replications = sum(!is.na(values)),
      mean = mean(values), mcse = sd(values) / sqrt(sum(!is.na(values))))
  }))
}))
write.csv(summary, paste0(args[2], "-summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
saveRDS(list(results = results, summary = summary, session = sessionInfo()), paste0(args[2], ".rds"))
