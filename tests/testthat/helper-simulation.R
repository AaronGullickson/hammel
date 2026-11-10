# General demographic fixtures explicitly preserve unmatched births. Tests of
# cancellation and constructors separately exercise the production default TRUE.
make_model <- function(n = 1000, sex = "m", age = 20, death_rate = 0,
           birth_rate = 0, birth_interval = 0.75) {
 persons <- data.frame(person_id = seq_len(n), sex = rep(sex, n),
  birth_time = rep(-age, n), death_time = rep(NA_real_, n),
  mother_id = rep(NA_integer_, n), father_id = rep(NA_integer_, n))
 hammel_model() |>
  set_population(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")) |>
  set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
              age_end = Inf, rate = death_rate)) |>
  set_fertility(data.frame(age_start = 15, age_end = 50, rate = birth_rate),
         birth_interval = birth_interval) |>
  set_union_formation("on_demand",
    matching = matching_model(~ age_difference, c(age_difference = 0)),
    cancel_birth_if_no_match = FALSE)
}

make_rate_model <- function(n = 2, search_rate = 1, birth_rate = 0, birth_interval = 0) {
  model <- make_model(n, birth_rate = birth_rate, birth_interval = birth_interval)
  if (n) model$population$persons$sex[1] <- "f"
  rates <- model$processes$fertility$rates
  model <- set_fertility(model, rbind(
    transform(rates, union_type = "partnered"),
    transform(rates, union_type = "unpartnered", rate = 0)), birth_interval)
  set_union_formation(model, "rate",
    matching_model(~ age_difference, c(age_difference = 0)),
    rates = data.frame(age_start = 15, age_end = Inf, rate = search_rate),
    cancel_birth_if_no_match = FALSE)
}

# Simple synthetic parity schedules often give positive rates in young cells
# that no eligible woman can reach; silence only that advisory warning.
quiet_unreachable <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("no woman can be eligible", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

# Models that differ only in how equal rates are stratified create different
# cells, so their histories differ draw by draw while agreeing in distribution.
# Compare union and dissolution counts over replications.
expect_equivalent_unions <- function(a, b, years, seeds = 1:20) {
  counts <- function(model) vapply(seeds, function(seed) {
    unions <- run_simulation(model, years, seed)$unions
    c(nrow(unions), sum(unions$end_reason %in% "dissolution"))
  }, numeric(2))
  x <- counts(a)
  y <- counts(b)
  for (k in 1:2) {
    se <- sqrt(var(x[k, ]) / length(seeds) + var(y[k, ]) / length(seeds))
    expect_lt(abs(mean(x[k, ]) - mean(y[k, ])), 4 * se + 1e-9)
  }
}

# n copies of a woman aged `age` whose sons were born at `sons`, so a birth
# interval can keep every woman ineligible; unrealized births are then counted
# events with expected value n times the single-woman intensity.
spacing_population <- function(n, age = 30, sons = -0.1, groups = NULL, end_time = 0) {
  k <- length(sons)
  persons <- data.frame(person_id = seq_len(n * (k + 1)),
    sex = c(rep("f", n), rep("m", n * k)),
    birth_time = c(rep(end_time - age, n), end_time + rep(sons, times = n)),
    death_time = NA_real_,
    mother_id = c(rep(NA_real_, n), rep(seq_len(n), each = k)),
    father_id = NA_real_)
  if (!is.null(groups)) persons$group <- c(rep(groups[1], n), rep(groups[2], n * k))
  hammel_population(persons, end_time = end_time, groups = if ("group" %in% names(persons)) "group")
}
# Poisson check of a count against its expected value.
expect_count <- function(count, expected) {
  expect_lt(abs(sum(count) - expected), 4 * sqrt(expected) + 1e-9)
}

# Mortality and fertility counts in one table (the shape of the former internal
# collect_demographic_counts()): births are NA for men.
demographic_counts <- function(population, time_breaks, age_breaks, by = NULL) {
  out <- calculate_mortality_rates(population, time_breaks, age_breaks, by = by)
  births <- calculate_fertility_rates(population, time_breaks, age_breaks, by = by)
  keys <- setdiff(names(births), c("person_years", "births", "rate_fertility"))
  key <- function(x) do.call(paste, c(lapply(x[keys], as.character), sep = "\r"))
  out$births <- births$births[match(key(out), key(births))]
  out$rate_fertility <- as.numeric(ifelse(out$person_years > 0,
    out$births / out$person_years, NA_real_))
  out <- out[c(setdiff(names(out), c("deaths", "births", "rate_mortality", "rate_fertility")),
               "deaths", "births", "rate_mortality", "rate_fertility")]
  attr(out, "unassigned_births") <- attr(births, "unassigned_births")
  out
}

# An unordered pair label for two group labels, such as "A:B" for B and A.
unordered_pair <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  ifelse(a <= b, paste(a, b, sep = ":"), paste(b, a, sep = ":"))
}

# Sum directional wife/husband dissolution rows into unordered group pairs.
pair_rates <- function(x) {
  x$group_pair <- unordered_pair(x$group_w, x$group_h)
  keys <- intersect(c("time_start", "time_end", "dur_start", "dur_end", "group_pair", "union_type"),
                    names(x))
  aggregate(x[c("union_years", "dissolutions")], x[keys], sum)
}

# The recorded histories of a result, without its group declaration (which
# records the inheritance rule and so differs between equivalent models).
histories <- function(x) x[intersect(c("persons", "unions", "ancestry", "transitions"), names(x))]

# Change the inheritance rule of a model's (single) declared group dimension.
# `extra` appends categories (such as a custom rule's outputs) to the current
# ones, or to the population's values when none are listed.
set_inheritance <- function(model, rule, extra = NULL) {
  declaration <- hammel:::model_groups(model)
  if (is.null(declaration)) declaration <- list(group = group())
  current <- declaration$group
  categories <- current$categories
  if (!is.null(extra)) {
    if (is.null(categories)) categories <- unique(model$population$persons$group)
    categories <- c(categories, extra)
  }
  declaration$group <- group(categories = categories, ancestry = current$ancestry,
                             inheritance = rule)
  do.call(set_groups, c(list(model), declaration))
}

# Expand symmetric `group_pair` rate rows ("A:B") into directional rows with
# `group_w` and `group_h`, both directions sharing the row's rate.
directional_pairs <- function(rates) {
  out <- do.call(rbind, lapply(seq_len(nrow(rates)), function(i) {
    p <- strsplit(as.character(rates$group_pair[i]), ":", fixed = TRUE)[[1]]
    same <- p[1] == p[2]
    r <- rates[rep(i, if (same) 1L else 2L), , drop = FALSE]
    r$group_w <- if (same) p[1] else p
    r$group_h <- if (same) p[2] else rev(p)
    r$group_pair <- NULL
    r
  }))
  rownames(out) <- NULL
  out
}

# Translate a public matching model for a single dimension named `group` with
# categories `labels`, as run_simulation() does, for tests of internal helpers.
as_internal_matching <- function(model, labels) {
  if (is.null(labels)) return(internal_matching_model(model, NULL, NULL))
  dim <- if (length(matching_dimensions(model, "group"))) "group"
  internal_matching_model(model, dim, labels)
}

# Count recorded events strictly inside the latest run's interval. Input
# history at its starting instant is excluded; run-end events are not executed.
count_run_events <- function(population, event) {
  time <- switch(
    event,
    births = population$persons$birth_time,
    deaths = population$persons$death_time,
    dissolutions = population$unions$end_time[
      population$unions$end_reason %in% "dissolution"
    ],
    conversions = population$unions$end_time[
      population$unions$end_reason %in% "conversion"
    ],
    transitions = population$transitions$time
  )
  as.numeric(sum(
    time > population$simulation$start_time &
      time < population$simulation$end_time,
    na.rm = TRUE
  ))
}
