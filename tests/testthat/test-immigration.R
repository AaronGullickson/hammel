# External immigration, phase 1 (D184, D185).

# A template of unpartnered people with the given ages and sexes at end_time 0,
# optionally with group columns declared as groups.
immigration_template <- function(ages, sexes = "f", ..., groups = NULL,
                                 ancestry = NULL, unions = NULL,
                                 death_time = NA_real_) {
  n <- length(ages)
  persons <- data.frame(
    person_id = seq_len(n),
    sex = rep(sexes, length.out = n),
    birth_time = -ages,
    death_time = rep(death_time, length.out = n),
    mother_id = NA_integer_,
    father_id = NA_integer_,
    ...
  )
  hammel_population(
    persons,
    unions = unions,
    end_time = 0,
    groups = groups,
    ancestry = ancestry
  )
}

# The destination's race declaration, which templates must share.
race_group <- function() {
  group(categories = c("A", "B"))
}

# A destination where nothing happens except arrivals.
quiet_model <- function(n = 10, ...) {
  make_model(n, ...)
}

arrivals_of <- function(result) {
  result$persons[!is.na(result$persons$arrival_time), ]
}

test_that("immigration processes need named, unique streams", {
  stream <- immigration_stream(immigration_template(30), 1)
  process <- immigration(a = stream, b = stream)
  expect_s3_class(process, c("hammel_immigration", "hammel_process"))
  expect_identical(names(process$streams), c("a", "b"))
  expect_hammel_error(immigration(stream), "specification", "immigration")
  expect_hammel_error(
    immigration(a = stream, a = stream), "specification", "immigration"
  )
  expect_hammel_error(immigration(a = 1), "specification", "immigration")
  expect_hammel_error(
    set_immigration(hammel_model(), stream), "specification",
    "set_immigration"
  )
})

test_that("set_immigration replaces every stream and can remove them", {
  stream <- immigration_stream(immigration_template(30), 1)
  model <- set_immigration(hammel_model(), a = stream, b = stream)
  expect_identical(names(model$processes$immigration$streams), c("a", "b"))
  model <- set_immigration(model, c = stream)
  expect_identical(names(model$processes$immigration$streams), "c")
  model <- set_immigration(model)
  expect_null(model$processes$immigration)
  model <- set_process(model, immigration(d = stream))
  expect_identical(names(model$processes$immigration$streams), "d")
})

test_that("stream intensities are validated", {
  template <- immigration_template(30)
  expect_hammel_error(
    immigration_stream(template, -1), "argument", "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template, c(1, 2)), "argument", "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template), "argument", "immigration_stream"
  )
  table <- data.frame(time_start = c(0, 5), time_end = c(5, Inf),
                      intensity = c(1, 2))
  expect_s3_class(
    immigration_stream(template, table), "hammel_immigration_stream"
  )
  overlap <- transform(table, time_start = c(0, 4))
  error <- expect_hammel_error(
    immigration_stream(template, overlap), "rates", "immigration_stream"
  )
  expect_identical(error$rows, 1:2)
  expect_hammel_error(
    immigration_stream(template, transform(table, sex = "f")), "rates",
    "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template, transform(table, intensity = c(1, -1))),
    "table", "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template, table[0, ]), "rates", "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(data.frame(), 1), "population", "immigration_stream"
  )
})

test_that("stream weights and assignments are validated", {
  template <- immigration_template(30)
  expect_hammel_error(
    immigration_stream(
      template, 1,
      weights = weight_model(~parity, c(parity = 1))
    ),
    "specification", "immigration_stream"
  )
  schedule <- tibble::tibble(
    time_start = 0, time_end = Inf,
    model = list(weight_model(~age, c(age = 1)))
  )
  expect_hammel_error(
    immigration_stream(template, 1, weights = schedule), "specification",
    "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template, 1, assignment = migrant_assignment(
      "fixed", fixed_value = "a"
    )),
    "specification", "immigration_stream"
  )
  expect_hammel_error(
    immigration_stream(template, 1, assignment = list(region = 1)),
    "specification", "immigration_stream"
  )
})

test_that("migrant assignments are validated", {
  expect_s3_class(
    migrant_assignment("fixed", fixed_value = "a"), "hammel_migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("mother"), "argument", "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("fixed"), "argument", "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("fixed", fixed_value = c("a", "b")), "argument",
    "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("fixed", fixed_value = "a", probs = c(a = 1)),
    "argument", "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("probability", probs = c(a = 0.5, b = 0.4)),
    "argument", "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("probability", probs = c(0.5, 0.5)), "argument",
    "migrant_assignment"
  )
  expect_hammel_error(
    migrant_assignment("custom", fun = "f"), "argument", "migrant_assignment"
  )
})

test_that("arrival counts follow the stream intensity", {
  template <- immigration_template(c(20, 30, 40))
  model <- set_immigration(
    quiet_model(), a = immigration_stream(template, intensity = 30)
  )
  counts <- vapply(1:30, function(seed) {
    nrow(arrivals_of(run_simulation(model, 2, seed = seed)))
  }, numeric(1))
  # Poisson(60) per run; the mean of 30 runs has standard error sqrt(2).
  expect_lt(abs(mean(counts) - 60), 4 * sqrt(60 / 30))
  expect_gt(var(counts), 0)
})

test_that("calendar intensities apply piecewise and zero pauses a stream", {
  template <- immigration_template(c(20, 30))
  intensity <- data.frame(time_start = c(0, 2, 4), time_end = c(2, 4, Inf),
                          intensity = c(50, 0, 10))
  model <- set_immigration(
    quiet_model(), a = immigration_stream(template, intensity)
  )
  times <- unlist(lapply(1:10, function(seed) {
    arrivals_of(run_simulation(model, 6, seed = seed))$arrival_time
  }))
  expect_false(any(times >= 2 & times < 4))
  expect_count(sum(times < 2), 10 * 100)
  expect_count(sum(times >= 4), 10 * 20)
  gap <- data.frame(time_start = c(0, 3), time_end = c(2, Inf),
                    intensity = 1)
  expect_hammel_error(
    run_simulation(
      set_immigration(quiet_model(), a = immigration_stream(template, gap)), 5
    ),
    "coverage", "run_simulation"
  )
})

test_that("migrants keep age and sex and arrive without relatives", {
  template <- immigration_template(c(25, 61), c("f", "m"))
  model <- set_immigration(
    quiet_model(), origin = immigration_stream(template, 40)
  )
  result <- run_simulation(model, 3, seed = 4)
  arrived <- arrivals_of(result)
  expect_gt(nrow(arrived), 0)
  age <- arrived$arrival_time - arrived$birth_time
  expect_true(all(abs(age - ifelse(arrived$sex == "f", 25, 61)) < 1e-9))
  expect_true(all(is.na(arrived$mother_id) & is.na(arrived$father_id)))
  expect_identical(unique(arrived$birthplace), "origin")
  expect_false(any(duplicated(result$persons$person_id)))
  expect_true(all(arrived$person_id > 10))
  # The original people have no arrival record and are native.
  expect_true(all(is.na(result$persons$arrival_time[1:10])))
  expect_true(all(result$persons$birthplace[1:10] == "native"))
  # The template is untouched.
  expect_identical(
    model$processes$immigration$streams$origin$template, template
  )
})

test_that("cloning streams draw anyone alive without their histories", {
  unions <- data.frame(union_id = 1, wife_id = 1, husband_id = 2,
                       start_time = -5, end_time = NA_real_,
                       end_reason = NA_character_)
  template <- immigration_template(
    c(30, 32, 40, 50, 70), c("f", "m", "f", "m", "f"), unions = unions,
    death_time = c(NA, NA, NA, NA, -1)
  )
  template$persons$mother_id[3] <- 1
  expect_identical(migrant_parties(template, NULL)$rows, 1:4)
  # A template with histories warns that clones leave them behind.
  region <- list(region = migrant_assignment("fixed", fixed_value = "n"))
  expect_warning(
    stream <- immigration_stream(
      template, 50,
      assignment = region, household = NULL
    ),
    class = "hammel_warning_no_effect"
  )
  model <- make_model(10)
  model$population$persons$region <- "n"
  model <- model |>
    set_groups(region = group(categories = "n")) |>
    set_immigration(a = stream)
  result <- run_simulation(model, 2, seed = 2)
  arrived <- arrivals_of(result)
  age <- round(arrived$arrival_time - arrived$birth_time)
  # Partnered people are cloned too; the dead are not.
  expect_setequal(unique(age), c(30, 32, 40, 50))
  # Clones arrive never partnered, at parity zero, with unknown parents and
  # no imported records; each group is entered at arrival.
  expect_false(any(c(result$unions$wife_id, result$unions$husband_id) %in%
                     arrived$person_id))
  expect_true(all(is.na(arrived$mother_id)))
  expect_identical(NROW(result$external_relatives), 0L)
  entry <- result$transitions[
    result$transitions$person_id %in% arrived$person_id,
  ]
  # One entry per arrival in region and in the automatic birthplace.
  expect_identical(nrow(entry), 2L * nrow(arrived))
  expect_setequal(entry$group, c("region", "birthplace"))
  expect_true(all(is.na(entry$origin)))
  expect_equal(
    entry$time, arrived$arrival_time[match(entry$person_id, arrived$person_id)]
  )
  # Templates without histories do not warn.
  expect_no_warning(
    immigration_stream(immigration_template(30), 1, household = NULL)
  )
  # Nobody alive is an error before the run.
  nobody <- immigration_template(30, death_time = -1)
  expect_hammel_error(
    run_simulation(
      set_immigration(quiet_model(), a = immigration_stream(nobody, 1)), 1
    ),
    "population", "run_simulation"
  )
})

test_that("weights select template people in proportion", {
  members <- tibble::tibble(person_id = 1:2, age = c(20, 40), sex = "f")
  model <- weight_model(~age, c(age = log(3) / 20))
  p <- migrant_probabilities(model, members, list(), "a")
  expect_equal(p, c(0.25, 0.75))
  # Scaling every weight leaves the selection unchanged.
  shifted <- weight_model(fun = function(members, bin) {
    members$age * log(3) / 20 + 7
  })
  expect_equal(migrant_probabilities(shifted, members, list(), "a"), p)
  expect_equal(migrant_probabilities(NULL, members, list(), "a"), c(0.5, 0.5))
  bad <- weight_model(fun = function(members, bin) 1)
  expect_error(
    migrant_probabilities(bad, members, list(), "a"),
    class = "hammel_error_runtime"
  )
  # In simulation, draws follow the probabilities.
  template <- immigration_template(c(20, 40))
  sim <- set_immigration(
    quiet_model(), a = immigration_stream(template, 200, weights = model)
  )
  arrived <- arrivals_of(run_simulation(sim, 5, seed = 3))
  old <- sum(round(arrived$arrival_time - arrived$birth_time) == 40)
  n <- nrow(arrived)
  expect_lt(abs(old - 0.75 * n), 4 * sqrt(n * 0.75 * 0.25))
})

test_that("weights may use sex and shared groups", {
  members <- tibble::tibble(
    person_id = 1:4, age = 30, sex = c("f", "m", "f", "m"),
    region = c("n", "n", "s", "s")
  )
  model <- weight_model(~ sex + region, c(sexm = log(2), regions = log(3)))
  p <- migrant_probabilities(
    model, members, list(region = c("n", "s")), "a"
  )
  expect_equal(p, c(1, 2, 3, 6) / 12)
})

# A destination declaring `region` (ancestry tracked) and `race`.
grouped_destination <- function() {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "m", "f", "m"), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    region = c("n", "n", "s", "s"), race = c("A", "B", "A", "B")
  )
  population <- hammel_population(persons, end_time = 0)
  make_model(4) |>
    set_population(population) |>
    set_groups(
      region = group(categories = c("n", "s", "w"), ancestry = TRUE),
      race = group(categories = c("A", "B"))
    )
}

test_that("shared groups come from the template; missing ones are assigned", {
  template <- immigration_template(
    c(20, 30), race = c("A", "B"), education = c("hs", "college"),
    groups = c("race", "education")
  )
  model <- grouped_destination()
  stream <- immigration_stream(
    template, 100,
    assignment = list(region = migrant_assignment("fixed", fixed_value = "w"))
  )
  result <- run_simulation(set_immigration(model, a = stream), 1, seed = 9)
  arrived <- arrivals_of(result)
  expect_true(all(arrived$region == "w"))
  age <- round(arrived$arrival_time - arrived$birth_time)
  expect_identical(arrived$race, ifelse(age == 20, "A", "B"))
  expect_false("education" %in% names(result$persons))
  # Ancestry for an assigned dimension is all in the assigned category.
  shares <- result$ancestry$region
  shares <- shares[match(arrived$person_id, shares$person_id), ]
  expect_true(all(shares$w == 1 & shares$n == 0 & shares$unknown == 0))
})

test_that("assignment probabilities are drawn per migrant", {
  template <- immigration_template(
    c(20, 60), race = "A", groups = list(race = race_group())
  )
  model <- grouped_destination()
  probability <- immigration_stream(
    template, 300,
    assignment = list(region = migrant_assignment(
      "probability", probs = c(n = 0.25, s = 0.75)
    ))
  )
  arrived <- arrivals_of(
    run_simulation(set_immigration(model, a = probability), 2, seed = 1)
  )
  n <- nrow(arrived)
  expect_false(any(arrived$region == "w"))
  expect_lt(abs(sum(arrived$region == "s") - 0.75 * n), 4 * sqrt(n * 0.1875))
  custom <- immigration_stream(
    template, 300,
    assignment = list(region = migrant_assignment(
      "custom", fun = function(members) {
        ifelse(members$age < 40, "n", "w")
      }
    ))
  )
  arrived <- arrivals_of(
    run_simulation(set_immigration(model, a = custom), 1, seed = 1)
  )
  age <- round(arrived$arrival_time - arrived$birth_time)
  expect_identical(arrived$region, ifelse(age == 20, "n", "w"))
})

test_that("template ancestry is copied where the template tracks it", {
  ancestry <- data.frame(person_id = 1:2, n = c(0.5, 0), s = c(0.25, 1),
                         unknown = c(0.25, 0))
  template <- immigration_template(
    c(20, 40), region = c("n", "s"), race = "A",
    groups = list(region = group(ancestry = TRUE), race = race_group()),
    ancestry = ancestry
  )
  model <- grouped_destination()
  model <- set_groups(
    model,
    region = group(categories = c("n", "s"), ancestry = TRUE),
    race = group(categories = c("A", "B"))
  )
  result <- run_simulation(
    set_immigration(model, a = immigration_stream(template, 100)), 1, seed = 2
  )
  arrived <- arrivals_of(result)
  shares <- result$ancestry$region
  shares <- shares[match(arrived$person_id, shares$person_id), ]
  young <- round(arrived$arrival_time - arrived$birth_time) == 20
  expect_true(all(shares$n[young] == 0.5 & shares$unknown[young] == 0.25))
  expect_true(all(shares$s[!young] == 1))
  # Without template ancestry, all ancestry is in the migrant's category.
  plain <- immigration_template(
    c(20, 40), region = c("n", "s"), race = "A",
    groups = list(region = group(), race = race_group())
  )
  result <- run_simulation(
    set_immigration(model, a = immigration_stream(plain, 100)), 1, seed = 2
  )
  arrived <- arrivals_of(result)
  shares <- result$ancestry$region
  shares <- shares[match(arrived$person_id, shares$person_id), ]
  expect_identical(
    shares$n, as.numeric(arrived$region == "n")
  )
})

test_that("group compatibility errors name the problem", {
  model <- grouped_destination()
  lacking <- immigration_template(30, race = "A",
                                  groups = list(race = race_group()))
  error <- expect_hammel_error(
    run_simulation(
      set_immigration(model, a = immigration_stream(lacking, 1)), 1
    ),
    "specification", "run_simulation"
  )
  expect_identical(error$group, "region")
  overriding <- immigration_stream(
    immigration_template(30, race = "A", region = "n",
                         groups = list(
                           race = race_group(),
                           region = group(categories = c("n", "s", "w"))
                         )),
    1,
    assignment = list(region = migrant_assignment("fixed", fixed_value = "n"))
  )
  error <- expect_hammel_error(
    run_simulation(set_immigration(model, a = overriding), 1),
    "specification", "run_simulation"
  )
  expect_identical(error$group, "region")
  mismatch <- immigration_template(
    30, race = "A", groups = list(race = group(categories = c("A", "B", "C"))),
    region = "n"
  )
  stream <- immigration_stream(
    mismatch, 1,
    assignment = list(region = migrant_assignment("fixed", fixed_value = "n"))
  )
  error <- expect_hammel_error(
    run_simulation(set_immigration(model, a = stream), 1),
    "groups", "run_simulation"
  )
  expect_identical(error$group, "race")
  unknown_category <- immigration_stream(
    immigration_template(30, race = "A", groups = list(race = race_group())),
    1,
    assignment = list(region = migrant_assignment("fixed", fixed_value = "e"))
  )
  expect_hammel_error(
    run_simulation(set_immigration(model, a = unknown_category), 1),
    "groups", "run_simulation"
  )
})

test_that("observation starts at arrival", {
  template <- immigration_template(30, "f")
  model <- set_immigration(
    quiet_model(death_rate = 0), a = immigration_stream(template, 20)
  )
  result <- run_simulation(model, 4, seed = 5)
  arrived <- arrivals_of(result)
  expect_gt(nrow(arrived), 0)
  # Population counts include migrants only from their arrival.
  sizes <- calculate_population_size(result, times = c(0, 2, 4))
  expected <- 10 + vapply(c(0, 2, 4), function(t) {
    sum(arrived$arrival_time <= t)
  }, numeric(1))
  expect_identical(as.numeric(sizes$population), expected)
  # Exposure begins at arrival, and migrants' own births are not births.
  rates <- calculate_mortality_rates(
    result, time_breaks = c(0, 4), age_breaks = c(0, Inf)
  )
  exposure <- 10 * 4 + sum(4 - arrived$arrival_time)
  expect_equal(sum(rates$person_years), exposure)
  lexis <- person_lexis_cpp(
    result$persons, result$unions, rep(0L, nrow(result$persons)), 1L,
    c(0, 4), c(0, Inf), numeric(1),
    list(typed = FALSE, exclusive = FALSE, counts_cohab = FALSE, axis = 0L,
         utype_filter = 0L, previous_mask = c(TRUE, TRUE, TRUE), sex = FALSE,
         group = FALSE, union_type = FALSE, previous_union = FALSE,
         previous_map = 0:2)
  )
  expect_identical(sum(lexis$unassigned_births), 0)
})

test_that("population change is arrivals plus births minus deaths", {
  model <- make_model(300, birth_rate = 0.2, death_rate = 0.05)
  model$population$persons$sex[1:150] <- "f"
  template <- immigration_template(c(20, 25, 70), c("f", "m", "m"))
  model <- set_immigration(model, a = immigration_stream(template, 25))
  result <- run_simulation(model, 8, seed = 21)
  persons <- result$persons
  living <- sum(is.na(persons$death_time))
  arrivals <- sum(!is.na(persons$arrival_time))
  births <- sum(persons$birth_time > 0 & is.na(persons$arrival_time))
  deaths <- sum(!is.na(persons$death_time))
  expect_gt(arrivals, 0)
  expect_gt(births, 0)
  expect_identical(living, 300L + arrivals + births - deaths)
})

test_that("a constant death rate is recovered from migrants' exposure", {
  # Old migrants die at rate 0.2 from arrival; counting exposure from birth
  # would understate it badly.
  template <- immigration_template(60, "m")
  model <- set_immigration(
    quiet_model(death_rate = 0.2), a = immigration_stream(template, 100)
  )
  result <- run_simulation(model, 5, seed = 8)
  rates <- calculate_mortality_rates(
    result, time_breaks = c(0, 5), age_breaks = c(60, Inf)
  )
  deaths <- sum(rates$deaths)
  expect_gt(deaths, 50)
  expect_lt(abs(deaths / sum(rates$person_years) - 0.2),
            4 * sqrt(deaths) / sum(rates$person_years))
})

test_that("arrival columns appear only with immigration and survive restarts", {
  model <- quiet_model()
  plain <- run_simulation(model, 1, seed = 1)
  expect_false(
    any(c("arrival_time", "birthplace") %in% names(plain$persons))
  )
  migrating <- set_immigration(
    model, a = immigration_stream(immigration_template(30), 20)
  )
  none <- run_simulation(migrating, 0)
  expect_true(all(c("arrival_time", "birthplace") %in% names(none$persons)))
  expect_false("stream" %in% names(none$persons))
  result <- run_simulation(migrating, 1, seed = 1)
  restarted <- run_simulation(set_population(model, result), 1, seed = 2)
  kept <- restarted$persons[seq_len(nrow(result$persons)), ]
  expect_identical(kept$arrival_time, result$persons$arrival_time)
  expect_identical(kept$birthplace, result$persons$birthplace)
  expect_true(restarted$groups$birthplace$automatic)
})

test_that("models without streams are unchanged", {
  model <- make_model(200, birth_rate = 0.1, death_rate = 0.01)
  model$population$persons$sex[1:100] <- "f"
  a <- run_simulation(model, 10, seed = 11)
  removed <- set_immigration(
    set_immigration(model, a = immigration_stream(immigration_template(30), 5))
  )
  b <- run_simulation(removed, 10, seed = 11)
  expect_identical(a, b)
})

test_that("arrival records are validated in populations", {
  result <- run_simulation(
    set_immigration(
      quiet_model(), a = immigration_stream(immigration_template(30), 20)
    ),
    1, seed = 1
  )
  persons <- result$persons
  bad <- persons
  row <- which(!is.na(bad$arrival_time))[1]
  bad$birthplace[row] <- "native"
  expect_hammel_error(
    hammel_population(bad, end_time = 1), "population", "hammel_population"
  )
  bad <- persons
  bad$arrival_time[row] <- bad$birth_time[row] - 1
  expect_hammel_error(
    hammel_population(bad, end_time = 1), "population", "hammel_population"
  )
  expect_hammel_error(
    hammel_population(persons, end_time = 0.5), "population",
    "hammel_population"
  )
})

test_that("immigration summaries describe each stream", {
  template <- immigration_template(c(20, 30))
  process <- immigration(
    a = immigration_stream(template, 5),
    b = immigration_stream(
      template,
      data.frame(time_start = 0, time_end = Inf, intensity = 2),
      weights = weight_model(~age, c(age = 0.1))
    )
  )
  x <- summary(process)
  expect_s3_class(x, "summary.hammel_immigration")
  expect_identical(x$streams$a$eligible, 2L)
  expect_identical(x$streams$a$intensity, 5)
  expect_s3_class(x$streams$b$intensity, "tbl_df")
  text <- cli::cli_fmt(value <- print(x))
  expect_identical(value, x)
  expect_true(any(grepl("Immigration process", text, fixed = TRUE)))
  expect_true(any(grepl("Stream: b", text, fixed = TRUE)))
  expect_true(any(grepl("2 per year in [0, Inf)", text, fixed = TRUE)))
  model <- set_process(quiet_model(), process)
  text <- cli::cli_fmt(print(summary(model)))
  expect_true(any(grepl("Immigration process", text, fixed = TRUE)))
})

test_that("population summaries count immigrants by stream", {
  template <- immigration_template(c(30, 40))
  model <- set_immigration(
    quiet_model(),
    b = immigration_stream(template, 10),
    a = immigration_stream(template, 20)
  )
  result <- run_simulation(model, 2, seed = 3)
  s <- summary(result)
  arrived <- !is.na(result$persons$arrival_time)
  streams <- result$persons$birthplace[arrived]
  expect_identical(s$n_immigrants, sum(arrived))
  expect_identical(
    s$immigrants_by_stream,
    c(a = sum(streams %in% "a"), b = sum(streams %in% "b"))
  )
  expect_match(
    paste(cli::cli_fmt(print(s)), collapse = "\n"),
    sprintf(
      "Total immigrants: %d \\(a %d, b %d\\)", s$n_immigrants,
      s$immigrants_by_stream[["a"]], s$immigrants_by_stream[["b"]]
    )
  )
  # Populations without arrival records report no immigrants.
  plain <- summary(quiet_model()$population)
  expect_true(is.na(plain$n_immigrants))
  expect_null(plain$immigrants_by_stream)
  expect_no_match(
    paste(cli::cli_fmt(print(plain)), collapse = "\n"), "immigrants"
  )
})
