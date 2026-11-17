# Template cuts every `step` years (D196).

# A hand-built template simulated from 0 to 20: woman 1 (born -30) married
# man 2 (born -32) at -5; he dies at 12; she has child 3 at 14 by an unknown
# father. Man 4 (born -20) lives alone throughout.
cut_template <- function() {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "m", "f", "m"),
    birth_time = c(-30, -32, 14, -20), death_time = c(NA, 12, NA, NA),
    mother_id = c(NA, NA, 1, NA), father_id = NA_integer_
  )
  unions <- data.frame(
    union_id = 1, wife_id = 1, husband_id = 2, start_time = -5,
    end_time = 12, end_reason = "death_husband"
  )
  template <- hammel_population(persons, unions = unions, end_time = 20)
  template$unrealized_since <- 0
  template
}

# A destination of 10 men aged 20 that starts at 10.
cut_destination <- function() {
  persons <- data.frame(
    person_id = 1:10, sex = "m", birth_time = -10, death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_
  )
  make_model(10) |>
    set_population(hammel_population(persons, end_time = 10))
}

test_that("a population at a time keeps what had happened by then", {
  template <- cut_template()
  at <- population_at(template, 10)
  expect_identical(at$end_time, 10)
  # The child is not yet born; the husband is still alive and married.
  expect_identical(at$persons$person_id, c(1L, 2L, 4L))
  expect_true(all(is.na(at$persons$death_time)))
  expect_true(is.na(at$unions$end_time))
  expect_true(is.na(at$unions$end_reason))
  expect_no_error(validate_population(at))
  # At 15 the widow has her child; the marriage ended by death.
  later <- population_at(template, 15)
  expect_identical(later$persons$person_id, 1:4)
  expect_identical(later$unions$end_reason, "death_husband")
  expect_identical(population_at(template, 20), template)
})

test_that("a population at a time has its groups and history then", {
  model <- make_model(30, death_rate = 0.05)
  model$population$persons$region <- "n"
  model <- model |>
    set_groups(region = group(categories = c("n", "s"))) |>
    set_group_transition(
      data.frame(
        origin = c("n", "s"), destination = c("s", "n"), age_start = 0,
        age_end = Inf, rate = 0.3
      ),
      group = "region"
    )
  result <- run_simulation(model, 20, seed = 1)
  at <- population_at(result, 10)
  expect_no_error(validate_population(at))
  expect_true(all(at$transitions$time <= 10))
  expect_true(all(is.na(at$persons$death_time) | at$persons$death_time <= 10))
  # Each person's region at 10 is where their history put them then.
  tx <- result$transitions
  for (id in at$persons$person_id) {
    moves <- tx[tx$person_id == id, ]
    before <- moves[moves$time <= 10, ]
    expected <- if (nrow(before)) {
      tail(before$destination[order(before$time)], 1)
    } else if (nrow(moves)) {
      moves$origin[which.min(moves$time)]
    } else {
      result$persons$region[result$persons$person_id == id]
    }
    expect_identical(at$persons$region[at$persons$person_id == id], expected)
  }
})

test_that("arrivals come from the template as it stood at each cut", {
  model <- cut_destination() |>
    set_immigration(
      a = immigration_stream(per_year(40), cut_template(), step = 5)
    )
  result <- run_simulation(model, 10, seed = 2)
  people <- result$persons[!is.na(result$persons$arrival_time), ]
  early <- people[people$arrival_time < 15, ]
  late <- people[people$arrival_time >= 15, ]
  expect_gt(nrow(early), 0)
  expect_gt(nrow(late), 0)
  age <- function(x) round(x$arrival_time - x$birth_time, 6)
  # From the cut at 10: the couple (40 and 42) and the man (30).
  expect_setequal(unique(age(early)), c(40, 42, 30))
  # From the cut at 15: the widow (45) with her child (1), and the man (35).
  expect_setequal(unique(age(late)), c(45, 1, 35))
  # Couples from the first cut arrive married; the widow arrives widowed.
  live <- result$unions[is.na(result$unions$end_time), ]
  wives <- early$person_id[age(early) == 40]
  expect_true(all(wives %in% live$wife_id))
  widows <- late$person_id[age(late) == 45]
  expect_false(any(widows %in% live$wife_id))
  ended <- result$unions[result$unions$wife_id %in% widows, ]
  expect_true(all(ended$end_reason == "death_husband"))
  # Her marriage ended 3 years before the cut, shifted to her arrival.
  arrived <- late$arrival_time[match(ended$wife_id, late$person_id)]
  expect_equal(ended$end_time, arrived - 3)
  expect_no_error(validate_population(result))
})

test_that("cuts must lie within the template's simulated history", {
  # Cuts at 10, 15, 20 and 25; the template ends at 20.
  model <- cut_destination() |>
    set_immigration(
      a = immigration_stream(per_year(5), cut_template(), step = 5)
    )
  expect_hammel_error(run_simulation(model, 20), "coverage", "run_simulation")
  # A template without a simulated history can't be cut.
  plain <- cut_template()
  plain$unrealized_since <- NULL
  model <- cut_destination() |>
    set_immigration(a = immigration_stream(per_year(5), plain, step = 5))
  expect_hammel_error(run_simulation(model, 5), "population", "run_simulation")
  for (bad in list(0, -1, c(1, 2), "1", Inf)) {
    expect_hammel_error(
      immigration_stream(per_year(1), cut_template(), step = bad),
      "argument", "immigration_stream"
    )
  }
  text <- cli::cli_fmt(print(summary(immigration(
    a = immigration_stream(per_year(1), cut_template(), step = 5)
  ))))
  expect_true(any(grepl("Template cuts: every 5 years", text, fixed = TRUE)))
})

# Two unpartnered people aged 30 (id 1) and 40 (id 2) at the end, 20.
age_template <- function() {
  template <- hammel_population(
    data.frame(
      person_id = 1:2, sex = "m", birth_time = c(-10, -20),
      death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
    ),
    end_time = 20
  )
  template$unrealized_since <- 0
  template
}

only <- function(id) {
  weight_model(fun = function(members, bin) {
    ifelse(members$person_id == id, 0, -1000)
  })
}

# Equal weights for everyone.
uniform <- weight_model(fun = function(members, bin) rep(0, nrow(members)))

test_that("weight schedules switch selection at their boundaries", {
  schedule <- tibble::tibble(
    time_start = c(0, 15, 18), time_end = c(15, 18, Inf),
    model = list(only(1), uniform, only(2))
  )
  model <- cut_destination() |>
    set_immigration(
      a = immigration_stream(
        per_year(40), age_template(),
        weights = schedule, household = NULL
      )
    )
  result <- run_simulation(model, 10, seed = 3)
  people <- result$persons[!is.na(result$persons$arrival_time), ]
  age <- round(people$arrival_time - people$birth_time, 6)
  early <- people$arrival_time < 15
  middle <- people$arrival_time >= 15 & people$arrival_time < 18
  late <- people$arrival_time >= 18
  expect_true(all(age[early] == 30))
  expect_setequal(unique(age[middle]), c(30, 40))
  expect_true(all(age[late] == 40))
})

test_that("weight schedules combine with cuts and keep household intensity", {
  seen <- list()
  record <- function(members, bin) {
    seen[[length(seen) + 1L]] <<- bin
    ifelse(members$partnered == 1, 0, -1000)
  }
  singles <- function(members, bin) {
    seen[[length(seen) + 1L]] <<- bin
    ifelse(members$partnered == 1, -1000, 0)
  }
  schedule <- tibble::tibble(
    time_start = c(0, 13), time_end = c(13, Inf),
    model = list(weight_model(fun = record), weight_model(fun = singles))
  )
  model <- cut_destination() |>
    set_immigration(
      a = immigration_stream(
        per_year(60), cut_template(),
        weights = schedule, step = 5
      )
    )
  result <- run_simulation(model, 10, seed = 4)
  # Slices start at the cuts (10, 15) and the boundary at 13.
  bins <- do.call(rbind, lapply(seen, function(b) {
    data.frame(cut = b$cut, time_start = b$time_start)
  }))
  expect_equal(unique(bins[c("cut", "time_start")]),
               data.frame(cut = c(10, 10, 15), time_start = c(0, 13, 13)),
               ignore_attr = TRUE)
  # Couples (2 people) before 13 and singles after: arrivals still count
  # people at 60 a year in both periods.
  arrival <- result$persons$arrival_time[!is.na(result$persons$arrival_time)]
  before <- sum(arrival < 13)
  after <- sum(arrival >= 13)
  expect_lt(abs(before - 180), 4 * sqrt(2 * 180))
  expect_lt(abs(after - 420), 4 * sqrt(420))
})

test_that("weight schedules must cover the run", {
  schedule <- tibble::tibble(
    time_start = 12, time_end = Inf, model = list(only(1))
  )
  model <- cut_destination() |>
    set_immigration(
      a = immigration_stream(
        per_year(5), age_template(),
        weights = schedule, household = NULL
      )
    )
  expect_hammel_error(run_simulation(model, 5), "coverage", "run_simulation")
  household_schedule <- tibble::tibble(
    time_start = 0, time_end = Inf,
    model = list(weight_model(~partnered, c(partnered = 1)))
  )
  expect_hammel_error(
    immigration_stream(
      per_year(5), age_template(),
      weights = household_schedule, household = NULL
    ),
    "specification", "immigration_stream"
  )
})
