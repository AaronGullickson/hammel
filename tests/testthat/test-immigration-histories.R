# External immigration, phase 2: imported histories (D186).

# A template ending at 0 whose only eligible person drawn is woman 1, aged
# 40: widowed at -5 from man 2 (married -20), with children 3 (born -15,
# alive), 4 (born -12, died -10) by him, and 5 (born -2) by man 6, who was
# never her partner. `extra` adds person columns.
history_template <- function(groups = NULL, ..., transitions = NULL) {
  persons <- data.frame(
    person_id = 1:6,
    sex = c("f", "m", "f", "m", "f", "m"),
    birth_time = c(-40, -45, -15, -12, -2, -42),
    death_time = c(NA, -5, NA, -10, NA, NA),
    mother_id = c(NA, NA, 1, 1, 1, NA),
    father_id = c(NA, NA, 2, 2, 6, NA),
    ...
  )
  unions <- data.frame(
    union_id = 1, wife_id = 1, husband_id = 2, start_time = -20,
    end_time = -5, end_reason = "death_husband"
  )
  hammel_population(
    persons, unions = unions, end_time = 0, groups = groups,
    transitions = transitions
  )
}

# Households without dependents: each template person arrives alone, with
# their history (D188).
alone <- function() {
  household(age_independence = 0)
}

# Draw only template person 1.
only_first <- function() {
  weight_model(fun = function(members, bin) {
    ifelse(members$person_id == 1, 0, -1000)
  })
}

history_stream <- function(template = history_template(), intensity = 5, ...) {
  immigration_stream(
    per_year(intensity), template,
    weights = only_first(), household = alone(), ...
  )
}

first_arrival <- function(result) {
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  arrived[which.min(arrived$arrival_time), ]
}

test_that("an arrival imports unions, partners and children", {
  model <- set_immigration(make_model(10), a = history_stream())
  result <- run_simulation(model, 1, seed = 1)
  p <- first_arrival(result)
  offset <- p$arrival_time
  expect_equal(p$birth_time, offset - 40)
  external <- result$external_relatives
  union <- result$unions[result$unions$wife_id == p$person_id, ]
  expect_identical(nrow(union), 1L)
  expect_true(union$pre_migration)
  expect_equal(c(union$start_time, union$end_time), offset + c(-20, -5))
  expect_identical(union$end_reason, "death_husband")
  husband <- external[external$person_id == union$husband_id, ]
  expect_equal(husband$death_time, offset - 5)
  expect_equal(husband$censor_time, offset - 5)
  expect_true(is.na(husband$mother_id))
  kids <- external[external$mother_id %in% p$person_id, ]
  kids <- kids[order(kids$birth_time), ]
  expect_equal(kids$birth_time, offset + c(-15, -12, -2))
  # The other parent is known only when they were her partner.
  expect_identical(
    kids$father_id, c(union$husband_id, union$husband_id, NA_real_)
  )
  expect_equal(kids$censor_time, c(offset, offset - 10, offset))
  # Every arrival gets its own copies with fresh IDs.
  all_ids <- c(result$persons$person_id, external$person_id)
  expect_false(anyDuplicated(all_ids) > 0)
  expect_false(any(external$person_id %in% 1:6))
  expect_identical(
    nrow(external), 4L * sum(!is.na(result$persons$arrival_time))
  )
  expect_false(any(result$persons$person_id %in% external$person_id))
})

test_that("imported parity selects parity-specific fertility", {
  # Only women at parity 3 or more can give birth.
  rates <- tibble::tibble(
    age_start = 15, age_end = 50, parity_start = c(0, 3),
    parity_end = c(3, Inf), rate = c(0, 2)
  )
  model <- quiet_unreachable(
    make_model(10) |>
      set_fertility(rates, birth_interval = 0) |>
      set_immigration(a = history_stream(intensity = 3))
  )
  result <- quiet_unreachable(run_simulation(model, 3, seed = 2))
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  births <- result$persons[result$persons$mother_id %in% arrived$person_id, ]
  expect_gt(nrow(births), 0)
  # Their exposure is at parity 3 or more from arrival.
  observed <- calculate_fertility_rates(
    result, c(0, 3), c(15, 50), parity_breaks = c(0, 3, Inf)
  )
  expect_identical(sum(observed$births[observed$parity_start == 0]), 0)
  expect_gt(sum(observed$person_years[observed$parity_start == 3]), 0)
  expect_identical(
    sum(observed$person_years[observed$parity_start == 0]), 0
  )
})

test_that("imported last births keep birth spacing", {
  # Her last birth was 2 years before arrival; spacing is 5 years.
  model <- quiet_unreachable(
    make_model(10) |>
      set_fertility(
        data.frame(age_start = 15, age_end = 50, rate = 50),
        birth_interval = 5
      ) |>
      set_immigration(a = history_stream(intensity = 2))
  )
  result <- run_simulation(model, 6, seed = 3)
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  births <- result$persons[result$persons$mother_id %in% arrived$person_id, ]
  expect_gt(nrow(births), 0)
  mother <- arrived[match(births$mother_id, arrived$person_id), ]
  expect_true(all(births$birth_time >= mother$arrival_time + 3 - 1e-9))
})

test_that("imported previous unions select remarriage schedules", {
  # Only widows search.
  model <- make_rate_model(10)
  model$population$persons$sex <- "m"
  model <- set_union_formation(
    model, "rate",
    matching_model(~age_difference, c(age_difference = 0)),
    rates = data.frame(
      age_start = 15, age_end = Inf, rate = c(0, 0, 50),
      previous_union = c("none", "separated", "widowed")
    ),
    cancel_birth_if_no_match = FALSE
  )
  template <- history_template()
  # Person 5 is a never-married woman who could also be drawn.
  weights <- weight_model(fun = function(members, bin) {
    ifelse(members$person_id %in% c(1, 5), 0, -1000)
  })
  model <- set_immigration(
    model,
    a = immigration_stream(
      per_year(5), template,
      weights = weights,
      household = alone()
    )
  )
  result <- run_simulation(model, 4, seed = 4)
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  formed <- result$unions[!result$unions$pre_migration, ]
  wives <- arrived[arrived$person_id %in% formed$wife_id, ]
  widowed <- round(arrived$arrival_time - arrived$birth_time) == 40
  expect_gt(nrow(wives), 0)
  expect_true(all(round(wives$arrival_time - wives$birth_time) == 40))
  expect_true(any(!widowed))
})

test_that("group duration continues from the template or starts at arrival", {
  template <- history_template(
    groups = list(race = group(categories = c("A", "B"))),
    race = c("B", "A", "A", "A", "A", "A"),
    transitions = data.frame(
      person_id = 1, time = -1, group = "race", origin = "A",
      destination = "B"
    )
  )
  destination <- hammel_population(
    data.frame(
      person_id = 1:4, sex = c("f", "m", "f", "m"), birth_time = -30,
      death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
      race = "A", region = "n"
    ),
    end_time = 0
  )
  model <- make_model(4) |>
    set_population(destination) |>
    set_groups(
      race = group(categories = c("A", "B")),
      region = group(categories = c("n", "s"))
    ) |>
    set_immigration(
      a = history_stream(
        template, 20,
        assignment = list(
          region = migrant_assignment("fixed", fixed_value = "s")
        )
      )
    )
  # Moves only after 4 years in a category.
  moves <- function(origin, destination) {
    data.frame(
      origin = origin, destination = destination, dur_start = c(0, 4),
      dur_end = c(4, Inf), rate = c(0, 100)
    )
  }
  model <- model |>
    set_group_transition(moves("B", "A"), group = "race") |>
    set_group_transition(moves("s", "n"), group = "region")
  result <- run_simulation(model, 6, seed = 5)
  arrived <- result$persons[!is.na(result$persons$arrival_time), ]
  tx <- result$transitions
  tx <- tx[tx$person_id %in% arrived$person_id, ]
  # Entry rows record the assigned region at arrival.
  entries <- tx[is.na(tx$origin) & tx$group != "birthplace", ]
  expect_identical(nrow(entries), nrow(arrived))
  expect_true(all(entries$group == "region" & entries$destination == "s"))
  expect_false(any(entries$pre_migration))
  at <- arrived$arrival_time[match(entries$person_id, arrived$person_id)]
  expect_equal(entries$time, at)
  # Imported race moves are pre-migration.
  imported <- tx[tx$pre_migration, ]
  expect_true(all(imported$group == "race" & imported$origin == "A"))
  # Race (entered 1 year before arrival) can change after 3 years; region
  # (entered at arrival) only after 4.
  later <- tx[!tx$pre_migration & !is.na(tx$origin), ]
  start <- arrived$arrival_time[match(later$person_id, arrived$person_id)]
  race <- later$group == "race"
  expect_gt(sum(race), 0)
  expect_gt(sum(!race), 0)
  expect_true(all(later$time[race] - start[race] >= 3 - 1e-9))
  expect_true(all(later$time[!race] - start[!race] >= 4 - 1e-9))
  expect_true(any(later$time[race] - start[race] < 4))
  # Duration-based observation starts region duration at arrival and does not
  # count entries as moves.
  for (by in list(NULL, "race")) {
    rates <- calculate_transition_rates(
      result, c(0, 6), dur_breaks = c(0, 4, Inf), by = by, group = "region"
    )
    early <- rates$origin == "s" & rates$dur_start == 0
    expect_identical(sum(rates$transitions[early]), 0)
    expect_equal(sum(rates$transitions), sum(!race))
    # Migrants spend their first four years in "s" in the first band.
    expect_gt(sum(rates$person_years[early]), 0)
  }
  validate_population(result)
})

test_that("observation never counts external relatives", {
  model <- set_immigration(make_model(10), a = history_stream(intensity = 10))
  result <- run_simulation(model, 2, seed = 6)
  arrivals <- sum(!is.na(result$persons$arrival_time))
  sizes <- calculate_population_size(result, 2)
  expect_identical(sizes$population, 10L + arrivals)
  rates <- calculate_mortality_rates(result, c(-50, 2), c(0, Inf))
  expect_identical(sum(rates$deaths), 0)
  dissolution <- calculate_dissolution_rates(result, c(-50, 2), c(0, Inf))
  expect_identical(sum(dissolution$union_years), 0)
  augmented <- augment_unions(result)
  expect_false(anyNA(augmented$birth_time_h))
  expect_s3_class(extract_unpartnered_births(result), "tbl_df")
})

test_that("kin functions include external relatives only on request", {
  model <- set_immigration(make_model(10), a = history_stream(intensity = 3))
  result <- run_simulation(model, 1, seed = 7)
  p <- first_arrival(result)
  external <- result$external_relatives
  kids <- external$person_id[external$mother_id %in% p$person_id]
  without <- extract_kin(result, p$person_id, 1)
  expect_false(any(without$relative_id %in% kids))
  with <- extract_kin(result, p$person_id, 1, include_external = TRUE)
  expect_setequal(with$relative_id[with$relationship == "child"], kids)
  # Two children are alive and known at her arrival; none afterwards.
  at_arrival <- count_kin(
    result, p$arrival_time, 1, relationships = "child",
    include_external = TRUE
  )
  expect_gte(at_arrival$kin, 2L)
  later <- count_kin(
    result, 1, 1, relationships = "child", include_external = TRUE
  )
  expect_identical(later$kin, 0L)
  expect_identical(
    count_kin(result, 1, 1, relationships = "child")$persons, later$persons
  )
})

test_that("restarts keep imported histories and derived state", {
  rates <- tibble::tibble(
    age_start = 15, age_end = 50, parity_start = c(0, 3),
    parity_end = c(3, Inf), rate = c(0, 2)
  )
  model <- quiet_unreachable(
    make_model(10) |>
      set_fertility(rates, birth_interval = 0) |>
      set_immigration(a = history_stream(intensity = 3))
  )
  first <- quiet_unreachable(run_simulation(model, 1, seed = 8))
  plain <- quiet_unreachable(set_population(model, first) |> set_immigration())
  second <- quiet_unreachable(run_simulation(plain, 3, seed = 9))
  validate_population(second)
  expect_identical(
    second$external_relatives, first$external_relatives
  )
  ids <- c(second$persons$person_id, second$external_relatives$person_id)
  expect_false(anyDuplicated(ids) > 0)
  expect_false(anyDuplicated(second$unions$union_id) > 0)
  # Arrived women keep parity 3+ after the restart and give birth.
  arrived <- first$persons$person_id[!is.na(first$persons$arrival_time)]
  births <- second$persons[
    second$persons$mother_id %in% arrived & second$persons$birth_time > 1,
  ]
  expect_gt(nrow(births), 0)
})

test_that("template external relatives are flattened with their censor times", {
  model <- set_immigration(make_model(10), a = history_stream(intensity = 3))
  first <- run_simulation(model, 1, seed = 10)
  p <- first_arrival(first)
  # The first run's migrant becomes the only person drawn from it.
  weights <- weight_model(fun = function(members, bin) {
    ifelse(members$person_id == p$person_id, 0, -1000)
  })
  again <- set_immigration(
    make_model(10),
    b = immigration_stream(
      per_year(3), first,
      weights = weights,
      household = alone()
    )
  )
  second <- run_simulation(again, 1, seed = 11)
  q <- first_arrival(second)
  offset <- q$arrival_time - first$end_time
  kids <- second$external_relatives[
    second$external_relatives$mother_id %in% q$person_id,
  ]
  kids <- kids[order(kids$birth_time), ]
  original <- first$external_relatives[
    first$external_relatives$mother_id %in% p$person_id,
  ]
  original <- original[order(original$birth_time), ]
  expect_equal(kids$birth_time, original$birth_time + offset)
  expect_equal(kids$censor_time, original$censor_time + offset)
  # The template's arrival records are not imported.
  expect_equal(
    q$arrival_time - q$birth_time, first$end_time - p$birth_time
  )
})

test_that("populations with external relatives are validated", {
  model <- set_immigration(make_model(10), a = history_stream(intensity = 3))
  result <- run_simulation(model, 1, seed = 12)
  rebuild <- function(x) {
    hammel_population(
      x$persons, unions = x$unions, end_time = x$end_time,
      external_relatives = x$external_relatives
    )
  }
  expect_s3_class(rebuild(result), "hammel_population")
  bad <- result
  bad$external_relatives$censor_time[1] <- 99
  expect_hammel_error(rebuild(bad), "population", "hammel_population")
  bad <- result
  bad$external_relatives$person_id[1] <- result$persons$person_id[1]
  expect_hammel_error(rebuild(bad), "population", "hammel_population")
  bad <- result
  bad$external_relatives <- rbind(
    bad$external_relatives,
    transform(bad$external_relatives[1, ], person_id = 9999, mother_id = NA,
              father_id = NA)
  )
  expect_hammel_error(rebuild(bad), "population", "hammel_population")
  bad <- result
  bad$external_relatives$censor_time <- NULL
  expect_hammel_error(rebuild(bad), "table", "hammel_population")
  bad <- result
  bad$unions$pre_migration <- NULL
  expect_hammel_error(rebuild(bad), "population", "hammel_population")
  bad <- result
  husband <- bad$unions$husband_id[bad$unions$pre_migration][1]
  bad$unions$end_time[bad$unions$husband_id == husband] <- NA
  bad$unions$end_reason[bad$unions$husband_id == husband] <- NA
  expect_hammel_error(rebuild(bad), "population", "hammel_population")
})

test_that("cohabitations in templates need cohabitation in the simulation", {
  template <- history_template()
  template$unions$union_type <- "cohab"
  model <- set_immigration(make_model(10), a = history_stream(template))
  expect_hammel_error(
    run_simulation(model, 1), "specification", "run_simulation"
  )
})
