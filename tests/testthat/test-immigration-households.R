# External immigration, phase 3: household parties (D187).

# The specification's worked example, ending at 0. Woman 1 (30) married man
# 2 (34) at -10. His first wife 3 died at -13; their son 4 (15) is attached
# to him. Daughter 5 of 1 and 2 died at -7; children 6 and 7 are alive.
# Man 8 (50) lives alone; woman 9 (25) has child 10 by an unknown father;
# man 11 (10) is an orphan whose mother 12 died at -2. Households: {1, 2, 4,
# 6, 7}, {8}, {9, 10} and {11}. `cohab` makes the marriage a conversion of a
# cohabitation that began at -12. `at` moves every date to end at `at`.
household_template <- function(cohab = FALSE, at = 0) {
  persons <- data.frame(
    person_id = 1:12,
    sex = c("f", "m", "f", "m", "f", "f", "m", "m", "f", "f", "m", "f"),
    birth_time = c(-30, -34, -40, -15, -8, -5, -1, -50, -25, -3, -10, -35),
    death_time = c(NA, NA, -13, NA, -7, NA, NA, NA, NA, NA, NA, -2),
    mother_id = c(NA, NA, NA, 3, 1, 1, 1, NA, NA, 9, 12, NA),
    father_id = c(NA, NA, NA, 2, 2, 2, 2, NA, NA, NA, NA, NA)
  )
  unions <- data.frame(
    union_id = 1:2, wife_id = c(3, 1), husband_id = c(2, 2),
    start_time = c(-20, -10), end_time = c(-13, NA),
    end_reason = c("death_wife", NA)
  )
  if (cohab) {
    unions <- rbind(
      unions,
      data.frame(
        union_id = 3, wife_id = 1, husband_id = 2, start_time = -12,
        end_time = -10, end_reason = "conversion"
      )
    )
    unions$union_type <- c("marriage", "marriage", "cohab")
    unions$previous_cohab_id <- c(NA, 3, NA)
  }
  persons$birth_time <- persons$birth_time + at
  persons$death_time <- persons$death_time + at
  unions$start_time <- unions$start_time + at
  unions$end_time <- unions$end_time + at
  hammel_population(persons, unions = unions, end_time = at)
}

# Draw only the household carried by template person `id`.
only_carrier <- function(id = 1) {
  weight_model(fun = function(members, bin) {
    ifelse(members$person_id == id, 0, -1000)
  })
}

household_stream <- function(
  template = household_template(),
  intensity = 5,
  carrier = 1,
  ...
) {
  immigration_stream(
    per_year(intensity), template,
    weights = if (!is.null(carrier)) only_carrier(carrier),
    household = household(),
    ...
  )
}

arrived <- function(result) {
  result$persons[!is.na(result$persons$arrival_time), ]
}

# Members of the earliest party, keyed by template person.
first_party <- function(result) {
  people <- arrived(result)
  people <- people[people$arrival_time == min(people$arrival_time), ]
  people[order(people$person_id), ]
}

test_that("templates divide into households under D140's rule", {
  parties <- migrant_parties(household_template(), household())
  starts <- parties$party_start
  members <- lapply(
    seq_len(length(starts) - 1L),
    function(p) parties$rows[seq(starts[p] + 1L, starts[p + 1L])]
  )
  expect_setequal(
    lapply(members, sort), list(c(1L, 2L, 4L, 6L, 7L), 8L, 9:10, 11L)
  )
  # Every living person belongs to exactly one household.
  expect_setequal(parties$rows, c(1L, 2L, 4L, 6L, 7L, 8L, 9L, 10L, 11L))
  expect_false(anyDuplicated(parties$rows) > 0)
  # The carrier comes first, then her partner; parents precede children.
  family <- members[[which(vapply(members, length, 1L) == 5L)]]
  expect_identical(family[1:2], c(1L, 2L))
  at <- starts[which(vapply(members, length, 1L) == 5L)]
  expect_identical(parties$partner[at + 1L], 1L)
  expect_identical(parties$union_row[at + 1L], 2L)
  son <- at + match(4L, family)
  expect_identical(parties$father[son], 1L)
  expect_identical(parties$mother[son], -1L)
  kids <- at + match(c(6L, 7L), family)
  expect_identical(parties$mother[kids], c(0L, 0L))
  expect_identical(parties$father[kids], c(1L, 1L))
  # With a lower age of independence, the son (15) lives alone.
  older <- migrant_parties(
    household_template(), household(age_independence = 12)
  )
  expect_identical(length(older$party_start) - 1L, 5L)
})

test_that("a dependent whose mother is an external relative joins the father", {
  persons <- data.frame(
    person_id = c(1, 2, 3), sex = c("m", "f", "f"),
    birth_time = c(-40, -10, -8), death_time = NA_real_,
    mother_id = c(NA, 100, 100), father_id = c(NA, 1, NA)
  )
  external <- data.frame(
    person_id = 100, sex = "f", birth_time = -38, death_time = NA_real_,
    mother_id = NA_real_, father_id = NA_real_, censor_time = -5
  )
  template <- hammel_population(
    persons, end_time = 0, external_relatives = external
  )
  parties <- migrant_parties(template, household())
  # Child 2 joins father 1; child 3 has no parent in the template.
  expect_identical(parties$rows, c(1L, 2L, 3L))
  expect_identical(parties$party_start, c(0L, 2L, 3L))
  expect_identical(parties$father, c(-1L, 0L, -1L))
  expect_identical(parties$mother, c(-1L, -1L, -1L))
})

test_that("a household arrives together with its links and relatives", {
  model <- set_immigration(make_model(10), a = household_stream())
  result <- run_simulation(model, 1, seed = 1)
  party <- first_party(result)
  expect_identical(nrow(party), 5L)
  offset <- party$arrival_time[1]
  expect_true(all(party$birthplace == "a"))
  expect_equal(sort(party$birth_time), offset + c(-34, -30, -15, -5, -1))
  wife <- party[party$birth_time == offset - 30, ]
  husband <- party[party$birth_time == offset - 34, ]
  son <- party[party$birth_time == offset - 15, ]
  kids <- party[party$birth_time > offset - 10, ]
  # The marriage arrives live, at its imported duration.
  live <- result$unions[
    result$unions$wife_id == wife$person_id & is.na(result$unions$end_time),
  ]
  expect_identical(nrow(live), 1L)
  expect_identical(live$husband_id, husband$person_id)
  expect_equal(live$start_time, offset - 10)
  expect_true(live$pre_migration)
  # Children in the household link to members; the son's mother is his
  # father's first wife, an external relative created once.
  expect_identical(kids$mother_id, rep(wife$person_id, 2))
  expect_identical(kids$father_id, rep(husband$person_id, 2))
  expect_identical(son$father_id, husband$person_id)
  external <- result$external_relatives
  first_wife <- external[external$person_id == son$mother_id, ]
  expect_identical(nrow(first_wife), 1L)
  expect_equal(first_wife$death_time, offset - 13)
  first <- result$unions[result$unions$wife_id == son$mother_id, ]
  expect_identical(first$husband_id, husband$person_id)
  expect_true(first$pre_migration)
  # The daughter who died is the only other relative.
  mine <- external$person_id %in% c(son$mother_id) |
    external$mother_id %in% wife$person_id
  expect_identical(sum(mine), 2L)
  daughter <- external[external$mother_id %in% wife$person_id, ]
  expect_equal(daughter$death_time, offset - 7)
  expect_identical(daughter$father_id, husband$person_id)
  # Two relatives per household arrival, each with a fresh ID.
  parties <- length(unique(arrived(result)$arrival_time))
  expect_identical(nrow(external), 2L * parties)
  all_ids <- c(result$persons$person_id, external$person_id)
  expect_false(anyDuplicated(all_ids) > 0)
})

test_that("an orphan's parents stay unknown", {
  model <- set_immigration(make_model(10), a = household_stream(carrier = 11))
  result <- run_simulation(model, 1, seed = 2)
  people <- arrived(result)
  expect_gt(nrow(people), 0)
  expect_true(all(is.na(people$mother_id)))
  expect_identical(NROW(result$external_relatives), 0L)
})

test_that("imported state carries over for household members", {
  # Only women at parity 3 or more can give birth; the carrier has 3.
  rates <- tibble::tibble(
    age_start = 15, age_end = 50, parity_start = c(0, 3),
    parity_end = c(3, Inf), rate = c(0, 2)
  )
  model <- quiet_unreachable(
    make_model(10) |>
      set_fertility(rates, birth_interval = 0) |>
      set_immigration(a = household_stream(intensity = 15))
  )
  result <- quiet_unreachable(run_simulation(model, 3, seed = 3))
  people <- arrived(result)
  births <- result$persons[
    result$persons$mother_id %in% people$person_id &
      result$persons$birth_time > 0,
  ]
  expect_gt(nrow(births), 0)
  # Births in the destination have the carrier's husband as father.
  wives <- people[people$birth_time == people$arrival_time - 30, ]
  live <- result$unions[is.na(result$unions$end_time), ]
  expect_identical(
    births$father_id, live$husband_id[match(births$mother_id, live$wife_id)]
  )
  expect_true(all(births$mother_id %in% wives$person_id))
})

test_that("members share their household's assigned categories", {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "m", "f", "m"), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    region = c("n", "n", "s", "s")
  )
  model <- make_model(4) |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_groups(region = group(categories = c("n", "s"), ancestry = TRUE)) |>
    set_immigration(
      a = household_stream(
        intensity = 100,
        carrier = NULL,
        assignment = list(
          region = migrant_assignment("probability", probs = c(n = .5, s = .5))
        )
      )
    )
  result <- run_simulation(model, 2, seed = 4)
  people <- arrived(result)
  by_party <- tapply(people$region, people$arrival_time, function(x) {
    length(unique(x))
  })
  expect_true(all(by_party == 1L))
  expect_true(all(c("n", "s") %in% people$region))
  # Each member has all their ancestry in the assigned category.
  shares <- result$ancestry$region
  shares <- shares[match(people$person_id, shares$person_id), ]
  expect_equal(shares$s, as.numeric(people$region == "s"))
  # Each member enters the assigned group at arrival.
  entry <- result$transitions[
    result$transitions$person_id %in% people$person_id &
      result$transitions$group == "region",
  ]
  expect_identical(nrow(entry), nrow(people))
  expect_true(all(is.na(entry$origin)))
})

test_that("intensity counts people and households arrive in proportion", {
  template <- household_template()
  # Carrier 1's household is three times as likely as the others.
  weights <- weight_model(fun = function(members, bin) {
    ifelse(members$person_id == 1, log(3), 0)
  })
  stream <- immigration_stream(
    per_year(90), template,
    weights = weights, household = household()
  )
  prepared <- prepare_immigration(
    immigration(a = stream), list(), list(), list(), integer(), 0, 20
  )$streams[[1]]
  slice <- prepared$slices[[1]]
  sizes <- diff(slice$party_start)
  probs <- diff(c(0, slice$cumulative))
  expected_size <- sum(probs * sizes)
  expect_equal(expected_size, (3 * 5 + 1 + 2 + 1) / 6)
  # The engine runs the clock at intensity / expected size.
  expect_identical(prepared$intensity, 90)
  expect_equal(slice$expected_size, expected_size)
  model <- set_immigration(make_model(10), a = stream)
  result <- run_simulation(model, 20, seed = 5)
  people <- arrived(result)
  parties <- unique(people[c("arrival_time")])
  n_parties <- nrow(parties)
  # Parties are Poisson with mean 1800 / E[size]; people are compound
  # Poisson with variance (90 / E[size]) * 20 * E[size^2].
  mean_parties <- 1800 / expected_size
  expect_lt(abs(n_parties - mean_parties), 4 * sqrt(mean_parties))
  sd_people <- sqrt(90 / expected_size * 20 * sum(probs * sizes^2))
  expect_lt(abs(nrow(people) - 1800), 4 * sd_people)
  # Household frequencies follow the carriers' weights.
  size <- table(people$arrival_time)
  big <- sum(size == 5)
  expect_lt(abs(big / n_parties - 0.5), 4 * sqrt(0.25 / n_parties))
})

test_that("live unions dissolve at their imported duration", {
  # No dissolution before 12 years; the couple has been married 10.
  model <- make_model(10) |>
    set_union_dissolution(
      data.frame(dur_start = c(0, 12), dur_end = c(12, Inf), rate = c(0, 5))
    ) |>
    set_immigration(a = household_stream(intensity = 15))
  result <- run_simulation(model, 6, seed = 6)
  # The imported marriages, which only dissolution can end.
  live <- result$unions[
    result$unions$pre_migration & result$unions$end_reason %in% c(
      NA, "dissolution"
    ),
  ]
  expect_gt(nrow(live), 0)
  ended <- live[!is.na(live$end_time), ]
  expect_gt(nrow(ended), 0)
  expect_true(all(ended$end_time - ended$start_time >= 12))
  # Rates count them from arrival only, at their imported durations.
  rates <- calculate_dissolution_rates(result, c(0, 6), c(0, 10, 12, Inf))
  expect_lt(rates$union_years[rates$dur_start == 0], 1e-9)
  expect_gt(rates$union_years[rates$dur_start == 10], 0)
  expect_identical(rates$dissolutions[rates$dur_start < 12], c(0, 0))
  expect_identical(
    rates$dissolutions[rates$dur_start == 12], as.numeric(nrow(ended))
  )
  arrival <- result$persons$arrival_time[
    match(live$wife_id, result$persons$person_id)
  ]
  end <- ifelse(is.na(live$end_time), 6, live$end_time)
  expect_equal(sum(rates$union_years), sum(end - arrival))
})

test_that("imported unions are not destination formations", {
  model <- make_model(10) |>
    set_immigration(a = household_stream(intensity = 15))
  result <- run_simulation(model, 4, seed = 7)
  expect_gt(sum(result$unions$pre_migration), 0)
  rates <- calculate_marriage_rates(result, c(0, 4), c(0, 100))
  expect_identical(sum(rates$marriages), 0)
})

test_that("partnered arrivals do not search", {
  model <- make_rate_model(0, search_rate = 50) |>
    set_immigration(a = household_stream(intensity = 15))
  result <- run_simulation(model, 3, seed = 8)
  people <- arrived(result)
  expect_gt(nrow(people), 0)
  # No destination union forms while the imported marriages last.
  expect_true(all(result$unions$pre_migration))
})

test_that("cohabiting histories link to live unions and convert", {
  model <- make_rate_model(0, search_rate = 0) |>
    set_union_formation(
      "rate",
      rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 0.5),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = 0)
    ) |>
    set_immigration(
      a = household_stream(household_template(cohab = TRUE), intensity = 15)
    )
  result <- run_simulation(model, 2, seed = 9)
  unions <- result$unions
  live <- unions[is.na(unions$end_time), ]
  expect_gt(nrow(live), 0)
  cohab <- unions[match(live$previous_cohab_id, unions$union_id), ]
  expect_identical(cohab$union_type, rep("cohab", nrow(live)))
  expect_identical(cohab$wife_id, live$wife_id)
  expect_equal(live$start_time - cohab$start_time, rep(2, nrow(live)))
  expect_true(all(cohab$pre_migration))

  # On another calendar, translated dates must still meet exactly at the
  # conversion.
  model <- set_immigration(
    model,
    a = household_stream(
      household_template(cohab = TRUE, at = 2000),
      intensity = 15
    )
  )
  result <- run_simulation(model, 2, seed = 9)
  unions <- result$unions
  live <- unions[is.na(unions$end_time), ]
  expect_gt(nrow(live), 0)
  cohab <- unions[match(live$previous_cohab_id, unions$union_id), ]
  expect_identical(cohab$end_time, live$start_time)
})

test_that("events at the template's end are not destination events", {
  # A newborn and a union dated at the template's end arrive at the
  # arrival time itself; they are pre-migration.
  persons <- data.frame(
    person_id = 1:3, sex = c("f", "m", "f"), birth_time = c(-30, -32, 0),
    death_time = NA_real_, mother_id = c(NA, NA, 1), father_id = c(NA, NA, 2)
  )
  unions <- data.frame(
    union_id = 1, wife_id = 1, husband_id = 2, start_time = 0,
    end_time = NA_real_, end_reason = NA_character_
  )
  template <- hammel_population(persons, unions = unions, end_time = 0)
  model <- make_rate_model(0, search_rate = 0) |>
    set_immigration(a = immigration_stream(per_year(3), template))
  result <- run_simulation(model, 3, seed = 1)
  expect_gt(sum(!is.na(result$persons$arrival_time)), 0)
  expect_true(all(result$unions$pre_migration))
  births <- calculate_fertility_rates(result, c(0, 3), c(0, 100))
  expect_identical(sum(births$births), 0)
  unions <- calculate_union_rates(result, c(0, 3), c(0, 100))
  expect_identical(sum(unions$unions), 0)
})

test_that("a live cohabitation converts after arrival", {
  template <- household_template(cohab = TRUE)
  template$unions <- template$unions[template$unions$union_id != 2, ]
  template$unions$end_time[template$unions$union_id == 3] <- NA
  template$unions$end_reason[template$unions$union_id == 3] <- NA
  template$unions$previous_cohab_id <- NA_real_
  model <- make_rate_model(0, search_rate = 0) |>
    set_union_formation(
      "rate",
      rates = data.frame(age_start = 15, age_end = Inf, rate = 0),
      cohab_prob = data.frame(age_start = 0, age_end = Inf, prob = 0.5),
      union_convert_rate = data.frame(dur_start = 0, dur_end = Inf, rate = 5)
    ) |>
    set_immigration(a = household_stream(template, intensity = 15))
  result <- run_simulation(model, 3, seed = 10)
  unions <- result$unions
  converted <- unions[unions$end_reason %in% "conversion", ]
  expect_gt(nrow(converted), 0)
  expect_true(all(converted$pre_migration))
  marriages <- unions[unions$previous_cohab_id %in% converted$union_id, ]
  expect_identical(nrow(marriages), nrow(converted))
  expect_false(any(marriages$pre_migration))
  rates <- calculate_cohab_convert_rate(result, c(0, 3), dur_breaks = c(0, Inf))
  expect_identical(rates$conversions, as.numeric(nrow(converted)))
  arrival <- result$persons$arrival_time[
    match(converted$wife_id, result$persons$person_id)
  ]
  expect_true(all(converted$end_time > arrival))
})

test_that("households persist through restarts", {
  model <- make_model(10) |>
    set_immigration(a = household_stream(intensity = 15))
  first <- run_simulation(model, 2, seed = 11)
  second <- run_simulation(set_population(model, first), 2, seed = 12)
  expect_no_error(validate_population(second))
  ids <- c(second$persons$person_id, second$external_relatives$person_id)
  expect_false(anyDuplicated(ids) > 0)
  expect_false(anyDuplicated(second$unions$union_id) > 0)
  # Live imported unions from the first run keep their records.
  old <- first$unions[first$unions$pre_migration, ]
  kept <- second$unions[match(old$union_id, second$unions$union_id), ]
  expect_true(all(kept$pre_migration))
  expect_equal(kept$start_time, old$start_time)
  # Children whose mother is an external relative keep the link.
  people <- second$persons[
    second$persons$person_id %in% first$persons$person_id,
  ]
  before <- first$persons[match(people$person_id, first$persons$person_id), ]
  expect_identical(people$mother_id, before$mother_id)
})

test_that("household rules are validated", {
  expect_hammel_error(
    immigration_stream(per_year(1), household_template(), household = 18),
    "specification", "immigration_stream"
  )
  process <- immigration(a = household_stream(carrier = NULL))
  stream <- summary(process)$streams$a
  expect_identical(stream$eligible, 4L)
  expect_identical(stream$household, 18)
})

# Selection probabilities of a household stream's parties, in carrier order.
party_probabilities <- function(stream) {
  prepared <- prepare_immigration(
    immigration(a = stream), list(), list(), list(), integer(), 0, 1
  )$streams[[1]]
  diff(c(0, prepared$slices[[1]]$cumulative))
}

test_that("household variables describe each carrier's household", {
  seen <- NULL
  weights <- weight_model(fun = function(members, bin) {
    seen <<- members
    rep(0, nrow(members))
  })
  stream <- immigration_stream(
    per_year(10), household_template(),
    weights = weights
  )
  party_probabilities(stream)
  # Carriers 1 (married to 2, three dependents), 8, 9 (one child) and 11.
  expect_identical(seen$person_id, c(1L, 8L, 9L, 11L))
  expect_identical(seen$partnered, c(1L, 0L, 0L, 0L))
  expect_identical(seen$household_size, c(5L, 1L, 2L, 1L))
  expect_identical(seen$n_dependents, c(3L, 0L, 1L, 0L))
  expect_equal(seen$partner_age, c(34, NA, NA, NA))
})

test_that("formula weights use household variables", {
  template <- household_template()
  stream <- function(formula, coefficients) {
    immigration_stream(
      per_year(10), template,
      weights = weight_model(formula, coefficients)
    )
  }
  p <- party_probabilities(
    stream(
      ~ partnered + n_dependents,
      c(partnered = log(2), n_dependents = log(3))
    )
  )
  w <- c(2 * 3^3, 1, 3, 1)
  expect_equal(p, w / sum(w))
  p <- party_probabilities(stream(~household_size, c(household_size = log(2))))
  w <- 2^c(5, 1, 2, 1)
  expect_equal(p, w / sum(w))
  # The partner's age counts only for couples.
  p <- party_probabilities(stream(~partner_age, c(partner_age = log(1.1))))
  w <- c(1.1^34, 1, 1, 1)
  expect_equal(p, w / sum(w))
})

test_that("custom assignment sees household variables", {
  persons <- data.frame(
    person_id = 1:2, sex = c("f", "m"), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    region = c("n", "s")
  )
  # Couples go north and everyone else south.
  rule <- migrant_assignment("custom", fun = function(members) {
    ifelse(members$partnered == 1, "n", "s")
  })
  model <- make_model(2) |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_groups(region = group(categories = c("n", "s"))) |>
    set_immigration(
      a = immigration_stream(
        per_year(60), household_template(),
        assignment = list(region = rule)
      )
    )
  result <- run_simulation(model, 2, seed = 13)
  people <- arrived(result)
  size <- table(people$arrival_time)[as.character(people$arrival_time)]
  expect_true(all(people$region[size == 5] == "n"))
  expect_true(all(people$region[size < 5] == "s"))
  expect_true(any(size == 5) && any(size < 5))
})

test_that("cloning streams can't weight by household variables", {
  template <- generate_uniform_population(20)
  expect_hammel_error(
    immigration_stream(
      per_year(1), template,
      weights = weight_model(~partnered, c(partnered = 1)),
      household = NULL
    ),
    "specification", "immigration_stream"
  )
})

test_that("population summaries leave out immigrants' births and unions", {
  # Imported marriages dissolve only after 2 years in the destination.
  model <- make_model(10) |>
    set_union_dissolution(
      data.frame(dur_start = c(0, 12), dur_end = c(12, Inf), rate = c(0, 5))
    ) |>
    set_immigration(a = household_stream(intensity = 15))
  result <- run_simulation(model, 6, seed = 14)
  s <- summary(result)
  people <- result$persons
  expect_identical(s$n_immigrants, sum(!is.na(people$arrival_time)))
  expect_identical(s$n_births, sum(is.na(people$arrival_time)))
  unions <- result$unions
  expect_gt(sum(unions$pre_migration), 0)
  expect_identical(s$n_unions, sum(!unions$pre_migration))
  # Dissolutions after arrival count; none happened before it.
  dissolved <- unions$end_reason %in% "dissolution"
  expect_gt(sum(dissolved), 0)
  expect_identical(s$n_dissolutions, sum(dissolved))
})

# Ten unpartnered adults, built without random numbers.
clone_template <- function() {
  hammel_population(
    data.frame(
      person_id = 1:10, sex = rep(c("f", "m"), 5), birth_time = -(21:30),
      death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
    ),
    end_time = 0
  )
}

# A destination declaring `generation`, which templates lack.
generation_model <- function(rule, template = household_template(), ...) {
  persons <- data.frame(
    person_id = 1:2, sex = c("f", "m"), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    generation = "3+"
  )
  make_model(2) |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_groups(generation = group(categories = c("1", "1.5", "2", "3+"))) |>
    set_immigration(
      a = immigration_stream(
        per_year(60), template,
        assignment = list(generation = rule), ...
      )
    )
}

test_that("custom assignment gives each arriving member a category", {
  rule <- migrant_assignment("custom", fun = function(members) {
    ifelse(members$age < 13, "1.5", "1")
  })
  result <- run_simulation(generation_model(rule), 2, seed = 15)
  people <- arrived(result)
  expect_gt(nrow(people), 0)
  age <- people$arrival_time - people$birth_time
  expect_identical(people$generation, ifelse(age < 13, "1.5", "1"))
  # The worked example's household brings both.
  size <- table(people$arrival_time)[as.character(people$arrival_time)]
  expect_setequal(people$generation[size == 5], c("1", "1.5"))
})

test_that("custom assignment sees each member and its household", {
  seen <- list()
  rule <- migrant_assignment("custom", fun = function(members) {
    seen[[length(seen) + 1L]] <<- members
    rep("1", nrow(members))
  })
  result <- run_simulation(generation_model(rule), 1, seed = 16)
  family <- Filter(function(x) nrow(x) == 5L, seen)[[1]]
  expect_identical(family$person_id, c(1L, 2L, 4L, 6L, 7L))
  expect_identical(
    family$role, c("carrier", "partner", rep("dependent", 3))
  )
  expect_identical(family$household_size, rep(5L, 5))
  expect_identical(family$partner_age, rep(34, 5))
  expect_length(unique(family$household), 1L)
  households <- vapply(seen, function(x) x$household[1], numeric(1))
  expect_false(anyDuplicated(households) > 0)
  # Called at arrival, with the arrival time.
  times <- vapply(seen, function(x) x$arrival_time[1], numeric(1))
  expect_setequal(times, unique(arrived(result)$arrival_time))
  # Cloning streams describe people alone.
  clone <- NULL
  rule <- migrant_assignment("custom", fun = function(members) {
    clone <<- members
    "1"
  })
  run_simulation(
    generation_model(rule, clone_template(), household = NULL),
    1, seed = 17
  )
  expect_identical(nrow(clone), 1L)
  expect_false(any(c("role", "household") %in% names(clone)))
})

test_that("custom assignment can sample once per household", {
  rule <- migrant_assignment("custom", fun = function(members) {
    rep(sample(c("1", "1.5"), 1), nrow(members))
  })
  model <- generation_model(rule)
  result <- run_simulation(model, 3, seed = 18)
  people <- arrived(result)
  by_party <- tapply(people$generation, people$arrival_time, function(x) {
    length(unique(x))
  })
  expect_true(all(by_party == 1L))
  expect_setequal(people$generation, c("1", "1.5"))
  # Sampling inside the function follows the seed.
  expect_identical(run_simulation(model, 3, seed = 18), result)
})

test_that("custom assignment results are validated", {
  short <- migrant_assignment("custom", fun = function(members) "1")
  expect_hammel_error(
    run_simulation(generation_model(short), 2, seed = 19),
    "runtime", "run_simulation"
  )
  unknown <- migrant_assignment("custom", fun = function(members) {
    rep("9", nrow(members))
  })
  expect_hammel_error(
    run_simulation(generation_model(unknown), 2, seed = 19),
    "runtime", "run_simulation"
  )
})
