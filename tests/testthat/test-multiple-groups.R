# Populations with two group dimensions, race and region.
two_dimension_model <- function(n = 400, birth_rate = 0, union = FALSE) {
  persons <- data.frame(person_id = seq_len(n), sex = rep(c("f", "m"), length.out = n),
    birth_time = -30, death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    race = rep(c("A", "A", "B", "B"), length.out = n),
    region = rep(c("N", "N", "N", "N", "S", "S", "S", "S"), length.out = n))
  unions <- if (union) data.frame(union_id = seq_len(n / 2), wife_id = seq(1, n, 2),
    husband_id = seq(2, n, 2), start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  hammel_model() |>
    set_population(hammel_population(persons, unions = unions, end_time = 0,
                                     groups = c("race", "region"))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = birth_rate), birth_interval = 0) |>
    set_union_formation("on_demand", kinship_generations = 0)
}

test_that("rates stratify by any combination of declared dimensions", {
  model <- two_dimension_model()
  persons <- model$population$persons
  # Only race A in region S dies.
  rates <- stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1e8),
    list(race = c(A = 1, B = 0), region = c(N = 0, S = 1)))
  out <- simulate(set_mortality(model, rates), .1, seed = 1)
  died <- !is.na(out$persons$death_time)
  expect_identical(died, persons$race == "A" & persons$region == "S")
  # Rates by one dimension pool the other.
  by_region <- stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1e8),
    list(region = c(N = 1, S = 0)))
  out <- simulate(set_mortality(model, by_region), .1, seed = 1)
  expect_identical(!is.na(out$persons$death_time), persons$region == "N")
  # Every combination must be covered, and categories must be known.
  partial <- rates[!(rates$race == "B" & rates$region == "N"), ]
  expect_error(simulate(set_mortality(model, partial), 0), "missing schedules.*race=B, region=N")
  declared <- set_groups(model, race = group(categories = c("A", "B")), region = group())
  expect_error(simulate(set_mortality(declared, rbind(rates, transform(rates[1:2, ], race = "C"))), 0),
               "not categories of group `race`: C")
  # A rate table may use a dimension's column only once declared.
  undeclared <- set_groups(model, race = group())
  expect_error(simulate(set_mortality(undeclared, rates), 0), "`region`.*not a declared group")
})

test_that("union rates stratify by either partner's categories in any dimension", {
  model <- two_dimension_model(union = TRUE)
  persons <- model$population$persons
  # Persons alternate f/m, so couples are (1, 2), (3, 4), ...
  wives <- persons[seq(1, nrow(persons), 2), ]
  husbands <- persons[seq(2, nrow(persons), 2), ]
  rates <- stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = 1e8),
    list(race_w = c(A = 1, B = 0), region_h = c(N = 0, S = 1)))
  out <- simulate(set_union_dissolution(model, rates), .1, seed = 2)
  expect_identical(out$unions$end_reason %in% "dissolution",
                   wives$race == "A" & husbands$region == "S")
  expect_error(simulate(set_union_dissolution(model, transform(rates, race = "A")), 0), "_w`")
})

test_that("newborns get a category in every dimension from its own rule", {
  model <- two_dimension_model(400, birth_rate = 2)
  model <- set_groups(model,
    race = group(categories = c("A", "B", "AB"), ancestry = TRUE,
                 inheritance = create_group_inheritance("custom", fun = function(mother_group, father_group, ...)
                   if (is.na(father_group) || mother_group == father_group) mother_group else "AB")),
    region = group(inheritance = "father"))
  out <- simulate(model, 1, seed = 3)
  kids <- out$persons[out$persons$birth_time >= 0, ]
  expect_gt(nrow(kids), 50)
  mothers <- match(kids$mother_id, out$persons$person_id)
  fathers <- match(kids$father_id, out$persons$person_id)
  expect_identical(kids$region, out$persons$region[fathers])
  mixed <- out$persons$race[mothers] != out$persons$race[fathers]
  expect_true(all(kids$race[mixed] == "AB"))
  expect_true(all(kids$race[!mixed] == out$persons$race[mothers][!mixed]))
  # Ancestry is tracked for race only, named by category.
  expect_named(out$ancestry, "race")
  expect_named(out$ancestry$race, c("person_id", "A", "B", "unknown"))
  kid_rows <- match(kids$person_id, out$ancestry$race$person_id)
  expect_equal(out$ancestry$race$A[kid_rows][mixed], rep(0.5, sum(mixed)))
  # Restarts continue from the carried declaration and ancestry.
  again <- simulate(set_population(model, out), .5, seed = 4)
  expect_identical(again$ancestry$race[seq_len(nrow(out$persons)), ], out$ancestry$race)
})

test_that("unrealized births report every dimension", {
  # Spacing longer than the run: after one birth each, every cell lacks eligible women.
  model <- two_dimension_model(400, birth_rate = 3)
  model <- set_fertility(model, data.frame(age_start = 15, age_end = 50, rate = 3), birth_interval = 5)
  out <- simulate(model, 1, seed = 5)
  u <- out$simulation$unrealized_births
  expect_true(all(c("race", "region") %in% names(u)))
  expect_gt(sum(u$unrealized_births), 0)
  expect_setequal(unique(u$race), c("A", "B"))
  # Split by both dimensions, each record lands in its own row.
  rates <- calculate_fertility_rates(out, c(0, 1), c(15, 50), by = c("region", "race"), unrealized = TRUE)
  rates <- rates[rates$age_start == 15, ]
  for (k in seq_len(nrow(rates))) {
    mine <- u$race == rates$race[k] & u$region == rates$region[k]
    expect_equal(rates$unrealized_births[k], sum(u$unrealized_births[mine]))
  }
  expect_equal(sum(rates$unrealized_births), sum(u$unrealized_births))
})

# Stage 5: matching effects in one named dimension while several are declared.
matching_pair_model <- function(matching, father_matching = NULL) {
  persons <- data.frame(person_id = 1:600, sex = rep(c("f", "m"), 300), birth_time = -25,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    race = rep(c("A", "A", "B", "B"), 150), region = rep(rep(c("N", "S"), each = 6), 50))
  hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = c("race", "region"))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .5)) |>
    set_union_formation("on_demand", matching = matching, father_matching = father_matching,
                        kinship_generations = 0)
}

pair_categories <- function(out, d) {
  persons <- out$persons
  u <- out$unions
  data.frame(w = persons[[d]][match(u$wife_id, persons$person_id)],
             h = persons[[d]][match(u$husband_id, persons$person_id)])
}

test_that("matching group effects use one named dimension among several", {
  # Restrict region only: every union joins one region, races mix.
  out <- simulate(matching_pair_model(create_age_gap_matching_model(0, 5,
    prevent_exogamy = list(region = TRUE))), 2, seed = 41)
  expect_gt(nrow(out$unions), 100)
  region <- pair_categories(out, "region"); race <- pair_categories(out, "race")
  expect_true(all(region$w == region$h))
  expect_gt(mean(race$w != race$h), .2)
  # A strong race pair effect, given by name, separates races instead.
  out <- simulate(matching_pair_model(create_age_gap_matching_model(0, 5,
    group_exogamy = list(race = c("A:B" = -30)))), 2, seed = 41)
  race <- pair_categories(out, "race"); region <- pair_categories(out, "region")
  expect_true(all(race$w == race$h))
  expect_gt(mean(region$w != region$h), .2)
  # Formula predictors name the dimension; the first category is the reference.
  same <- create_matching_model(~ age_difference + same_region, c(age_difference = 0, same_region = 30))
  out <- simulate(matching_pair_model(same), 2, seed = 41)
  region <- pair_categories(out, "region")
  expect_true(all(region$w == region$h))
  female <- create_matching_model(~ female_race:age_difference + same_race,
    c("female_raceA:age_difference" = .1, "female_raceB:age_difference" = .2, same_race = 1))
  expect_s3_class(simulate(matching_pair_model(female), .1, seed = 41), "hammel_population")
  wrong <- create_matching_model(~ female_race:age_difference, c("female_raceA:age_difference" = 1))
  expect_error(simulate(matching_pair_model(wrong), 0), "coefficient names.*female_raceB:age_difference")
})

test_that("matching settings must name one declared dimension", {
  bare <- create_age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1))
  expect_error(simulate(matching_pair_model(bare), 0),
               "Name the dimension.*group_exogamy = list\\(race")
  both <- create_matching_model(~ same_race, c(same_race = 1),
                                prevent_exogamy = list(region = TRUE))
  expect_error(simulate(matching_pair_model(both), 0), "one group dimension for now.*`race`, `region`")
  race <- create_matching_model(~ same_race, c(same_race = 1))
  region <- create_matching_model(~ same_region, c(same_region = 1))
  expect_error(simulate(matching_pair_model(race, father_matching = region), 0), "one group dimension")
  expect_error(simulate(matching_pair_model(create_matching_model(~ same_education,
    c(same_education = 1))), 0), "`education`, which is not a declared group dimension")
  expect_error(create_matching_model(~ age_difference, c(age_difference = 0),
    group_exogamy = list(c("A:B" = 1))), "named by group dimension")
  expect_error(create_matching_model(~ age_difference, c(age_difference = 0),
    delta = list(race = c(A = NA))), "delta")
  expect_error(create_matching_model(~ partner_race, c(partner_race = 1)), "one-sided matching formula")
})

test_that("a single dimension may have any name, including in analysis", {
  persons <- data.frame(person_id = 1:200, sex = rep(c("f", "m"), 100), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    region = rep(c("N", "N", "S", "S"), 50))
  model <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = "region")) |>
    set_mortality(stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .5),
                                 list(region = c(N = 1, S = 0)))) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .2)) |>
    set_union_formation("on_demand", kinship_generations = 0)
  out <- simulate(model, 2, seed = 7)
  rates <- calculate_mortality_rates(out, c(0, 2), c(0, Inf), by = "region")
  expect_true("region" %in% names(rates))
  expect_equal(rates$deaths[rates$region == "S"], c(0, 0))
  expect_gt(sum(rates$deaths[rates$region == "N"]), 0)
  expect_error(calculate_mortality_rates(out, c(0, 2), c(0, Inf), by = "group"), "`region`")
  sizes <- calculate_population_size(out, 2, by = "region")
  expect_named(sizes, c("time", "region", "population"))
  # Also when the dimension tracks ancestry.
  tracked <- simulate(set_groups(model, region = group(ancestry = TRUE)), 1, seed = 7)
  expect_named(calculate_population_size(tracked, 1, by = "region"), c("time", "region", "population"))
  expect_true(all(c("region_w", "region_h") %in% names(get_detailed_unions(tracked))))
})

test_that("remaining single-dimension features and ancestry inputs are checked", {
  model <- two_dimension_model()
  w <- create_weight_model(~ group, c(groupB = 1))
  expect_error(simulate(set_mortality(model, model$processes$mortality$rates, weights = w), 0),
               "not available for this process: group")
  # Analysis without group splits works with several dimensions.
  expect_s3_class(calculate_mortality_rates(simulate(model, 1, seed = 8), c(0, 1), c(0, Inf)), "tbl_df")
  # A bare ancestry table needs exactly one tracked dimension.
  persons <- model$population$persons
  shares <- data.frame(person_id = persons$person_id, A = 1, unknown = 0)
  expect_error(hammel_population(persons, end_time = 0, ancestry = shares,
    groups = list(race = group(ancestry = TRUE), region = group(ancestry = TRUE))), "named list")
})

test_that("unrealized births follow a woman's current group after a transition", {
  # Fertility pools groups; women in spacing move from A to B, then their
  # cell's unrealized births are reported under B.
  persons <- data.frame(person_id = 1:20, sex = "f", birth_time = -20, death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_, group = "A")
  model <- hammel_model() |>
    set_population(hammel_population(persons, end_time = 0, groups = "group")) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 20), birth_interval = 5) |>
    set_union_formation("on_demand", kinship_generations = 0, cancel_birth_if_no_match = FALSE) |>
    set_group_transition(data.frame(origin = "A", destination = "B", age_start = 0, age_end = Inf,
                                    rate = 50))
  withr::local_options(hammel.check_cells = TRUE)
  out <- simulate(model, .5, seed = 9)
  expect_true(all(out$persons$group[1:20] == "B"))
  u <- out$simulation$unrealized_births
  expect_gt(sum(u$unrealized_births[u$group == "B"]), sum(u$unrealized_births[u$group == "A"]))
})

test_that("detailed unions add each dimension's partner categories", {
  out <- simulate(two_dimension_model(union = TRUE), .5, seed = 10)
  u <- get_detailed_unions(out)
  expect_true(all(c("race_w", "race_h", "region_w", "region_h") %in% names(u)))
  expect_identical(u$race_h, out$persons$race[match(u$husband_id, out$persons$person_id)])
  # One dimension of any name, replayed through transitions to formation.
  p <- create_baseline_population(200, groups = list(race = c(A = .5, B = .5)))
  m <- hammel_model() |> set_population(p) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .2)) |>
    set_union_formation("on_demand", kinship_generations = 0) |>
    set_group_transition(data.frame(origin = "A", destination = "B", age_start = 0,
                                    age_end = Inf, rate = .5))
  out <- simulate(m, 3, seed = 11)
  u <- get_detailed_unions(out)
  expect_false(any(c("group_w", "group_h") %in% names(u)))
  moved <- out$transitions$person_id[out$transitions$time > 0]
  late <- u[u$wife_id %in% moved & u$start_time > 0, ]
  expect_true(all(late$race_w %in% c("A", "B")))
})

# Transitions in each of two dimensions (stage 3).
transition_pair_model <- function(n = 400) {
  model <- two_dimension_model(n)
  set_mortality(model, data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0))
}

test_that("each dimension has its own transition process, history and duration", {
  model <- transition_pair_model()
  race <- data.frame(origin = c("A", "B"), destination = c("B", "A"), age_start = 0, age_end = Inf,
                     rate = c(.3, .3))
  region <- data.frame(origin = c("N", "S"), destination = c("S", "N"), dur_start = c(0, 1),
                       dur_end = c(Inf, Inf), rate = c(.4, 1e8))
  expect_error(set_group_transition(model, race), "Name the dimension.*\"race\", \"region\"")
  model <- set_group_transition(model, race, group = "race")
  model <- set_group_transition(model, region, group = "region")
  expect_named(model$processes$group_transition, c("race", "region"))
  # Replacing one dimension's process keeps the other's.
  replaced <- set_group_transition(model, transform(race, rate = 0), group = "race")
  expect_identical(replaced$processes$group_transition$region, model$processes$group_transition$region)
  out <- simulate(model, 4, seed = 21)
  tx <- out$transitions
  expect_named(tx, c("person_id", "time", "group", "origin", "destination"))
  expect_setequal(unique(tx$group), c("race", "region"))
  expect_true(all(tx$origin[tx$group == "race"] %in% c("A", "B")))
  expect_true(all(tx$origin[tx$group == "region"] %in% c("N", "S")))
  # Duration counts from entering the current region only: a return S -> N
  # happens one year after the move to S, whatever happens to race. People
  # who start in S have no move to S, so their duration is their age (30).
  back <- tx[tx$group == "region" & tx$origin == "S", ]
  arrived <- vapply(seq_len(nrow(back)), function(k) {
    h <- tx$time[tx$person_id == back$person_id[k] & tx$group == "region" & tx$time < back$time[k]]
    if (length(h)) max(h) else NA_real_
  }, 0)
  expect_gt(sum(!is.na(arrived)), 20)
  expect_equal(back$time[!is.na(arrived)], arrived[!is.na(arrived)] + 1, tolerance = 1e-6)
  expect_true(all(back$time[is.na(arrived)] < 1e-6))
  expect_true(any(tx$group == "race" & tx$person_id %in% back$person_id[!is.na(arrived)]))
  # Each person's history chains to their current category in each dimension.
  final <- out$persons
  for (d in c("race", "region")) {
    last <- tx[tx$group == d, ]
    last <- last[order(last$person_id, last$time), ]
    last <- last[!duplicated(last$person_id, fromLast = TRUE), ]
    expect_identical(last$destination, final[[d]][match(last$person_id, final$person_id)])
  }
  # Restarts validate and continue every dimension's history.
  again <- simulate(set_population(model, out), 1, seed = 22)
  expect_identical(again$transitions[seq_len(nrow(tx)), ], tx)
  broken <- out
  broken$transitions$destination[which(tx$group == "region")[1]] <- "A"
  expect_error(set_population(model, broken), "chain")
  # Transition rates are counted per dimension.
  expect_error(calculate_transition_rates(out, c(0, 4), age_breaks = c(0, Inf)), "Name the dimension")
  counts <- calculate_transition_rates(out, c(0, 4), dur_breaks = c(0, 1, Inf), group = "region")
  expect_named(counts, c("time_start", "time_end", "dur_start", "dur_end", "origin", "destination",
                         "person_years", "transitions", "rate_transition"))
  expect_equal(sum(counts$transitions), sum(tx$group == "region"))
  race_counts <- calculate_transition_rates(out, c(0, 4), age_breaks = c(0, Inf), group = "race")
  expect_equal(sum(race_counts$transitions), sum(tx$group == "race"))
  expect_equal(race_counts$rate_transition, c(.3, .3), tolerance = .1)
})

test_that("transition rates can be stratified by another dimension", {
  model <- transition_pair_model()
  # Everyone starts in N; race A -> B happens only in region S.
  model$population$persons$region <- "N"
  model <- set_groups(model, race = group(), region = group(categories = c("N", "S")))
  model <- set_group_transition(model, data.frame(origin = "N", destination = "S",
    age_start = 0, age_end = Inf, rate = .5), group = "region")
  model <- set_group_transition(model, data.frame(origin = "A", destination = "B", region = "S",
    age_start = 0, age_end = Inf, rate = 2), group = "race")
  out <- simulate(model, 3, seed = 23)
  tx <- out$transitions
  moves <- tx[tx$group == "race", ]
  expect_gt(nrow(moves), 20)
  arrived <- tx$time[tx$group == "region"][match(moves$person_id, tx$person_id[tx$group == "region"])]
  expect_true(all(!is.na(arrived) & arrived <= moves$time))
  # The transition's own dimension is not a stratum; other names must be declared.
  bad <- set_group_transition(model, data.frame(origin = "A", destination = "B", race = "A",
    age_start = 0, age_end = Inf, rate = 2), group = "race")
  expect_error(simulate(bad, 0), "cannot be stratified by `race` itself")
  undeclared <- set_group_transition(model, data.frame(origin = "A", destination = "B",
    age_start = 0, age_end = Inf, rate = 2), group = "education")
  expect_error(simulate(undeclared, 0), "`education`, which is not a declared group dimension")
  expect_error(group_transition(data.frame(group_origin = "A", group_dest = "B", age_start = 0,
    age_end = Inf, rate = 1)), "`origin` and `destination`")
})

test_that("transition weights and custom weight functions are per dimension", {
  model <- transition_pair_model()
  seen <- new.env()
  fun <- function(members, bin) {
    seen$bins <- c(seen$bins, paste(bin$origin, bin$destination))
    seen$dest <- c(seen$dest, unique(members$destination))
    ifelse(members$sex == "f", 2, 0)
  }
  model <- set_group_transition(model, data.frame(origin = "N", destination = "S",
    age_start = 0, age_end = Inf, rate = .5), group = "region",
    weights = create_weight_model(fun = fun))
  model <- set_group_transition(model, data.frame(origin = "A", destination = "B",
    age_start = 0, age_end = Inf, rate = .5), group = "race",
    weights = create_weight_model(~ sex, c(sexm = 3)))
  out <- simulate(model, 1, seed = 24)
  expect_identical(unique(seen$bins), "N S")
  expect_identical(unique(seen$dest), "S")
  persons <- out$persons
  tx <- out$transitions
  sex <- persons$sex[match(tx$person_id, persons$person_id)]
  # Region moves favor women; race moves favor men.
  expect_gt(mean(sex[tx$group == "region"] == "f"), .7)
  expect_gt(mean(sex[tx$group == "race"] == "m"), .7)
})

test_that("each dimension's inheritance can follow its own calendar schedule", {
  model <- two_dimension_model(400, birth_rate = 2)
  schedule <- tibble::tibble(time_start = c(-100, .5), time_end = c(.5, Inf),
    model = list(create_group_inheritance("mother"), create_group_inheritance("father")))
  model <- set_groups(model, race = group(inheritance = "father"), region = group(inheritance = schedule))
  out <- simulate(model, 1, seed = 25)
  kids <- out$persons[out$persons$birth_time >= 0, ]
  mothers <- match(kids$mother_id, out$persons$person_id)
  fathers <- match(kids$father_id, out$persons$person_id)
  early <- kids$birth_time < .5
  expect_gt(sum(early), 20)
  expect_gt(sum(!early), 20)
  expect_identical(kids$region[early], out$persons$region[mothers][early])
  expect_identical(kids$region[!early], out$persons$region[fathers][!early])
  expect_identical(kids$race, out$persons$race[fathers])
  gap <- tibble::tibble(time_start = .5, time_end = Inf, model = list(create_group_inheritance("mother")))
  expect_error(simulate(set_groups(model, race = group(), region = group(inheritance = gap)), 1),
               "`region` inheritance model schedule has a coverage gap")
})

test_that("detailed unions replay transitions in every dimension", {
  model <- set_group_transition(transition_pair_model(),
    data.frame(origin = "N", destination = "S", age_start = 0, age_end = Inf, rate = 1), group = "region")
  model <- set_group_transition(model,
    data.frame(origin = "A", destination = "B", age_start = 0, age_end = Inf, rate = 1), group = "race")
  model <- set_union_formation(model, "rate", rates = data.frame(age_start = 15, age_end = 50, rate = .5))
  out <- suppressWarnings(simulate(model, 2, seed = 26))
  u <- get_detailed_unions(out)
  expect_gt(nrow(u), 20)
  tx <- out$transitions
  # A wife's category at formation: her latest move before it, or her
  # origin before her first move.
  at <- function(id, time, d) {
    h <- tx[tx$person_id == id & tx$group == d, ]
    h <- h[order(h$time), ]
    if (!nrow(h)) return(out$persons[[d]][out$persons$person_id == id])
    if (any(h$time <= time)) tail(h$destination[h$time <= time], 1) else h$origin[1]
  }
  for (d in c("race", "region"))
    expect_identical(u[[paste0(d, "_w")]], mapply(at, u$wife_id, u$start_time, d, USE.NAMES = FALSE))
})

test_that("transition settings and histories are checked against the declaration", {
  model <- transition_pair_model()
  moves <- data.frame(origin = "A", destination = "B", age_start = 0, age_end = Inf, rate = 1)
  expect_error(group_transition(moves, group = c("race", "region")), "`group` must be NULL")
  # A process stored before the declaration resolves at simulation.
  early <- set_group_transition(hammel_model(), moves)
  expect_named(early$processes$group_transition, ".default")
  early$population <- model$population
  early$processes[c("mortality", "fertility", "union_formation")] <-
    model$processes[c("mortality", "fertility", "union_formation")]
  expect_error(simulate(early, 0), "without `group` needs exactly one")
  early$processes$group_transition$race <- group_transition(moves, group = "race")
  early$processes$group_transition$.default$group <- "race"
  expect_error(simulate(early, 0), "more than one transition process")
  # Histories.
  persons <- model$population$persons
  tx <- data.frame(person_id = 1, time = -1, origin = "B", destination = "A")
  expect_error(hammel_population(persons, end_time = 0, transitions = tx, groups = c("race", "region")),
               "`group` column")
  expect_no_error(hammel_population(persons, end_time = 0, transitions = tx, groups = "race"))
  expect_error(hammel_population(persons, end_time = 0, transitions = transform(tx, group = "education"),
                                 groups = c("race", "region")), "not declared: education")
  expect_error(hammel_population(persons, end_time = 0,
    transitions = data.frame(person_id = 1, time = -1, group_origin = "B", group_dest = "A"),
    groups = "race"), "`origin` and `destination`")
  expect_error(hammel_population(persons, end_time = 0, transitions = transform(tx, group = "race")),
               "require a declared group")
  expect_error(calculate_transition_rates(simulate(model, 0), c(0, 1), age_breaks = c(0, Inf),
                                          group = "education"), "must name one declared")
})

# Stage 4: splitting analysis by several dimensions. The reference is the same
# population with one combined column, `rr` ("A.N"), whose history merges the
# race and region moves (computed independently here).
combined_reference <- function(out) {
  persons <- out$persons
  persons$rr <- paste(persons$race, persons$region, sep = ".")
  tx <- out$transitions
  rows <- list()
  for (id in unique(tx$person_id)) {
    h <- tx[tx$person_id == id, ]
    h <- h[order(h$time, match(h$group, c("race", "region"))), ]
    state <- c(race = NA, region = NA)
    for (d in names(state)) {
      first <- h[h$group == d, ]
      state[[d]] <- if (nrow(first)) first$origin[1] else persons[[d]][persons$person_id == id]
    }
    for (k in seq_len(nrow(h))) {
      from <- paste(state[["race"]], state[["region"]], sep = ".")
      state[[h$group[k]]] <- h$destination[k]
      rows[[length(rows) + 1L]] <- data.frame(person_id = id, time = h$time[k], group = "rr",
        origin = from, destination = paste(state[["race"]], state[["region"]], sep = "."))
    }
  }
  merged <- if (length(rows)) do.call(rbind, rows) else NULL
  # Same-instant moves merge into one move in the single combined column.
  if (!is.null(merged)) {
    key <- paste(merged$person_id, merged$time)
    first <- !duplicated(key); last <- !duplicated(key, fromLast = TRUE)
    merged <- data.frame(person_id = merged$person_id[first], time = merged$time[first], group = "rr",
      origin = merged$origin[first], destination = merged$destination[last])
  }
  hammel_population(persons, unions = out$unions, end_time = out$end_time,
                    transitions = merged, groups = "rr")
}

split_rr <- function(x, suffix = "") {
  column <- paste0("rr", suffix)
  parts <- do.call(rbind, strsplit(x[[column]], ".", fixed = TRUE))
  x[[paste0("race", suffix)]] <- parts[, 1]; x[[paste0("region", suffix)]] <- parts[, 2]
  x[[column]] <- NULL
  x
}

# Sum a table's counts and exposures within the given key columns.
sum_by <- function(x, keys, values) {
  key <- do.call(paste, c(lapply(x[keys], as.character), sep = "|"))
  out <- x[!duplicated(key), keys]
  for (v in values) out[[v]] <- as.vector(tapply(x[[v]], factor(key, unique(key)), sum))
  out
}

compare_split <- function(ours, reference, keys, values) {
  ours <- as.data.frame(ours); reference <- as.data.frame(reference)
  k1 <- do.call(paste, c(lapply(ours[keys], as.character), sep = "|"))
  k2 <- do.call(paste, c(lapply(reference[keys], as.character), sep = "|"))
  # Our grid is the full cross of categories; the reference has observed
  # combinations only. Rows missing there have no exposure or events.
  expect_true(all(k2 %in% k1))
  missing <- !k1 %in% k2
  for (v in values) {
    expect_equal(ours[[v]][match(k2, k1)], reference[[v]], tolerance = 1e-9)
    expect_true(all(ours[[v]][missing] == 0))
  }
}

split_model <- function() {
  model <- two_dimension_model(600, birth_rate = .3)
  model <- set_mortality(model, data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .03))
  model <- set_group_transition(model, data.frame(origin = c("A", "B"), destination = c("B", "A"),
    age_start = 0, age_end = Inf, rate = c(.2, .1)), group = "race")
  set_group_transition(model, data.frame(origin = c("N", "S"), destination = c("S", "N"),
    age_start = 0, age_end = Inf, rate = c(.15, .25)), group = "region")
}

test_that("person-level rates split by several dimensions match a combined column", {
  out <- simulate(split_model(), 4, seed = 31)
  # A same-instant move in both dimensions.
  p <- out$persons[out$persons$birth_time < 0 & is.na(out$persons$death_time), ][1, ]
  extra <- data.frame(person_id = p$person_id, time = 3.5, group = c("race", "region"),
    origin = c(p$race, p$region), destination = c(setdiff(c("A", "B"), p$race), setdiff(c("N", "S"), p$region)))
  out$transitions <- rbind(out$transitions, extra)
  out$persons$race[out$persons$person_id == p$person_id] <- extra$destination[1]
  out$persons$region[out$persons$person_id == p$person_id] <- extra$destination[2]
  out <- hammel_population(out$persons, unions = out$unions, end_time = out$end_time,
                           transitions = out$transitions, groups = c("race", "region"))
  ref <- combined_reference(out)
  times <- c(0, 1.5, 3.5, 4)
  ours <- calculate_population_size(out, times, by = c("race", "sex", "region"))
  expect_named(ours, c("time", "race", "sex", "region", "population"))
  # The first `by` variable varies fastest after time; categories keep the
  # order of their first appearance, as with one dimension.
  expect_identical(ours$race[1:8], rep(unique(ours$race), each = 4))
  expect_identical(ours$sex[1:16], rep(unique(ours$sex), each = 8))
  expect_identical(ours$region, rep(unique(ours$region), each = 16))
  compare_split(ours, split_rr(calculate_population_size(ref, times, by = c("rr", "sex"))),
                c("time", "race", "sex", "region"), "population")
  ours <- calculate_mortality_rates(out, c(0, 2, 4), c(0, 40, Inf), by = c("region", "race"))
  expect_named(ours, c("time_start", "time_end", "sex", "age_start", "age_end", "region", "race",
                       "person_years", "deaths", "rate_mortality"))
  compare_split(ours, split_rr(calculate_mortality_rates(ref, c(0, 2, 4), c(0, 40, Inf), by = "rr")),
                c("time_start", "age_start", "sex", "region", "race"), c("person_years", "deaths"))
  ours <- calculate_fertility_rates(out, c(0, 4), c(15, 30, 50), by = c("race", "region"))
  compare_split(ours, split_rr(calculate_fertility_rates(ref, c(0, 4), c(15, 30, 50), by = "rr")),
                c("time_start", "age_start", "race", "region"), c("person_years", "births"))
  ours <- calculate_union_rates(out, c(0, 4), c(15, 50), by = c("region", "race"))
  compare_split(ours, split_rr(calculate_union_rates(ref, c(0, 4), c(15, 50), by = "rr")),
                c("time_start", "age_start", "sex", "region", "race"), c("person_years", "unions"))
  ours <- calculate_population_size(out, times, by = "region")
  reference <- sum_by(split_rr(calculate_population_size(ref, times, by = "rr")), c("time", "region"), "population")
  compare_split(ours, reference, c("time", "region"), "population")
})

test_that("union rates split by each partner's dimensions match a combined column", {
  model <- split_model()
  model <- set_union_dissolution(model, data.frame(dur_start = 0, dur_end = Inf, rate = .2))
  out <- simulate(model, 4, seed = 32)
  ref <- combined_reference(out)
  reference <- split_rr(split_rr(calculate_dissolution_rates(ref, c(0, 2, 4), c(0, Inf),
    by = c("rr_w", "rr_h")), "_w"), "_h")
  ours <- calculate_dissolution_rates(out, c(0, 2, 4), c(0, Inf), by = c("race_w", "region_h"))
  expect_named(ours, c("time_start", "time_end", "dur_start", "dur_end", "race_w", "region_h",
                       "union_years", "dissolutions", "rate_dissolution"))
  expect_equal(nrow(ours), 2 * 2 * 2)
  compare_split(ours, sum_by(reference, c("time_start", "race_w", "region_h"), c("union_years", "dissolutions")),
                c("time_start", "race_w", "region_h"), c("union_years", "dissolutions"))
  expect_equal(ours$rate_dissolution, ours$dissolutions / ours$union_years)
  ours <- calculate_dissolution_rates(out, c(0, 4), c(0, Inf), by = c("region_w", "race_w", "race_h"),
                                      conf.level = .95)
  compare_split(ours, sum_by(reference, c("region_w", "race_w", "race_h"), c("union_years", "dissolutions")),
                c("region_w", "race_w", "race_h"), c("union_years", "dissolutions"))
  expect_true(all(ours$conf_low <= ours$rate_dissolution))
})

test_that("transition rates split by another dimension count only their own moves", {
  out <- simulate(split_model(), 4, seed = 33)
  ref <- combined_reference(out)
  ours <- calculate_transition_rates(out, c(0, 4), age_breaks = c(0, 40, Inf), by = c("region", "sex"),
                                     group = "race")
  expect_named(ours, c("time_start", "time_end", "age_start", "age_end", "origin", "destination",
                       "region", "sex", "person_years", "transitions", "rate_transition"))
  # Reference rows: origin/destination are combined labels; race moves keep region.
  r <- as.data.frame(calculate_transition_rates(ref, c(0, 4), age_breaks = c(0, 40, Inf), by = "sex"))
  o <- do.call(rbind, strsplit(r$origin, ".", fixed = TRUE))
  d <- do.call(rbind, strsplit(r$destination, ".", fixed = TRUE))
  keep <- o[, 2] == d[, 2] & o[, 1] != d[, 1]
  r <- data.frame(age_start = r$age_start[keep], origin = o[keep, 1], destination = d[keep, 1],
                  region = o[keep, 2], sex = r$sex[keep], person_years = r$person_years[keep],
                  transitions = r$transitions[keep])
  compare_split(ours, r, c("age_start", "origin", "destination", "region", "sex"),
                c("person_years", "transitions"))
  expect_equal(sum(ours$transitions), sum(out$transitions$group == "race"))
  expect_error(calculate_transition_rates(out, c(0, 4), age_breaks = c(0, Inf), by = "race", group = "race"),
               "not by `race` in `by`")
})

test_that("durations by another dimension count from entering this dimension's category", {
  persons <- data.frame(person_id = 1, sex = "f", birth_time = -10, death_time = NA_real_,
    mother_id = NA_real_, father_id = NA_real_, race = "B", region = "S")
  # Race A -> B at 1; region N -> S at 2 (does not reset race duration).
  tx <- data.frame(person_id = 1, time = c(1, 2), group = c("race", "region"),
                   origin = c("A", "N"), destination = c("B", "S"))
  pop <- hammel_population(persons, end_time = 4, transitions = tx, groups = c("race", "region"))
  out <- calculate_transition_rates(pop, c(0, 4), dur_breaks = c(0, 1, 2, Inf), by = "region", group = "race")
  b <- out[out$origin == "B", ]
  # In B from 1: duration [0,1) in N (time 1-2), [1,2) in S (2-3), [2,Inf) in S (3-4).
  expect_equal(b$person_years[b$region == "N"], c(1, 0, 0))
  expect_equal(b$person_years[b$region == "S"], c(0, 1, 1))
  a <- out[out$origin == "A", ]
  expect_equal(sum(a$transitions[a$region == "N"]), 1)
  expect_equal(a$person_years[a$region == "N"], c(0, 0, 1))   # age 10-11 in A before the move
})

test_that("summaries keep group dimension columns and errors name the dimensions", {
  out <- simulate(split_model(), 3, seed = 34)
  asfr <- calculate_fertility_rates(out, c(0, 3), seq(15, 50, 5), by = c("race", "region"))
  tfr <- calculate_tfr(asfr)
  expect_named(tfr, c("time_start", "time_end", "race", "region", "age_start", "age_end", "tfr"))
  expect_equal(nrow(tfr), 4)
  one <- asfr[asfr$race == tfr$race[1] & asfr$region == tfr$region[1], ]
  expect_equal(tfr$tfr[1], sum(one$rate_fertility * 5))
  mx <- calculate_mortality_rates(out, c(0, 3), c(0, 1, seq(5, 85, 5), Inf), by = "region")
  life <- calculate_life_table(mx)
  expect_true(all(c("sex", "region") %in% names(life)))
  expect_equal(sum(life$age_start == 0), 4)
  expect_error(calculate_mortality_rates(out, c(0, 3), c(0, Inf), by = "education"),
               "`by` accepts \"race\", \"region\", \"union_type\"")
  expect_error(calculate_dissolution_rates(out, c(0, 3), c(0, Inf), by = "race"),
               "accepts \"race_w\", \"region_w\", \"race_h\", \"region_h\"")
  expect_error(calculate_mortality_rates(out, c(0, 3), c(0, Inf), by = "group"),
               "not `group`")
})

test_that("weight models use any dimension's variables", {
  model <- two_dimension_model(400)
  persons <- model$population$persons
  # Deaths in a pooled cell go almost only to region S.
  rates <- data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .05)
  weighted <- set_mortality(model, rates, weights = create_weight_model(~ region, c(regionS = 12)))
  out <- simulate(weighted, 3, seed = 51)
  dead <- out$persons[!is.na(out$persons$death_time), ]
  expect_gt(nrow(dead), 30)
  expect_true(all(dead$region == "S"))
  # Race moves weighted by region; ancestry of a tracked dimension.
  tracked <- set_groups(model, race = group(ancestry = TRUE), region = group())
  moves <- set_group_transition(tracked, data.frame(origin = c("A", "B"), destination = c("B", "A"),
    age_start = 0, age_end = Inf, rate = .3), group = "race",
    weights = create_weight_model(~ region + ancestry_race_A, c(regionS = 12, ancestry_race_A = 0)))
  out <- simulate(moves, 2, seed = 52)
  movers <- out$transitions$person_id[out$transitions$group == "race"]
  expect_gt(length(movers), 20)
  expect_true(all(out$persons$region[match(movers, out$persons$person_id)] == "S"))
  expect_error(simulate(set_mortality(model, rates, weights = create_weight_model(~ ancestry_race_A,
    c(ancestry_race_A = 1))), 0), "not available for this process: ancestry_race_A")
  # Union variables by dimension: dissolutions only where regions differ.
  persons <- data.frame(person_id = 1:400, sex = rep(c("f", "m"), 200), birth_time = -30,
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_,
    race = "A", region = c(rep(c("N", "N"), 100), rep(c("N", "S"), 100)))
  unions <- data.frame(union_id = 1:200, wife_id = seq(1, 400, 2), husband_id = seq(2, 400, 2),
    start_time = -1, end_time = NA_real_, end_reason = NA_character_)
  couples <- hammel_model() |>
    set_population(hammel_population(persons, unions = unions, end_time = 0, groups = c("race", "region"))) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0)) |>
    set_union_formation("on_demand") |>
    set_union_dissolution(data.frame(dur_start = 0, dur_end = Inf, rate = .1),
      weights = create_weight_model(~ same_region + region_h, c(same_region = -12, region_hS = 0)))
  out <- simulate(couples, 3, seed = 53)
  ended <- !is.na(out$unions$end_time)
  expect_gt(sum(ended), 20)
  expect_true(all(out$unions$union_id[ended] > 100))
  expect_error(set_union_dissolution(couples, data.frame(dur_start = 0, dur_end = Inf, rate = .1),
    weights = create_weight_model(~ partnered, c(partnered = 1))), "not available")
})

test_that("custom weight functions see every dimension in members and cells", {
  model <- two_dimension_model(200, union = TRUE)
  seen <- new.env()
  record <- function(name) function(members, bin) {
    if (is.null(seen[[name]])) assign(name, list(members = members, bin = bin), envir = seen)
    rep(0, nrow(members))
  }
  model <- set_mortality(model, stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf,
    rate = .2), list(region = c(N = 1, S = 2))), weights = create_weight_model(fun = record("mortality")))
  model <- set_union_dissolution(model, stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = .3),
    list(race_w = c(A = 1, B = 2))), weights = create_weight_model(fun = record("dissolution")))
  model <- set_group_transition(model, data.frame(origin = "N", destination = "S", race = "B",
    age_start = 0, age_end = Inf, rate = .5), group = "region",
    weights = create_weight_model(fun = record("transition")))
  simulate(model, 2, seed = 54)
  m <- seen$mortality
  expect_true(all(c("race", "region", "partner_race", "partner_region") %in% names(m$members)))
  expect_true(m$bin$region %in% c("N", "S"))
  expect_true(is.na(m$bin$race))
  d <- seen$dissolution
  expect_true(all(c("race_w", "race_h", "region_w", "region_h", "same_race", "same_region") %in%
                  names(d$members)))
  expect_true(d$bin$race_w %in% c("A", "B"))
  expect_true(is.na(d$bin$race_h) && is.na(d$bin$region_w))
  expect_null(d$bin$group_pair)
  t <- seen$transition
  expect_identical(c(t$bin$origin, t$bin$destination, t$bin$race), c("N", "S", "B"))
  expect_null(t$bin$region)
  expect_identical(unique(t$members$destination), "S")
})

test_that("calendar bands switch rates stratified by several dimensions", {
  model <- two_dimension_model(400)
  persons <- model$population$persons
  # Before time 1 only region S dies; from time 1 only race B.
  early <- stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1e8,
    time_start = -100, time_end = 1), list(race = c(A = 1, B = 1), region = c(N = 0, S = 1)))
  late <- stratify_rates(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 1e8,
    time_start = 1, time_end = Inf), list(race = c(A = 0, B = 1), region = c(N = 1, S = 1)))
  out <- simulate(set_mortality(model, rbind(early, late)), 2, seed = 61)
  died <- out$persons$death_time[seq_len(nrow(persons))]
  expect_true(all(died[persons$region == "S"] < 1e-6))
  expect_true(all(abs(died[persons$region == "N" & persons$race == "B"] - 1) < 1e-6))
  expect_true(all(is.na(died[persons$region == "N" & persons$race == "A"])))
  # Every active slice must cover each combination.
  gap <- rbind(early, late[!(late$race == "A" & late$region == "S"), ])
  expect_error(simulate(set_mortality(model, gap), 2), "missing schedules.*race=A, region=S")
})

test_that("cohabitation probabilities and conversion use declared dimensions", {
  model <- matching_pair_model(create_age_gap_matching_model(0, 5))
  model <- set_union_formation(model, "rate",
    rates = data.frame(age_start = 0, age_end = Inf, rate = 2),
    cohab_prob = data.frame(race = c("A", "B"), age_start = 0, age_end = Inf, prob = c(1, 0)),
    union_convert_rate = stratify_rates(data.frame(dur_start = 0, dur_end = Inf, rate = 1e8),
                                        list(region_h = c(N = 1, S = 0))),
    kinship_generations = 0)
  model <- set_fertility(model, data.frame(union_type = c("unpartnered", "cohab", "marriage"),
    age_start = 15, age_end = 50, rate = 0))
  out <- simulate(model, 1, seed = 62)
  u <- get_detailed_unions(out)
  first <- u[is.na(u$previous_cohab_id), ]
  expect_gt(nrow(first), 50)
  # Women of race A cohabit, women of race B marry.
  expect_true(all(first$union_type[first$race_w == "A"] == "cohab"))
  expect_true(all(first$union_type[first$race_w == "B"] == "marriage"))
  # Cohabitations with a northern husband convert at once; southern ones never.
  cohab <- first[first$union_type == "cohab", ]
  expect_true(all(cohab$end_reason[cohab$region_h == "N"] %in% "conversion"))
  expect_true(all(is.na(cohab$end_reason[cohab$region_h == "S"])))
  rates <- calculate_cohab_convert_rate(out, c(0, 1), dur_breaks = c(0, Inf), by = c("region_h", "race_w"))
  expect_named(rates, c("time_start", "time_end", "dur_start", "dur_end", "region_h", "race_w",
                        "union_years", "conversions", "rate_cohab_convert"))
  expect_equal(sum(rates$conversions), sum(cohab$region_h == "N"))
})
