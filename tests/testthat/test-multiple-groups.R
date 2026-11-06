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
})

test_that("features of a single dimension explain that they are not yet available", {
  model <- two_dimension_model(union = TRUE)
  exo <- set_union_formation(model, "on_demand", kinship_generations = 0,
    matching = create_age_gap_matching_model(0, 5, group_exogamy = c("A:B" = -1)))
  expect_error(simulate(exo, 0), "Group effects in matching")
  expect_error(calculate_mortality_rates(simulate(model, 1, seed = 6), c(0, 1), c(0, Inf), by = "race"),
               "exactly one declared group dimension")
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
  w <- create_weight_model(~ group, c(groupB = 1), group_levels = c("A", "B"))
  expect_error(simulate(set_mortality(model, model$processes$mortality$rates, weights = w), 0),
               "Group variables in weight models")
  moves <- data.frame(origin = "N", destination = "S", age_start = 0, age_end = Inf, rate = .1)
  w <- create_weight_model(~ destination, c(destinationS = 1), group_levels = c("N", "S"))
  expect_warning(weighted <- set_group_transition(model, moves, group = "region", weights = w), "no effect")
  expect_error(simulate(weighted, 0), "Group variables in weight models")
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
