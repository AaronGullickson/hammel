declaration_model <- function(n = 20, rate = 5) {
  model <- make_model(n, birth_rate = rate, birth_interval = 0)
  model$population$persons$sex <- rep(c("f", "m"), length.out = n)
  model$population$persons$group <- rep(c("A", "B"), length.out = n)
  model
}

test_that("group() validates its settings and normalizes inheritance", {
  g <- group(inheritance = "mother")
  expect_s3_class(g, "hammel_group")
  expect_identical(g$inheritance, group_inheritance("mother"))
  expect_false(g$ancestry)
  expect_identical(group()$inheritance$method, "random_parent")
  expect_error(group(categories = c("A", "A")), "unique")
  expect_error(group(categories = character()), "nonempty")
  expect_error(group(categories = c("A", NA)), "group")
  expect_error(group(ancestry = NA), "TRUE or FALSE")
  expect_error(group(inheritance = "bogus"), "method")
})

test_that("declarations name dimensions, which must not be reserved", {
  model <- declaration_model()
  expect_error(set_groups(model, sex = group()), "reserved")
  expect_error(set_groups(model, race_w = group()), "reserved")
  expect_error(set_groups(model, partner_race = group()), "reserved")
  expect_error(set_groups(model, race = group()), "not a persons column")
  expect_error(set_groups(model, group = group(), group = group()), "unique")
  model$population$persons$region <- rep(c("N", "S"), length.out = nrow(model$population$persons))
  m <- set_groups(model, group = group(), region = group())
  expect_named(m$population$groups, c("group", "region"))
  expect_error(set_groups(model, group()), "names")
  expect_error(set_groups(model, group = "mother"), "group\\(\\)")
  expect_error(hammel_population(model$population$persons, end_time = 0, groups = 1), "groups")
})

test_that("set_groups works before or after the population and replaces its declaration", {
  model <- declaration_model()
  before <- hammel_model() |> set_groups(group = group(inheritance = "mother"))
  before <- set_population(before, model$population)
  before$processes <- model$processes
  after <- set_groups(model, group = group(inheritance = "mother"))
  expect_identical(histories(run_simulation(before, 1, seed = 3)), histories(
    run_simulation(after, 1, seed = 3)))
  out <- run_simulation(after, 1, seed = 3)
  expect_identical(out$groups$group$inheritance$method, "mother")
  kids <- out$persons[out$persons$birth_time >= 0, ]
  moms <- out$persons$group[match(kids$mother_id, out$persons$person_id)]
  expect_identical(kids$group, moms)
  # The result carries the declaration into restarts and analysis.
  restarted <- set_population(hammel_model(), out)
  expect_identical(restarted$population$groups, out$groups)
  expect_true("group" %in% names(calculate_population_size(out, 1, by = "group")))
  # A later declaration replaces a population's own (here the baseline's).
  base <- create_baseline_population(10, groups = list(group = c(A = 1, B = 1)))
  expect_false(base$groups$group$ancestry)
  m <- set_groups(set_population(hammel_model(), base), group = group(ancestry = TRUE))
  expect_true(m$population$groups$group$ancestry)
  expect_null(set_groups(m)$population$groups)
})

test_that("an undeclared group column is an inert custom column", {
  model <- declaration_model()
  out <- run_simulation(model, 1, seed = 4)
  expect_null(out$groups)
  expect_null(out$ancestry)
  expect_null(out$transitions)
  n <- nrow(model$population$persons)
  expect_identical(out$persons$group[seq_len(n)], model$population$persons$group)
  expect_true(all(is.na(out$persons$group[-seq_len(n)])))
  # It is not a group: rates and analysis cannot use it.
  stratified <- model
  stratified$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
    list(group = c(A = 1, B = 1)))
  expect_error(run_simulation(stratified, 1),
    "not a declared group dimension.*set_groups")
  expect_error(calculate_mortality_rates(out, c(0, 1), c(0, Inf), by = "group"), "declared group")
})

test_that("categories fix the complete ordered list", {
  model <- declaration_model()
  expect_error(set_groups(model, group = group(categories = "A")), "missing from its categories.*B")
  # Extra categories may start empty; group-specific rates must cover them.
  m <- set_groups(model, group = group(categories = c("B", "A", "C")))
  m$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
    list(group = c(A = 1, B = 1)))
  expect_error(run_simulation(m, 0), "missing schedules.*C")
  m$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
    list(group = c(A = 1, B = 1, C = 1)))
  out <- run_simulation(m, 1, seed = 5)
  expect_identical(out$groups$group$categories, c("B", "A", "C"))
  # A custom rule must return a listed category.
  rule <- group_inheritance("custom", fun = function(...) "D")
  expect_error(run_simulation(set_groups(model, group = group(inheritance =
    rule)), 1, seed = 5),
               "undeclared group identity")
})

test_that("ancestry is tracked only when requested", {
  model <- declaration_model()
  plain <- run_simulation(set_groups(model, group = group()), 1, seed = 6)
  expect_null(plain$ancestry)
  tracked <- run_simulation(set_groups(model, group = group(ancestry = TRUE)),
    1, seed = 6)
  expect_named(tracked$ancestry, "group")
  expect_s3_class(tracked$ancestry$group, "tbl_df")
  expect_identical(tracked$ancestry$group$person_id, tracked$persons$person_id)
  expect_identical(histories(plain)[c("persons", "unions")], histories(tracked)[c("persons", "unions")])
  expect_error(run_simulation(set_groups(model, group = group(inheritance =
    "ancestry")), 1),
               "ancestry = TRUE")
  # Turning tracking off drops an existing table.
  expect_null(set_groups(set_population(hammel_model(), tracked), group = group())$population$ancestry)
})

test_that("baseline populations sample independent or joint group distributions", {
  set.seed(7)
  independent <- create_baseline_population(4000, groups = list(group = c(A = 3, B = 1)))
  expect_lt(abs(mean(independent$persons$group == "A") - .75), .03)
  expect_named(independent$groups, "group")
  set.seed(7)
  joint <- create_baseline_population(4000,
    groups = data.frame(group = c("A", "B"), weight = c(3, 1)))
  expect_lt(abs(mean(joint$persons$group == "A") - .75), .03)
  expect_named(joint$groups, "group")
  for (bad in list(data.frame(group = "A"), data.frame(weight = 1), data.frame(group = "A", weight = -1),
                   data.frame(group = NA, weight = 1), list(), list(group = c(1, 2)),
                   list(sex = c(A = 1)), c(A = 1)))
    expect_error(create_baseline_population(10, groups = bad))
  # An empty baseline keeps its typed column and declaration.
  empty <- create_baseline_population(0, groups = data.frame(group = c("A", "B"), weight = 1))
  expect_identical(empty$persons$group, character())
})

test_that("baselines declare categories in the order given, and set_groups keeps them", {
  set.seed(3)
  # B may well appear first in the persons table; the declaration still lists A first.
  pop <- create_baseline_population(200, groups = list(race = c(A = .1, B = .9, AB = 0)))
  expect_identical(pop$groups$race$categories, c("A", "B", "AB"))
  expect_false("AB" %in% pop$persons$race)
  joint <- create_baseline_population(50, groups = data.frame(race = c("W", "B", "W", "B"),
    region = c("S", "S", "N", "N"), weight = c(1, 1, 1, 0)))
  expect_identical(joint$groups$race$categories, c("W", "B"))
  expect_identical(joint$groups$region$categories, c("S", "N"))
  codes <- create_baseline_population(20, groups = data.frame(code = c(20L, 10L), weight = 1))
  expect_identical(codes$groups$code$categories, c(20L, 10L))
  # Omitted categories are kept, before or after the population is set.
  model <- hammel_model() |> set_population(pop) |> set_groups(race = group(ancestry = TRUE))
  expect_identical(model$population$groups$race$categories, c("A", "B", "AB"))
  expect_true(model$population$groups$race$ancestry)
  early <- hammel_model() |> set_groups(race = group(inheritance = "mother")) |> set_population(pop)
  expect_identical(early$population$groups$race$categories, c("A", "B", "AB"))
  # Explicit categories replace them, and must cover the population.
  model <- set_groups(model, race = group(categories = c("B", "A", "AB")))
  expect_identical(model$population$groups$race$categories, c("B", "A", "AB"))
  expect_error(set_groups(model, race = group(categories = "A")), "missing from its categories")
  # The kept order sets the reference level at simulation.
  model <- model |>
    set_groups(race = group()) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = .01),
                  weights = weight_model(~ race, c(raceA = 0, raceAB = 0))) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = .1)) |>
    set_union_formation("on_demand")
  expect_s3_class(run_simulation(model, 1, seed = 4), "hammel_population")
})
