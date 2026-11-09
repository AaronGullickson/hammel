inheritance_test_model <- function(rule = group_inheritance(), birth_rate = 20,
                                   categories = NULL) {
  model <- make_model(2, birth_rate = birth_rate, birth_interval = 0)
  model$population$persons$sex <- c("f", "m")
  model$population$persons$group <- c("A", "B")
  set_groups(model, group = group(categories = categories, ancestry = TRUE, inheritance = rule))
}

inheritance_children <- function(result) result$persons[result$persons$birth_time >= 0, ]

test_that("inheritance method strings are equivalent to constructed models", {
  for (method in c("mother", "father", "random_parent", "same_sex", "opposite_sex", "ancestry")) {
    model <- inheritance_test_model(method)
    explicit <- inheritance_test_model(group_inheritance(method))
    expect_identical(model$groups$group$inheritance, explicit$groups$group$inheritance)
    expect_identical(run_simulation(model, 1, seed = 52), run_simulation(
      explicit, 1, seed = 52))
  }
  for (method in list("other", NA_character_, c("mother", "father"), character()))
    expect_error(inheritance_test_model(method), "method")
  expect_error(inheritance_test_model("hierarchy"), "order")
  expect_error(inheritance_test_model("custom"), "fun")
})

test_that("custom inheritance receives parents' union status at birth", {
  rule <- group_inheritance("custom",
    fun = function(parents_in_union, ...) {
      stopifnot(is.logical(parents_in_union), length(parents_in_union) == 1L,
                !is.na(parents_in_union))
      if (parents_in_union) "partnered" else "nonmarital"
    })
  model <- inheritance_test_model(rule, categories = c("A", "B", "partnered", "nonmarital"))
  result <- run_simulation(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 0)
  expect_true(all(children$group == "partnered"))
  expect_gt(nrow(result$unions), 0)

  model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0)
  model <- set_union_formation(model, "rate", rates = data.frame(
    age_start = 0, age_end = Inf, rate = 0), cancel_birth_if_no_match = FALSE)
  for (known_father in c(TRUE, FALSE)) {
    unpartnered <- model
    if (!known_father) unpartnered$population$persons <- model$population$persons[1, ]
    result <- run_simulation(unpartnered, 1, seed = 52)
    children <- inheritance_children(result)
    expect_gt(nrow(children), 0)
    expect_true(all(children$group == "nonmarital"))
    expect_true(all(is.na(children$father_id) == !known_father))
    expect_equal(nrow(result$unions), 0)
  }
})

test_that("built-in inheritance uses parents and newborn sex", {
  for (method in c("mother", "father", "same_sex", "opposite_sex")) {
    result <- run_simulation(inheritance_test_model(group_inheritance(method)),
      1, seed = 52)
    children <- inheritance_children(result)
    expect_gt(nrow(children), 5)
    expect_true(all(c("m", "f") %in% children$sex))
    expected <- switch(method, mother = rep("A", nrow(children)), father = rep("B", nrow(children)),
      same_sex = ifelse(children$sex == "f", "A", "B"),
      opposite_sex = ifelse(children$sex == "f", "B", "A"))
    expect_identical(children$group, expected)
    shares <- result$ancestry$group[result$persons$birth_time >= 0, ]
    expect_true(all(shares$A == 0.5 & shares$B == 0.5))
  }
})

test_that("default, explicit, and missing inheritance specifications use random parent", {
  model <- inheritance_test_model()
  set.seed(654)
  original <- run_simulation(model, 1)
  after <- .Random.seed
  expect_identical(group_inheritance()$method, "random_parent")
  expect_identical(model$groups$group$inheritance$method, "random_parent")
  explicit <- inheritance_test_model("random_parent")
  set.seed(654)
  expect_identical(run_simulation(explicit, 1), original)
  expect_identical(.Random.seed, after)
  model <- set_inheritance(model, NULL)
  set.seed(654)
  expect_identical(run_simulation(model, 1), original)
  expect_identical(.Random.seed, after)
  # Without a declared group, inheritance does not apply.
  plain <- make_model(2)
  expect_no_error(run_simulation(plain, 0))
  expect_null(run_simulation(plain, 1, seed = 654)$groups)
})

test_that("missing fathers fall back to maternal identity without extra parent-choice draws", {
  original <- inheritance_test_model()
  original$population$persons <- original$population$persons[1, ]
  reference <- run_simulation(original, 1, seed = 52)
  for (method in c("father", "random_parent", "same_sex", "opposite_sex", "hierarchy")) {
    rule <- if (method == "hierarchy") group_inheritance(method, order = c("B", "A"), passing = 0.5) else
      group_inheritance(method)
    model <- original
    model <- set_inheritance(model, rule)
    histories <- c("persons", "unions", "ancestry", "transitions")
    expect_identical(run_simulation(model, 1, seed = 52)[histories], reference[
      histories])
  }
})

test_that("random-parent and hierarchy probabilities agree with analytical cases", {
  model <- inheritance_test_model(group_inheritance("random_parent"), birth_rate = 2000)
  result <- run_simulation(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 1000)
  expect_lt(abs(mean(children$group == "A") - 0.5), 5 * sqrt(0.25/nrow(children)))
  for (passing in c(0, 1, 0.05)) {
    model <- set_inheritance(model, group_inheritance("hierarchy", order = c("A", "B"), passing = passing))
    children <- inheritance_children(run_simulation(model, 1, seed = 52))
    if (passing %in% c(0, 1)) expect_true(all(children$group == if (passing == 0) "A" else "B")) else
      expect_lt(abs(mean(children$group == "B") - passing), 5 * sqrt(passing*(1-passing)/nrow(children)))
  }
  model <- inheritance_test_model(group_inheritance("hierarchy", order = "B", passing = 1))
  histories <- c("persons", "unions", "ancestry", "transitions")
  expect_identical(run_simulation(model, 1, seed = 52)[histories],
                   run_simulation(inheritance_test_model("mother"), 1, seed =
                     52)[histories])
})

test_that("ancestry assignment conditions on known shares and handles all-unknown ancestry", {
  model <- inheritance_test_model(group_inheritance("ancestry"), birth_rate = 2000)
  # Both parents identify A, but carry a known B reference component.
  p <- model$population$persons
  p$group <- "A"
  a <- data.frame(person_id = p$person_id, A = 0.35, B = 0.15, unknown = 0.5)
  model <- set_population(model, hammel_population(p, end_time = 0, ancestry = a,
    groups = list(group = group(ancestry = TRUE))))
  result <- run_simulation(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_lt(abs(mean(children$group == "A") - 0.7), 5 * sqrt(0.7*0.3/nrow(children)))
  expect_true(all(result$ancestry$group$unknown == 0.5))
  model$population$ancestry$group$A <- 0
  model$population$ancestry$group$B <- 0
  model$population$ancestry$group$unknown <- 1
  model$processes$fertility$rates$rate <- 20
  maternal <- model; maternal$groups$group$inheritance <- group_inheritance()
  expect_identical(histories(run_simulation(model, 1, seed = 52)), histories(
    run_simulation(maternal, 1, seed = 52)))
})

test_that("custom callbacks receive copies of all inputs and use the unified R RNG stream", {
  custom <- group_inheritance("custom", fun = function(mother_group, father_group,
      mother_ancestry, father_ancestry, child_ancestry, child_sex) {
    stopifnot(mother_group == "A", father_group == "B", child_sex %in% c("m", "f"),
      identical(names(child_ancestry), c("A", "B", "unknown")),
      mother_ancestry[["A"]] == 1, father_ancestry[["B"]] == 1,
      child_ancestry[["A"]] == 0.5)
    mother_ancestry[] <- 0 # Copied values cannot alter engine ancestry.
    if (runif(1) < 0.5) mother_group else father_group
  })
  model <- inheritance_test_model(custom)
  set.seed(654); state <- .Random.seed
  actual <- run_simulation(model, 1, seed = 52)
  expect_identical(.Random.seed, state)
  native <- inheritance_test_model(group_inheritance("random_parent"))
  expect_identical(histories(actual), histories(run_simulation(native, 1, seed =
    52)))
  set.seed(52)
  actual <- run_simulation(model, 1)
  after <- .Random.seed
  set.seed(52)
  expect_identical(histories(run_simulation(native, 1)), histories(actual))
  expect_identical(.Random.seed, after)
  # Multinomial-logit-style ancestry probabilities can be implemented in R.
  rule <- group_inheritance("custom", fun = function(child_ancestry, ...) {
    weights <- exp(c(A = child_ancestry[["A"]], B = child_ancestry[["B"]]))
    sample(names(weights), 1, prob = weights)
  })
  expect_no_error(run_simulation(inheritance_test_model(rule), 1, seed = 52))
})

test_that("custom mixed identities are declared, use their own schedules, and retain A/B ancestry", {
  rule <- group_inheritance("custom", fun = function(mother_group, father_group, ...) {
    if (is.na(father_group) || mother_group == father_group) mother_group else "AB"
  })
  model <- inheritance_test_model(rule, categories = c("A", "B", "AB"))
  # AB infants die immediately; A/B founders survive. No AB ancestry component.
  rates <- data.frame(sex = rep(c("m", "f"), each = 2), age_start = c(0, 1),
    age_end = c(1, Inf), rate = 0)
  rates <- rbind(transform(rates, group = "A"), transform(rates, group = "B"),
                 transform(rates, group = "AB", rate = ifelse(age_start == 0, 1e8, 0)))
  model <- set_mortality(model, rates)
  result <- run_simulation(model, 0.5, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 0)
  expect_true(all(children$group == "AB"))
  expect_true(all(!is.na(children$death_time)))
  expect_false("AB" %in% names(result$ancestry$group))
  expect_true(all(result$ancestry$group$A[result$persons$birth_time >= 0] == 0.5))
  continued <- run_simulation(set_population(model, result), 0.5, seed = 53)
  expect_identical(continued$ancestry$group[seq_len(nrow(result$ancestry$group)), ], result$ancestry$group)
})

test_that("future identities are checked against rates and categorical matching before execution", {
  model <- inheritance_test_model(group_inheritance("custom", fun = function(...) "AB"),
                                  categories = c("A", "B", "AB"))
  for (name in c("mortality", "fertility")) {
    edited <- model
    rates <- edited$processes[[name]]$rates
    edited$processes[[name]]$rates <- rbind(transform(rates, group = "A"), transform(rates, group = "B"))
    expect_error(run_simulation(edited, 0), "missing schedules.*AB")
  }
  expect_warning(matching <- matching_model(~ male_group, c(male_groupB = 0)), "male_group` are absorbed")
  model$processes$union_formation$matching <- matching
  # Declared categories include AB, so its coefficient is required.
  expect_error(run_simulation(model, 0), "coefficient names.*male_groupAB")
  # Positive reference shares also infer identities; zero-only columns do not.
  model <- inheritance_test_model()
  model <- set_population(model, hammel_population(model$population$persons, end_time = 0, groups = if ("group" %in% names(model$population$persons)) "group"))
  model$population$ancestry$group$C <- 0
  model$population$ancestry$group$D <- 0
  expect_identical(group_registries(model$population)$group, c("A", "B"))
  model$population$ancestry$group$A[1] <- 0.5
  model$population$ancestry$group$C[1] <- 0.5
  expect_identical(group_registries(model$population)$group, c("A", "B", "C"))
})

test_that("custom errors and invalid returns restore explicit RNG state", {
  for (fun in list(function(...) NA_character_, function(...) NULL, function(...) c("A", "B"),
                   function(...) "undeclared", function(...) stop("custom failed"))) {
    model <- inheritance_test_model(group_inheritance("custom", fun = fun))
    set.seed(654); state <- .Random.seed
    expect_error(run_simulation(model, 1, seed = 52))
    expect_identical(.Random.seed, state)
  }
  rule <- group_inheritance("custom", fun = function(mother_group, father_group, father_ancestry, ...) {
    stopifnot(is.na(father_group), is.null(father_ancestry)); mother_group
  })
  model <- inheritance_test_model(rule)
  model$population$persons <- model$population$persons[1, ]
  expect_true(all(inheritance_children(run_simulation(model, 1, seed =
    52))$group == "A"))
})

test_that("inheritance validates edited specifications and configuration", {
  expect_error(group_inheritance("other"), "method")
  expect_error(group_inheritance("custom"), "fun")
  expect_error(group_inheritance("hierarchy"), "order")
  expect_error(group_inheritance("hierarchy", order = c("A", "A")), "unique")
  expect_error(group_inheritance("mother", passing = 0.1), "hierarchy")
  for (passing in list(-1, 2, NA, "0", c(0, 1)))
    expect_error(group_inheritance("hierarchy", order = c("A", "B"), passing = passing), "passing")
  expect_error(group(categories = c("A", NA)), "group")
  model <- inheritance_test_model()
  model$groups$group$inheritance$method <- "bad"
  expect_error(run_simulation(model, 0), "method")
})

test_that("numeric identities remain numeric, and empty populations stay typed", {
  model <- inheritance_test_model(group_inheritance("father"))
  model$population$persons$group <- c(20L, 10L)
  model$population$groups <- list(group = group())
  result <- run_simulation(model, 0.5, seed = 52)
  expect_type(result$persons$group, "integer")
  expect_true(all(inheritance_children(result)$group == 10L))
  model$population$persons <- model$population$persons[FALSE, ]
  result <- run_simulation(model, 1, seed = 52)
  expect_equal(nrow(result$persons), 0)
  expect_type(result$persons$group, "integer")
  expect_s3_class(result$ancestry$group, "tbl_df")
})

test_that("new AB identities are available to categorical partner matching during the run", {
  rule <- group_inheritance("custom", fun = function(...) "AB")
  model <- inheritance_test_model(rule, categories = c("A", "B", "AB"))
  # All newborns are male: until AB women enter the search ages, AB men's
  # mirrored baseline uses the overall female rate rather than a zero rate
  # averaged across newborn girls. The high delta makes their selection clear.
  model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0, ffab = 0)
  expect_warning(matching <- matching_model(~ male_group,
    c(male_groupB = 0, male_groupAB = .5),
    delta = c(AB = 20)), "male_group` are absorbed")
  model <- set_union_formation(model, "rate", matching, rates = data.frame(
    age_start = 0, age_end = Inf, rate = 0), min_male_age = 0, kinship_generations = 0,
    cancel_birth_if_no_match = FALSE)
  result <- run_simulation(model, 1, seed = 52)
  children <- inheritance_children(result)
  father_groups <- result$persons$group[match(children$father_id, result$persons$person_id)]
  expect_true(any(father_groups == "AB"))
  expect_true(all(children$group == "AB"))
})

test_that("future search schedules, stratified rates, and numeric declarations are supported", {
  rule <- group_inheritance("custom", fun = function(...) "AB")
  model <- inheritance_test_model(rule, categories = c("A", "B", "AB"))
  model <- set_fertility(model, data.frame(union_type = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0)
  model <- set_union_formation(model, "rate", rates = data.frame(group = c("A", "B"),
    age_start = 0, age_end = Inf, rate = 0))
  expect_error(run_simulation(model, 0), "missing schedules.*AB")
  model <- inheritance_test_model(rule, categories = c("A", "B", "AB"))
  model$processes$mortality$rates$rate <- 1e8
  model$processes$mortality$rates <- stratify_rates(model$processes$mortality$rates,
    list(group = c(A = 0, B = 0, AB = 1)))
  result <- run_simulation(model, 0.5, seed = 52)
  expect_true(all(!is.na(inheritance_children(result)$death_time)))
  model <- inheritance_test_model(group_inheritance("custom", fun = function(...) 30L))
  model$population$persons$group <- c(20L, 10L)
  model <- set_groups(model, group = group(categories = c(20L, 10L, 30L),
    inheritance = group_inheritance("custom", fun = function(...) 30L)))
  result <- run_simulation(model, 0.5, seed = 52)
  expect_type(result$persons$group, "integer")
  expect_true(all(inheritance_children(result)$group == 30L))
})

test_that("fixed inheritance validates a single explicit category", {
  for (value in list(NULL, character(), c("A", "B"), NA_character_, "", " ", Inf,
                     TRUE, matrix("A"), 1 + 1i))
    expect_error(group_inheritance("fixed", fixed_value = value), "value")
  expect_error(group_inheritance("mother", fixed_value = "A"), "only to fixed")
  expect_error(inheritance_test_model("fixed"), "value")
  model <- inheritance_test_model(group_inheritance("fixed", fixed_value = "none"))
  expect_error(run_simulation(model, 0), "undeclared")
})

test_that("fixed newborn categories preserve founders, ancestry and deterministic RNG", {
  fixed <- group_inheritance("fixed", fixed_value = "none")
  custom <- group_inheritance("custom", fun = function(...) "none")
  categories <- c("A", "B", "none")
  model <- inheritance_test_model(fixed, categories = categories)
  reference <- inheritance_test_model(custom, categories = categories)
  set.seed(726)
  result <- run_simulation(model, 1)
  after <- .Random.seed
  set.seed(726)
  expected <- run_simulation(reference, 1)
  expect_identical(.Random.seed, after)
  expect_identical(result[c("persons", "unions", "ancestry", "transitions")],
                   expected[c("persons", "unions", "ancestry", "transitions")])
  children <- inheritance_children(result)
  expect_gt(nrow(children), 5)
  expect_true(all(c("m", "f") %in% children$sex))
  expect_true(all(children$group == "none"))
  expect_identical(result$persons$group[1:2], c("A", "B"))
  shares <- result$ancestry$group[result$persons$birth_time >= 0, ]
  expect_true(all(shares$A == .5 & shares$B == .5))
  expect_identical(run_simulation(model, 1, seed = 726), run_simulation(model,
    1, seed = 726))

  model$population$persons <- model$population$persons[1, ]
  model <- set_union_formation(model, "rate", rates = data.frame(
    age_start = 0, age_end = Inf, rate = 0), cancel_birth_if_no_match = FALSE)
  result <- suppressWarnings(run_simulation(model, 1, seed = 726))
  children <- inheritance_children(result)
  expect_gt(nrow(children), 0)
  expect_true(all(is.na(children$father_id)))
  expect_true(all(children$group == "none"))
})

test_that("fixed inheritance supports numeric categories and calendar switches", {
  model <- inheritance_test_model()
  model$population$persons$group <- c(1, 2)
  model <- set_groups(model, group = group(categories = c(0, 1, 2),
    inheritance = group_inheritance("fixed", fixed_value = 0)))
  children <- inheritance_children(run_simulation(model, 1, seed = 726))
  expect_gt(nrow(children), 0)
  expect_identical(children$group, rep(0, nrow(children)))

  schedule <- tibble::tibble(time_start = c(0, .5), time_end = c(.5, Inf),
    model = list(group_inheritance("fixed", fixed_value = "none"),
                 group_inheritance("fixed", fixed_value = "primary")))
  model <- inheritance_test_model(schedule, birth_rate = 100,
                                  categories = c("A", "B", "none", "primary"))
  children <- inheritance_children(run_simulation(model, 1, seed = 726))
  expect_true(any(children$birth_time < .5))
  expect_true(any(children$birth_time >= .5))
  expect_identical(children$group, ifelse(children$birth_time < .5, "none", "primary"))
  model$groups$group$inheritance$model[[2]] <- group_inheritance("fixed", fixed_value = "other")
  expect_error(run_simulation(model, .1), "undeclared")
})
