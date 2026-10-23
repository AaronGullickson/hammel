inheritance_test_model <- function(rule = create_group_inheritance(), birth_rate = 20) {
  model <- make_model(2, birth_rate = birth_rate, birth_interval = 0)
  model$population$persons$sex <- c("f", "m")
  model$population$persons$group <- c("A", "B")
  set_fertility(model, model$processes$fertility$rates, 0, inheritance = rule)
}

inheritance_children <- function(result) result$persons[result$persons$birth_time >= 0, ]

test_that("inheritance method strings are equivalent to constructed models", {
  for (method in c("mother", "father", "random_parent", "same_sex", "opposite_sex", "ancestry")) {
    model <- inheritance_test_model(method)
    explicit <- inheritance_test_model(create_group_inheritance(method))
    expect_identical(model$processes$fertility$inheritance,
                     explicit$processes$fertility$inheritance)
    expect_identical(simulate(model, 1, seed = 52), simulate(explicit, 1, seed = 52))
    model$processes$fertility$inheritance <- method
    expect_identical(set_process(model, model$processes$fertility)$processes$fertility$inheritance,
                     explicit$processes$fertility$inheritance)
    expect_identical(simulate(model, 1, seed = 52), simulate(explicit, 1, seed = 52))
  }
  for (method in list("other", NA_character_, c("mother", "father"), character()))
    expect_error(inheritance_test_model(method), "method")
  expect_error(inheritance_test_model("hierarchy"), "order")
  expect_error(inheritance_test_model("custom"), "fun")
})

test_that("custom inheritance receives parents' union status at birth", {
  rule <- create_group_inheritance("custom", groups = c("partnered", "nonmarital"),
    fun = function(parents_in_union, ...) {
      stopifnot(is.logical(parents_in_union), length(parents_in_union) == 1L,
                !is.na(parents_in_union))
      if (parents_in_union) "partnered" else "nonmarital"
    })
  model <- inheritance_test_model(rule)
  result <- simulate(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 0)
  expect_true(all(children$group == "partnered"))
  expect_gt(nrow(result$unions), 0)

  model <- set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0, inheritance = rule)
  model <- set_union_formation(model, "rate", rates = data.frame(
    age_start = 0, age_end = Inf, rate = 0), cancel_birth_if_no_match = FALSE)
  for (known_father in c(TRUE, FALSE)) {
    unpartnered <- model
    if (!known_father) unpartnered$population$persons <- model$population$persons[1, ]
    result <- simulate(unpartnered, 1, seed = 52)
    children <- inheritance_children(result)
    expect_gt(nrow(children), 0)
    expect_true(all(children$group == "nonmarital"))
    expect_true(all(is.na(children$father_id) == !known_father))
    expect_equal(nrow(result$unions), 0)
  }
})

test_that("built-in inheritance uses parents and newborn sex", {
  for (method in c("mother", "father", "same_sex", "opposite_sex")) {
    result <- simulate(inheritance_test_model(create_group_inheritance(method)), 1, seed = 52)
    children <- inheritance_children(result)
    expect_gt(nrow(children), 5)
    expect_true(all(c("m", "f") %in% children$sex))
    expected <- switch(method, mother = rep("A", nrow(children)), father = rep("B", nrow(children)),
      same_sex = ifelse(children$sex == "f", "A", "B"),
      opposite_sex = ifelse(children$sex == "f", "B", "A"))
    expect_identical(children$group, expected)
    shares <- result$ancestry[result$persons$birth_time >= 0, ]
    expect_true(all(shares$group_A == 0.5 & shares$group_B == 0.5))
  }
})

test_that("default, explicit, and older maternal specifications preserve histories and RNG", {
  model <- inheritance_test_model()
  set.seed(654)
  original <- simulate(model, 1)
  after <- .Random.seed
  expect_identical(model$processes$fertility$inheritance$method, "mother")
  model$processes$fertility$inheritance <- NULL
  set.seed(654)
  expect_identical(simulate(model, 1), original)
  expect_identical(.Random.seed, after)
  plain <- make_model(2)
  expect_no_error(simulate(plain, 0))
  plain$processes$fertility$inheritance <- create_group_inheritance("father")
  expect_error(simulate(plain, 0), "group identities")
})

test_that("missing fathers fall back to maternal identity without extra parent-choice draws", {
  original <- inheritance_test_model()
  original$population$persons <- original$population$persons[1, ]
  reference <- simulate(original, 1, seed = 52)
  for (method in c("father", "random_parent", "same_sex", "opposite_sex", "hierarchy")) {
    rule <- if (method == "hierarchy") create_group_inheritance(method, order = c("B", "A"), passing = 0.5) else
      create_group_inheritance(method)
    model <- original
    model$processes$fertility$inheritance <- rule
    expect_identical(simulate(model, 1, seed = 52), reference)
  }
})

test_that("random-parent and hierarchy probabilities agree with analytical cases", {
  model <- inheritance_test_model(create_group_inheritance("random_parent"), birth_rate = 2000)
  result <- simulate(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 1000)
  expect_lt(abs(mean(children$group == "A") - 0.5), 5 * sqrt(0.25/nrow(children)))
  for (passing in c(0, 1, 0.05)) {
    model$processes$fertility$inheritance <- create_group_inheritance("hierarchy", order = c("A", "B"), passing = passing)
    children <- inheritance_children(simulate(model, 1, seed = 52))
    if (passing %in% c(0, 1)) expect_true(all(children$group == if (passing == 0) "A" else "B")) else
      expect_lt(abs(mean(children$group == "B") - passing), 5 * sqrt(passing*(1-passing)/nrow(children)))
  }
  model <- inheritance_test_model(create_group_inheritance("hierarchy", order = "B", passing = 1))
  expect_identical(simulate(model, 1, seed = 52), simulate(inheritance_test_model(), 1, seed = 52))
})

test_that("ancestry assignment conditions on known shares and handles all-unknown ancestry", {
  model <- inheritance_test_model(create_group_inheritance("ancestry"), birth_rate = 2000)
  # Both parents identify A, but carry a known B reference component.
  p <- model$population$persons
  p$group <- "A"
  a <- data.frame(person_id = p$person_id, group_A = 0.35, group_B = 0.15, unknown = 0.5)
  model <- set_population(model, hammel_population(p, end_time = 0, ancestry = a))
  result <- simulate(model, 1, seed = 52)
  children <- inheritance_children(result)
  expect_lt(abs(mean(children$group == "A") - 0.7), 5 * sqrt(0.7*0.3/nrow(children)))
  expect_true(all(result$ancestry$unknown == 0.5))
  model$population$ancestry$group_A <- 0
  model$population$ancestry$group_B <- 0
  model$population$ancestry$unknown <- 1
  model$processes$fertility$rates$rate <- 20
  maternal <- model; maternal$processes$fertility$inheritance <- create_group_inheritance()
  expect_identical(simulate(model, 1, seed = 52), simulate(maternal, 1, seed = 52))
})

test_that("custom callbacks receive copies of all inputs and use the unified R RNG stream", {
  custom <- create_group_inheritance("custom", fun = function(mother_group, father_group,
      mother_ancestry, father_ancestry, child_ancestry, child_sex) {
    stopifnot(mother_group == "A", father_group == "B", child_sex %in% c("m", "f"),
      identical(names(child_ancestry), c("group_A", "group_B", "unknown")),
      mother_ancestry[["group_A"]] == 1, father_ancestry[["group_B"]] == 1,
      child_ancestry[["group_A"]] == 0.5)
    mother_ancestry[] <- 0 # Copied values cannot alter engine ancestry.
    if (runif(1) < 0.5) mother_group else father_group
  })
  model <- inheritance_test_model(custom)
  set.seed(654); state <- .Random.seed
  actual <- simulate(model, 1, seed = 52)
  expect_identical(.Random.seed, state)
  native <- inheritance_test_model(create_group_inheritance("random_parent"))
  expect_identical(actual, simulate(native, 1, seed = 52))
  set.seed(52)
  actual <- simulate(model, 1)
  after <- .Random.seed
  set.seed(52)
  expect_identical(simulate(native, 1), actual)
  expect_identical(.Random.seed, after)
  # Multinomial-logit-style ancestry probabilities can be implemented in R.
  rule <- create_group_inheritance("custom", fun = function(child_ancestry, ...) {
    weights <- exp(c(A = child_ancestry[["group_A"]], B = child_ancestry[["group_B"]]))
    sample(names(weights), 1, prob = weights)
  })
  expect_no_error(simulate(inheritance_test_model(rule), 1, seed = 52))
})

test_that("custom mixed identities are declared, use their own schedules, and retain A/B ancestry", {
  rule <- create_group_inheritance("custom", groups = "AB", fun = function(mother_group, father_group, ...) {
    if (is.na(father_group) || mother_group == father_group) mother_group else "AB"
  })
  model <- inheritance_test_model(rule)
  # AB infants die immediately; A/B founders survive. No AB ancestry component.
  rates <- data.frame(sex = rep(c("m", "f"), each = 2), age_start = c(0, 1),
    age_end = c(1, Inf), rate = 0)
  rates <- rbind(transform(rates, group = "A"), transform(rates, group = "B"),
                 transform(rates, group = "AB", rate = ifelse(age_start == 0, 1e8, 0)))
  model <- set_mortality(model, rates)
  result <- simulate(model, 0.5, seed = 52)
  children <- inheritance_children(result)
  expect_gt(nrow(children), 0)
  expect_true(all(children$group == "AB"))
  expect_true(all(!is.na(children$death_time)))
  expect_false("group_AB" %in% names(result$ancestry))
  expect_true(all(result$ancestry$group_A[result$persons$birth_time >= 0] == 0.5))
  continued <- simulate(set_population(model, result), 0.5, seed = 53)
  expect_identical(continued$ancestry[seq_len(nrow(result$ancestry)), ], result$ancestry)
})

test_that("future identities are checked against rates and categorical matching before execution", {
  model <- inheritance_test_model(create_group_inheritance("custom", groups = "AB", fun = function(...) "AB"))
  for (name in c("mortality", "fertility")) {
    edited <- model
    rates <- edited$processes[[name]]$rates
    edited$processes[[name]]$rates <- rbind(transform(rates, group = "A"), transform(rates, group = "B"))
    expect_error(simulate(edited, 0), "missing schedules.*AB")
  }
  matching <- create_matching_model(~ male_group, c(male_groupB = 0), group_levels = c("A", "B"), matching_method = "female_choice")
  model$processes$union_formation$matching <- matching
  expect_error(simulate(model, 0), "covered.*group_levels")
  # Positive reference shares also infer identities; zero-only columns do not.
  model <- inheritance_test_model()
  model <- set_population(model, hammel_population(model$population$persons, end_time = 0))
  model$population$ancestry$group_C <- 0
  model$population$ancestry$group_D <- 0
  prepared <- prepare_group_inheritance(model$processes$fertility$inheritance, model$population)
  expect_identical(prepared$values, c("A", "B"))
  model$population$ancestry$group_A[1] <- 0.5
  model$population$ancestry$group_C[1] <- 0.5
  expect_identical(prepare_group_inheritance(model$processes$fertility$inheritance, model$population)$values, c("A", "B", "C"))
})

test_that("custom errors and invalid returns restore explicit RNG state", {
  for (fun in list(function(...) NA_character_, function(...) NULL, function(...) c("A", "B"),
                   function(...) "undeclared", function(...) stop("custom failed"))) {
    model <- inheritance_test_model(create_group_inheritance("custom", fun = fun))
    set.seed(654); state <- .Random.seed
    expect_error(simulate(model, 1, seed = 52))
    expect_identical(.Random.seed, state)
  }
  rule <- create_group_inheritance("custom", fun = function(mother_group, father_group, father_ancestry, ...) {
    stopifnot(is.na(father_group), is.null(father_ancestry)); mother_group
  })
  model <- inheritance_test_model(rule)
  model$population$persons <- model$population$persons[1, ]
  expect_true(all(inheritance_children(simulate(model, 1, seed = 52))$group == "A"))
})

test_that("inheritance validates edited specifications and configuration", {
  expect_error(create_group_inheritance("other"), "method")
  expect_error(create_group_inheritance("custom"), "fun")
  expect_error(create_group_inheritance("hierarchy"), "order")
  expect_error(create_group_inheritance("hierarchy", order = c("A", "A")), "unique")
  expect_error(create_group_inheritance("mother", passing = 0.1), "hierarchy")
  for (passing in list(-1, 2, NA, "0", c(0, 1)))
    expect_error(create_group_inheritance("hierarchy", order = c("A", "B"), passing = passing), "passing")
  expect_error(create_group_inheritance(groups = c("A", NA)), "group")
  model <- inheritance_test_model()
  model$processes$fertility$inheritance$method <- "bad"
  expect_error(set_process(model, model$processes$fertility), "method")
  expect_error(simulate(model, 0), "method")
})

test_that("numeric identities remain numeric, and empty populations stay typed", {
  model <- inheritance_test_model(create_group_inheritance("father"))
  model$population$persons$group <- c(20L, 10L)
  result <- simulate(model, 0.5, seed = 52)
  expect_type(result$persons$group, "integer")
  expect_true(all(inheritance_children(result)$group == 10L))
  model$population$persons <- model$population$persons[FALSE, ]
  result <- simulate(model, 1, seed = 52)
  expect_equal(nrow(result$persons), 0)
  expect_type(result$persons$group, "integer")
  expect_s3_class(result$ancestry, "tbl_df")
})

test_that("new AB identities are available to categorical partner matching during the run", {
  rule <- create_group_inheritance("custom", groups = "AB", fun = function(...) "AB")
  model <- inheritance_test_model(rule)
  model <- set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0, inheritance = rule)
  matching <- create_matching_model(~ male_group,
    c(male_groupB = 0, male_groupAB = 1000), group_levels = c("A", "B", "AB"), matching_method = "female_choice")
  model <- set_union_formation(model, "rate", matching, rates = data.frame(
    age_start = 0, age_end = Inf, rate = 0), min_male_age = 0, kinship_generations = 0,
    cancel_birth_if_no_match = FALSE)
  result <- simulate(model, 1, seed = 52)
  children <- inheritance_children(result)
  father_groups <- result$persons$group[match(children$father_id, result$persons$person_id)]
  expect_true(any(father_groups == "AB"))
  expect_true(all(children$group == "AB"))
})

test_that("future search schedules, group multipliers, and numeric declarations are supported", {
  rule <- create_group_inheritance("custom", groups = "AB", fun = function(...) "AB")
  model <- inheritance_test_model(rule)
  model <- set_fertility(model, data.frame(union_status = c("partnered", "unpartnered"),
    age_start = 15, age_end = 50, rate = c(0, 20)), 0, inheritance = rule)
  model <- set_union_formation(model, "rate", rates = data.frame(group = c("A", "B"),
    age_start = 0, age_end = Inf, rate = 0))
  expect_error(simulate(model, 0), "missing schedules.*AB")
  model <- inheritance_test_model(rule)
  model$processes$mortality$rates$rate <- 1e8
  model$processes$mortality$group_prop_hazard <- c(A = 0, B = 0, AB = 1)
  result <- simulate(model, 0.5, seed = 52)
  expect_true(all(!is.na(inheritance_children(result)$death_time)))
  model <- inheritance_test_model(create_group_inheritance("custom", groups = 30L, fun = function(...) 30L))
  model$population$persons$group <- c(20L, 10L)
  result <- simulate(model, 0.5, seed = 52)
  expect_type(result$persons$group, "integer")
  expect_true(all(inheritance_children(result)$group == 30L))
})
