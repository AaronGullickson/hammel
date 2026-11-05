ancestry_population <- function() {
  p <- make_model(2)$population
  p$persons$person_id <- c(10, 30)
  p$persons$sex <- c("f", "m")
  p$persons$group <- c("A", "B")
  p$groups <- list(group = group(ancestry = TRUE))
  hammel_population(p$persons, end_time = 0, groups = list(group = group(ancestry = TRUE)))
}

ancestry_birth_model <- function(population = ancestry_population()) {
  make_model(0, birth_rate = 20, birth_interval = 0) |>
    set_population(population)
}

test_that("grouped constructors initialize all starting people from identity", {
  p <- ancestry_population()
  expect_s3_class(p$ancestry, "tbl_df")
  expect_named(p$ancestry, c("person_id", "group_A", "group_B", "unknown"))
  expect_equal(p$ancestry$person_id, c(10, 30))
  expect_equal(p$ancestry$group_A, c(1, 0))
  expect_equal(p$ancestry$group_B, c(0, 1))
  expect_equal(p$ancestry$unknown, c(0, 0))
  # Initial identities are the reference, even where earlier parents are known.
  p$persons$mother_id[2] <- 10
  rebuilt <- hammel_population(p$persons, end_time = 0, groups = list(group = group(ancestry = TRUE)))
  expect_identical(rebuilt$ancestry, p$ancestry)
  ungrouped <- make_model(2)$population
  expect_null(ungrouped$ancestry)
  expect_null(simulate(make_model(2), 1, seed = 123)$ancestry)
})

test_that("supplied shares and ancestral groups survive normalization and restarts", {
  p <- ancestry_population()
  a <- data.frame(person_id = c(30, 10), unknown = c(0, 0.25),
                  group_B = c(1, 0.5), group_C = c(0, 0.25))
  p <- hammel_population(p$persons, end_time = 0, ancestry = a, groups = list(group = group(ancestry = TRUE)))
  expect_s3_class(p$ancestry, "tbl_df")
  expect_named(p$ancestry, c("person_id", "group_B", "group_C", "unknown"))
  expect_equal(p$ancestry$person_id, c(10, 30))
  expect_equal(p$ancestry$group_C, c(0.25, 0))
  result <- simulate(ancestry_birth_model(p), 0.3, seed = 654)
  expect_identical(result$ancestry[1:2, ], p$ancestry)
  expect_gt(nrow(result$ancestry), 2)
  expect_equal(result$ancestry$group_C[-(1:2)], rep(0.125, nrow(result$ancestry) - 2))
  expect_equal(result$ancestry$unknown[-(1:2)], rep(0.125, nrow(result$ancestry) - 2))
  continued <- simulate(set_population(ancestry_birth_model(p), result), 0.3, seed = 123)
  expect_identical(continued$ancestry[seq_len(nrow(result$ancestry)), ], result$ancestry)
  expect_equal(rowSums(as.matrix(continued$ancestry[-1])), rep(1, nrow(continued$ancestry)))
})

test_that("cross-group newborns average parental ancestry under explicit maternal inheritance", {
  model <- ancestry_birth_model()
  model <- set_inheritance(model, create_group_inheritance("mother"))
  result <- simulate(model, 0.3, seed = 654)
  children <- result$persons$birth_time >= 0
  expect_true(any(children))
  expect_true(all(result$persons$group[children] == "A"))
  expect_equal(result$ancestry$group_A[children], rep(0.5, sum(children)))
  expect_equal(result$ancestry$group_B[children], rep(0.5, sum(children)))
  expect_equal(result$ancestry$unknown[children], rep(0, sum(children)))
  expect_identical(result$ancestry$person_id, result$persons$person_id)
})

test_that("missing fathers contribute an unknown half and known unknown shares transmit", {
  p <- ancestry_population()
  p$persons <- p$persons[1, ]
  p$ancestry <- p$ancestry[1, ]
  result <- simulate(ancestry_birth_model(p), 0.3, seed = 654)
  children <- result$persons$birth_time >= 0
  expect_true(any(children))
  expect_true(all(is.na(result$persons$father_id[children])))
  expect_equal(result$ancestry$group_A[children], rep(0.5, sum(children)))
  expect_equal(result$ancestry$unknown[children], rep(0.5, sum(children)))
  p <- ancestry_population()
  p$ancestry$group_A[1] <- 0.5
  p$ancestry$unknown[1] <- 0.5
  result <- simulate(ancestry_birth_model(p), 0.3, seed = 654)
  children <- result$persons$birth_time >= 0
  expect_equal(result$ancestry$group_A[children], rep(0.25, sum(children)))
  expect_equal(result$ancestry$group_B[children], rep(0.5, sum(children)))
  expect_equal(result$ancestry$unknown[children], rep(0.25, sum(children)))
})

test_that("multiple generations inherit the live ancestry state at birth", {
  p <- ancestry_population()
  model <- ancestry_birth_model(p) |>
    set_fertility(data.frame(union_type = c("partnered", "unpartnered"),
      age_start = 0, age_end = 100, rate = c(0, 2)), birth_interval = 0) |>
    set_union_formation("rate", rates = data.frame(age_start = 0, age_end = Inf, rate = 0),
      kinship_generations = 0, cancel_birth_if_no_match = FALSE)
  result <- simulate(model, 5, seed = 654)
  a <- as.matrix(result$ancestry[-1])
  children <- which(result$persons$birth_time >= 0)
  expect_true(any(result$persons$mother_id[children] > 30))
  for (child in children) {
    mother <- match(result$persons$mother_id[child], result$persons$person_id)
    father <- match(result$persons$father_id[child], result$persons$person_id)
    expect_false(is.na(father))
    expect_equal(a[child, ], 0.5 * a[mother, ] + 0.5 * a[father, ])
  }
  expect_true(any(result$ancestry$group_A == 0.25 & result$ancestry$group_B == 0.75))
  expect_equal(rowSums(a), rep(1, nrow(a)))
})

test_that("supplied ancestry needs complete IDs and bounded unit-sum shares", {
  p <- ancestry_population()
  a <- p$ancestry
  bad <- a[-1, ]
  expect_error(hammel_population(p$persons, end_time = 0, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "one row")
  bad <- a; bad$person_id[2] <- bad$person_id[1]
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "one row")
  bad <- a; bad$person_id[2] <- 99
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "one row")
  for (value in c(NA_real_, Inf, -0.1, 1.1)) {
    bad <- a; bad$group_A[1] <- value
    expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "proportions")
  }
  bad <- a; bad$group_A[1] <- 0.9
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "sum to 1")
  bad <- a; bad$unknown <- NULL
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "unknown")
  bad <- a; names(bad)[2] <- "A"
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "group_<label>")
  bad <- a; names(bad)[2] <- "group_"
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "nonempty labels")
  bad <- a; bad$group_A <- as.character(bad$group_A)
  expect_error(hammel_population(p$persons, ancestry = bad, groups = list(group = group(ancestry = TRUE))), "proportions")
  p$ancestry$unknown[1] <- NA_real_
  expect_error(set_population(hammel_model(), p), "proportions")
  model <- ancestry_birth_model()
  model$population$ancestry$unknown[1] <- 0.5
  expect_error(simulate(model, 0), "sum to 1")
})

test_that("automatic ancestry supports numeric labels, collisions, and empty groups", {
  p <- ancestry_population()
  p$persons$group <- c(20L, 10L)
  p$groups <- list(group = group(ancestry = TRUE))
  p <- hammel_population(p$persons, end_time = 0, groups = list(group = group(ancestry = TRUE)))
  expect_named(p$ancestry, c("person_id", "group_20", "group_10", "unknown"))
  p$persons$group <- c("unknown", "person_id")
  p$groups <- list(group = group(ancestry = TRUE))
  p <- hammel_population(p$persons, end_time = 0, groups = list(group = group(ancestry = TRUE)))
  expect_named(p$ancestry, c("person_id", "group_unknown", "group_person_id", "unknown"))
  empty <- hammel_population(p$persons[FALSE, ], end_time = 0, groups = list(group = group(ancestry = TRUE)))
  expect_s3_class(empty$ancestry, "tbl_df")
  expect_named(empty$ancestry, c("person_id", "unknown"))
  expect_equal(nrow(empty$ancestry), 0)
  expect_s3_class(simulate(ancestry_birth_model(empty), 1)$ancestry, "tbl_df")
  expect_identical(simulate(ancestry_birth_model(p), 0)$ancestry, p$ancestry)
})

test_that("ancestry belongs to a declared group dimension", {
  p <- ancestry_population()
  # Ancestry requires its dimension's column and an ancestry = TRUE declaration.
  expect_error(hammel_population(p$persons[setdiff(names(p$persons), "group")], end_time = 0,
    ancestry = p$ancestry, groups = list(group = group(ancestry = TRUE))), "not a persons column")
  expect_error(hammel_population(p$persons, end_time = 0, ancestry = p$ancestry, groups = "group"),
               "ancestry = TRUE")
  expect_error(hammel_population(p$persons, end_time = 0, ancestry = p$ancestry), "declared group")
  unknown <- data.frame(person_id = p$persons$person_id, unknown = 1)
  p <- hammel_population(p$persons, end_time = 0, ancestry = unknown, groups = list(group = group(ancestry = TRUE)))
  expect_true(all(simulate(ancestry_birth_model(p), 0.3, seed = 654)$ancestry$unknown == 1))
})

test_that("ancestry is deterministic, independent, and retained after deaths", {
  model <- ancestry_birth_model()
  set.seed(321)
  state <- .Random.seed
  a <- simulate(model, 0.3, seed = 654)
  expect_identical(.Random.seed, state)
  expect_identical(a, simulate(model, 0.3, seed = 654))
  model$population$ancestry$group_A[1] <- 0
  model$population$ancestry$unknown[1] <- 1
  b <- simulate(model, 0.3, seed = 654)
  expect_identical(a$persons, b$persons)
  expect_identical(a$unions, b$unions)
  expect_true(all(a$ancestry$unknown == 0))
  model <- set_mortality(ancestry_birth_model(), data.frame(sex = c("m", "f"),
    age_start = 0, age_end = Inf, rate = 1e8))
  result <- simulate(model, 1, seed = 654)
  expect_true(all(!is.na(result$persons$death_time)))
  expect_identical(result$ancestry, model$population$ancestry)
})

test_that("integer shares and cancelled births retain valid ancestry rows", {
  p <- ancestry_population()
  for (column in setdiff(names(p$ancestry), "person_id")) p$ancestry[[column]] <- as.integer(p$ancestry[[column]])
  result <- simulate(ancestry_birth_model(p), 0.3, seed = 654)
  expect_equal(as.matrix(result$ancestry[1:2, -1]), as.matrix(p$ancestry[-1]))
  p <- ancestry_population()
  p$persons <- p$persons[1, ]
  p$ancestry <- p$ancestry[1, ]
  model <- set_union_formation(ancestry_birth_model(p), "on_demand")
  result <- simulate(model, 0.3, seed = 654)
  expect_equal(result$simulation$births, 0)
  expect_identical(result$ancestry, p$ancestry)
})

test_that("directly edited grouped populations initialize ancestry at attachment and execution", {
  model <- make_model(2)
  model$population$persons$group <- c("A", "B")
  model$population$groups <- list(group = group(ancestry = TRUE))
  expect_null(model$population$ancestry)
  attached <- set_population(model, model$population)
  expect_named(attached$population$ancestry, c("person_id", "group_A", "group_B", "unknown"))
  result <- simulate(model, 0)
  expect_identical(result$ancestry, attached$population$ancestry)
  expect_null(model$population$ancestry)
})
