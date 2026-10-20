test_that("shared distances identify direct, sibling, half, and collateral relatives", {
  population <- make_kinship_population()
  cases <- data.frame(relative = c(50, 70, 20, 85, 86, 60, 90, 95, 110),
                      ego_distance = c(1, 1, 2, 1, 1, 2, 2, 2, 0),
                      relative_distance = c(0, 0, 0, 1, 1, 1, 2, 2, 1))
  for (j in seq_len(nrow(cases))) {
    result <- extract_kin(population, 80, relative_ids = cases$relative[j])
    expect_gt(nrow(result), 0)
    expect_true(all(result$ego_distance == cases$ego_distance[j]))
    expect_true(all(result$relative_distance == cases$relative_distance[j]))
  }
  full <- extract_kin(population, 80, relative_ids = 85)
  expect_equal(full$common_ancestor_id, c(50, 70))
  half <- extract_kin(population, 80, relative_ids = 86)
  expect_equal(half$common_ancestor_id, 50)
  cousin <- extract_kin(population, 80, relative_ids = 90)
  expect_equal(cousin$common_ancestor_id, c(10, 20))
  expect_true(all(cousin$ego_lineage == "maternal"))
  expect_true(all(cousin$relative_lineage == "paternal"))
  expect_true(all(cousin$ego_ancestry_complete))
  expect_true(all(!cousin$relative_ancestry_complete))
  population$persons$father_id[population$persons$person_id == 60] <- NA_real_
  half_cousin <- extract_kin(population, 80, relative_ids = 90)
  expect_equal(half_cousin$common_ancestor_id, 10)
  expect_equal(half_cousin$ego_distance, 2L)
  expect_equal(half_cousin$relative_distance, 2L)
})

test_that("the cutoff bounds both distances, including direct ancestry", {
  population <- make_kinship_population()
  expect_equal(nrow(extract_kin(population, 80, relative_ids = 120)), 0)
  removed <- extract_kin(population, 80, generations = 3, relative_ids = 120)
  expect_true(all(removed$ego_distance == 2 & removed$relative_distance == 3))
  expect_gt(nrow(removed), 0)
  expect_equal(nrow(extract_kin(population, 110, relative_ids = 120)), 0)
  second <- extract_kin(population, 110, generations = 3, relative_ids = 120)
  expect_true(all(second$ego_distance == 3 & second$relative_distance == 3))
  expect_gt(nrow(second), 0)
  expect_equal(nrow(extract_kin(population, 110, relative_ids = 20)), 0)
  great <- extract_kin(population, 110, generations = 3, relative_ids = 20)
  expect_equal(great$ego_distance, 3L)
  expect_equal(great$relative_distance, 0L)
  expect_equal(nrow(extract_kin(population, 80, generations = 0)), 0)
})

test_that("missing parents do not hide paternal kin or create artificial shared ancestors", {
  population <- make_kinship_population()
  paternal <- extract_kin(population, 80, relative_ids = 95)
  expect_equal(paternal$common_ancestor_id, c(10, 20))
  expect_true(all(paternal$relative_lineage == "paternal"))
  expect_true(all(!paternal$relative_ancestry_complete))
  expect_equal(nrow(extract_kin(population, 80, relative_ids = 150)), 0)
  # Unrelated co-parents share children, not a common ancestor.
  expect_equal(nrow(extract_kin(population, 50, relative_ids = 70)), 0)
  expect_equal(nrow(extract_kin(population, 10, relative_ids = 30)), 0)
  direct <- extract_kin(population, 80, relative_ids = 70)
  expect_equal(direct$relative_lineage, "self")
  expect_equal(direct$ego_lineage, "paternal")
})

test_that("kinship extraction preserves IDs, query order, deaths, and empty schemas", {
  population <- make_kinship_population()
  population$persons$death_time[population$persons$person_id %in% c(10, 20)] <- -1
  expected <- extract_kin(population, c(80, 90), relative_ids = c(95, 85, 20))
  shuffled <- population
  shuffled$persons <- shuffled$persons[rev(seq_len(nrow(shuffled$persons))), ]
  expect_identical(extract_kin(shuffled, c(80, 90), relative_ids = c(95, 85, 20)), expected)
  expect_true(any(expected$common_ancestor_id == 10))
  expect_true(any(expected$relative_id == 20))
  expect_s3_class(expected, "tbl_df")
  expect_identical(extract_kin(as.data.frame(population$persons), c(80, 90),
    relative_ids = c(95, 85, 20)), expected)
  empty <- extract_kin(population, integer())
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 0)
  expect_identical(vapply(empty, typeof, character(1)), vapply(expected, typeof, character(1)))
  expect_equal(nrow(extract_kin(population, 80, relative_ids = numeric())), 0)
  expect_equal(nrow(extract_kin(create_baseline_population(0), numeric())), 0)
  expect_equal(nrow(extract_kin(population, 80, relative_ids = 80)), 0)
  expect_error(extract_kin(population, 999), "present")
  expect_error(extract_kin(population, c(80, 80)), "unique")
  expect_error(extract_kin(population, 80, relative_ids = 999), "present")
  for (value in list(-1, NA_real_, Inf, 2.5, "2", c(1, 2)))
    expect_error(extract_kin(population, 80, generations = value), "generations")
})

test_that("repeated ancestry is deduplicated with shortest-route lineage", {
  persons <- data.frame(person_id = c(10, 20, 30, 40), sex = c("f", "f", "m", "f"),
    birth_time = c(-80, -50, -50, -20), death_time = NA_real_,
    mother_id = c(NA, 10, 10, 20), father_id = c(NA, NA, NA, 30))
  population <- hammel_population(persons, end_time = 0)
  result <- extract_kin(population, 40, generations = 100, relative_ids = 10)
  expect_equal(nrow(result), 1)
  expect_equal(result$ego_distance, 2L)
  expect_equal(result$ego_lineage, "both")
  population$persons$mother_id[4] <- 10
  shorter <- extract_kin(population, 40, generations = 100, relative_ids = 10)
  expect_equal(shorter$ego_distance, 1L)
  expect_equal(shorter$ego_lineage, "maternal")
})

test_that("ancestry cycles fail and deep traversal is iterative", {
  model <- make_model(2)
  model$population$persons$mother_id <- c(2, 1)
  expect_error(extract_kin(model$population, 1), "cycle")
  expect_error(extract_kin(model$population, 1, generations = 0), "cycle")
  expect_error(simulate(model, 1), "cycle")
  n <- 5000
  ids <- seq_len(n) * 10
  persons <- data.frame(person_id = ids, sex = "f", birth_time = -10,
    death_time = NA_real_, mother_id = c(ids[-1], NA_real_), father_id = NA_real_)
  result <- extract_kin(hammel_population(persons, end_time = 0), ids[1],
                         generations = n, relative_ids = ids[n])
  expect_equal(result$ego_distance, n - 1L)
  expect_equal(result$relative_distance, 0L)
})
