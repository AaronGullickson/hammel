test_that("shared distances identify direct, sibling, half, and collateral relatives", {
  population <- make_kinship_population()
  cases <- data.frame(relative = c(50, 70, 20, 85, 86, 60, 90, 95, 110),
                      ego_distance = c(1, 1, 2, 1, 1, 2, 2, 2, 0),
                      relative_distance = c(0, 0, 0, 1, 1, 1, 2, 2, 1))
  for (j in seq_len(nrow(cases))) {
    result <- extract_kin(population, 80, relative_ids = cases$relative[j])
    expect_equal(nrow(result), 1)
    expect_true(all(result$ego_distance == cases$ego_distance[j]))
    expect_true(all(result$relative_distance == cases$relative_distance[j]))
  }
  full <- extract_kin(population, 80, relative_ids = 85)
  expect_equal(full$n_common_ancestors, 2L)
  expect_equal(full$relationship, "sibling")
  expect_equal(full$ego_lineage, "both")
  half <- extract_kin(population, 80, relative_ids = 86)
  expect_equal(half$n_common_ancestors, 1L)
  expect_equal(half$relationship, "sibling")
  cousin <- extract_kin(population, 80, relative_ids = 90)
  expect_equal(cousin$n_common_ancestors, 2L)
  expect_equal(cousin$relationship, "first_cousin")
  expect_true(all(cousin$ego_lineage == "maternal"))
  expect_true(all(cousin$relative_lineage == "paternal"))
  expect_true(all(cousin$ego_ancestry_complete))
  expect_true(all(!cousin$relative_ancestry_complete))
  population$persons$father_id[population$persons$person_id == 60] <- NA_real_
  half_cousin <- extract_kin(population, 80, relative_ids = 90)
  expect_equal(half_cousin$n_common_ancestors, 1L)
  expect_equal(half_cousin$ego_distance, 2L)
  expect_equal(half_cousin$relative_distance, 2L)
})

test_that("the cutoff bounds both distances, including direct ancestry", {
  population <- make_kinship_population()
  expect_equal(nrow(extract_kin(population, 80, relative_ids = 120)), 0)
  removed <- extract_kin(population, 80, generations = 3, relative_ids = 120)
  expect_true(all(removed$ego_distance == 2 & removed$relative_distance == 3))
  expect_equal(removed$relationship, "first_cousin_once_removed")
  expect_equal(nrow(extract_kin(population, 110, relative_ids = 120)), 0)
  second <- extract_kin(population, 110, generations = 3, relative_ids = 120)
  expect_true(all(second$ego_distance == 3 & second$relative_distance == 3))
  expect_equal(second$relationship, "second_cousin")
  expect_equal(nrow(extract_kin(population, 110, relative_ids = 20)), 0)
  great <- extract_kin(population, 110, generations = 3, relative_ids = 20)
  expect_equal(great$ego_distance, 3L)
  expect_equal(great$relative_distance, 0L)
  expect_equal(great$relationship, "great_grandparent")
  expect_equal(nrow(extract_kin(population, 80, generations = 0)), 0)
})

test_that("missing parents do not hide paternal kin or create artificial shared ancestors", {
  population <- make_kinship_population()
  paternal <- extract_kin(population, 80, relative_ids = 95)
  expect_equal(paternal$n_common_ancestors, 2L)
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
  expect_false(anyDuplicated(expected[c("ego_id", "relative_id")]) > 0)
  expect_true(any(expected$relative_id == 20))
  expect_s3_class(expected, "tbl_df")
  expect_error(
    extract_kin(as.data.frame(population$persons), c(80, 90)),
    "hammel_population"
  )
  empty <- extract_kin(population, integer())
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 0)
  expect_identical(vapply(empty, typeof, character(1)), vapply(expected, typeof, character(1)))
  expect_equal(nrow(extract_kin(population, 80, relative_ids = numeric())), 0)
  expect_equal(nrow(extract_kin(generate_uniform_population(0), numeric())), 0)
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
  population <- hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")
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
  expect_error(run_simulation(model, 1), "cycle")
  n <- 5000
  ids <- seq_len(n) * 10
  persons <- data.frame(person_id = ids, sex = "f", birth_time = -10,
    death_time = NA_real_, mother_id = c(ids[-1], NA_real_), father_id = NA_real_)
  result <- extract_kin(hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group"), ids[1],
                         generations = n, relative_ids = ids[n])
  expect_equal(result$ego_distance, n - 1L)
  expect_equal(result$relative_distance, 0L)
})

living_kin_population <- function(groups = NULL) {
  persons <- data.frame(
    person_id = 1:9, sex = c("f", "m", "f", "m", "f", "m", "f", "f", "m"),
    birth_time = c(-60, -60, -30, -30, -30, -30, -5, -5, -5),
    death_time = c(-10, rep(NA, 8)),
    mother_id = c(NA, NA, 1, 1, NA, NA, 3, 5, 3),
    father_id = c(NA, NA, 2, 2, NA, NA, 6, 4, NA)
  )
  if (!is.null(groups)) {
    persons$race <- groups
  }
  hammel_population(persons, end_time = 0,
    groups = if (!is.null(groups)) "race"
  )
}

test_that("relationship labels follow the generation distances", {
  expect_identical(
    kin_relationship(c(1, 0, 2, 1, 0, 2, 1, 2, 3, 3, 2, 3, 4, 5, 2, 4),
      c(0, 1, 0, 1, 2, 1, 2, 2, 0, 1, 3, 3, 4, 2, 5, 0)),
    c("parent", "child", "grandparent", "sibling", "grandchild",
      "aunt_uncle", "niece_nephew", "first_cousin", "great_grandparent",
      "great_aunt_uncle", "first_cousin_once_removed", "second_cousin",
      "third_cousin", "first_cousin_thrice_removed",
      "first_cousin_thrice_removed", "great_great_grandparent")
  )
  expect_identical(kin_relationship_types(2)$relationship, c(
    "parent", "child", "grandparent", "sibling", "grandchild",
    "aunt_uncle", "niece_nephew", "first_cousin"
  ))
})

test_that("living kin follow the worked example", {
  population <- living_kin_population()
  ego <- extract_kin(population, 7)
  expect_identical(
    ego$relationship[order(ego$relative_id)],
    c("grandparent", "grandparent", "parent", "aunt_uncle", "parent",
      "first_cousin", "sibling")
  )
  expect_identical(
    ego$n_common_ancestors[order(ego$relative_id)],
    c(1L, 1L, 1L, 2L, 1L, 2L, 1L)
  )
  x <- count_kin(population, c(-20, 0))
  expect_s3_class(x, "tbl_df")
  expect_named(x, c("time", "relationship", "persons", "kin", "mean_kin",
    "share_with_kin"))
  early <- x[x$time == -20, ]
  expect_equal(unique(early$persons), 6L)
  expect_equal(early$kin, c(4, 4, 0, 2, 0, 0, 0, 0))
  expect_equal(early$mean_kin, early$kin / 6)
  # By hand at time 0, for each living ego (1 has died).
  late <- x[x$time == 0, ]
  expect_equal(unique(late$persons), 8L)
  expected <- c(parent = 7, child = 7, grandparent = 3, sibling = 4,
    grandchild = 3, aunt_uncle = 3, niece_nephew = 3, first_cousin = 4)
  expect_equal(late$kin, unname(expected[late$relationship]))
  expect_equal(late$share_with_kin[late$relationship == "first_cousin"], 3 / 8)
  # Rows vary by time fastest.
  expect_identical(x$time[1:2], c(-20, 0))
})

test_that("living kin split by the ego's sex and group at each time", {
  # Person 7 moves from A to B at -2, so B is their current category.
  population <- living_kin_population(
    groups = c("A", "B", "A", "B", "A", "B", "B", "B", "A")
  )
  population$transitions <- tibble::tibble(
    person_id = 7, time = -2, origin = "A", destination = "B"
  )
  population <- hammel_population(population$persons,
    end_time = 0, transitions = population$transitions, groups = "race"
  )
  pooled <- count_kin(population, c(-3, 0))
  split <- count_kin(population, c(-3, 0), by = c("sex", "race"))
  expect_named(split, c("time", "sex", "race", "relationship", "persons",
    "kin", "mean_kin", "share_with_kin"))
  sums <- aggregate(cbind(persons, kin) ~ time + relationship, split, sum)
  merged <- merge(sums, pooled, by = c("time", "relationship"))
  expect_equal(merged$persons.x, merged$persons.y)
  expect_equal(merged$kin.x, merged$kin.y)
  # Ego 7, a woman, is in A at -3 and in B at 0: women in A are 3, 5 and 7
  # at -3, then 3 and 5; women in B are 8, then 7 and 8.
  women <- split[split$relationship == "parent" & split$sex == "f", ]
  count <- function(time, race) {
    women$persons[women$time == time & women$race == race]
  }
  expect_equal(c(count(-3, "A"), count(0, "A")), c(3L, 2L))
  expect_equal(c(count(-3, "B"), count(0, "B")), c(1L, 2L))
  empty <- split[split$persons == 0, ]
  expect_true(all(is.na(empty$mean_kin) & is.na(empty$share_with_kin)))
})

test_that("living kin validate inputs", {
  population <- living_kin_population()
  expect_error(count_kin(population$persons, 0), "hammel_population")
  expect_error(count_kin(population, 1), "endpoint")
  expect_error(count_kin(population, numeric()), "times")
  expect_error(count_kin(population, NA_real_), "times")
  expect_error(count_kin(population, 0, generations = -1),
    "generations")
  expect_error(count_kin(population, 0, by = "race"), "by")
  three <- count_kin(population, 0, generations = 3)
  expect_equal(sum(three$relationship == "first_cousin_once_removed"), 1)
  expect_equal(nrow(count_kin(population, 0, generations = 0)), 0)
})

test_that("relationship selection retains full-result rows and denominators", {
  population <- living_kin_population()
  times <- c(0, -20, 0)
  all <- count_kin(population, times)
  for (selected in list("first_cousin", c("first_cousin", "parent"))) {
    actual <- count_kin(population, times, relationships = selected)
    expected <- all[all$relationship %in% selected, ]
    expect_identical(actual, expected)
  }
  cousins <- count_kin(population, 0, relationships = "first_cousin")
  # Half-siblings 7 and 9 share grandparents, but remain siblings, not cousins.
  expect_identical(cousins$kin, 4L)
  expect_identical(cousins$persons, 8L)
  expect_equal(cousins$share_with_kin, 3 / 8)
  set.seed(42)
  rng <- .Random.seed
  count_kin(population, times, relationships = "first_cousin")
  expect_identical(.Random.seed, rng)
})

test_that("selective counts observe group changes and empty strata", {
  population <- living_kin_population(
    groups = c("A", "B", "A", "B", "A", "B", "B", "B", "A")
  )
  population <- hammel_population(
    population$persons,
    end_time = 0,
    groups = "race",
    transitions = data.frame(
      person_id = 7, time = -2, origin = "A", destination = "B"
    )
  )
  all <- count_kin(population, c(-3, -2, 0), by = c("sex", "race"))
  selected <- count_kin(
    population, c(-3, -2, 0),
    by = c("sex", "race"),
    relationships = "first_cousin"
  )
  expect_identical(selected, all[all$relationship == "first_cousin", ])
  women_b <- selected[
    selected$time == -2 & selected$sex == "f" & selected$race == "B",
  ]
  expect_identical(women_b$persons, 2L)
  expect_identical(women_b$kin, 3L)
  expect_identical(women_b$share_with_kin, 1)
  empty <- generate_uniform_population(0, groups = list(race = c(A = 1, B = 1)))
  out <- count_kin(
    empty, 0,
    by = c("sex", "race"),
    relationships = "first_cousin"
  )
  expect_equal(nrow(out), 4L)
  expect_true(all(out$persons == 0L & out$kin == 0L))
  expect_true(all(is.na(out$mean_kin) & is.na(out$share_with_kin)))
  founders <- count_kin(
    generate_uniform_population(10), 0, relationships = "first_cousin"
  )
  expect_identical(founders$kin, 0L)
  expect_identical(founders$share_with_kin, 0)
})

test_that("closest-relationship ties retain the full ancestry bound", {
  persons <- data.frame(
    person_id = 1:4, sex = c("f", "f", "m", "f"),
    birth_time = c(-90, -60, -40, -20), death_time = NA_real_,
    mother_id = c(NA, 1, 2, 1), father_id = c(NA, NA, NA, 3)
  )
  population <- hammel_population(persons, end_time = 0)
  # Person 4 is both sibling and grandchild of 2. At depth two, the (0,2)
  # evidence wins over (1,1) from ego 2; (1,1) wins over (2,0) from ego 4.
  expect_identical(
    count_kin(population, 0, relationships = "sibling")$kin, 1L
  )
  expect_identical(
    count_kin(population, 0, generations = 1, relationships = "sibling")$kin,
    2L
  )
})

test_that("counts agree with extracted histories on a synthetic pedigree", {
  set.seed(17)
  size <- 20L
  cohorts <- 6L
  persons <- data.frame(
    person_id = seq_len(size * cohorts),
    sex = rep(c("f", "m"), size * cohorts / 2),
    birth_time = rep(seq(-60, -10, by = 10), each = size),
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
  )
  for (cohort in 2:cohorts) {
    rows <- seq_len(size) + size * (cohort - 1L)
    parents <- rows - size
    persons$mother_id[rows] <- sample(parents[persons$sex[parents] == "f"],
      size, replace = TRUE
    )
    persons$father_id[rows] <- sample(parents[persons$sex[parents] == "m"],
      size, replace = TRUE
    )
  }
  persons$death_time[c(1, 4, 25, 52)] <- -5
  population <- hammel_population(persons, end_time = 0)
  for (depth in 1:3) {
    kin <- extract_kin(population, persons$person_id, generations = depth)
    counts <- count_kin(population, c(0, -20, -5, -20),
      generations = depth, by = "sex"
    )
    for (i in seq_len(nrow(counts))) {
      row <- counts[i, ]
      alive <- persons$birth_time <= row$time &
        (is.na(persons$death_time) | persons$death_time > row$time)
      egos <- persons$person_id[alive & persons$sex == row$sex]
      relatives <- persons$person_id[alive]
      pairs <- kin[
        kin$ego_id %in% egos & kin$relative_id %in% relatives &
          kin$relationship == row$relationship,
      ]
      expect_equal(row$persons, length(egos))
      expect_equal(row$kin, nrow(pairs))
      expect_equal(row$mean_kin, nrow(pairs) / length(egos))
      expect_equal(
        row$share_with_kin, length(unique(pairs$ego_id)) / length(egos)
      )
    }
    chosen <- unique(kin$relationship[grepl("cousin", kin$relationship)])
    if (length(chosen)) {
      selected <- count_kin(population, c(0, -20, -5, -20),
        generations = depth, by = "sex", relationships = chosen
      )
      expect_identical(selected, counts[counts$relationship %in% chosen, ])
    }
  }
})

test_that("relationship selection validates labels and retains cycle checks", {
  population <- living_kin_population()
  bad <- list(
    character(), NA_character_, "", "cousin", "second_cousin", 1, TRUE,
    c("parent", "parent"), matrix("parent")
  )
  for (selected in bad) {
    expect_error(count_kin(population, 0, relationships = selected),
      "relationships"
    )
  }
  expect_error(
    count_kin(population, 0, generations = 0, relationships = "parent"),
    "relationships"
  )
  expect_identical(nrow(count_kin(population, 0, generations = 0)), 0L)
  dead <- make_model(2)$population
  dead$persons$mother_id <- c(2, 1)
  dead$persons$death_time <- -10
  expect_error(count_kin(dead, 0, relationships = "first_cousin"), "cycle")
  expect_error(count_kin(dead, 0, generations = 0), "cycle")
})
