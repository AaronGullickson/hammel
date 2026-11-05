detailed_union_population <- function() {
  persons <- data.frame(person_id = c(30, 10, 70), sex = c("m", "f", "m"),
    birth_time = c(-30, -20, -25), death_time = c(1, NA, NA),
    mother_id = NA_integer_, father_id = NA_integer_, group = c(20L, 10L, 30L))
  unions <- data.frame(union_id = c(9, 4), wife_id = 10, husband_id = c(70, 30),
    start_time = c(1.5, 0), end_time = c(NA, 1), end_reason = c(NA, "death_husband"),
    note = c("current", "earlier"))
  hammel_population(persons, unions, end_time = 2, groups = if ("group" %in% names(persons)) "group")
}

test_that("detailed unions attach spouses by ID and preserve episodes and input", {
  population <- detailed_union_population()
  original <- population
  result <- get_detailed_unions(population)
  expect_s3_class(result, "tbl_df")
  expect_identical(result[names(population$unions)], population$unions)
  expect_equal(result$birth_time_w, c(-20, -20))
  expect_true(all(is.na(result$death_time_w)))
  expect_equal(result$birth_time_h, c(-25, -30))
  expect_equal(result$death_time_h, c(NA, 1))
  expect_identical(result$group_w, c(10L, 10L))
  expect_identical(result$group_h, c(30L, 20L))
  expect_identical(population, original)
  population["end_time"] <- list(NULL)
  expect_identical(get_detailed_unions(population), result)
})

test_that("groups are optional and their labels retain their input types", {
  population <- detailed_union_population()
  population$persons$group <- c("B", "A", "C")
  population$groups <- list(group = group())
  expect_identical(get_detailed_unions(population)$group_h, c("C", "B"))
  population$persons$group <- NULL
  population$groups <- NULL
  expect_false(any(c("group_w", "group_h") %in% names(get_detailed_unions(population))))
})

test_that("empty union histories have the full typed enriched schema", {
  population <- detailed_union_population()
  population$unions <- population$unions[FALSE, ]
  result <- get_detailed_unions(population)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_type(result$group_w, "integer")
  expect_type(result$birth_time_h, "double")
  expect_named(result, c(names(population$unions), "birth_time_w", "death_time_w", "group_w",
                         "birth_time_h", "death_time_h", "group_h"))
  population <- create_baseline_population(0)
  expect_equal(nrow(get_detailed_unions(population)), 0)
  expect_named(get_detailed_unions(population), c(names(population$unions),
    "birth_time_w", "death_time_w", "birth_time_h", "death_time_h"))
})

test_that("detailed unions reject invalid histories and conflicting columns", {
  population <- detailed_union_population()
  expect_error(get_detailed_unions(population$persons), "hammel_population")
  population$unions$birth_time_w <- 100
  expect_error(get_detailed_unions(population), "conflict.*birth_time_w")
  population <- detailed_union_population()
  population$unions$husband_id[1] <- 99
  expect_error(get_detailed_unions(population), "refer to two distinct persons")
})

test_that("partner groups describe formation rather than identity after later moves", {
  population <- detailed_union_population()
  # Current identities: husband 30 = 20, wife 10 = 10, husband 70 = 30.
  # Unsorted moves include a return move and a move exactly at formation.
  population$transitions <- tibble::tibble(person_id = c(10, 70, 10, 30),
    time = c(1.75, 1.5, .5, .75), group_origin = c(40L, 50L, 10L, 60L),
    group_dest = c(10L, 30L, 40L, 20L))
  original <- population
  out <- get_detailed_unions(population)
  expect_identical(out$group_w, c(40L, 10L))
  expect_identical(out$group_h, c(50L, 60L))
  expect_identical(out[names(population$unions)], population$unions)
  expect_identical(population, original)
  # A subsequent restart/additional move cannot alter these earlier pairings.
  population$persons$group[population$persons$person_id == 10] <- 80L
  population$transitions <- rbind(population$transitions,
    tibble::tibble(person_id = 10, time = 1.9, group_origin = 10L, group_dest = 80L))
  expect_identical(get_detailed_unions(population)$group_w, out$group_w)
  population$unions <- population$unions[0, ]
  empty <- get_detailed_unions(population)
  expect_identical(empty$group_w, integer())
  expect_identical(empty$group_h, integer())
})

test_that("formation reconstruction preserves character labels and works without an endpoint", {
  population <- detailed_union_population()
  population$persons$group <- c("B", "A", "C")
  population$groups <- list(group = group())
  population$transitions <- tibble::tibble(person_id = c(10, 70), time = c(1, 1.75),
    group_origin = c("X", "Y"), group_dest = c("A", "C"))
  population["end_time"] <- list(NULL)
  out <- get_detailed_unions(population)
  expect_identical(out$group_w, c("A", "X"))
  expect_identical(out$group_h, c("Y", "B"))
})
