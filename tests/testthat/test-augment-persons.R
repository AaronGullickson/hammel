parent_info_population <- function() {
  persons <- data.frame(person_id = c(30, 10, 70, 50), sex = c("f", "m", "f", "m"),
    birth_time = c(2000.5, 1970, 1975.25, 2002),
    death_time = c(NA, 2003, 2010, NA),
    mother_id = c(70, NA, NA, 70), father_id = c(10, NA, NA, NA),
    group = c(30L, 10L, 20L, 20L), note = c("child", "father", "mother", "half sibling"))
  hammel_population(persons, end_time = 2020, groups = if ("group" %in% names(persons)) "group")
}

test_that("parent information matches IDs and preserves person histories", {
  population <- parent_info_population()
  original <- population
  result <- augment_persons(population)
  expect_s3_class(result, "tbl_df")
  expect_identical(result[names(population$persons)], population$persons)
  expect_identical(result$age_mother, c(25.25, NA_real_, NA_real_, 26.75))
  expect_identical(result$age_father, c(30.5, NA_real_, NA_real_, NA_real_))
  expect_identical(result$group_mother, c(20L, NA_integer_, NA_integer_, 20L))
  expect_identical(result$group_father, c(10L, NA_integer_, NA_integer_, NA_integer_))
  expect_identical(population, original)
  population["end_time"] <- list(NULL)
  expect_error(augment_persons(population), "end_time")
  population$end_time <- original$end_time
  population$persons <- as.data.frame(population$persons)
  expect_identical(augment_persons(population), result)
})

test_that("parent groups are optional and retain their types", {
  population <- parent_info_population()
  population$persons$group <- c("AB", "A", "B", "B")
  population$groups <- list(group = group())
  result <- augment_persons(population)
  expect_identical(result$group_mother, c("B", NA_character_, NA_character_, "B"))
  expect_identical(result$group_father, c("A", NA_character_, NA_character_, NA_character_))
  population$persons$group <- NULL
  population$groups <- NULL
  result <- augment_persons(population)
  expect_named(result, c(names(population$persons), "age_mother", "age_father"))
  expect_identical(result$age_mother, c(25.25, NA_real_, NA_real_, 26.75))
})

test_that("empty populations and unknown parents keep an enriched typed schema", {
  population <- parent_info_population()
  population$persons <- population$persons[FALSE, ]
  population$ancestry <- population$ancestry[FALSE, ]
  result <- augment_persons(population)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_named(result, c(names(population$persons), "age_mother", "age_father",
    "group_mother", "group_father"))
  expect_type(result$age_mother, "double")
  expect_type(result$group_father, "integer")
  empty <- generate_uniform_population(0)
  expect_named(augment_persons(empty), c(names(empty$persons), "age_mother", "age_father"))
  founders <- augment_persons(generate_uniform_population(10))
  expect_true(all(is.na(founders$age_mother)))
  expect_true(all(is.na(founders$age_father)))
})

test_that("parent information rejects invalid input and conflicting columns", {
  population <- parent_info_population()
  expect_error(augment_persons(population$persons), "hammel_population")
  for (column in c("age_mother", "age_father", "group_mother", "group_father")) {
    bad <- population
    bad$persons[[column]] <- 1
    expect_error(augment_persons(bad), paste0("conflict.*", column))
  }
  population$persons$mother_id[1] <- 99
  expect_error(augment_persons(population), "mother_id")
})
