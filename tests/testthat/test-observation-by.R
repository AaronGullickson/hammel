# One woman: cohabits 20-22, converts to marriage, dissolves it at 30 while
# giving birth, gives birth again at 31, cohabits again from 33 until her
# partner dies at 38, and moves from group A to B at 35.
by_history <- function() {
  persons <- tibble::tibble(person_id = 1:5, sex = c("f", "m", "m", "f", "m"),
    birth_time = c(0, -2, -5, 30, 31), death_time = c(NA, NA, 38, NA, NA),
    mother_id = c(NA, NA, NA, 1L, 1L), father_id = c(NA, NA, NA, 2L, NA),
    group = c("B", "A", "A", "A", "A"))
  unions <- tibble::tibble(union_id = 1:3, wife_id = 1L, husband_id = c(2L, 2L, 3L),
    start_time = c(20, 22, 33), end_time = c(22, 30, 38),
    end_reason = c("conversion", "dissolution", "death_husband"),
    union_type = c("cohab", "marriage", "cohab"), previous_cohab_id = c(NA, 1L, NA))
  transitions <- tibble::tibble(person_id = 1L, time = 35, origin = "A", destination = "B")
  hammel_population(persons, unions, end_time = 40, transitions = transitions, groups = if ("group" %in% names(persons)) "group")
}

woman <- function(x) x[x$sex == "f", ]

test_that("mortality exposure splits by union type and previous union", {
  pop <- by_history()
  # Female exposure: the woman (0-40) and her daughter (born 30, unpartnered).
  x <- woman(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = "union_type"))
  expect_identical(x$union_type, c("unpartnered", "cohab", "marriage"))
  expect_equal(x$person_years, c(20 + 3 + 2 + 10, 2 + 5, 8))
  x <- woman(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = "previous_union"))
  expect_identical(x$previous_union, c("none", "separated", "widowed"))
  # The woman is none to 30 and separated after: a cohabitation ended by death
  # keeps the legal category under the default policy.
  expect_equal(x$person_years, c(30 + 10, 10, 0))
})

test_that("previous union follows the recorded cohabitation policy", {
  pop <- by_history()
  pop$simulation$cohab_as_previous_union <- TRUE
  x <- woman(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = "previous_union"))
  # Now the cohabitation ended by death at 38 makes her widowed.
  expect_equal(x$person_years, c(40, 8, 2))
  pooled <- woman(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = "previous_union",
                                            use_detailed_previous_union = FALSE))
  expect_identical(pooled$previous_union, c("none", "dissolved"))
  expect_equal(pooled$person_years, c(40, 10))
})

test_that("births use the state before them, with birth before dissolution", {
  pop <- by_history()
  x <- calculate_fertility_rates(pop, c(0, 40), c(0, Inf), by = c("union_type", "previous_union"))
  births <- x[x$births > 0, c("union_type", "previous_union", "births")]
  # The birth at 30 precedes the dissolution; the birth at 31 is after it.
  expect_equal(as.data.frame(births), data.frame(union_type = c("marriage", "unpartnered"),
    previous_union = c("none", "separated"), births = c(1, 1)), ignore_attr = TRUE)
  expect_equal(sum(x$person_years), 40 + 10)
})

test_that("population size counts states after events at each time", {
  pop <- by_history()
  x <- calculate_population_size(pop, c(25, 30, 38), by = c("sex", "union_type", "previous_union"))
  w <- x[x$sex == "f" & x$population > 0, ]
  # At 30 the dissolution has happened and the daughter is born; at 38 the
  # partner's death ends the cohabitation.
  # The daughter (from 30) and the woman; the partner's death at 38 ends her
  # cohabitation. Rows vary by time, then the by variables in order.
  expect_equal(as.data.frame(w[c("time", "union_type", "previous_union", "population")]),
    data.frame(time = c(30, 38, 25, 30, 38),
      union_type = c("unpartnered", "unpartnered", "marriage", "unpartnered", "unpartnered"),
      previous_union = c("none", "none", "none", "separated", "separated"),
      population = 1L), ignore_attr = TRUE)
  expect_equal(sum(x$population[x$time == 30]), 4)
})

test_that("transitions split by previous union and union type", {
  pop <- by_history()
  x <- calculate_transition_rates(pop, c(0, 40), age_breaks = c(0, Inf),
                                  by = c("previous_union", "union_type"))
  moved <- x[x$transitions > 0, ]
  expect_equal(nrow(moved), 1)
  expect_identical(c(moved$origin, moved$destination, moved$previous_union, moved$union_type),
                   c("A", "B", "separated", "cohab"))
})

test_that("by is validated and rows follow the documented order", {
  pop <- by_history()
  expect_error(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = "sex"),
               "cannot split by \"sex\".*accepts \"group\", \"union_type\", \"previous_union\"")
  expect_error(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), by = c("group", "group")), "distinct")
  expect_error(
    calculate_mortality_rates(pop$persons, c(0, 40), c(0, Inf),
      by = "union_type"
    ),
    "must be a hammel_population"
  )
  expect_error(calculate_mortality_rates(pop, c(0, 40), c(0, Inf), use_detailed_previous_union = FALSE), "previous_union")
  x <- calculate_mortality_rates(pop, c(0, 20, 40), c(0, 30, Inf), by = c("union_type", "group"))
  n <- 2 * 2 * 2
  expect_identical(x$time_start[1:2], c(0, 20))
  expect_identical(x$age_start[1:4], c(0, 0, 30, 30))
  expect_identical(x$sex[1:n], rep(c("m", "f"), each = 4))
  expect_identical(unique(x$union_type), c("unpartnered", "cohab", "marriage"))
  expect_identical(x$union_type[c(1, n + 1)], c("unpartnered", "cohab"))
  expect_identical(x$group[c(1, 3 * n + 1)], c("B", "A"))
  # Strata sum to the pooled table.
  pooled <- calculate_mortality_rates(pop, c(0, 20, 40), c(0, 30, Inf))
  expect_equal(as.vector(rowsum(x$person_years, rep(seq_len(n), 6))), pooled$person_years)
})

test_that("observation inputs are validated before splitting", {
  expect_error(
    calculate_mortality_rates(list(1), c(0, 1), c(0, Inf)),
    "must be a hammel_population"
  )
  pop <- by_history()
  pop$persons$group <- NULL
  pop$groups <- NULL
  pop$transitions <- NULL
  expect_error(calculate_dissolution_rates(pop, c(0, 40), c(0, Inf), by = "group_w"),
               "declared group dimension")
})
