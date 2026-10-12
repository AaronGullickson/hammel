library(hammel)

set.seed(42)
population <- create_baseline_population(50000, end_time = 2000)
set.seed(42)
stopifnot(identical(population, create_baseline_population(50000, end_time = 2000)))
validate_population(population)
p <- population$persons
ages <- population$end_time - p$birth_time
stopifnot(
  identical(p$person_id, seq_len(50000)),
  all(ages >= 0 & ages <= 70),
  all(is.na(p$mother_id)), all(is.na(p$father_id)), all(is.na(p$death_time)),
  nrow(population$unions) == 0L,
  identical(set_population(hammel_model(), population)$population, population)
)
# Broad six-standard-error checks against the specified distributions.
stopifnot(abs(mean(p$sex == "f") - 0.5) < 6 * sqrt(0.25 / nrow(p)),
          abs(mean(ages) - 35) < 6 * sqrt(70^2 / (12 * nrow(p))))

newborns <- create_baseline_population(3, max_age = 0, end_time = 2000,
                                female_probability = 1)
stopifnot(identical(newborns$persons$sex, rep("f", 3)),
          identical(newborns$persons$birth_time, rep(2000, 3)),
          all(create_baseline_population(3, female_probability = 0)$persons$sex == "m"))
empty <- create_baseline_population(0)
validate_population(empty)
stopifnot(nrow(empty$persons) == 0L, nrow(empty$unions) == 0L)
fractional <- create_baseline_population(10, max_age = 0.5, end_time = -1.25)
stopifnot(all(fractional$persons$birth_time <= -1.25),
          all(fractional$persons$birth_time >= -1.75))

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
for (size in list(-1, 1.5, NA_real_, Inf, "10", c(1, 2), 2^31)) {
  expect_error(create_baseline_population(size))
}
for (age in list(-1, NA_real_, Inf, "70", c(1, 2))) {
  expect_error(create_baseline_population(2, max_age = age))
}
for (time in list(NULL, NA_real_, Inf, "2000", c(1, 2))) {
  expect_error(create_baseline_population(2, end_time = time))
}
for (probability in list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0, 1))) {
  expect_error(create_baseline_population(2, female_probability = probability))
}
