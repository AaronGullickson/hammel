library(hammel)

persons <- tibble::tibble(person_id = c(1L, 2L), sex = c("m", "f"),
                      birth_time = c(-25, -23), death_time = NA_real_,
                      mother_id = NA_integer_, father_id = NA_integer_)
unions <- hammel_population(persons)$unions
mortality <- tibble::tibble(sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01)
fertility <- tibble::tibble(age_start = 15, age_end = 50, rate = 0.1)

empty <- hammel_model()
configured <- empty |>
  set_population(hammel_population(persons, unions, end_time = 0)) |>
  set_mortality(mortality) |>
  set_fertility(fertility)
stopifnot(
  inherits(configured, "hammel_model"),
  is.null(empty$population),
  identical(empty$processes, list()),
  identical(configured$population$persons, persons),
  identical(configured$population$unions, unions),
  identical(configured$population$end_time, 0),
  identical(configured$processes$mortality$rates, mortality),
  identical(configured$processes$fertility$rates, fertility),
  identical(set_mortality(empty, mortality),
            set_process(empty, hammel::mortality(mortality))),
  identical(set_fertility(empty, fertility),
            set_process(empty, hammel::fertility(fertility)))
)

replacement <- mortality; replacement$rate <- 0.02
updated <- set_mortality(configured, replacement)
stopifnot(
  identical(configured$processes$mortality$rates, mortality),
  identical(updated$processes$mortality$rates, replacement),
  identical(updated$processes$fertility$rates, fertility),
  length(updated$processes) == 2L
)

# Restart preserves historical tables and uses the endpoint, not event dates.
results <- hammel_population(persons, unions, end_time = 120)
restarted <- set_population(configured, results)
stopifnot(
  inherits(results, "hammel_population"),
  identical(restarted$population$persons, persons),
  identical(restarted$population$unions, unions),
  identical(restarted$population$end_time, 120),
  identical(restarted$processes, configured$processes),
  identical(configured$population$end_time, 0)
)
unknown <- set_population(restarted, hammel_population(persons, unions))
stopifnot("end_time" %in% names(unknown$population), is.null(unknown$population$end_time))

# Input subclasses are normalized to tibbles while retaining column contents.
custom_persons <- structure(persons, class = c("custom_persons", "data.frame"))
stopifnot(identical(
  set_population(empty, hammel_population(custom_persons, unions))$population$persons, tibble::as_tibble(custom_persons)
))

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
expect_error(hammel_population(NULL, unions))
expect_error(hammel_population(persons, list()))
expect_error(set_population(empty, persons, list()))
expect_error(set_population(empty, results, unions))
expect_error(set_population(empty, results, start_time = 0))
expect_error(set_process(list(), hammel::mortality(mortality)))
expect_error(set_process(empty, mortality))
expect_error(set_process(empty, NULL))
for (month in list(NA_real_, Inf, c(0, 1), "120")) {
  expect_error(hammel_population(persons, unions, end_time = month))
  expect_error(set_population(empty, persons, unions, start_time = month))
}

# An unmarried population can omit unions or supply NULL, clearing old unions.
unmarried <- set_population(configured, hammel_population(persons))
stopifnot(
  identical(unmarried, set_population(configured, hammel_population(persons, unions = NULL))),
  identical(unmarried$population$persons, persons),
  is.data.frame(unmarried$population$unions),
  nrow(unmarried$population$unions) == 0L,
  identical(unmarried$processes, configured$processes),
  identical(configured$population$unions, unions)
)

# Fractional years are accepted and carried through restarts.
stopifnot(identical(set_population(empty,
  hammel_population(persons, end_time = 1.5))$population$end_time, 1.5))

# Both entry points reject malformed schemas and inconsistent histories.
check_rejected <- function(p = persons, u = unions) {
  expect_error(hammel_population(p, u))
  invalid <- structure(list(persons = p, unions = u, end_time = NULL),
                       class = "hammel_population")
  expect_error(set_population(empty, invalid))
}
check_rejected(persons[setdiff(names(persons), "sex")])
p <- persons; p$person_id[2] <- p$person_id[1]; check_rejected(p)
p <- persons; p$person_id[1] <- 1.5; check_rejected(p)
p <- persons; p$sex[1] <- 2; check_rejected(p)
p <- persons; p$birth_time[1] <- Inf; check_rejected(p)
p <- persons; p$death_time[1] <- -30; check_rejected(p)
p <- persons; p$mother_id[1] <- 999L; check_rejected(p)
p <- persons; p$father_id[1] <- 1L; check_rejected(p)
u <- data.frame(union_id = 1L, partner1_id = 1L, partner2_id = 2L,
                start_time = -1, end_time = NA_real_, end_reason = NA_character_)
stopifnot(identical(hammel_population(persons, u, end_time = 0)$unions, tibble::as_tibble(u)))
v <- u; v$partner2_id <- 999L; check_rejected(u = v)
v <- u; v$partner2_id <- 1L; check_rejected(u = v)
v <- u; v$end_time <- -2; v$end_reason <- "dissolution"; check_rejected(u = v)
v <- u; v$end_reason <- "death"; check_rejected(u = v)
v <- u; v$start_time <- -30; check_rejected(u = v)
check_rejected(u = rbind(u, u))
expect_error(hammel_population(persons, u, end_time = -2))
expect_error(set_population(empty, persons, u, start_time = -2))

# Extra columns and empty populations retain their intended structure.
p <- persons; p$custom_attribute <- c("a", "b")
stopifnot(identical(hammel_population(p)$persons, p),
          nrow(hammel_population(persons[FALSE, ])$persons) == 0L)

# Object validation is read-only; malformed results are not normalized.
checked <- withVisible(validate_population(results))
stopifnot(identical(checked$value, results), !checked$visible)
invalid <- results; invalid$unions <- NULL
expect_error(validate_population(invalid))
expect_error(set_population(empty, invalid))
invalid <- results; invalid$end_time <- Inf
expect_error(validate_population(invalid))
invalid <- results; invalid$persons$person_id[2] <- 1L
expect_error(validate_population(invalid))
expect_error(validate_population(persons))
stopifnot(identical(hammel_population(persons, data.frame())$unions, unions))

# Store the entire result, including optional observation metadata, unchanged.
annotated <- results
annotated$observations <- list(description = "preserve supplied components")
stored <- set_population(empty, annotated)
stopifnot(identical(stored$population, annotated),
          identical(names(stored), c("population", "processes")),
          identical(set_population(empty, hammel_population(persons, unions, end_time = 120))$population,
                    results))

expect_error(set_population(empty, persons))

# Sex uses explicit character labels; numeric codes and missing labels fail.
for (sex in list(c(0L, 1L), c("m", NA_character_), c("male", "female"))) {
  p <- persons; p$sex <- sex; check_rejected(p)
}

stopifnot(tibble::is_tibble(results$persons), tibble::is_tibble(results$unions),
          tibble::is_tibble(hammel_population(as.data.frame(persons))$persons),
          tibble::is_tibble(hammel_population(persons[FALSE, ])$persons))
