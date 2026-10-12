library(hammel)

persons <- data.frame(person_id = c(1L, 2L), sex = c(0L, 1L),
                      birth_time = c(-25, -23), death_time = NA_real_,
                      mother_id = NA_integer_, father_id = NA_integer_)
unions <- hammel_results(persons)$unions
mortality <- data.frame(age = 0, rate = 0.01)
fertility <- data.frame(age = 20, rate = 0.1)

empty <- hammel_model()
configured <- empty |>
  set_population(persons, unions, start_time = 0) |>
  set_mortality(mortality) |>
  set_fertility(fertility)
stopifnot(
  inherits(configured, "hammel_model"),
  is.null(empty$persons),
  identical(empty$processes, list()),
  identical(configured$persons, persons),
  identical(configured$unions, unions),
  identical(configured$start_time, 0),
  identical(configured$processes$mortality, mortality),
  identical(configured$processes$fertility, fertility),
  identical(set_mortality(empty, mortality),
            set_process(empty, "mortality", mortality)),
  identical(set_fertility(empty, fertility),
            set_process(empty, "fertility", fertility))
)

replacement <- data.frame(age = 0, rate = 0.02)
updated <- set_mortality(configured, replacement)
stopifnot(
  identical(configured$processes$mortality, mortality),
  identical(updated$processes$mortality, replacement),
  identical(updated$processes$fertility, fertility),
  length(updated$processes) == 2L
)

# Restart preserves historical tables and uses the endpoint, not event dates.
results <- hammel_results(persons, unions, end_time = 120)
restarted <- set_population(configured, results)
stopifnot(
  inherits(results, "hammel_results"),
  identical(restarted$persons, persons),
  identical(restarted$unions, unions),
  identical(restarted$start_time, 120),
  identical(restarted$processes, configured$processes),
  identical(configured$start_time, 0)
)
unknown <- set_population(restarted, hammel_results(persons, unions))
stopifnot("start_time" %in% names(unknown), is.null(unknown$start_time))

# Data-frame subclasses and custom columns survive construction unchanged.
custom_persons <- structure(persons, class = c("custom_persons", "data.frame"))
stopifnot(identical(
  set_population(empty, custom_persons, unions)$persons, custom_persons
))

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
expect_error(hammel_results(NULL, unions))
expect_error(hammel_results(persons, list()))
expect_error(set_population(empty, persons, list()))
expect_error(set_population(empty, results, unions))
expect_error(set_population(empty, results, start_time = 0))
expect_error(set_process(list(), "mortality", mortality))
expect_error(set_process(empty, "", mortality))
expect_error(set_process(empty, NA_character_, mortality))
expect_error(set_process(empty, "mortality", NULL))
for (month in list(NA_real_, Inf, c(0, 1), "120")) {
  expect_error(hammel_results(persons, unions, end_time = month))
  expect_error(set_population(empty, persons, unions, start_time = month))
}

# An unmarried population can omit unions or supply NULL, clearing old unions.
unmarried <- set_population(configured, persons)
stopifnot(
  identical(unmarried, set_population(configured, persons, unions = NULL)),
  identical(unmarried$persons, persons),
  is.data.frame(unmarried$unions),
  nrow(unmarried$unions) == 0L,
  identical(unmarried$processes, configured$processes),
  identical(configured$unions, unions)
)

# Fractional years are accepted and carried through restarts.
stopifnot(identical(set_population(empty,
  hammel_results(persons, end_time = 1.5))$start_time, 1.5))

# Both entry points reject malformed schemas and inconsistent histories.
check_rejected <- function(p = persons, u = unions) {
  expect_error(hammel_results(p, u))
  expect_error(set_population(empty, p, u))
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
stopifnot(identical(hammel_results(persons, u, end_time = 0)$unions, u))
v <- u; v$partner2_id <- 999L; check_rejected(u = v)
v <- u; v$partner2_id <- 1L; check_rejected(u = v)
v <- u; v$end_time <- -2; v$end_reason <- "dissolution"; check_rejected(u = v)
v <- u; v$end_reason <- "death"; check_rejected(u = v)
v <- u; v$start_time <- -30; check_rejected(u = v)
check_rejected(u = rbind(u, u))
expect_error(hammel_results(persons, u, end_time = -2))
expect_error(set_population(empty, persons, u, start_time = -2))

# Extra columns and empty populations retain their intended structure.
p <- persons; p$custom_attribute <- c("a", "b")
stopifnot(identical(hammel_results(p)$persons, p),
          nrow(hammel_results(persons[FALSE, ])$persons) == 0L)
