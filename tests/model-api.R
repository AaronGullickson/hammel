library(hammel)

persons <- data.frame(person_id = c(1L, 2L))
unions <- data.frame(union_id = integer())
mortality <- data.frame(age = 0, rate = 0.01)
fertility <- data.frame(age = 20, rate = 0.1)

empty <- hammel_model()
configured <- empty |>
  set_population(persons, unions, start_month = 0) |>
  set_mortality(mortality) |>
  set_fertility(fertility)
stopifnot(
  inherits(configured, "hammel_model"),
  is.null(empty$persons),
  identical(empty$processes, list()),
  identical(configured$persons, persons),
  identical(configured$unions, unions),
  identical(configured$start_month, 0),
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
results <- hammel_results(persons, unions, end_month = 120)
restarted <- set_population(configured, results)
stopifnot(
  inherits(results, "hammel_results"),
  identical(restarted$persons, persons),
  identical(restarted$unions, unions),
  identical(restarted$start_month, 120),
  identical(restarted$processes, configured$processes),
  identical(configured$start_month, 0)
)
unknown <- set_population(restarted, hammel_results(persons, unions))
stopifnot("start_month" %in% names(unknown), is.null(unknown$start_month))

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
expect_error(set_population(empty, results, start_month = 0))
expect_error(set_process(list(), "mortality", mortality))
expect_error(set_process(empty, "", mortality))
expect_error(set_process(empty, NA_character_, mortality))
expect_error(set_process(empty, "mortality", NULL))
for (month in list(NA_real_, Inf, 1.5, c(0, 1), "120")) {
  expect_error(hammel_results(persons, unions, end_month = month))
  expect_error(set_population(empty, persons, unions, start_month = month))
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
