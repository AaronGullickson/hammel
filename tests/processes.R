library(hammel)

rates <- data.frame(sex = c("m", "m", "f", "f"),
                    age_start = c(0, 1, 0, 1),
                    age_end = c(1, Inf, 1, Inf), rate = c(0.02, 0.01, 0.01, 0))
births <- data.frame(age_start = c(15, 20), age_end = c(20, 50), rate = c(0, 1.2))
death_process <- mortality(rates)
birth_process <- fertility(births)
stopifnot(inherits(death_process, "hammel_process"),
          inherits(death_process, "hammel_mortality"),
          inherits(birth_process, "hammel_fertility"),
          identical(death_process$rates, tibble::as_tibble(rates)),
          identical(birth_process$rate_unit, "person_year"))
model <- hammel_model() |> set_process(death_process) |> set_process(birth_process)
stopifnot(identical(model, hammel_model() |> set_mortality(rates) |> set_fertility(births)))
replacement <- rates; replacement$rate <- 0.03
updated <- set_mortality(model, replacement)
stopifnot(identical(model$processes$mortality, death_process),
          identical(updated$processes$fertility, birth_process),
          identical(updated$processes$mortality$rates, tibble::as_tibble(replacement)))

# Validation does not reorder or strip custom attributes.
custom <- rates[c(4, 2, 1, 3), ]; custom$source <- "example"
class(custom) <- c("custom_rates", "data.frame")
stopifnot(identical(mortality(custom)$rates, tibble::as_tibble(custom)))

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
expect_error(mortality(rates[rates$sex == "f", ]))
expect_error(mortality(rates[FALSE, ]))
expect_error(fertility(births[setdiff(names(births), "rate")]))
for (value in list(-1, Inf, NA_real_, "0.1")) {
  bad <- births; bad$rate <- value; expect_error(fertility(bad))
}
for (value in list(-1, Inf, NA_real_)) {
  bad <- births; bad$age_start[1] <- value; expect_error(fertility(bad))
}
bad <- births; bad$age_end[1] <- 15; expect_error(fertility(bad))
bad <- births; bad$age_end[1] <- 21; expect_error(fertility(bad))
expect_error(fertility(rbind(births, births[1, ])))
expect_error(fertility(data.frame(age_start = c(0, 1, 3),
                                age_end = c(10, 2, 4), rate = 0.1)))
bad <- rates; bad$sex <- c(0, 0, 1, 1); expect_error(mortality(bad))
bad <- rates; bad$time_start <- 2000; expect_error(mortality(bad))
bad <- births; bad$time_end <- 2100; expect_error(fertility(bad))
edited <- death_process; edited$rates$rate[1] <- -1
expect_error(set_process(model, edited))
edited <- death_process; edited$rate_unit <- "month"
expect_error(set_process(model, edited))
expect_error(set_process(model, rates))

# Constructors do not invent gap behavior or require full lifetime coverage.
gap <- births; gap$age_start[2] <- 25
stopifnot(identical(fertility(gap)$rates, tibble::as_tibble(gap)))

stopifnot(tibble::is_tibble(death_process$rates),
          tibble::is_tibble(birth_process$rates))
edited <- death_process; edited$rates <- as.data.frame(edited$rates)
stopifnot(tibble::is_tibble(set_process(model, edited)$processes$mortality$rates))
