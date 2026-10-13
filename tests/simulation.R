library(hammel)

make_model <- function(n = 1000, sex = "m", age = 20, death_rate = 0,
                       birth_rate = 0) {
  persons <- data.frame(person_id = seq_len(n), sex = rep(sex, n),
    birth_time = rep(-age, n), death_time = rep(NA_real_, n),
    mother_id = rep(NA_integer_, n), father_id = rep(NA_integer_, n))
  hammel_model() |>
    set_population(hammel_population(persons, end_time = 0)) |>
    set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
                             age_end = Inf, rate = death_rate)) |>
    set_fertility(data.frame(age_start = 15, age_end = 50, rate = birth_rate))
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
expect_error(simulate(hammel_model(), 1))
stopifnot(simulate(make_model(), 0.1)$end_time == 0.1)
expect_error(simulate(make_model(), 1, birth_interval = 0.1))
model <- make_model()
zero <- simulate(model, 1, seed = 42)
stopifnot(zero$simulation$births == 0, zero$simulation$deaths == 0,
          all(is.na(zero$persons$death_time)), zero$end_time == 1,
          identical(model$population$end_time, 0),
          identical(simulate(model, 0)$persons, model$population$persons))
empty <- simulate(make_model(0), 1, seed = 42)
stopifnot(nrow(empty$persons) == 0, empty$end_time == 1)

# Analytical continuous-time mortality survival.
n <- 50000
model <- make_model(n, death_rate = 0.2)
set.seed(99); before <- .Random.seed
result <- simulate(model, 2, seed = 42)
stopifnot(identical(.Random.seed, before),
          identical(result, simulate(model, 2, seed = 42)))
expected <- exp(-0.4)
stopifnot(abs(mean(is.na(result$persons$death_time)) - expected) <
            6 * sqrt(expected * (1 - expected) / n),
          all(result$persons$death_time < 2, na.rm = TRUE))
restart <- hammel_model() |> set_population(result) |>
  set_mortality(model$processes$mortality$rates) |>
  set_fertility(model$processes$fertility$rates)
later <- simulate(restart, 1, seed = 43)
dead <- !is.na(result$persons$death_time)
stopifnot(identical(later$persons$death_time[dead], result$persons$death_time[dead]),
          later$end_time == 3)

# Exact mortality times match the sampled exponential waits, with no rounding.
set.seed(42)
waits <- rexp(100, rate = 1e8)
exact <- simulate(make_model(100, death_rate = 1e8), 1/12, seed = 42)
stopifnot(isTRUE(all.equal(exact$persons$death_time, waits, tolerance = 1e-14)),
          all(exact$persons$death_time > 0))
# A candidate exactly at the half-open endpoint is excluded.
set.seed(42)
endpoint <- rexp(1)
stopifnot(is.na(simulate(make_model(1, death_rate = 1), endpoint,
                        seed = 42)$persons$death_time))
# Birth precedes death within the same former monthly bin. Choose the rate
# from the known draws so that the first birth occurs at 0.02 and death at 0.06.
set.seed(42)
draws <- rexp(2)
ordered <- simulate(make_model(1, sex = "f", death_rate = draws[1]/0.06,
                               birth_rate = draws[2]/0.02),
                    0.03, seed = 42, birth_interval = 0)
stopifnot(ordered$simulation$births == 1,
          abs(ordered$persons$birth_time[2] - 0.02) < 1e-14)

# Female-only populations still give birth, with unknown fathers.
model <- make_model(1000, sex = "f", birth_rate = 0.5)
model$population$persons$custom <- seq_len(1000)
result <- simulate(model, 1, seed = 42)
new <- result$persons$person_id > 1000
stopifnot(any(new), all(result$persons$mother_id[new] %in% seq_len(1000)),
          all(is.na(result$persons$custom[new])),
          identical(result$persons$custom[!new], seq_len(1000)),
          nrow(result$unions) == 0L)
# A newborn male can become a father immediately after birth in this basic
# model. The first birth, before any such male exists, has an unknown father.
stopifnot(is.na(result$persons$father_id[which(new)[1]]))

# Birth accounting, living father references, spacing, and ID continuation.
model <- make_model(2000, birth_rate = 0.2)
model$population$persons$sex[seq_len(1000)] <- "f"
result <- simulate(model, 5, seed = 42)
new <- result$persons$person_id > 2000
stopifnot(sum(new) == result$simulation$births,
          nrow(result$persons) == 2000 + result$simulation$births,
          !anyDuplicated(result$persons$person_id),
          all(result$persons$father_id[new] %in%
                result$persons$person_id[result$persons$sex == "m"]))
father_rows <- match(result$persons$father_id[new], result$persons$person_id)
child_times <- result$persons$birth_time[new]
stopifnot(all(result$persons$birth_time[father_rows] <= child_times),
          all(is.na(result$persons$death_time[father_rows]) |
                result$persons$death_time[father_rows] >= child_times))
for (times in split(result$persons$birth_time[new], result$persons$mother_id[new])) {
  stopifnot(all(diff(sort(times)) >= 0.75 - 1e-9))
}
validate_population(result)
restarted <- hammel_model() |> set_population(result) |>
  set_mortality(model$processes$mortality$rates) |>
  set_fertility(model$processes$fertility$rates)
continued <- simulate(restarted, 1, seed = 43)
for (times in split(continued$persons$birth_time[!is.na(continued$persons$mother_id)],
                    continued$persons$mother_id[!is.na(continued$persons$mother_id)])) {
  stopifnot(all(diff(sort(times)) >= 0.75 - 1e-9))
}

# Continuous waits with retained adjustment have a known renewal mean.
model <- make_model(10000, sex = "f", birth_rate = 0.2)
result <- simulate(model, 20, seed = 123)
q <- -expm1(-0.2/12)
eligible_hazard <- 0.2 / (1 - 9*q)
expected_rate <- 1 / (0.75 + 1 / eligible_hazard)
observed <- sum(result$persons$mother_id <= 10000 &
                  result$persons$birth_time >= 5, na.rm = TRUE) / (10000 * 15)
stopifnot(abs(observed - expected_rate) < 0.005,
          abs(expected_rate - 0.2) < 0.005)

# Death closes an existing union; no additional union rows appear.
model <- make_model(2, death_rate = 1e8)
model$population$persons$sex[2] <- "f"
model$population$unions <- data.frame(union_id = 1L, partner1_id = 1L,
  partner2_id = 2L, start_time = -1, end_time = NA_real_, end_reason = NA_character_)
result <- simulate(model, 1/12, seed = 42)
stopifnot(nrow(result$unions) == 1L, result$unions$end_time == min(result$persons$death_time),
          result$unions$end_reason == "death")

# No execution step; fractional durations are supported independently of compensation.
result <- simulate(make_model(), 0.123, seed = 42)
stopifnot(result$end_time == 0.123,
          identical(names(result$simulation),
                    c("start_time", "end_time", "birth_interval", "births", "deaths")))
expect_error(simulate(make_model(), -1))
expect_error(simulate(make_model(), 1, step_duration = 1))
stopifnot(tibble::is_tibble(result$persons), tibble::is_tibble(result$unions),
          tibble::is_tibble(simulate(make_model(0), 1)$persons))
