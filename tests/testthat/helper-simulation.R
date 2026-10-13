make_model <- function(n = 1000, sex = "m", age = 20, death_rate = 0,
           birth_rate = 0, birth_interval = 0.75) {
 persons <- data.frame(person_id = seq_len(n), sex = rep(sex, n),
  birth_time = rep(-age, n), death_time = rep(NA_real_, n),
  mother_id = rep(NA_integer_, n), father_id = rep(NA_integer_, n))
 hammel_model() |>
  set_population(hammel_population(persons, end_time = 0)) |>
  set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
              age_end = Inf, rate = death_rate)) |>
  set_fertility(data.frame(age_start = 15, age_end = 50, rate = birth_rate),
         birth_interval = birth_interval)
}
