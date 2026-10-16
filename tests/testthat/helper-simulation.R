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
         birth_interval = birth_interval) |>
  set_union_formation("on_demand",
    create_matching_model(~ age_difference, c(age_difference = 0)))
}

make_rate_model <- function(n = 2, search_rate = 1, birth_rate = 0, birth_interval = 0) {
  model <- make_model(n, birth_rate = birth_rate, birth_interval = birth_interval)
  if (n) model$population$persons$sex[1] <- "f"
  rates <- model$processes$fertility$rates
  model <- set_fertility(model, rbind(
    transform(rates, union_status = "partnered"),
    transform(rates, union_status = "unpartnered", rate = 0)), birth_interval)
  set_union_formation(model, "rate",
    create_matching_model(~ age_difference, c(age_difference = 0)),
    rates = data.frame(age_start = 15, age_end = Inf, rate = search_rate),
    sample_size = 0)
}
