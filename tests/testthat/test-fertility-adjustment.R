# General fertility execution checks. Midpoint compensation (D028/D099) was
# removed when cell-accounted fertility (D101) became the only method; its
# analytical tests are in Git history.

test_that("explicit zero fertility tails preserve histories and positive tails execute", {
  model <- make_model(50,sex="f",age=45,birth_rate=.1)
  rates <- model$processes$fertility$rates
  zero <- set_fertility(model,rbind(rates,data.frame(age_start=50,age_end=Inf,rate=0)))
  expect_identical(run_simulation(model,10,42),run_simulation(zero,10,42))
  positive <- make_model(100,sex="f",age=60,birth_rate=0)
  positive <- set_fertility(positive,data.frame(age_start=55,age_end=Inf,rate=.2))
  result <- run_simulation(positive,2,42)
  expect_gt(count_run_events(result, "births"),0)
})

test_that("independent mortality preserves constant fertility after burn-in", {
  n <- 10000
  model <- make_model(n, sex = "f", death_rate = 0.1, birth_rate = 0.2,
                      birth_interval = 0.37)
  result <- run_simulation(model, 20, seed = 456)
  deaths <- result$persons$death_time[seq_len(n)]
  deaths[is.na(deaths)] <- 20
  exposure <- sum(pmax(0, pmin(deaths, 20) - 5))
  births <- sum(result$persons$mother_id <= n &
                result$persons$birth_time >= 5, na.rm = TRUE)
  expect_lt(abs(births / exposure - 0.2), 0.01)
})
