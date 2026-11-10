test_that("diagnostics validate and pass through the setter", {
  matching <- matching_model(~ age_difference, c(age_difference = 0))
  expect_false(union_formation("on_demand", matching = matching)$diagnostics)
  expect_true(set_union_formation(hammel_model(), "on_demand", matching = matching,
    diagnostics = TRUE)$processes$union_formation$diagnostics)
  for (bad in list(NULL, NA, 1, c(TRUE, FALSE), "yes")) {
    expect_error(union_formation("on_demand", matching = matching, diagnostics = bad), "diagnostics")
    model <- make_model(0)
    model$processes$union_formation$diagnostics <- bad
    expect_error(run_simulation(model, 0), "diagnostics")
  }
})

test_that("diagnostics add no candidate-set columns and do not change draws", {
  for (mode in c("on_demand", "rate")) {
    model <- if (mode == "rate") make_rate_model(5, search_rate = 100) else
      make_model(5, birth_rate = 100, birth_interval = 0)
    model$population$persons$sex[1] <- "f"
    model$population$persons$birth_time <- c(-25, -20, -25, -30, -5)
    plain <- run_simulation(model, 0.1, seed = 92)
    model$processes$union_formation$diagnostics <- TRUE
    observed <- run_simulation(model, 0.1, seed = 92)
    expect_identical(observed, plain)
    expect_false(any(c(".n_eligible", ".utility_rank") %in% names(observed$unions)))
  }
})
