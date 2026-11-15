# Expect a Hammel error of `kind` naming the test's stand-in function `fn`.
expect_engine_error <- function(object, kind, fn) {
  error <- expect_error(object, class = paste0("hammel_error_", kind))
  expect_identical(rlang::call_name(error$call), fn)
  invisible(error)
}

native_error <- function(message) {
  structure(
    class = c("Rcpp::exception", "C++Error", "error", "condition"),
    list(message = message, call = NULL)
  )
}

test_that("native messages are split into known tags and text", {
  expect_identical(
    parse_native_message("[limits] Too large."),
    list(tag = "limits", text = "Too large.")
  )
  expect_null(parse_native_message("[unknown] Text.")$tag)
  expect_null(parse_native_message("Invalid exact-draw inputs.")$tag)
})

test_that("call_native translates native errors and warnings", {
  engine <- function(message) {
    call_native(stop(native_error(message)))
  }
  error <- expect_engine_error(
    engine("[limits] Population exceeds supported ID or row limits."),
    "runtime",
    "engine"
  )
  expect_identical(error$tag, "limits")
  expect_match(conditionMessage(error), "Use a smaller population")
  error <- expect_engine_error(
    engine("[population] Kin {cycle}."),
    "population",
    "engine"
  )
  expect_match(conditionMessage(error), "Kin {cycle}.", fixed = TRUE)
  error <- expect_engine_error(
    engine("Internal error: cell totals differ."),
    "internal",
    "engine"
  )
  expect_match(conditionMessage(error), "Internal error: cell totals differ.",
    fixed = TRUE)
  expect_match(conditionMessage(error), "likely a bug", fixed = TRUE)
  expect_warning(
    call_native(warning("[weight_missing] Missing values were zero.")),
    "Missing values were zero.",
    class = "hammel_warning_no_effect"
  )
  expect_warning(call_native(warning("Other.")), "Other.")
})

test_that("callback errors keep their class or become runtime errors", {
  engine <- function(expr) {
    call_native(expr)
  }
  expect_error(
    engine(abort_hammel("Bad category.", "runtime")),
    class = "hammel_error_runtime"
  )
  error <- expect_engine_error(
    engine(stop("my function broke")),
    "runtime",
    "engine"
  )
  expect_match(conditionMessage(error$parent), "my function broke")
})

test_that("engine errors reach users with classes and function names", {
  persons <- tibble::tibble(
    person_id = c(1L, 2L),
    sex = "f",
    birth_time = c(-30, -40),
    death_time = NA_real_,
    mother_id = c(2L, 1L),
    father_id = NA_integer_
  )
  population <- hammel_population(persons, end_time = 0)
  error <- expect_hammel_error(
    extract_kin(population, 1),
    "population",
    "extract_kin"
  )
  expect_identical(error$tag, "population")

  population <- generate_uniform_population(20)
  women <- population$persons$person_id[population$persons$sex == "f"]
  men <- population$persons$person_id[population$persons$sex == "m"]
  unions <- tibble::tibble(
    union_id = 1:2,
    wife_id = women[1],
    husband_id = men[1:2],
    start_time = -1,
    end_time = NA_real_,
    end_reason = NA_character_
  )
  population <- hammel_population(population$persons, unions, end_time = 0)
  model <- hammel_model() |>
    set_population(population) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(usa_2000_fertility) |>
    set_union_formation("on_demand")
  expect_hammel_error(run_simulation(model, 1), "population", "run_simulation")

  weights <- weight_model(
    fun = function(members, bin) {
      rep(0, nrow(members) + 1)
    }
  )
  model <- hammel_model() |>
    set_population(generate_uniform_population(50)) |>
    set_mortality(usa_2000_mortality) |>
    set_fertility(
      data.frame(age_start = 0, age_end = Inf, rate = 2),
      weights = weights
    ) |>
    set_union_formation("on_demand")
  error <- expect_hammel_error(
    run_simulation(model, 1, seed = 1),
    "runtime",
    "run_simulation"
  )
  expect_identical(error$tag, "weight_function")
})

test_that("native messages keep their layout", {
  engine <- function(message) {
    call_native(stop(native_error(message)))
  }
  expect_snapshot(
    engine("[weight_function] A custom fertility weight function failed."),
    error = TRUE
  )
  expect_snapshot(engine("Invalid exact-draw inputs."), error = TRUE)
})
