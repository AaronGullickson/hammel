# Run explicitly with HAMMEL_RUN_ONLINE_TESTS=true. Ordinary package checks
# never contact the databases or require credentials.
if (identical(Sys.getenv("HAMMEL_RUN_ONLINE_TESTS"), "true")) {
  library(hammel)
  if (!requireNamespace("HMDHFDplus", quietly = TRUE)) {
    stop("The online test requires HMDHFDplus.", call. = FALSE)
  }
  credentials <- c("HMD_USERNAME", "HMD_PASSWORD", "HFD_USERNAME", "HFD_PASSWORD")
  if (any(!nzchar(Sys.getenv(credentials)))) {
    stop("The online test requires HMD and HFD credential environment variables.",
         call. = FALSE)
  }

  # Credentials are consumed by the helpers; never print them or downloaded rows.
  model <- hammel_model() |>
    set_mortality_hmd("USA", year = 2000) |>
    set_fertility_hfd("USA", year = 2000)
  deaths <- model$processes$mortality
  births <- model$processes$fertility
  stopifnot(
    inherits(deaths, "hammel_mortality"),
    inherits(births, "hammel_fertility"),
    identical(deaths$rate_unit, "person_year"),
    identical(births$rate_unit, "person_year"),
    setequal(deaths$rates$sex, c("m", "f")),
    any(is.infinite(deaths$rates$age_end[deaths$rates$sex == "m"])),
    any(is.infinite(deaths$rates$age_end[deaths$rates$sex == "f"])),
    all(is.finite(deaths$rates$rate)), all(deaths$rates$rate >= 0),
    all(is.finite(births$rates$rate)), all(births$rates$rate >= 0),
    any(deaths$rates$rate > 0), any(births$rates$rate > 0),
    all(births$rates$age_end - births$rates$age_start == 1),
    deaths$source$country == "USA", births$source$country == "USA",
    deaths$source$year == 2000, births$source$year == 2000,
    deaths$source$item == "Mx_1x1", births$source$item == "asfrRR",
    births$source$open_ages == "omit",
    !any(c("username", "password") %in% names(deaths$source)),
    !any(c("username", "password") %in% names(births$source))
  )
  cat("Online HMD and HFD download tests passed.\n")
} else {
  cat("Skipping online tests; set HAMMEL_RUN_ONLINE_TESTS=true to enable.\n")
}
