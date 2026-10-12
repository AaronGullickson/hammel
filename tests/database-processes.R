library(hammel)

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
hmd <- data.frame(Year = c(1999, 2000, 2000), Age = c("0", "0", "110+"),
                  Female = c(0.5, 0.01, 0.6), Male = c(0.5, 0.02, 0.7))
process <- hammel:::mortality_from_hmd(hmd, 2000)
stopifnot(identical(process$rates$sex, c("m", "m", "f", "f")),
          identical(process$rates$age_start, c(0, 110, 0, 110)),
          identical(process$rates$age_end, c(1, Inf, 1, Inf)),
          identical(process$rates$rate, c(0.02, 0.7, 0.01, 0.6)))
expect_error(hammel:::mortality_from_hmd(hmd, 2001))
bad <- hmd; bad$Male[2] <- NA_real_
expect_error(hammel:::mortality_from_hmd(bad, 2000))
bad <- hmd; bad$Age[2] <- "0-4"
expect_error(hammel:::mortality_from_hmd(bad, 2000))
expect_error(hammel:::mortality_from_hmd(rbind(hmd, hmd[2, ]), 2000))
expect_error(mortality_hmd("USA", 2000, username = "", password = ""))
expect_error(mortality_hmd("USA", c(2000, 2001), "test", "test"))
expect_error(set_mortality_hmd(list(), "USA", 2000, "test", "test"))

# Verify local parsing with the optional package, without authentication.
if (requireNamespace("HMDHFDplus", quietly = TRUE)) {
  path <- tempfile(fileext = ".txt")
  writeLines(c("Synthetic HMD fixture", "", "Year Age Female Male Total",
               "2000 0 0.01 0.02 0.015", "2000 110+ 0.6 0.7 0.65"), path)
  read <- HMDHFDplus::readHMD(path, fixup = FALSE)
  stopifnot(identical(hammel:::mortality_from_hmd(read, 2000), process))
  unlink(path)
}

hfd <- data.frame(Year = c(1999, rep(2000, 4)),
                  Age = c("20", "12-", "13", "14", "55+"),
                  ASFR = c(0.9, 0.003, 0.01, 0.02, 0.001))
births <- hammel:::fertility_from_hfd(hfd, 2000)
stopifnot(identical(births$rates$age_start, c(13, 14)),
          identical(births$rates$age_end, c(14, 15)),
          identical(births$rates$rate, c(0.01, 0.02)))
expect_error(hammel:::fertility_from_hfd(hfd, 2001))
bad <- hfd; bad$ASFR[3] <- NA_real_
expect_error(hammel:::fertility_from_hfd(bad, 2000))
expect_error(hammel:::fertility_from_hfd(hfd[c(2, 5), ], 2000))
bad <- hfd; bad$Age[3] <- "13-17"
expect_error(hammel:::fertility_from_hfd(bad, 2000))

if (requireNamespace("HMDHFDplus", quietly = TRUE)) {
  path <- tempfile(fileext = ".txt")
  writeLines(c("Synthetic HFD fixture", "", "Year Age ASFR",
               "2000 12- 0.003", "2000 13 0.01", "2000 14 0.02",
               "2000 55+ 0.001"), path)
  read <- HMDHFDplus::readHFD(path, fixup = FALSE)
  stopifnot(identical(hammel:::fertility_from_hfd(read, 2000), births))
  unlink(path)

  # Verify public download plumbing with temporary reader replacements.
  # Restore each binding even on failure; no network calls or real credentials.
  with_reader <- function(name, replacement, code) {
    ns <- asNamespace("HMDHFDplus")
    original <- get(name, ns)
    unlockBinding(name, ns)
    assign(name, replacement, ns)
    lockBinding(name, ns)
    on.exit({
      unlockBinding(name, ns)
      assign(name, original, ns)
      lockBinding(name, ns)
    })
    force(code)
  }
  with_reader("readHMDweb", function(CNTRY, item, username, password, fixup) {
    stopifnot(CNTRY == "USA", item == "Mx_1x1", !fixup)
    hmd
  }, {
    model <- set_mortality_hmd(hammel_model(), "USA", 2000, "test", "test")
    imported <- model$processes$mortality
    stopifnot(identical(imported$rates, process$rates),
              imported$source$year == 2000,
              !any(c("username", "password") %in% names(imported$source)))
  })
  with_reader("readHFDweb", function(CNTRY, item, username, password, fixup, Update) {
    stopifnot(CNTRY == "USA", item == "asfrRR", !fixup, Update == "20230328")
    hfd
  }, {
    model <- set_fertility_hfd(hammel_model(), "USA", 2000, "test", "test",
                              update = "20230328")
    imported <- model$processes$fertility
    stopifnot(identical(imported$rates, births$rates),
              imported$source$open_ages == "omit",
              imported$source$update == "20230328")
  })
}
