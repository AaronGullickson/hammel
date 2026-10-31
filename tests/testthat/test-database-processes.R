test_that("HMD conversion selects the year and preserves rates and open ages", {
  hmd <- data.frame(Year = c(1999, 2000, 2000), Age = c("0", "0", "110+"), Female = c(0.5,
    0.01, 0.6), Male = c(0.5, 0.02, 0.7))
  process <- hammel:::mortality_from_hmd(hmd, 2000)
  expect_identical(process$rates$sex, c("m", "m", "f", "f"))
  expect_identical(process$rates$age_start, c(0, 110, 0, 110))
  expect_identical(process$rates$age_end, c(1, Inf, 1, Inf))
  expect_identical(process$rates$rate, c(0.02, 0.7, 0.01, 0.6))
  expect_error(hammel:::mortality_from_hmd(hmd, 2001))
  bad <- hmd
  bad$Male[2] <- NA_real_
  expect_error(hammel:::mortality_from_hmd(bad, 2000))
  bad <- hmd
  bad$Age[2] <- "0-4"
  expect_error(hammel:::mortality_from_hmd(bad, 2000))
  expect_error(hammel:::mortality_from_hmd(rbind(hmd, hmd[2, ]), 2000))
  expect_error(get_mortality_hmd("USA", 2000, username = "", password = ""))
  expect_error(get_mortality_hmd("USA", c(2000, 2000), "test", "test"))
})

test_that("HFD conversion represents both endpoints and rejects invalid rates", {
  hfd <- data.frame(Year = c(1999, rep(2000, 4)), Age = c("20", "12-", "13", "14", "55+"),
    ASFR = c(0.9, 0.003, 0.01, 0.02, 0.001))
  births <- hammel:::fertility_from_hfd(hfd, 2000)
  expect_identical(births$rates$age_start, c(12, 13, 14, 55))
  expect_identical(births$rates$age_end, c(13, 14, 15, Inf))
  expect_identical(births$rates$rate, c(0.003, 0.01, 0.02, 0.001))
  expect_error(hammel:::fertility_from_hfd(hfd, 2001))
  bad <- hfd
  bad$ASFR[3] <- NA_real_
  expect_error(hammel:::fertility_from_hfd(bad, 2000))
  expect_equal(nrow(hammel:::fertility_from_hfd(hfd[c(2, 5), ], 2000)$rates), 2)
  bad <- hfd
  bad$Age[3] <- "13-17"
  expect_error(hammel:::fertility_from_hfd(bad, 2000))
})

test_that("local HMD readers preserve the upper open interval", {
  skip_if_not_installed("HMDHFDplus")
  hmd <- data.frame(Year = c(1999, 2000, 2000), Age = c("0", "0", "110+"), Female = c(0.5,
    0.01, 0.6), Male = c(0.5, 0.02, 0.7))
  process <- hammel:::mortality_from_hmd(hmd, 2000)
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("Synthetic HMD fixture", "", "Year Age Female Male Total", "2000 0 0.01 0.02 0.015",
    "2000 110+ 0.6 0.7 0.65"), path)
  read <- HMDHFDplus::readHMD(path, fixup = FALSE)
  expect_identical(hammel:::mortality_from_hmd(read, 2000), process)
  parsed <- HMDHFDplus::readHMD(path, fixup = TRUE)
  expect_identical(hammel:::mortality_from_hmd(parsed, 2000), process)
  unlink(path)
})

test_that("local HFD readers preserve retained age bands", {
  skip_if_not_installed("HMDHFDplus")
  hfd <- data.frame(Year = c(1999, rep(2000, 4)), Age = c("20", "12-", "13", "14", "55+"),
    ASFR = c(0.9, 0.003, 0.01, 0.02, 0.001))
  births <- hammel:::fertility_from_hfd(hfd, 2000)
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("Synthetic HFD fixture", "", "Year Age ASFR", "2000 12- 0.003", "2000 13 0.01",
    "2000 14 0.02", "2000 55+ 0.001"), path)
  read <- HMDHFDplus::readHFD(path, fixup = FALSE)
  expect_identical(hammel:::fertility_from_hfd(read, 2000), births)
  unlink(path)
})

test_that("HMD downloads preserve provenance and rates can be reused without downloads",
  {
    skip_if_not_installed("HMDHFDplus")
    hmd <- data.frame(Year = c(1999, 2000, 2000), Age = c("0", "0", "110+"), Female = c(0.5,
      0.01, 0.6), Male = c(0.5, 0.02, 0.7))
    process <- hammel:::mortality_from_hmd(hmd, 2000)
    with_database_reader("readHMDweb", function(CNTRY, item, username, password, fixup) {
      expect_equal(CNTRY, "USA", tolerance = 0)
      expect_equal(item, "Mx_1x1", tolerance = 0)
      expect_true(!fixup)
      hmd
    }, {
      rates <- get_mortality_hmd("USA", 2000, "test", "test")
      path <- tempfile()
      saveRDS(rates, path)
      restored <- readRDS(path)
      unlink(path)
      expect_true(inherits(rates, "tbl_df"))
      expect_identical(rates, restored)
      model <- set_mortality(hammel_model(), restored)
      with_database_reader("readHMDweb", function(...) stop("Unexpected download"),
        {
         reused <- set_mortality(hammel_model(), restored)
         expect_identical(reused, model)
        })
      imported <- model$processes$mortality
      expect_identical(as.data.frame(imported$rates)[names(process$rates)], as.data.frame(process$rates))
      expect_equal(imported$source$year, 2000, tolerance = 0)
      expect_true(!any(c("username", "password") %in% names(imported$source)))
    })
  })

test_that("HFD downloads preserve revision metadata and rates can be reused without downloads",
  {
    skip_if_not_installed("HMDHFDplus")
    hfd <- data.frame(Year = c(1999, rep(2000, 4)), Age = c("20", "12-", "13", "14",
      "55+"), ASFR = c(0.9, 0.003, 0.01, 0.02, 0.001))
    births <- hammel:::fertility_from_hfd(hfd, 2000)
    with_database_reader("readHFDweb", function(CNTRY, item, username, password, fixup,
      Update) {
      expect_equal(CNTRY, "USA", tolerance = 0)
      expect_equal(item, "asfrRR", tolerance = 0)
      expect_true(!fixup)
      expect_equal(Update, "20230328", tolerance = 0)
      hfd
    }, {
      rates <- get_fertility_hfd("USA", 2000, "test", "test", update = "20230328")
      path <- tempfile()
      saveRDS(rates, path)
      restored <- readRDS(path)
      unlink(path)
      expect_true(inherits(rates, "tbl_df"))
      expect_identical(rates, restored)
      model <- set_fertility(hammel_model(), rates)
      with_database_reader("readHFDweb", function(...) stop("Unexpected download"),
        {
         reused <- set_fertility(hammel_model(), restored)
         expect_identical(reused, model)
        })
      imported <- model$processes$fertility
      expect_identical(imported$rates$rate, births$rates$rate)
      expect_equal(imported$source$open_ages, "lower_as_single_year_upper_open", tolerance = 0)
      expect_equal(imported$source$update, "20230328", tolerance = 0)
    })
  })

