hmd_year_fixture <- function() {
  data.frame(Year = rep(2000:2002, each = 3), Age = rep(c("0", "1", "110+"), 3),
    Male = seq_len(9) / 100, Female = seq_len(9) / 200)
}
hmd_group_fixture <- function(periods = c("2000-2004", "2005-2009")) {
  data.frame(Year = rep(periods, each = 4),
    Age = rep(c("0", "1-4", "5-9", "110+"), length(periods)),
    Male = .02, Female = .01)
}

test_that("year requests reject invalid vectors before any download", {
  for (years in list(numeric(), c(2000, 2000), c(2000, NA), c(2000, Inf),
      c(2000, 2001.5), "2000", matrix(2000))) {
    expect_error(get_mortality_hmd("SWE", years, "test", "test"), "year")
    expect_error(get_fertility_hfd("SWE", years, "test", "test"), "year")
  }
  for (shape in list(NULL, "5x2", "5", c("1x1", "5x1"), NA, 5, matrix("1x1")))
    expect_error(get_mortality_hmd("SWE", 2000, "test", "test", shape = shape), "shape")
})

test_that("HMD annual years have calendar columns only for vector requests", {
  raw <- hmd_year_fixture()
  single <- mortality_from_hmd(raw, 2001)$rates
  expect_false(any(c("time_start", "time_end") %in% names(single)))
  out <- mortality_from_hmd(raw, c(2002, 2000))$rates
  expect_s3_class(out, "tbl_df")
  expect_equal(unique(out$time_start), c(2000, 2002))
  expect_equal(unique(out$time_end), c(2001, 2003))
  expect_equal(out$rate, c(raw$Male[1:3], raw$Female[1:3], raw$Male[7:9], raw$Female[7:9]))
  expect_error(mortality_from_hmd(raw, c(2000, 2004)), "requested year.*2004")
  bad <- raw
  bad$Male[7] <- NA_real_
  expect_error(mortality_from_hmd(bad, c(2000, 2002)))
  expect_s3_class(mortality_from_hmd(bad, 2000)$rates, "tbl_df")
})

test_that("all native HMD shapes retain source age and period widths", {
  for (shape in c("1x1", "1x5", "1x10", "5x1", "5x5", "5x10")) {
    age_width <- as.numeric(strsplit(shape, "x")[[1]][1])
    time_width <- as.numeric(strsplit(shape, "x")[[1]][2])
    starts <- c(2000, 2000 + time_width)
    periods <- if (time_width == 1) starts else
      paste0(starts, "-", starts + time_width - 1)
    raw <- if (age_width == 5) hmd_group_fixture(periods) else {
      data.frame(Year = rep(periods, each = 3), Age = rep(c("0", "1", "110+"), 2),
        Male = .02, Female = .01)
    }
    out <- mortality_from_hmd(raw, starts, shape)$rates
    expect_equal(unique(out$time_start), starts)
    expect_equal(unique(out$time_end), starts + time_width)
    if (age_width == 5) {
      expect_equal(head(out$age_start, 4), c(0, 1, 5, 110))
      expect_equal(head(out$age_end, 4), c(1, 5, 10, Inf))
    } else expect_equal(head(out$age_end, 3), c(1, 2, Inf))
    expect_true(all(out$rate[out$sex == "m"] == .02))
    expect_false("time_start" %in% names(mortality_from_hmd(raw, starts[1], shape)$rates))
  }
  expect_error(mortality_from_hmd(hmd_group_fixture(), 2001, "5x5"), "period starts")
  expect_error(mortality_from_hmd(hmd_group_fixture(), 2000, "1x5"), "age labels")
  bad <- hmd_group_fixture()
  bad$Age[2] <- "1-5"
  expect_error(mortality_from_hmd(bad, 2000, "5x5"), "inconsistent")
})

test_that("raw HMD periods preserve shortened endpoints and reject conflicts", {
  raw <- hmd_group_fixture(c("1998-1999", "2000-2003"))
  out <- mortality_from_hmd(raw, c(1998, 2000), "5x5")$rates
  expect_equal(unique(out$time_end), c(2000, 2004))
  for (periods in list(c("2000-2005", "2006-2009"),
      c("2000-2004", "2004-2008"), c("2000-1999", "2005-2009"),
      c("2000-2004", "nonsense"), c("2000-2004", "2000-2003")))
    expect_error(mortality_from_hmd(hmd_group_fixture(periods), 2000, "5x5"))
})

test_that("parsed HMD bounds use annual coverage and preserve open ages", {
  raw <- hmd_group_fixture(c("1998-1999", "2000-2003"))
  parsed <- raw
  parsed$Year <- rep(c(1998, 2000), each = 4)
  parsed$Age <- rep(c(0, 1, 5, 110), 2)
  parsed$OpenInterval <- rep(c(FALSE, FALSE, FALSE, TRUE), 2)
  expect_identical(mortality_from_hmd(parsed, c(1998, 2000), "5x5", 1998:2003),
    mortality_from_hmd(raw, c(1998, 2000), "5x5"))
  expect_error(mortality_from_hmd(parsed, 2000, "5x5"), "annual coverage")
  expect_error(mortality_from_hmd(parsed, 2000, "5x5", c(1998, 2000:2003)), "annual coverage")
  bad <- parsed
  bad$OpenInterval[4] <- NA
  expect_error(mortality_from_hmd(bad, 1998, "5x5", 1998:2003), "OpenInterval")
})

test_that("HMD web getters download each required product once and keep provenance", {
  skip_if_not_installed("HMDHFDplus")
  calls <- character()
  with_database_reader("readHMDweb", function(CNTRY, item, username, password, fixup) {
    calls <<- c(calls, item)
    expect_equal(CNTRY, "SWE")
    expect_false(fixup)
    hmd_year_fixture()
  }, {
    out <- get_mortality_hmd("SWE", c(2002, 2000), "test", "test")
    expect_identical(calls, "Mx_1x1")
    source <- attr(out, "source")
    expect_equal(source$year, c(2000, 2002))
    expect_equal(source$shape, "1x1")
    expect_equal(source$periods$time_end, c(2001, 2003))
    expect_false(any(c("username", "password") %in% names(source)))
    expect_identical(set_mortality(hammel_model(), out)$processes$mortality$source, source)
  })
  calls <- character()
  raw <- hmd_group_fixture(c("1998-1999", "2000-2003"))
  parsed <- raw
  parsed$Year <- rep(c(1998, 2000), each = 4)
  parsed$Age <- rep(c(0, 1, 5, 110), 2)
  parsed$OpenInterval <- rep(c(FALSE, FALSE, FALSE, TRUE), 2)
  with_database_reader("readHMDweb", function(CNTRY, item, ...) {
    calls <<- c(calls, item)
    if (item == "Mx_1x1") data.frame(Year = 1998:2003) else parsed
  }, {
    out <- get_mortality_hmd("SWE", c(1998, 2000), "test", "test", shape = "5x5")
    expect_identical(calls, c("Mx_5x5", "Mx_1x1"))
    expect_equal(unique(out$time_end), c(2000, 2004))
    expect_equal(attr(out, "source")$coverage_item, "Mx_1x1")
    expect_equal(attr(out, "source")$item, "Mx_5x5")
  })
  calls <- character()
  with_database_reader("readHMDweb", function(CNTRY, item, ...) {
    calls <<- c(calls, item)
    raw
  }, {
    out <- get_mortality_hmd("SWE", 1998, "test", "test", shape = "5x5")
    expect_identical(calls, "Mx_5x5")
    expect_false("time_start" %in% names(out))
    expect_equal(attr(out, "source")$periods$time_end, 2000)
  })
})

test_that("local grouped HMD parsing produces the same intervals as raw labels", {
  skip_if_not_installed("HMDHFDplus")
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("Synthetic HMD fixture", "", "Year Age Female Male Total",
    "2000-2003 0 .01 .02 .015", "2000-2003 1-4 .01 .02 .015",
    "2000-2003 5-9 .01 .02 .015", "2000-2003 110+ .01 .02 .015"), path)
  raw <- HMDHFDplus::readHMD(path, fixup = FALSE)
  parsed <- HMDHFDplus::readHMD(path, fixup = TRUE)
  expect_identical(mortality_from_hmd(raw, 2000, "5x5"),
    mortality_from_hmd(parsed, 2000, "5x5", 2000:2003))
})

test_that("HFD annual and conditional products select multiple years locally", {
  raw <- data.frame(Year = rep(2000:2002, each = 3), Age = rep(c("12-", "20", "55+"), 3),
    ASFR = seq_len(9) / 100)
  out <- fertility_from_hfd(raw, c(2002, 2000))$rates
  expect_equal(unique(out$time_start), c(2000, 2002))
  expect_equal(unique(out$time_end), c(2001, 2003))
  expect_equal(out$rate, raw$ASFR[c(1:3, 7:9)])
  expect_equal(out$age_end, rep(c(13, 21, Inf), 2))
  expect_error(fertility_from_hfd(raw, c(2000, 1990)), "requested year.*1990")
  for (column in c(paste0("m", 1:4, "x"), "m5px")) raw[[column]] <- raw$ASFR
  raw$m2x[1] <- NA_real_
  raw$m3x[7] <- NA_real_
  conditional <- fertility_by_parity_from_hfd(raw, c(2002, 2000))
  expect_equal(unique(conditional$time_start), c(2000, 2002))
  expect_equal(conditional$parity_start, rep(rep(0:4, each = 3), 2))
  expect_equal(sum(conditional$rate == 0), 2)
  expect_false(anyNA(conditional$rate))
  expect_s3_class(set_fertility(hammel_model(), conditional), "hammel_model")
  skip_if_not_installed("HMDHFDplus")
  calls <- 0
  with_database_reader("readHFDweb", function(CNTRY, item, username, password, fixup, Update) {
    calls <<- calls + 1
    expect_equal(item, "mi")
    expect_equal(Update, "20260702")
    expect_false(fixup)
    raw
  }, {
    rates <- get_fertility_hfd("SWE", c(2002, 2000), "test", "test",
      update = "20260702", by_parity = TRUE)
    expect_equal(calls, 1)
    expect_identical(rates$rate, conditional$rate)
    expect_equal(attr(rates, "source")$missing_rates_replaced, 2)
    expect_equal(attr(rates, "source")$year, c(2000, 2002))
    expect_equal(attr(rates, "source")$periods$time_end, c(2001, 2003))
  })
  calls <- 0
  with_database_reader("readHFDweb", function(...) {
    calls <<- calls + 1
    raw
  }, {
    rates <- get_fertility_hfd("SWE", 2000:2002, "test", "test")
    expect_equal(calls, 1)
    expect_equal(unique(rates$time_start), 2000:2002)
    expect_equal(rates$rate, raw$ASFR)
  })
})

test_that("prepared multi-year rates execute directly with native bands", {
  hmd <- data.frame(Year = rep(2000:2001, each = 2), Age = rep(c("0", "1+"), 2),
    Male = 0, Female = 0)
  hfd <- data.frame(Year = 2000:2001, Age = "20", ASFR = c(0, 10))
  m <- make_model(1, sex = "f", birth_interval = 2)
  m$population$persons$birth_time <- 1981
  m$population$end_time <- 2000
  m <- set_mortality(m, mortality_from_hmd(hmd, 2000:2001)$rates)
  m <- set_fertility(m, fertility_from_hfd(hfd, 2000:2001)$rates, birth_interval = 2)
  out <- run_simulation(m, 2, seed = 91)
  expect_equal(count_run_events(out, "births"), 1)
  expect_gte(tail(out$persons$birth_time, 1), 2001)
  expect_lt(tail(out$persons$birth_time, 1), 2002)
  # Selecting separated years leaves a real gap; preparation never fills it.
  hmd$Year[hmd$Year == 2001] <- 2002
  m <- set_mortality(m, mortality_from_hmd(hmd, c(2000, 2002))$rates)
  expect_error(run_simulation(m, 2), "calendar coverage gap")
})
