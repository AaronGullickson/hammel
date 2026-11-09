hmd_open_fixture <- function() {
  mx <- data.frame(Year = rep(2000:2001, each = 5),
    Age = rep(c("0", "1", "2", "3", "4+"), 2),
    Male = rep(c(.01, .02, NA, NA, NA), 2),
    Female = rep(c(.005, .01, NA, NA, NA), 2))
  deaths <- exposures <- mx
  deaths$Male <- c(1, 1, 2, 3, 0, 1, 1, 4, 6, 0)
  deaths$Female <- c(1, 1, 1, 2, 0, 1, 1, 2, 1, 0)
  exposures$Male <- c(100, 50, 100, 20, 0, 100, 50, 50, 10, 0)
  exposures$Female <- c(200, 100, 80, 20, 0, 200, 100, 40, 10, 0)
  list(mx = mx, deaths = deaths, exposures = exposures)
}
pooled_hmd_fixture <- function(x = hmd_open_fixture(), years = 2000:2001, age = 2,
                               shape = "1x1", coverage = NULL) {
  mortality_from_hmd(x$mx, years, shape, coverage, age, x$deaths, x$exposures)
}

test_that("HMD open-age pooling uses ratios of sums separately by sex and period", {
  x <- hmd_open_fixture()
  out <- pooled_hmd_fixture(x)
  r <- out$rates
  expect_s3_class(r, "tbl_df")
  expect_equal(r$age_start, rep(c(0, 1, 2), 4))
  expect_equal(r$age_end, rep(c(1, 2, Inf), 4))
  expect_equal(r$rate, c(.01, .02, 5 / 120, .005, .01, 3 / 100,
    .01, .02, 10 / 60, .005, .01, 3 / 50))
  expect_equal(r$time_start, rep(c(2000, 2001), each = 6))
  expect_equal(r$time_end, rep(c(2001, 2002), each = 6))
  expect_equal(out$pooled_tail$deaths, c(5, 3, 10, 3))
  expect_equal(out$pooled_tail$person_years, c(120, 100, 60, 50))
  expect_s3_class(out$pooled_tail, "tbl_df")
  single <- pooled_hmd_fixture(x, years = 2000)
  expect_false(any(c("time_start", "time_end") %in% names(single$rates)))
  expect_equal(single$pooled_tail$time_start, c(2000, 2000))
  expect_error(mortality_from_hmd(x$mx, 2000), "year 2000, sex m, age \\[2, 3\\)")
  # Pooling does not change the input or hide missing younger rates.
  expect_true(all(is.na(x$mx$Male[x$mx$Age == "4+"])))
  bad <- x
  bad$mx$Female[1] <- NA_real_
  expect_error(pooled_hmd_fixture(bad), "year 2000, sex f, age \\[0, 1\\)")
})

test_that("cutoffs require valid source age boundaries", {
  for (age in list(NA_real_, Inf, -1, 2.5, "2", numeric(), c(2, 3), matrix(2)))
    expect_error(get_mortality_hmd("SWE", 2000, "test", "test", open_age = age), "open_age")
  expect_error(get_mortality_hmd("SWE", 2000, "test", "test", shape = "5x1",
    open_age = 97), "align")
  expect_error(pooled_hmd_fixture(age = 5), "complete source tail")
  x <- hmd_open_fixture()
  # A zero-age cutoff pools the entire life course without a younger block.
  all <- pooled_hmd_fixture(x, years = 2000, age = 0)
  expect_equal(all$rates$age_start, c(0, 0))
  expect_equal(all$rates$age_end, c(Inf, Inf))
  expect_equal(all$rates$rate, c(7 / 270, 5 / 400))
  # The existing open source band can itself be the cutoff, if defined.
  x$exposures$Male[x$exposures$Age == "4+"] <- 2
  x$exposures$Female[x$exposures$Age == "4+"] <- 3
  x$mx$Male[is.na(x$mx$Male)] <- .1
  x$mx$Female[is.na(x$mx$Female)] <- .1
  out <- pooled_hmd_fixture(x, years = 2000, age = 4)
  expect_equal(out$rates$rate[out$rates$age_start == 4], c(0, 0))
})

test_that("pooling requires finite nonnegative counts and positive aggregate exposure", {
  for (source in c("deaths", "exposures")) for (value in c(NA_real_, Inf, -1)) {
    x <- hmd_open_fixture()
    x[[source]]$Male[3] <- value
    expect_error(pooled_hmd_fixture(x), "year 2000, sex m, age \\[2, 3\\)")
  }
  x <- hmd_open_fixture()
  x$exposures$Female[3:5] <- 0
  expect_error(pooled_hmd_fixture(x), "positive finite exposure.*year 2000, sex f")
  x <- hmd_open_fixture()
  x$deaths$Male[3:5] <- 0
  out <- pooled_hmd_fixture(x)
  expect_equal(out$rates$rate[out$rates$time_start == 2000 & out$rates$sex == "m" &
    out$rates$age_start == 2], 0)
  # Counts/exposures below the cutoff and in unselected years do not enter sums.
  x <- hmd_open_fixture()
  x$deaths$Female[c(1, 8)] <- NA_real_
  x$exposures$Male[c(1, 8)] <- NA_real_
  expect_s3_class(pooled_hmd_fixture(x, years = 2000)$rates, "tbl_df")
  x <- hmd_open_fixture()
  x$deaths$Male <- as.character(x$deaths$Male)
  expect_error(pooled_hmd_fixture(x), "numeric")
})

test_that("pooling rejects missing, duplicate, mismatched and incomplete tail bands", {
  for (source in c("deaths", "exposures")) {
    x <- hmd_open_fixture()
    x[[source]] <- x[[source]][-3, ]
    expect_error(pooled_hmd_fixture(x), "bands do not match")
    x <- hmd_open_fixture()
    x[[source]] <- rbind(x[[source]], x[[source]][3, ])
    expect_error(pooled_hmd_fixture(x), "bands do not match")
    x <- hmd_open_fixture()
    x[[source]]$Age[5] <- "5+"
    expect_error(pooled_hmd_fixture(x), "bands do not match")
    x <- hmd_open_fixture()
    x[[source]] <- x[[source]][x[[source]]$Year != 2001, ]
    expect_error(pooled_hmd_fixture(x), "missing requested source periods")
  }
  x <- hmd_open_fixture()
  x$mx <- x$mx[-4, ]
  expect_error(pooled_hmd_fixture(x), "complete source tail")
})

test_that("pooling supports native grouped ages and raw or parsed calendar periods", {
  x <- hmd_open_fixture()
  for (source in names(x)) {
    x[[source]]$Age <- rep(c("0", "1-4", "95-99", "100-104", "105+"), 2)
    x[[source]]$Year <- rep(c("2000-2004", "2005-2008"), each = 5)
  }
  out <- pooled_hmd_fixture(x, c(2000, 2005), 95, "5x5")
  expect_equal(out$rates$age_start, rep(c(0, 1, 95), 4))
  expect_equal(out$rates$age_end, rep(c(1, 5, Inf), 4))
  expect_equal(unique(out$rates$time_end), c(2005, 2009))
  parsed <- lapply(x, function(d) {
    d$Age <- rep(c(0, 1, 95, 100, 105), 2)
    d$OpenInterval <- rep(c(FALSE, FALSE, FALSE, FALSE, TRUE), 2)
    d$Year <- rep(c(2000, 2005), each = 5)
    d
  })
  expect_identical(pooled_hmd_fixture(parsed, c(2000, 2005), 95, "5x5", 2000:2008), out)
  bad <- x
  bad$exposures$Year[6:10] <- "2005-2009"
  expect_error(pooled_hmd_fixture(bad, c(2000, 2005), 95, "5x5"), "bands do not match")
})

test_that("getters fetch pooling products only on request and retain safe provenance", {
  skip_if_not_installed("HMDHFDplus")
  x <- hmd_open_fixture()
  calls <- character()
  with_database_reader("readHMDweb", function(CNTRY, item, username, password, fixup) {
    calls <<- c(calls, item)
    expect_equal(CNTRY, "SWE")
    expect_false(fixup)
    switch(item, Mx_1x1 = x$mx, Deaths_1x1 = x$deaths, Exposures_1x1 = x$exposures,
      stop("Unexpected source product"))
  }, {
    r <- get_mortality_hmd("SWE", 2000:2001, "test", "test", open_age = 2)
    expect_identical(calls, c("Mx_1x1", "Deaths_1x1", "Exposures_1x1"))
    source <- attr(r, "source")
    expect_equal(source$open_age, 2)
    expect_equal(source$open_age_method, "pooled_deaths_exposure")
    expect_equal(unname(source$pooling_items), c("Deaths_1x1", "Exposures_1x1"))
    expect_equal(source$pooled_tail$deaths, c(5, 3, 10, 3))
    expect_false(any(c("username", "password") %in% names(source)))
    expect_identical(set_mortality(hammel_model(), r)$processes$mortality$source, source)
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    saveRDS(r, path)
    expect_identical(readRDS(path), r)
    calls <- character()
    expect_error(get_mortality_hmd("SWE", 2000, "test", "test"), "year 2000")
    expect_identical(calls, "Mx_1x1")
  })
})

test_that("pooled complete schedules run directly in the continuous engine", {
  r <- pooled_hmd_fixture()$rates
  m <- make_model(1000, age = 2, death_rate = 0)
  m$population$persons$birth_time <- 1998
  m$population$end_time <- 2000
  out <- run_simulation(set_mortality(m, r), 2, seed = 911)
  expected <- exp(-5 / 120 - 10 / 60)
  expect_lt(abs(mean(is.na(out$persons$death_time)) - expected),
    5 * sqrt(expected * (1 - expected) / 1000))
})

test_that("pooling aligns source rows by intervals rather than row positions", {
  x <- hmd_open_fixture()
  x$deaths <- x$deaths[10:1, ]
  x$exposures <- x$exposures[c(5:1, 10:6), ]
  expect_identical(pooled_hmd_fixture(x), pooled_hmd_fixture())
})

test_that("grouped-time pooling fetches each product once and preserves short periods", {
  skip_if_not_installed("HMDHFDplus")
  x <- lapply(hmd_open_fixture(), function(d) {
    d$Year <- rep(c(2000, 2005), each = 5)
    d$Age <- rep(0:4, 2)
    d$OpenInterval <- rep(c(FALSE, FALSE, FALSE, FALSE, TRUE), 2)
    d
  })
  calls <- character()
  with_database_reader("readHMDweb", function(CNTRY, item, ...) {
    calls <<- c(calls, item)
    switch(item, Mx_1x5 = x$mx, Mx_1x1 = data.frame(Year = 2000:2008),
      Deaths_1x5 = x$deaths, Exposures_1x5 = x$exposures,
      stop("Unexpected product"))
  }, {
    out <- get_mortality_hmd("SWE", c(2000, 2005), "test", "test",
      shape = "1x5", open_age = 2)
    expect_identical(calls, c("Mx_1x5", "Mx_1x1", "Deaths_1x5", "Exposures_1x5"))
    expect_equal(unique(out$time_end), c(2005, 2009))
    source <- attr(out, "source")
    expect_equal(source$coverage_item, "Mx_1x1")
    expect_equal(source$pooled_tail$time_end, c(2005, 2005, 2009, 2009))
    expect_equal(source$pooled_tail$rate, c(5 / 120, 3 / 100, 10 / 60, 3 / 50))
  })
})
