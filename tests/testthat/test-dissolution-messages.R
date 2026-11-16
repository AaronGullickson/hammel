expect_dissolution_error <- function(rates, ...) {
  error <- expect_error(union_dissolution(rates), class = "hammel_error_rates")
  for (part in c(...)) {
    expect_match(conditionMessage(error), part, fixed = TRUE)
  }
  invisible(error)
}

test_that("missing dissolution types identify calendar coverage rather than pair jargon", {
  rates <- data.frame(time_start=c(1990,1992),time_end=2020,
    union_type=c("cohab","marriage"),dur_start=0,dur_end=Inf,rate=.1)
  expect_dissolution_error(
    rates, 'No "marriage" rates', "In calendar period [1990, 1992)."
  )
  rates$time_start <- c(1992,1990)
  expect_dissolution_error(
    rates, 'No "cohab" rates', "In calendar period [1990, 1992)."
  )
  rates$time_start <- 1992
  expect_s3_class(union_dissolution(rates),"hammel_union_dissolution")
})

test_that("missing dissolution types identify relevant group pairs and open periods", {
  rates <- data.frame(time_start=1992.5,time_end=Inf,
    union_type="marriage",group_w="B",group_h="A",dur_start=0,dur_end=Inf,rate=.1)
  expect_dissolution_error(
    rates,
    'No "cohab" rates for group stratum group_w:group_h = B:A.',
    "In calendar period [1992.5, Inf)."
  )
  rates$time_start <- rates$time_end <- NULL
  expect_dissolution_error(
    rates, 'No "cohab" rates for group stratum group_w:group_h = B:A.'
  )
  rates$group_w <- rates$group_h <- NULL
  expect_dissolution_error(rates, 'No "cohab" rates.')
  rates$union_type <- "unknown"
  error <- expect_error(union_dissolution(rates), class = "hammel_error_table")
  expect_identical(error$argument, "rates$union_type")
})
