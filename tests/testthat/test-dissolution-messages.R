test_that("missing dissolution types identify calendar coverage rather than pair jargon", {
  rates <- data.frame(time_start=c(1990,1992),time_end=2020,
    union_type=c("cohab","marriage"),dur_start=0,dur_end=Inf,rate=.1)
  expect_error(union_dissolution(rates),
    "missing a 'marriage' schedule during calendar interval \\[1990, 1992\\)")
  rates$time_start <- c(1992,1990)
  expect_error(union_dissolution(rates),
    "missing a 'cohab' schedule during calendar interval \\[1990, 1992\\)")
  rates$time_start <- 1992
  expect_s3_class(union_dissolution(rates),"hammel_union_dissolution")
})

test_that("missing dissolution types identify relevant group pairs and open periods", {
  rates <- data.frame(time_start=1992.5,time_end=Inf,
    union_type="marriage",group_w="B",group_h="A",dur_start=0,dur_end=Inf,rate=.1)
  expect_error(union_dissolution(rates),
    "missing a 'cohab' schedule during calendar interval \\[1992.5, Inf\\) for group stratum 'B:A'")
  rates$time_start <- rates$time_end <- NULL
  expect_error(union_dissolution(rates),"missing a 'cohab' schedule for group stratum 'B:A'")
  rates$group_w <- rates$group_h <- NULL
  expect_error(union_dissolution(rates),"missing a 'cohab' schedule\\. Supply both")
  rates$union_type <- "unknown"
  expect_error(union_dissolution(rates),"union_type.*cohab and marriage labels")
})
