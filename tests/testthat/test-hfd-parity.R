hfd_parity_fixture <- function() {
  data.frame(Year=c(1999,2000,2000,2000,2000),
    Age=c("20","12-","20","21","55+"),
    m1x=c(.9,.001,.12,.13,.001),m2x=c(.8,NA,.20,.21,.001),
    m3x=c(.7,NA,.15,.16,.001),m4x=c(.6,NA,.10,.11,.001),
    m5x=c(.5,NA,.08,.09,.001))
}

test_that("HFD conditional rates map birth order to half-open current parity", {
  rates <- hammel:::fertility_by_parity_from_hfd(hfd_parity_fixture(),2000)
  expect_s3_class(rates,"tbl_df")
  expect_named(rates,c("age_start","age_end","parity_start","parity_end","rate"))
  expect_equal(rates$age_start,rep(c(20,21),5))
  expect_equal(rates$age_end,rep(c(21,22),5))
  expect_equal(rates$parity_start,rep(0:4,each=2))
  expect_equal(rates$parity_end,rep(c(1,2,3,4,Inf),each=2))
  expect_equal(rates$rate,c(.12,.13,.20,.21,.15,.16,.10,.11,.08,.09))
  expect_error(fertility(rates),"not yet implemented")
})

test_that("conditional preparation retains missing cells and explicit zeros", {
  raw <- hfd_parity_fixture()
  raw$m4x[3] <- NA_real_; raw$m5x[3] <- 0
  expect_warning(rates <- hammel:::fertility_by_parity_from_hfd(raw,2000),"retained as NA")
  expect_true(is.na(rates$rate[7]))
  expect_equal(rates$rate[9],0)
  expect_equal(nrow(rates),10)
})

test_that("conditional preparation rejects incompatible products and invalid values", {
  raw <- hfd_parity_fixture()
  expect_error(hammel:::fertility_by_parity_from_hfd(raw,1900),"requested year")
  for (column in c("m1x","m5x")) {
    bad <- raw; bad[[column]] <- NULL
    expect_error(hammel:::fertility_by_parity_from_hfd(bad,2000),column)
  }
  for (value in c(-.1,Inf)) {
    bad <- raw; bad$m2x[3] <- value
    expect_error(hammel:::fertility_by_parity_from_hfd(bad,2000))
  }
  bad <- raw; bad$m1x <- as.character(bad$m1x)
  expect_error(hammel:::fertility_by_parity_from_hfd(bad,2000),"numeric")
  bad <- raw; bad$Age[3] <- "20-24"
  expect_error(hammel:::fertility_by_parity_from_hfd(bad,2000),"single-year")
  expect_error(hammel:::fertility_by_parity_from_hfd(raw[c(2,5), ],2000),"No closed")
  expect_error(hammel:::fertility_by_parity_from_hfd(rbind(raw,raw[3, ]),2000),"duplicate")
  for (bad in list(NULL,NA,1,c(TRUE,FALSE),matrix(TRUE)))
    expect_error(get_fertility_hfd("USA",2000,"test","test",by_parity=bad),"by_parity")
})

test_that("conditional HFD reader preserves values, revision and safe provenance", {
  skip_if_not_installed("HMDHFDplus")
  with_database_reader("readHFDweb",function(CNTRY,item,username,password,fixup,Update) {
    expect_equal(CNTRY,"USA")
    expect_equal(item,"mi")
    expect_false(fixup)
    expect_equal(Update,"20260702")
    hfd_parity_fixture()
  }, {
    rates <- get_fertility_hfd("USA",2000,"test","test",update="20260702",by_parity=TRUE)
    plain <- rates; attr(plain,"source") <- NULL
    expect_identical(plain,hammel:::fertility_by_parity_from_hfd(hfd_parity_fixture(),2000))
    source <- attr(rates,"source")
    expect_equal(source$item,"mi")
    expect_equal(source$update,"20260702")
    expect_equal(source$rate_type,"conditional_age_parity")
    expect_equal(source$missing_rates,"preserve_NA")
    expect_equal(source$open_ages,"omit")
    expect_false(any(c("username","password") %in% names(source)))
    path <- tempfile(); on.exit(unlink(path),add=TRUE)
    saveRDS(rates,path)
    expect_identical(readRDS(path),rates)
  })
  with_database_reader("readHFDweb",function(...) stop("Conditional product not available"), {
    expect_error(get_fertility_hfd("USA",2000,"test","test",by_parity=TRUE),"not available")
  })
})

test_that("local HFD parsing handles conditional headers and missing dots", {
  skip_if_not_installed("HMDHFDplus")
  path <- tempfile(fileext=".txt"); on.exit(unlink(path),add=TRUE)
  writeLines(c("Synthetic conditional HFD fixture","", "Year Age m1x m2x m3x m4x m5x",
    "2000 12- 0.001 . . . .", "2000 20 0.12 0.20 0.15 0.10 0.08",
    "2000 21 0.13 0.21 0.16 0.11 0.09", "2000 55+ 0.001 . . . ."),path)
  raw <- HMDHFDplus::readHFD(path,fixup=FALSE,item="mi")
  expect_equal(hammel:::fertility_by_parity_from_hfd(raw,2000),
    hammel:::fertility_by_parity_from_hfd(hfd_parity_fixture(),2000))
})
