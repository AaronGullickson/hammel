# Abridge the NCHS 1999-2001 U.S. decennial life tables to mortality rates in
# bands 0-1, 1-5, five-year bands to 100 and an open 100+ band, by sex. Run
# from the package root. Requires readxl. See
# usa_1999_2001_mortality.README.md.

root <- "validation/data/usa_2000_2001"
sources <- c(
  m = "usa_life_table_1999_2001_males.xls",
  f = "usa_life_table_1999_2001_females.xls"
)
titles <- c(
  m = "Table 2. Life table for males: United States, 1999-2001",
  f = "Table 3. Life table for females: United States, 1999-2001"
)
starts <- c(0, 1, seq(5, 100, 5))
ends <- c(starts[-1], Inf)

abridge <- function(sex) {
  x <- readxl::read_excel(
    file.path(root, sources[[sex]]),
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  stopifnot(
    identical(x[[1]][1], titles[[sex]]),
    identical(unname(unlist(x[3, 2:6])), c("nqx", "lx", "ndx", "nLx", "Tx"))
  )
  # Single-year rows follow the "Years" heading; earlier rows are in days.
  rows <- seq(which(x[[1]] == "Years") + 1L, length.out = 110L)
  ages <- strsplit(x[[1]][rows], "-", fixed = TRUE)
  age <- as.numeric(vapply(ages, `[`, character(1), 1L))
  stopifnot(identical(age, as.numeric(0:109)), identical(
    as.numeric(vapply(ages, `[`, character(1), 2L)), as.numeric(1:110)
  ))
  lx <- as.numeric(x[[3]][rows])
  dx <- as.numeric(x[[4]][rows])
  big_l <- as.numeric(x[[5]][rows])
  big_t <- as.numeric(x[[6]][rows])
  stopifnot(
    all(is.finite(c(lx, dx, big_l, big_t))),
    lx[1] == 1e5,
    isTRUE(all.equal(head(lx, -1) - dx[-110], lx[-1], tolerance = 1e-6))
  )
  # Closed bands: deaths over person-years lived in the band. Open band:
  # survivors at 100 over their remaining person-years (l100 / T100).
  rate <- vapply(seq_along(starts), function(i) {
    if (is.infinite(ends[i])) {
      return(lx[age == starts[i]] / big_t[age == starts[i]])
    }
    band <- age >= starts[i] & age < ends[i]
    sum(dx[band]) / sum(big_l[band])
  }, numeric(1))
  data.frame(sex = sex, age_start = starts, age_end = ends, rate = rate)
}

mortality <- rbind(abridge("f"), abridge("m"))
stopifnot(
  nrow(mortality) == 44L,
  all(is.finite(mortality$rate)),
  all(mortality$rate > 0)
)
write.csv(
  mortality,
  file.path(root, "usa_1999_2001_mortality.csv"),
  row.names = FALSE
)
