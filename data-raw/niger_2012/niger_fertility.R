# Build the bundled Niger fertility schedule from the retained UN extract.
# Run from the package root with tibble installed.

root <- "data-raw/niger_2012"
source_file <- file.path(root, "niger_2012_asfr_wpp2024.csv")
niger_fertility <- tibble::as_tibble(read.csv(source_file))
stopifnot(
  identical(names(niger_fertility), c("age_start", "age_end", "rate")),
  identical(niger_fertility$age_start, 15:49),
  identical(niger_fertility$age_end, 16:50),
  all(is.finite(niger_fertility$rate) & niger_fertility$rate >= 0),
  isTRUE(all.equal(sum(niger_fertility$rate), 7.468122))
)
save(
  niger_fertility,
  file = "data/niger_fertility.rda",
  compress = "xz",
  version = 2
)
