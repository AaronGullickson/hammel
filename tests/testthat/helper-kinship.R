make_kinship_population <- function() {
  # All IDs are deliberately unrelated to row positions; retain dead ancestors.
  persons <- data.frame(
    person_id = c(10, 20, 30, 40, 50, 60, 70, 100, 80, 85, 86, 90, 95, 110, 130, 120, 150),
    sex = c("f", "m", "f", "m", "f", "m", "m", "f", "f", "m", "m", "m", "m", "f", "f", "m", "m"),
    birth_time = -c(80, 80, 80, 80, 50, 50, 50, 50, 25, 25, 25, 25, 25, 1, 40, 1, 25),
    death_time = NA_real_,
    mother_id = c(NA, NA, NA, NA, 10, 10, 30, NA, 50, 50, 50, 100, NA, 80, NA, 130, NA),
    father_id = c(NA, NA, NA, NA, 20, 20, 40, NA, 70, 70, NA, 60, 60, NA, NA, 90, NA))
  hammel_population(persons, end_time = 0, groups = if ("group" %in% names(persons)) "group")
}
