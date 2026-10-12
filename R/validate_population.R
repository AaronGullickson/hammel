empty_unions <- function() {
  data.frame(
    union_id = integer(), partner1_id = integer(), partner2_id = integer(),
    start_time = numeric(), end_time = numeric(), end_reason = character()
  )
}

check_time <- function(time, argument) {
  if (!is.null(time) &&
      (!is.numeric(time) || length(time) != 1L || !is.finite(time))) {
    stop(sprintf("`%s` must be NULL or a finite time in years.", argument),
         call. = FALSE)
  }
}

check_columns <- function(table, columns, name) {
  if (!is.data.frame(table)) {
    stop(sprintf("`%s` must be a data frame.", name), call. = FALSE)
  }
  if (anyDuplicated(names(table))) {
    stop(sprintf("`%s` must have unique column names.", name), call. = FALSE)
  }
  absent <- setdiff(columns, names(table))
  if (length(absent)) {
    stop(sprintf("`%s` is missing columns: %s.", name,
                 paste(absent, collapse = ", ")), call. = FALSE)
  }
}

check_numeric_column <- function(x, name, missing = FALSE, id = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) ||
      any(is.nan(x)) || (!missing && anyNA(x)) ||
      any(!is.finite(x[!is.na(x)])) ||
      (id && any(x[!is.na(x)] <= 0 | x[!is.na(x)] != trunc(x[!is.na(x)])))) {
    stop(sprintf("`%s` must contain %s%s.", name,
                 if (id) "positive whole-number IDs" else "finite numeric values",
                 if (missing) " or NA" else " without missing values"),
         call. = FALSE)
  }
}

validate_population <- function(persons, unions, time = NULL) {
  check_columns(persons, c("person_id", "sex", "birth_time", "death_time",
                          "mother_id", "father_id"), "persons")
  if (is.null(unions) ||
      (is.data.frame(unions) && nrow(unions) == 0L && ncol(unions) == 0L)) {
    unions <- empty_unions()
  }
  check_columns(unions, names(empty_unions()), "unions")
  check_numeric_column(persons$person_id, "persons$person_id", id = TRUE)
  if (anyDuplicated(persons$person_id)) {
    stop("`persons$person_id` must be unique.", call. = FALSE)
  }
  check_numeric_column(persons$sex, "persons$sex")
  if (any(!persons$sex %in% c(0, 1))) {
    stop("`persons$sex` must be 0 (male) or 1 (female).", call. = FALSE)
  }
  for (field in c("birth_time", "death_time")) {
    check_numeric_column(persons[[field]], paste0("persons$", field),
                         missing = field == "death_time")
  }
  if (any(persons$death_time < persons$birth_time, na.rm = TRUE)) {
    stop("Death times must not precede birth times.", call. = FALSE)
  }
  for (field in c("mother_id", "father_id")) {
    x <- persons[[field]]
    check_numeric_column(x, paste0("persons$", field), missing = TRUE, id = TRUE)
    if (any(!is.na(x) & (!x %in% persons$person_id | x == persons$person_id))) {
      stop(sprintf("`persons$%s` must refer to another person or be NA.", field),
           call. = FALSE)
    }
  }
  for (field in c("union_id", "partner1_id", "partner2_id")) {
    check_numeric_column(unions[[field]], paste0("unions$", field), id = TRUE)
  }
  if (anyDuplicated(unions$union_id)) {
    stop("`unions$union_id` must be unique.", call. = FALSE)
  }
  if (any(!unions$partner1_id %in% persons$person_id |
          !unions$partner2_id %in% persons$person_id |
          unions$partner1_id == unions$partner2_id)) {
    stop("Union partners must refer to two distinct persons.", call. = FALSE)
  }
  check_numeric_column(unions$start_time, "unions$start_time")
  check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
  if (any(unions$end_time < unions$start_time, na.rm = TRUE)) {
    stop("Union end times must not precede start times.", call. = FALSE)
  }
  reason <- unions$end_reason
  if (!is.character(reason) || !is.null(dim(reason)) ||
      any(!is.na(reason) & !reason %in% c("death", "dissolution")) ||
      any(is.na(reason) != is.na(unions$end_time))) {
    stop("`unions$end_reason` must be death or dissolution for ended unions, and NA for ongoing unions.",
         call. = FALSE)
  }
  for (partner in c("partner1_id", "partner2_id")) {
    born <- persons$birth_time[match(unions[[partner]], persons$person_id)]
    if (any(unions$start_time < born)) {
      stop("Unions must not start before either partner's birth.", call. = FALSE)
    }
  }
  if (!is.null(time) &&
      any(c(persons$birth_time, persons$death_time, unions$start_time,
            unions$end_time) > time, na.rm = TRUE)) {
    stop("Recorded events must not occur after the population's reference time.",
         call. = FALSE)
  }
  unions
}
