# Ancestry shares are independent of current group identity. Missing ancestry
# belongs in the explicit unknown component, never in NA values or missing rows.
validate_ancestry <- function(ancestry, persons) {
  if (is.null(ancestry)) return(invisible(ancestry))
  check_columns(ancestry, c("person_id", "unknown"), "ancestry")
  check_numeric_column(ancestry$person_id, "ancestry$person_id", id = TRUE)
  if (anyDuplicated(ancestry$person_id) || nrow(ancestry) != nrow(persons) ||
      !setequal(ancestry$person_id, persons$person_id)) {
    stop("`ancestry` must have exactly one row per person, with matching person IDs.", call. = FALSE)
  }
  columns <- setdiff(names(ancestry), c("person_id", "unknown"))
  if (any(!startsWith(columns, "group_")) ||
      any(!nzchar(trimws(substring(columns, 7L))))) {
    stop("Known ancestry columns must be named `group_<label>` with nonempty labels.", call. = FALSE)
  }
  for (column in c(columns, "unknown")) {
    values <- ancestry[[column]]
    if (!is.numeric(values) || is.complex(values) || !is.null(dim(values)) ||
        anyNA(values) || any(!is.finite(values)) || any(values < 0 | values > 1)) {
      stop("Ancestry proportions must be finite numeric values between 0 and 1, without missing values.",
           call. = FALSE)
    }
  }
  sums <- rowSums(as.matrix(ancestry[c(columns, "unknown")]))
  if (any(abs(sums - 1) > 1e-8)) {
    stop("Ancestry proportions, including `unknown`, must sum to 1 for each person.", call. = FALSE)
  }
  invisible(ancestry)
}

normalize_population_ancestry <- function(population) {
  persons <- population$persons
  ancestry <- population[["ancestry"]]
  if (is.null(ancestry)) {
    if (!"group" %in% names(persons)) return(population)
    labels <- as.character(persons$group)
    ancestry <- tibble::tibble(person_id = persons$person_id)
    for (group in unique(labels)) ancestry[[paste0("group_", group)]] <- as.numeric(labels == group)
    ancestry$unknown <- rep(0, nrow(persons))
  } else {
    columns <- setdiff(names(ancestry), c("person_id", "unknown"))
    ancestry <- ancestry[match(persons$person_id, ancestry$person_id),
                         c("person_id", columns, "unknown"), drop = FALSE]
  }
  population$ancestry <- tibble::as_tibble(ancestry)
  population
}
