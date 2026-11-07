#' Add recorded partner information to union histories
#'
#' Match wives and husbands to the population's person table by ID and return
#' their recorded birth times, death times, and optional groups at formation alongside every
#' union record. Union row order and all existing columns are retained. The
#' population and its union table are not modified. No known endpoint is needed.
#'
#' Dates remain in years. Missing death times denote people alive at the recorded
#' endpoint. Each declared group dimension adds its wife's and husband's
#' categories (`<group>_w`, `<group>_h`); undeclared columns are not added;
#' otherwise their types are preserved. Groups are reconstructed at each union's
#' start from retained transitions; later moves do not relabel earlier unions.
#' A transition at the exact formation time follows formation, so its origin
#' group applies. Without transition history, use recorded person groups.
#' Existing columns whose
#' names conflict with added partner fields cause an error.
#'
#' @param population A validated `hammel_population` containing person and union
#'   histories, such as the result of [simulate()].
#' @return A separate tibble with all union columns plus `birth_time_w`,
#'   `death_time_w`, `birth_time_h`, and `death_time_h`, plus `<group>_w` and
#'   `<group>_h` for each declared group dimension. Suffixes `_w` and `_h` identify wives and
#'   husbands. Empty union histories retain the enriched column types.
#' @export
#' @examples
#' population <- create_baseline_population(10)
#' get_detailed_unions(population)
get_detailed_unions <- function(population) {
  validate_population(population)
  persons <- population$persons
  result <- tibble::as_tibble(population$unions)
  dims <- names(population$groups)
  columns <- c("birth_time", "death_time", dims)
  added <- c(paste0(columns, "_w"), paste0(columns, "_h"))
  conflicts <- intersect(names(result), added)
  if (length(conflicts)) {
    stop(sprintf("Union columns conflict with partner fields: %s.",
                 paste(conflicts, collapse = ", ")), call. = FALSE)
  }
  # Categories at each union's start, replaying each dimension's transitions.
  formation <- list()
  population <- normalize_population_transitions(population)
  tx <- population$transitions
  for (d in dims) {
    if (is.null(tx) || !any(tx$group == d)) next
    view <- dimension_view(population, d)
    indexed <- prepare_group_observations(view$persons, view, TRUE)
    groups <- observation_groups(indexed, TRUE)
    codes <- union_start_groups_cpp(indexed, population$unions, match(view$persons$group, groups) - 1L)
    labels <- groups
    # Histories hold labels as text; keep the persons column's type.
    if (is.numeric(persons[[d]])) storage.mode(labels) <- storage.mode(persons[[d]])
    formation[[d]] <- list(w = labels[codes[, 1L] + 1L], h = labels[codes[, 2L] + 1L])
  }
  for (suffix in c("w", "h")) {
    id_column <- if (suffix == "w") "wife_id" else "husband_id"
    rows <- match(result[[id_column]], persons$person_id)
    for (column in columns) {
      result[[paste0(column, "_", suffix)]] <- if (!is.null(formation[[column]]))
        formation[[column]][[suffix]] else persons[[column]][rows]
    }
  }
  result
}
