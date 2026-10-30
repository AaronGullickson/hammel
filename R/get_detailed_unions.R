#' Add recorded partner information to union histories
#'
#' Match wives and husbands to the population's person table by ID and return
#' their recorded birth times, death times, and optional groups at formation alongside every
#' union record. Union row order and all existing columns are retained. The
#' population and its union table are not modified. No known endpoint is needed.
#'
#' Dates remain in years. Missing death times denote people alive at the recorded
#' endpoint. Group columns are omitted when persons has no `group` column;
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
#'   `death_time_w`, `birth_time_h`, and `death_time_h`. When group is present,
#'   also `group_w` and `group_h`. Suffixes `_w` and `_h` identify wives and
#'   husbands. Empty union histories retain the enriched column types.
#' @export
#' @examples
#' population <- create_baseline_population(10)
#' get_detailed_unions(population)
get_detailed_unions <- function(population) {
  validate_population(population)
  persons <- population$persons
  result <- tibble::as_tibble(population$unions)
  columns <- c("birth_time", "death_time", if ("group" %in% names(persons)) "group")
  added <- c(paste0(columns, "_w"), paste0(columns, "_h"))
  conflicts <- intersect(names(result), added)
  if (length(conflicts)) {
    stop(sprintf("Union columns conflict with partner fields: %s.",
                 paste(conflicts, collapse = ", ")), call. = FALSE)
  }
  formation_groups <- NULL
  if ("group" %in% names(persons) && !is.null(population$transitions) && nrow(population$transitions)) {
    indexed <- prepare_group_observations(persons, population, TRUE)
    groups <- observation_groups(indexed, TRUE)
    formation_groups <- union_start_groups_cpp(indexed, population$unions,
      match(persons$group, groups) - 1L)
  }
  for (suffix in c("w", "h")) {
    id_column <- if (suffix == "w") "wife_id" else "husband_id"
    rows <- match(result[[id_column]], persons$person_id)
    for (column in columns) {
      result[[paste0(column, "_", suffix)]] <- if (column == "group" && !is.null(formation_groups))
        groups[formation_groups[, if (suffix == "w") 1L else 2L] + 1L] else persons[[column]][rows]
    }
  }
  result
}
