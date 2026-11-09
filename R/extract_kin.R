# Extract recorded kin and validate kinship generation limits.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Extract known kin within generation-distance bounds
#'
#' Find recorded common ancestors within `generations` parent steps of both an
#' ego and a relative, including each person at distance zero. Uses the same
#' native kinship module as union eligibility. Unknown parents terminate their
#' branches; missing evidence does not establish unrelatedness. At depth two,
#' known parents/children, grandparents/grandchildren, siblings, aunts/uncles,
#' nieces/nephews, and first cousins are included, including half relationships.
#' Great-grandparents and first cousins once removed require depth three.
#'
#' Returns the closest common ancestors by minimum sum of the two distances
#' within the requested bounds. Equally close ancestors have separate rows;
#' count distinct relative IDs when counting people. Does not assign
#' relationship categories or infer full versus half kin from missing links. The
#' ego itself is excluded. Includes deceased and living people in recorded
#' history, without automatic filtering to a census date or living kin. Ancestry
#' cycles are errors.
#'
#' @param population A `hammel_population` or a persons data frame, including
#'   deceased persons and recorded parent IDs. No known endpoint is required.
#' @param ego_ids Unique person IDs whose known kin should be extracted.
#' @param generations Nonnegative whole number of parent generations on EACH
#'   side, default `2`. Zero returns no kin between distinct people.
#' @param relative_ids Optional unique relative IDs to consider. `NULL` means
#'   all recorded persons. Unknown IDs are errors; empty vectors return no rows.
#' @return A tibble with `ego_id`, `relative_id`, `common_ancestor_id`,
#'   `ego_distance`, `relative_distance`, `ego_lineage`, `relative_lineage`,
#'   `ego_ancestry_complete`, and `relative_ancestry_complete`. Lineage is
#'   `"self"`, `"maternal"`, `"paternal"`, or `"both"` for shortest routes.
#'   Completeness refers only to required parent links within the chosen depth.
#'   Empty results retain the same column types; no rows does not prove
#'   unrelatedness or complete ancestry. Input ego/relative order is retained;
#'   tied ancestor IDs are sorted.
#' @export
#' @examples
#' persons <- data.frame(
#'   person_id = c(10, 20, 30), sex = c("f", "m", "f"),
#'   birth_time = c(-50, -50, -20), death_time = NA_real_,
#'   mother_id = c(NA, NA, 10), father_id = c(NA, NA, 20)
#' )
#' population <- hammel_population(persons, end_time = 0)
#' extract_kin(population, ego_ids = 30)
extract_kin <- function(
  population,
  ego_ids,
  generations = 2,
  relative_ids = NULL
) {
  if (is.data.frame(population)) {
    population <- hammel_population(population)
  }
  validate_population(population)
  validate_kinship_generations(generations, "generations")
  ids <- population$persons$person_id
  if (is.null(relative_ids)) {
    relative_ids <- ids
  }
  for (name in c("ego_ids", "relative_ids")) {
    values <- get(name)
    check_numeric_column(values, name, id = TRUE)
    if (anyDuplicated(values) || any(!values %in% ids)) {
      stop(sprintf(
        "`%s` must contain unique IDs present in the population.",
        name
      ), call. = FALSE)
    }
  }
  tibble::as_tibble(extract_kin_cpp(
    population$persons, ego_ids, relative_ids,
    as.integer(generations)
  ))
}

# HELPER FUNCTIONS ---------------------------------------------------------

validate_kinship_generations <- function(
  generations,
  name = "kinship_generations"
) {
  if (!is.numeric(generations) || length(generations) != 1L ||
    !is.null(dim(generations)) || !is.finite(generations) || generations <
    0 ||
    generations != floor(generations) || generations > .Machine$integer.max) {
    stop(
      sprintf(
        "`%s` must be a nonnegative whole number within integer range.", name
      ),
      call. = FALSE
    )
  }
  invisible(generations)
}
