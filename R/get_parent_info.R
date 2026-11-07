#' Add recorded parent information to person histories
#'
#' Match mothers and fathers by person ID and append their ages at each child's
#' birth and their optional recorded groups. All persons, including founders
#' and deceased persons, retain their original row order and columns. The input
#' population is not modified and no known endpoint is required.
#'
#' Ages are in continuous years: the child's birth time minus the parent's birth
#' time. Unknown parents yield missing ages and groups. Each declared group
#' dimension adds the parents' recorded categories, with their type preserved;
#' undeclared columns are not added. These
#' are recorded identities, not reconstructions of historical group transitions.
#' Existing columns that conflict with the added fields cause an error.
#'
#' @param population A validated `hammel_population`, such as the result of
#'   [simulate()].
#' @return A separate tibble containing all person columns plus `age_mother`
#'   and `age_father`, measured at that person's birth, and for each declared
#'   group dimension `<group>_mother` and `<group>_father` (such as
#'   `race_mother`). Empty populations retain the
#'   enriched column types.
#' @export
#' @examples
#' population <- create_baseline_population(10)
#' get_parent_info(population)
get_parent_info <- function(population) {
  validate_population(population)
  persons <- population$persons
  result <- tibble::as_tibble(persons)
  dims <- names(population$groups)
  added <- c("age_mother", "age_father", paste0(rep(dims, each = 2), c("_mother", "_father")))
  conflicts <- intersect(names(result), added)
  if (length(conflicts)) {
    stop(sprintf("Person columns conflict with parent fields: %s.",
                 paste(conflicts, collapse = ", ")), call. = FALSE)
  }
  for (parent in c("mother", "father")) {
    rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
    result[[paste0("age_", parent)]] <- persons$birth_time - persons$birth_time[rows]
  }
  for (d in dims)
    for (parent in c("mother", "father")) {
      rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
      result[[paste0(d, "_", parent)]] <- persons[[d]][rows]
    }
  result
}
