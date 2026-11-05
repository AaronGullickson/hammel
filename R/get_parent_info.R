#' Add recorded parent information to person histories
#'
#' Match mothers and fathers by person ID and append their ages at each child's
#' birth and their optional recorded groups. All persons, including founders
#' and deceased persons, retain their original row order and columns. The input
#' population is not modified and no known endpoint is required.
#'
#' Ages are in continuous years: the child's birth time minus the parent's birth
#' time. Unknown parents yield missing ages and groups. Group fields are omitted
#' unless a group dimension is declared; otherwise its type is preserved. These
#' are recorded identities, not reconstructions of historical group transitions.
#' Existing columns that conflict with the added fields cause an error.
#'
#' @param population A validated `hammel_population`, such as the result of
#'   [simulate()].
#' @return A separate tibble containing all person columns plus `age_mother`
#'   and `age_father`, measured at that person's birth. When groups are present,
#'   also `group_mother` and `group_father`. Empty populations retain the
#'   enriched column types.
#' @export
#' @examples
#' population <- create_baseline_population(10)
#' get_parent_info(population)
get_parent_info <- function(population) {
  validate_population(population)
  persons <- population$persons
  result <- tibble::as_tibble(persons)
  grouped <- declared_group(population)
  added <- c("age_mother", "age_father",
             if (grouped) c("group_mother", "group_father"))
  conflicts <- intersect(names(result), added)
  if (length(conflicts)) {
    stop(sprintf("Person columns conflict with parent fields: %s.",
                 paste(conflicts, collapse = ", ")), call. = FALSE)
  }
  for (parent in c("mother", "father")) {
    rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
    result[[paste0("age_", parent)]] <- persons$birth_time - persons$birth_time[rows]
  }
  if (grouped) {
    for (parent in c("mother", "father")) {
      rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
      result[[paste0("group_", parent)]] <- persons$group[rows]
    }
  }
  result
}
