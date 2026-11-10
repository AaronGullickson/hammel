# Extract parent, partner, and unpartnered-birth information from histories.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Add recorded partner information to union histories
#'
#' Match wives and husbands to the population's person table by ID and return
#' their recorded birth times, death times, and optional groups at formation
#' alongside every union record. Union row order and all existing columns are
#' retained. The population and its union table are not modified.
#'
#' Dates remain in years. Missing death times denote people alive at the
#' recorded endpoint. Each declared group dimension adds its wife's and
#' husband's categories (`<group>_w`, `<group>_h`); undeclared columns are not
#' added; otherwise their types are preserved. Groups are reconstructed at each
#' union's start from retained transitions; later moves do not relabel earlier
#' unions. A transition at the exact formation time follows formation, so its
#' origin group applies. Without transition history, use recorded person groups.
#' Existing columns whose names conflict with added partner fields cause an
#' error.
#'
#' @param population A validated `hammel_population` containing person and union
#'   histories, such as the result of [run_simulation()].
#' @return A separate tibble with all union columns plus `birth_time_w`,
#'   `death_time_w`, `birth_time_h`, and `death_time_h`, plus `<group>_w` and
#'   `<group>_h` for each declared group dimension. Suffixes `_w` and `_h`
#'   identify wives and husbands. Empty union histories retain the enriched
#'   column types.
#' @export
#' @examples
#' population <- generate_uniform_population(10)
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
    stop(
      sprintf(
        "Union columns conflict with partner fields: %s.",
        paste(conflicts, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # Categories at each union's start, replaying each dimension's transitions.
  formation <- list()
  population <- normalize_population_transitions(population)
  tx <- population$transitions
  for (d in dims) {
    if (is.null(tx) || !any(tx$group == d)) {
      next
    }
    view <- dimension_view(population, d)
    indexed <- prepare_group_observations(view$persons, view, TRUE)
    groups <- observation_groups(indexed, TRUE)
    codes <- union_start_groups_cpp(
      indexed, population$unions, match(view$persons$group, groups) - 1L
    )
    labels <- groups
    # Histories hold labels as text; keep the persons column's type.
    if (is.numeric(persons[[d]])) {
      storage.mode(labels) <- storage.mode(persons[[d]])
    }
    formation[[d]] <- list(
      w = labels[codes[, 1L] + 1L],
      h = labels[codes[, 2L] + 1L]
    )
  }
  for (suffix in c("w", "h")) {
    id_column <- if (suffix == "w") "wife_id" else "husband_id"
    rows <- match(result[[id_column]], persons$person_id)
    for (column in columns) {
      result[[paste0(column, "_", suffix)]] <- if (!is.null(
                                                      formation[[column]]
                                                    )) {
        formation[[column]][[suffix]]
      } else {
        persons[[column]][rows]
      }
    }
  }
  result
}

#' Add recorded parent information to person histories
#'
#' Match mothers and fathers by person ID and append their ages at each child's
#' birth and their optional recorded groups. All persons, including founders and
#' deceased persons, retain their original row order and columns. The input
#' population is not modified.
#'
#' Ages are in continuous years: the child's birth time minus the parent's birth
#' time. Unknown parents yield missing ages and groups. Each declared group
#' dimension adds the parents' recorded categories, with their type preserved;
#' undeclared columns are not added. These are recorded identities, not
#' reconstructions of historical group transitions. Existing columns that
#' conflict with the added fields cause an error.
#'
#' @param population A validated `hammel_population`, such as the result of
#'   [run_simulation()].
#' @return A separate tibble containing all person columns plus `age_mother` and
#'   `age_father`, measured at that person's birth, and for each declared group
#'   dimension `<group>_mother` and `<group>_father` (such as `race_mother`).
#'   Empty populations retain the enriched column types.
#' @export
#' @examples
#' population <- generate_uniform_population(10)
#' get_parent_info(population)
get_parent_info <- function(population) {
  validate_population(population)
  persons <- population$persons
  result <- tibble::as_tibble(persons)
  dims <- names(population$groups)
  added <- c(
    "age_mother", "age_father", paste0(
      rep(dims, each = 2), c("_mother", "_father")
    )
  )
  conflicts <- intersect(names(result), added)
  if (length(conflicts)) {
    stop(
      sprintf(
        "Person columns conflict with parent fields: %s.",
        paste(conflicts, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (parent in c("mother", "father")) {
    rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
    result[[paste0("age_", parent)]] <- persons$birth_time - persons$birth_time[
      rows
    ]
  }
  for (d in dims) {
    for (parent in c("mother", "father")) {
      rows <- match(persons[[paste0(parent, "_id")]], persons$person_id)
      result[[paste0(d, "_", parent)]] <- persons[[d]][rows]
    }
  }
  result
}

#' Extract births with known parents to unpartnered mothers
#'
#' Return births where the mother was outside both cohabitation and marriage,
#' with parent ages and optional groups supplied by [get_parent_info()]. Known
#' mother and father IDs are required. Each birth is retained separately,
#' including repeated parental pairs; missing-parent records and founders are
#' omitted. These rows can be used to examine nonpartnership father matching.
#'
#' Histories must be complete for the intended interpretation. Partnership
#' coverage is half-open: births at formation are partnered. Births exactly at
#' dissolution or conversion retain the ending partnership (birth-first ties).
#' Births at death endings are unpartnered unless another episode covers that
#' time. Overlapping episodes do not duplicate births. Untyped partnerships are
#' treated as partnered too. No time window is imposed and no endpoint is
#' needed.
#'
#' @param population A validated `hammel_population` with person and union
#'   histories, such as the result of [run_simulation()].
#' @return A separate tibble containing selected person rows and `age_mother`,
#'   `age_father`, plus `group_mother` and `group_father` when groups are
#'   present. Input order, column types and custom person columns are preserved.
#'   Empty results retain the enriched schema. The population is not modified.
#' @export
#' @examples
#' population <- generate_uniform_population(10)
#' get_unpartnered_births(population)
get_unpartnered_births <- function(population) {
  parents <- get_parent_info(population)
  rows <- unpartnered_birth_rows_cpp(population$persons, population$unions)
  parents[rows, , drop = FALSE]
}
