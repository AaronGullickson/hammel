# Augment person and union histories and extract unpartnered births.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Extract union data from simulation results with additional detail
#'
#' @description
#'
#' Extract the unions data frame from a `hammel_population` object, augmented
#' with additional information.
#'
#' @details
#'
#' The unions record will match wives and husbands to the population's person
#' table by ID and return their recorded birth times, death times, and group
#' membership at formation time alongside every union record.
#'
#' @param population A `hammel_population` object containing person and union
#'   histories, such as the result of [run_simulation()].
#' @return A separate tibble with all union columns plus `birth_time_w`,
#'   `death_time_w`, `birth_time_h`, and `death_time_h`, plus `<group>_w` and
#'   `<group>_h` for each declared group dimension. Suffixes `_w` and `_h`
#'   identify wives and husbands. Empty union histories retain the enriched
#'   column types.
#' @export
#' @examples
#' population <- generate_uniform_population(1000)
#' model <- hammel_model() |>
#'   set_population(population) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_union_formation(
#'     "on_demand",
#'     matching = age_gap_matching_model(3, 5)
#'   )
#' result <- run_simulation(model, years = 50, seed = 42)
#' result |>
#'   augment_unions()
augment_unions <- function(population) {
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

#' Extract persons data from simulation results with additional detail
#'
#' @description
#'
#' Extract the persons data frame from a `hammel_population` object, augmented
#' with additional information.
#'
#' @details
#'
#' Mothers and fathers will be matched with children by person ID and their
#' age at the person's birth will be appended along with group membership.
#'
#' @param population A `hammel_population`, such as the result of
#'   [run_simulation()].
#' @return A separate tibble containing all person columns plus `age_mother` and
#'   `age_father`, measured at that person's birth, and for each declared group
#'   dimension `<group>_mother` and `<group>_father` (such as `race_mother`).
#'   Empty populations retain the enriched column types.
#' @export
#' @examples
#' population <- generate_uniform_population(1000)
#' model <- hammel_model() |>
#'   set_population(population) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_union_formation(
#'     "on_demand",
#'     matching = age_gap_matching_model(3, 5)
#'   )
#' result <- run_simulation(model, years = 50, seed = 42)
#' result |>
#'   augment_persons() |>
#'   tail()
augment_persons <- function(population) {
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
    result[[paste0("age_", parent)]] <-
      persons$birth_time - persons$birth_time[rows]
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
#' @description
#'
#' Return births where the mother was outside both cohabitation and marriage,
#' with parent ages and optional groups supplied by [augment_persons()].
#'
#' @param population A `hammel_population` with person and union
#'   histories, such as the result of [run_simulation()].
#' @return A separate tibble containing selected person rows and `age_mother`,
#'   `age_father`, plus `group_mother` and `group_father` when groups are
#'   present. Input order, column types and custom person columns are preserved.
#'   Empty results retain the enriched schema. The population is not modified.
#' @export
#' @examples
#' population <- generate_uniform_population(1000)
#' model <- hammel_model() |>
#'   set_population(population) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility_union) |>
#'   set_union_formation(
#'     "rate",
#'     rates = usa_2000_marriage,
#'     matching = age_gap_matching_model(3, 5)
#'   )
#' result <- run_simulation(model, years = 50, seed = 42)
#' extract_unpartnered_births(result)
extract_unpartnered_births <- function(population) {
  parents <- augment_persons(population)
  rows <- unpartnered_birth_rows_cpp(population$persons, population$unions)
  parents[rows, , drop = FALSE]
}
