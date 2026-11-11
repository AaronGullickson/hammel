# Extract recorded kin, count living kin over time, and validate kinship
# generation limits.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Extract known kin within generation-distance bounds
#'
#' @description
#'
#' Extract kin of a set of ego individuals from simulated data. Kin relations
#' are governed by generational distance.
#'
#' @details
#'
#' Find recorded common ancestors within `generations` parent steps of both an
#' ego and a relative, including each person at distance zero. This function
#' uses the same native kinship module as union eligibility. Unknown parents
#' terminate their branches;. At depth two, known parents/children,
#' grandparents/grandchildren, siblings, aunts/uncles, nieces/nephews, and first
#' cousins are included, including half relationships. Great-grandparents and
#' first cousins once removed require depth three.
#'
#' Returns one row per relative, describing their closest common ancestors:
#' those with the smallest sum of the two distances within the requested bounds.
#' `n_common_ancestors` counts them, ordinarily 2 for full kin and 1 for half
#' kin; a count of 1 is only recorded evidence when a parent is unknown.
#' `relationship` names the distance pair, such as `"sibling"`, `"aunt_uncle"`
#' or `"first_cousin"`, with full and half kin together. The ego itself is
#' excluded. Includes deceased and living people in recorded history; use
#' [count_kin()] for mean counts of living kin at given times.
#'
#' @param population A `hammel_population` with a finite `end_time`.
#' @param ego_ids Unique person IDs whose known kin should be extracted.
#' @param generations Nonnegative whole number of parent generations on EACH
#'   side, default `2`. Zero returns no kin between distinct people.
#' @param relative_ids Optional unique relative IDs to consider. `NULL` means
#'   all recorded persons. Unknown IDs are errors; empty vectors return no rows.
#' @return A tibble with one row per ego and relative: `ego_id`, `relative_id`,
#'   `relationship`, `ego_distance`, `relative_distance`, `n_common_ancestors`,
#'   `ego_lineage`, `relative_lineage`, `ego_ancestry_complete`, and
#'   `relative_ancestry_complete`. Lineage is `"self"`, `"maternal"`,
#'   `"paternal"`, or `"both"` for shortest routes to the closest common
#'   ancestors.
#' @seealso [count_kin()]
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
#' extract_kin(result, max(result$person$person_id))
extract_kin <- function(
  population,
  ego_ids,
  generations = 2,
  relative_ids = NULL
) {
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
      stop(
        sprintf(
          "`%s` must contain unique IDs present in the population.",
          name
        ),
        call. = FALSE
      )
    }
  }
  kin <- extract_kin_cpp(
    population$persons, ego_ids, relative_ids,
    as.integer(generations)
  )
  kin$relationship <- kin_relationship(kin$ego_distance, kin$relative_distance)
  tibble::as_tibble(
    kin[c(
      "ego_id", "relative_id", "relationship", "ego_distance",
      "relative_distance", "n_common_ancestors", "ego_lineage",
      "relative_lineage", "ego_ancestry_complete", "relative_ancestry_complete"
    )]
  )
}

#' Calculate the average number of living kin over time
#'
#' @description
#'
#' Count, for each person living at each requested time, the living relatives
#' recorded within `generations` parent steps on each side, by relationship,
#' and average them over living persons.
#'
#' @details
#'
#' Counts use recorded parent links only. A starting population without
#' recorded parents has no recorded kin, so counts are too low until two or
#' three generations have been simulated. Allow a burn-in before interpreting
#' them. Relatives share kin, so egos are not independent.
#'
#' The following named relationships are available:
#'
#' - "parent"
#' - "child"
#' - "grandparent"
#' - "grandchild"
#'- "sibling"
#' - "aunt_uncle"
#' - "niece_nephew"
#' - "first_cousin"
#' - "great_grandparent" (requires `generations = 3`)
#' - "great_grandchild" (requires `generations = 3`)
#' - "great_aunt_uncle" (requires `generations = 3`)
#' - "great_niece_nephew" (requires `generations = 3`)
#' - "first_cousin_once_removed" (requires `generations = 3`)
#' - "second_cousin" (requires `generations = 3`)
#'
#' @param population A `hammel_population` with a known `end_time` no earlier
#'   than the requested times.
#' @param times A numeric vector of finite time points in years.
#' @param generations Nonnegative whole number of parent generations on each
#'   side, default `2`, which covers parents, children, grandparents,
#'   grandchildren, siblings, aunts and uncles, nieces and nephews, and first
#'   cousins.
#' @param by `NULL` (default) or a character vector of the ego's
#'   characteristics at each time to count separately: `"sex"` and declared
#'   group dimensions, such as `"race"`.
#' @param relationships `NULL` (default) for all relationships, or a nonempty
#'   character vector of unique relationship names within `generations`, such
#'   as `"first_cousin"`. Retains the usual row order and ego denominators.
#' @return A tibble with `time`, the `by` columns, `relationship`, `persons`
#'   (living egos), `kin` (their living kin of that type), `mean_kin`
#'   (`kin / persons`) and `share_with_kin` (share of egos with at least one).
#'   `mean_kin` and `share_with_kin` are `NA` with no egos. Every time,
#'   stratum and selected relationship within `generations` has a row. Rows vary
#'   by time fastest, then the `by` variables in order, then relationship.
#' @seealso [extract_kin()]
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
#' result <- run_simulation(model, years = 100, seed = 42)
#' count_kin(result, 100, relationships = "first_cousin")
count_kin <- function(
  population,
  times,
  generations = 2,
  by = NULL,
  relationships = NULL
) {
  if (!is.numeric(times) || !is.null(dim(times)) || !length(times) ||
        any(!is.finite(times))) {
    stop(
      "`times` must be a nonempty numeric vector of finite times in years.",
      call. = FALSE
    )
  }
  if (!inherits(population, "hammel_population")) {
    stop("`population` must be a hammel_population.", call. = FALSE)
  }
  validate_population(population)
  if (max(times) > population$end_time) {
    stop(
      "Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  }
  validate_kinship_generations(generations, "generations")
  dims <- names(population$groups)
  by <- check_by(by, c("sex", dims), "count_kin")
  distances <- kin_relationship_types(generations)
  types <- unique(distances$relationship)
  if (!is.null(relationships)) {
    if (!is.character(relationships) || !is.null(dim(relationships)) ||
          !length(relationships) ||
          anyNA(relationships) ||
          anyDuplicated(relationships) ||
          any(!relationships %in% types)) {
      stop(
        paste0(
          "`relationships` must contain unique relationship names ",
          "available within `generations`."
        ),
        call. = FALSE
      )
    }
    types <- types[types %in% relationships]
  }
  lookup <- matrix(-1L, nrow = generations + 1L, ncol = generations + 1L)
  type <- match(distances$relationship, types)
  selected <- !is.na(type)
  lookup[cbind(
    distances$ego_distance[selected] + 1L,
    distances$relative_distance[selected] + 1L
  )] <- type[selected] - 1L
  n_types <- length(types)
  persons <- population$persons
  n <- nrow(persons)
  levels <- lapply(
    by, function(name) {
      if (name == "sex") {
        c("m", "f")
      } else {
        living_kin_categories(population, name)
      }
    }
  )
  names(levels) <- by
  n_strata <- prod(lengths(levels))
  birth <- persons$birth_time
  death <- persons$death_time
  strata <- matrix(0L, nrow = n, ncol = length(times))
  for (i in seq_along(times)) {
    time <- times[i]
    alive <- birth <= time & (is.na(death) | death > time)
    stratum <- rep(1L, n)
    stride <- 1L
    for (name in by) {
      values <- if (name == "sex") {
        persons$sex
      } else {
        groups_at(dimension_view(population, name), seq_len(n), rep(time, n))
      }
      stratum <- stratum + stride * (match(values, levels[[name]]) - 1L)
      stride <- stride * length(levels[[name]])
    }
    stratum[!alive] <- 0L
    strata[, i] <- as.integer(stratum)
  }
  counts <- count_kin_cpp(
    persons, as.integer(generations), lookup, strata, n_strata, n_types
  )
  grid <- expand.grid(
    time = seq_along(times),
    stratum = seq_len(n_strata),
    type = seq_len(n_types)
  )
  cell <- grid$time + length(times) * (
    grid$stratum - 1L + n_strata * (grid$type - 1L)
  )
  persons_n <- counts$persons[
    grid$time + length(times) * (grid$stratum - 1L)
  ]
  kin_n <- counts$kin[cell]
  with_n <- counts$with[cell]
  out <- tibble::tibble(time = times[grid$time])
  stride <- 1L
  for (name in by) {
    size <- length(levels[[name]])
    out[[name]] <- levels[[name]][((grid$stratum - 1L) %/% stride) %% size + 1L]
    stride <- stride * size
  }
  out$relationship <- types[grid$type]
  out$persons <- persons_n
  out$kin <- kin_n
  out$mean_kin <- ifelse(persons_n > 0, kin_n / persons_n, NA_real_)
  out$share_with_kin <- ifelse(persons_n > 0, with_n / persons_n, NA_real_)
  out
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Name kin from the generation distances of the ego and relative to their
# closest shared ancestors (see design/kinship.md).
kin_relationship <- function(ego_distance, relative_distance) {
  # Name each distinct distance pair once; kin tables repeat a few pairs.
  key <- as.numeric(ego_distance) * 2^31 + as.numeric(relative_distance)
  first <- !duplicated(key)
  if (sum(first) < length(key)) {
    labels <- kin_relationship(ego_distance[first], relative_distance[first])
    return(labels[match(key, key[first])])
  }
  greats <- function(k) strrep("great_", pmax(k, 0L))
  ordinal <- function(k) {
    words <- c("first", "second", "third", "fourth", "fifth")
    ifelse(k <= 5L, words[pmin(k, 5L)], paste0(k, "th"))
  }
  removed <- function(k) {
    words <- c("once", "twice", "thrice")
    ifelse(
      k == 0L, "",
      paste0(
        "_", ifelse(
          k <= 3L, words[pmax(pmin(k, 3L), 1L)],
          paste0(k, "_times")
        ), "_removed"
      )
    )
  }
  a <- as.integer(ego_distance)
  b <- as.integer(relative_distance)
  out <- character(length(a))
  up <- b == 0L
  out[up] <- ifelse(
    a[up] == 1L, "parent",
    paste0(greats(a[up] - 2L), "grandparent")
  )
  down <- a == 0L
  out[down] <- ifelse(
    b[down] == 1L, "child",
    paste0(greats(b[down] - 2L), "grandchild")
  )
  side <- a >= 1L & b >= 1L
  sibling <- side & a == 1L & b == 1L
  out[sibling] <- "sibling"
  aunt <- side & b == 1L & a >= 2L
  out[aunt] <- paste0(greats(a[aunt] - 2L), "aunt_uncle")
  niece <- side & a == 1L & b >= 2L
  out[niece] <- paste0(greats(b[niece] - 2L), "niece_nephew")
  cousin <- a >= 2L & b >= 2L
  out[cousin] <- paste0(
    ordinal(pmin(a[cousin], b[cousin]) - 1L), "_cousin",
    removed(abs(a[cousin] - b[cousin]))
  )
  out
}

# Every relationship within `generations`, in increasing distance sum, with
# the ego's ancestors before descendants at equal sums.
kin_relationship_types <- function(generations) {
  g <- seq.int(0L, as.integer(generations))
  pairs <- expand.grid(ego_distance = g, relative_distance = g)
  pairs <- pairs[pairs$ego_distance + pairs$relative_distance > 0L, ]
  pairs <- pairs[
    order(pairs$ego_distance + pairs$relative_distance, -pairs$ego_distance),
  ]
  pairs$relationship <- kin_relationship(
    pairs$ego_distance, pairs$relative_distance
  )
  rownames(pairs) <- NULL
  pairs
}

# Categories of a group dimension for living-kin strata: the declared order,
# or the sorted recorded values when none are declared.
living_kin_categories <- function(population, name) {
  declared <- population$groups[[name]]$categories
  if (!is.null(declared)) {
    return(declared)
  }
  tx <- population$transitions
  values <- population$persons[[name]]
  if (!is.null(tx) && nrow(tx)) {
    moves <- if ("group" %in% names(tx)) tx[tx$group == name, ] else tx
    values <- c(values, moves$origin, moves$destination)
  }
  sort(unique(as.character(values)), method = "radix")
}

validate_kinship_generations <- function(
  generations,
  name = "kinship_generations"
) {
  if (!is.numeric(generations) || length(generations) != 1L ||
        !is.null(dim(generations)) ||
        !is.finite(generations) ||
        generations < 0 ||
        generations != floor(generations) ||
        generations > .Machine$integer.max) {
    stop(
      sprintf(
        "`%s` must be a nonnegative whole number within integer range.", name
      ),
      call. = FALSE
    )
  }
  invisible(generations)
}
