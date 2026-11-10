# Extract recorded kin, count living kin over time, and validate kinship
# generation limits.

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
#' Returns one row per relative, describing their closest common ancestors:
#' those with the smallest sum of the two distances within the requested
#' bounds. `n_common_ancestors` counts them, ordinarily 2 for full kin and 1
#' for half kin; a count of 1 is only recorded evidence when a parent is
#' unknown. `relationship` names the distance pair, such as `"sibling"`,
#' `"aunt_uncle"` or `"first_cousin"`, with full and half kin together. The ego
#' itself is excluded. Includes deceased and living people in recorded history;
#' use [count_kin()] for living kin at given times. Ancestry cycles
#' are errors.
#'
#' @param population A `hammel_population` or a persons data frame, including
#'   deceased persons and recorded parent IDs. No known endpoint is required.
#' @param ego_ids Unique person IDs whose known kin should be extracted.
#' @param generations Nonnegative whole number of parent generations on EACH
#'   side, default `2`. Zero returns no kin between distinct people.
#' @param relative_ids Optional unique relative IDs to consider. `NULL` means
#'   all recorded persons. Unknown IDs are errors; empty vectors return no rows.
#' @return A tibble with one row per ego and relative: `ego_id`,
#'   `relative_id`, `relationship`, `ego_distance`, `relative_distance`,
#'   `n_common_ancestors`, `ego_lineage`, `relative_lineage`,
#'   `ego_ancestry_complete`, and `relative_ancestry_complete`. Lineage is
#'   `"self"`, `"maternal"`, `"paternal"`, or `"both"` for shortest routes to
#'   the closest common ancestors. Completeness refers only to required parent
#'   links within the chosen depth. Empty results retain the same column
#'   types; no rows does not prove unrelatedness or complete ancestry. Input
#'   ego/relative order is retained.
#' @seealso [count_kin()]
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
#' A person is living at time `t` when `birth_time <= t` and `death_time` is
#' missing or after `t`, as in [calculate_population_size()]. Every living
#' person is an ego, and a relative counts when both are living at `t`.
#' Relatives are found and labelled as in [extract_kin()]: each counts once,
#' under the relationship of their closest shared ancestors, with full and half
#' kin together.
#'
#' Counts use recorded parent links only. A starting population without
#' recorded parents has no recorded kin, so counts are too low until two or
#' three generations have been simulated. Allow a burn-in before interpreting
#' them. Relatives share kin, so egos are not independent; use replicate
#' simulations to describe uncertainty.
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
#' @return A tibble with `time`, the `by` columns, `relationship`, `persons`
#'   (living egos), `kin` (their living kin of that type), `mean_kin`
#'   (`kin / persons`) and `share_with_kin` (share of egos with at least one).
#'   `mean_kin` and `share_with_kin` are `NA` with no egos. Every time,
#'   stratum and relationship within `generations` has a row. Rows vary by
#'   time fastest, then the `by` variables in order, then relationship.
#' @seealso [extract_kin()]
#' @export
#' @examples
#' persons <- data.frame(
#'   person_id = 1:4, sex = c("f", "m", "f", "m"),
#'   birth_time = c(-50, -50, -20, -18), death_time = c(-5, NA, NA, NA),
#'   mother_id = c(NA, NA, 1, 1), father_id = c(NA, NA, 2, 2)
#' )
#' population <- hammel_population(persons, end_time = 0)
#' count_kin(population, times = c(-10, 0))
count_kin <- function(
  population,
  times,
  generations = 2,
  by = NULL
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
  if (is.null(population$end_time) || max(times) > population$end_time) {
    stop(
      "Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  }
  validate_kinship_generations(generations, "generations")
  dims <- names(population$groups)
  by <- check_by(by, c("sex", dims), "count_kin")
  persons <- population$persons
  ids <- persons$person_id
  n <- length(ids)
  kin <- extract_kin_cpp(persons, ids, ids, as.integer(generations))
  ego <- match(kin$ego_id, ids)
  relative <- match(kin$relative_id, ids)
  # Labels pool distance pairs that share a name, such as (2,3) and (3,2).
  types <- unique(kin_relationship_types(generations)$relationship)
  type <- match(
    kin_relationship(kin$ego_distance, kin$relative_distance), types
  )
  n_types <- length(types)
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
  results <- lapply(
    times, function(time) {
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
      keep <- alive[ego] & alive[relative]
      cell <- (type[keep] - 1L) * n_strata + stratum[ego[keep]]
      # Egos with at least one living relative of each type.
      pairs <- unique((ego[keep] - 1L) * n_types + type[keep])
      pair_ego <- (pairs - 1L) %/% n_types + 1L
      pair_type <- (pairs - 1L) %% n_types + 1L
      list(
        persons = tabulate(stratum[alive], n_strata),
        kin = tabulate(cell, n_strata * n_types),
        with = tabulate(
          (pair_type - 1L) * n_strata + stratum[pair_ego],
          n_strata * n_types
        )
      )
    }
  )
  grid <- expand.grid(
    time = seq_along(times),
    stratum = seq_len(n_strata),
    type = seq_len(n_types)
  )
  cell <- (grid$type - 1L) * n_strata + grid$stratum
  persons_n <- vapply(
    seq_len(nrow(grid)), function(i) {
      results[[grid$time[i]]]$persons[grid$stratum[i]]
    }, integer(1)
  )
  kin_n <- vapply(
    seq_len(nrow(grid)), function(i) {
      results[[grid$time[i]]]$kin[cell[i]]
    }, integer(1)
  )
  with_n <- vapply(
    seq_len(nrow(grid)), function(i) {
      results[[grid$time[i]]]$with[cell[i]]
    }, integer(1)
  )
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
  pairs <- pairs[order(
    pairs$ego_distance + pairs$relative_distance, -pairs$ego_distance
    ), ]
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
