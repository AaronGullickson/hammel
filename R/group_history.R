# Normalize, validate, and query recorded group-transition histories.

# HELPER FUNCTIONS ---------------------------------------------------------

empty_transitions <- function() {
  tibble::tibble(
    person_id = numeric(),
    time = numeric(),
    group = character(),
    origin = character(),
    destination = character(),
    move_id = numeric(),
    initiator_id = numeric()
  )
}

# Grouped populations always carry a transitions tibble. A history without a
# `group` column belongs to the only declared dimension.
normalize_population_transitions <- function(population) {
  dims <- names(population$groups)
  if (!length(dims)) {
    return(population)
  }
  tx <- population$transitions
  if (is.null(tx)) {
    population$transitions <- empty_transitions()
    if (all(arrival_columns %in% names(population$persons))) {
      population$transitions$pre_migration <- logical()
    }
    return(population)
  }
  tx <- tibble::as_tibble(tx)
  counted <- counted_dimensions(population$groups)
  if (!"group" %in% names(tx) && length(counted) == 1L) {
    tx$group <- rep(counted, nrow(tx))
    tx <- tx[c(
      "person_id", "time", "group", setdiff(
        names(tx), c("person_id", "time", "group")
      )
    )]
  }
  for (column in intersect(c("group", "origin", "destination"), names(tx))) {
    if (is.factor(tx[[column]]) || is.numeric(tx[[column]])) {
      tx[[column]] <- as.character(tx[[column]])
    }
  }
  # Move identifiers (D140); unknown for supplied rows without them.
  for (column in c("move_id", "initiator_id")) {
    if (!column %in% names(tx)) {
      tx[[column]] <- rep(NA_real_, nrow(tx))
    }
  }
  population$transitions <- tx
  population
}

validate_transition_history <- function(
  population,
  call = rlang::caller_env()
) {
  tx <- population$transitions
  if (is.null(tx)) {
    return(invisible(population))
  }
  p <- population$persons
  old <- intersect(c("group_origin", "group_dest"), names(tx))
  if (length(old)) {
    abort_hammel(
      c(
        paste(
          "Transition histories use {.field group}, {.field origin} and",
          "{.field destination} columns."
        ),
        "x" = "Found {.field {old}}.",
        "i" = "Rename these columns."
      ),
      "population",
      argument = "transitions",
      columns = old,
      call = call
    )
  }
  if (nrow(tx) && !declared_group(population)) {
    abort_hammel(
      c(
        "Transition histories need a declared group dimension.",
        "i" = "Declare dimensions with {.arg groups} or {.fn set_groups}."
      ),
      "groups",
      argument = "transitions",
      call = call
    )
  }
  if (!"group" %in% names(tx) && length(population$groups) == 1L) {
    tx$group <- rep(names(population$groups), nrow(tx))
  }
  if (!"group" %in% names(tx) && nrow(tx)) {
    dimensions <- names(population$groups)
    abort_hammel(
      c(
        "Transition histories need a {.field group} column.",
        "x" = "Several dimensions are declared: {.val {dimensions}}.",
        "i" = "Name each row's dimension in {.field group}."
      ),
      "population",
      argument = "transitions",
      call = call
    )
  }
  check_columns(
    tx, c("person_id", "time", "origin", "destination"), "transitions",
    call = call
  )
  check_numeric_column(
    tx$person_id, "transitions$person_id",
    id = TRUE,
    call = call
  )
  check_numeric_column(tx$time, "transitions$time", call = call)
  for (column in intersect(c("move_id", "initiator_id"), names(tx))) {
    check_numeric_column(
      tx[[column]], paste0("transitions$", column),
      missing = TRUE,
      call = call
    )
  }
  # An NA origin marks entry into a group, such as an immigrant's assigned
  # category at arrival (D186); it is checked below as a first row.
  entry <- is.na(tx$origin)
  check_group_labels(tx$origin[!entry], "transitions$origin", call = call)
  check_group_labels(tx$destination, "transitions$destination", call = call)
  if (!nrow(tx)) {
    return(invisible(population))
  }
  unknown <- setdiff(unique(as.character(tx$group)), names(population$groups))
  if (anyNA(tx$group) || length(unknown)) {
    rows <- which(is.na(tx$group) | !tx$group %in% names(population$groups))
    abort_hammel(
      c(
        "Transition histories must name declared dimensions.",
        "x" = "Undeclared or missing {.field group} at {format_rows(rows)}."
      ),
      "groups",
      argument = "transitions",
      rows = rows,
      call = call
    )
  }
  i <- match(tx$person_id, p$person_id)
  abort_population_rows(
    which(is.na(i)),
    "Transition histories must refer to persons in {.field persons}.",
    "Unknown {.field person_id}",
    table = "transitions",
    call = call
  )
  abort_population_rows(
    which(
      tx$time < p$birth_time[i] |
        (!is.na(p$death_time[i]) & tx$time >= p$death_time[i]) |
        (!is.null(population$end_time) & tx$time > population$end_time)
    ),
    paste(
      "Transitions must happen while the person is alive, by the",
      "population's {.field end_time}."
    ),
    "Transitions outside a lifetime",
    table = "transitions",
    call = call
  )
  abort_population_rows(
    which(as.character(tx$origin) == as.character(tx$destination)),
    "Transitions must move between different categories.",
    "{.field origin} equals {.field destination}",
    table = "transitions",
    call = call
  )
  for (rows in split(
    seq_len(nrow(tx)), list(tx$person_id, tx$group),
    drop = TRUE
  )) {
    rows <- rows[order(tx$time[rows])]
    current <- as.character(p[[tx$group[rows[1]]]][i[rows[1]]])
    chained <- !any(entry[rows[-1L]]) && (length(rows) < 2L || all(
      as.character(tx$destination[head(rows, -1L)]) ==
        as.character(tx$origin[tail(rows, -1L)])
    ))
    if (anyDuplicated(tx$time[rows]) || !chained ||
          as.character(tail(tx$destination[rows], 1)) != current) {
      person <- tx$person_id[rows[1]]
      dimension <- tx$group[rows[1]]
      abort_population_rows(
        sort(rows),
        paste(
          "Each person's transitions must form a chronological chain ending",
          "in their current category."
        ),
        paste0(
          "The ", dimension, " transitions of person ", person,
          " don't chain"
        ),
        table = "transitions",
        call = call
      )
    }
  }
  invisible(population)
}

# A population seen through one dimension, as if it were the only one and
# named `group`: observation code works on this view.
dimension_view <- function(population, dim) {
  if (dim != "group") {
    population$persons$group <- population$persons[[dim]]
  }
  population$groups <- stats::setNames(population$groups[dim], "group")
  if (!is.null(population$ancestry)) {
    a <- population$ancestry[[dim]]
    population["ancestry"] <- list(if (is.null(a)) NULL else list(group = a))
  }
  tx <- population$transitions
  if (!is.null(tx)) {
    if ("group" %in% names(tx)) {
      tx <- tx[tx$group == dim, , drop = FALSE]
    }
    tx$group <- rep("group", nrow(tx))
    population$transitions <- tx
  }
  population
}

# Return identity at an observation time without changing the persons table.
groups_at <- function(population, rows, times, left = FALSE, histories = NULL) {
  p <- population$persons
  out <- p$group[rows]
  tx <- population$transitions
  if (is.null(tx) || !nrow(tx) || !length(rows)) {
    return(out)
  }
  if (is.null(histories)) {
    histories <- split(seq_len(nrow(tx)), tx$person_id)
  }
  for (j in seq_along(rows)) {
    if (is.na(rows[j])) {
      next
    }
    h <- histories[[as.character(p$person_id[rows[j]])]]
    if (!length(h)) {
      next
    }
    h <- h[order(tx$time[h])]
    prior <- if (left) tx$time[h] < times[j] else tx$time[h] <= times[j]
    out[j] <- if (any(prior)) {
      tail(tx$destination[h[prior]], 1)
    } else {
      tx$origin[h[1]]
    }
  }
  out
}
