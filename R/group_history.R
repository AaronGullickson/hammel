empty_transitions <- function() {
  tibble::tibble(person_id = numeric(), time = numeric(), group = character(),
                 origin = character(), destination = character(),
                 move_id = numeric(), initiator_id = numeric())
}

# Grouped populations always carry a transitions tibble. A history without a
# `group` column belongs to the only declared dimension.
normalize_population_transitions <- function(population) {
  dims <- names(population$groups)
  if (!length(dims)) return(population)
  tx <- population$transitions
  if (is.null(tx)) {
    population$transitions <- empty_transitions()
    return(population)
  }
  tx <- tibble::as_tibble(tx)
  if (!"group" %in% names(tx) && length(dims) == 1L) {
    tx$group <- rep(dims, nrow(tx))
    tx <- tx[c("person_id", "time", "group", setdiff(names(tx), c("person_id", "time", "group")))]
  }
  for (column in intersect(c("group", "origin", "destination"), names(tx)))
    if (is.factor(tx[[column]]) || is.numeric(tx[[column]])) tx[[column]] <- as.character(tx[[column]])
  # Move identifiers (D140); unknown for supplied rows without them.
  for (column in c("move_id", "initiator_id"))
    if (!column %in% names(tx)) tx[[column]] <- rep(NA_real_, nrow(tx))
  population$transitions <- tx
  population
}

validate_transition_history <- function(population) {
  tx <- population$transitions
  if (is.null(tx)) return(invisible(population))
  p <- population$persons
  old <- intersect(c("group_origin", "group_dest"), names(tx))
  if (length(old))
    stop("Transition histories use `group`, `origin` and `destination`, not ",
         paste0("`", old, "`", collapse = " and "), ".", call. = FALSE)
  if (nrow(tx) && !declared_group(population))
    stop("Transition histories require a declared group dimension.", call. = FALSE)
  if (!"group" %in% names(tx) && length(population$groups) == 1L)
    tx$group <- rep(names(population$groups), nrow(tx))
  if (!"group" %in% names(tx) && nrow(tx))
    stop("Transition histories need a `group` column naming each row's dimension when several are declared.",
         call. = FALSE)
  check_columns(tx, c("person_id", "time", "origin", "destination"), "transitions")
  check_numeric_column(tx$person_id, "transitions$person_id", id = TRUE)
  check_numeric_column(tx$time, "transitions$time")
  for (column in intersect(c("move_id", "initiator_id"), names(tx)))
    check_numeric_column(tx[[column]], paste0("transitions$", column), missing = TRUE)
  for (col in c("origin", "destination")) check_groups(data.frame(group = tx[[col]]), col)
  if (!nrow(tx)) return(invisible(population))
  unknown <- setdiff(unique(as.character(tx$group)), names(population$groups))
  if (anyNA(tx$group) || length(unknown))
    stop(sprintf("Transition histories name dimensions that are not declared: %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  i <- match(tx$person_id, p$person_id)
  if (anyNA(i) || any(tx$time < p$birth_time[i]) ||
      any(tx$time >= p$death_time[i], na.rm = TRUE) ||
      (!is.null(population$end_time) && any(tx$time > population$end_time)) ||
      any(as.character(tx$origin) == as.character(tx$destination)))
    stop("Transition history has invalid people, dates or self-transitions.", call. = FALSE)
  for (rows in split(seq_len(nrow(tx)), list(tx$person_id, tx$group), drop = TRUE)) {
    rows <- rows[order(tx$time[rows])]
    current <- as.character(p[[tx$group[rows[1]]]][i[rows[1]]])
    if (anyDuplicated(tx$time[rows]) ||
        (length(rows) > 1L && any(as.character(tx$destination[head(rows, -1L)]) !=
                                as.character(tx$origin[tail(rows, -1L)]))) ||
        as.character(tail(tx$destination[rows], 1)) != current)
      stop("Transition history must chain chronologically to each person's current category in each dimension.",
           call. = FALSE)
  }
  invisible(population)
}

# A population seen through one dimension, as if it were the only one and
# named `group`: observation code works on this view.
dimension_view <- function(population, dim) {
  if (dim != "group") population$persons$group <- population$persons[[dim]]
  population$groups <- stats::setNames(population$groups[dim], "group")
  if (!is.null(population$ancestry)) {
    a <- population$ancestry[[dim]]
    population["ancestry"] <- list(if (is.null(a)) NULL else list(group = a))
  }
  tx <- population$transitions
  if (!is.null(tx)) {
    if ("group" %in% names(tx)) tx <- tx[tx$group == dim, , drop = FALSE]
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
  if (is.null(tx) || !nrow(tx) || !length(rows)) return(out)
  if (is.null(histories)) histories <- split(seq_len(nrow(tx)), tx$person_id)
  for (j in seq_along(rows)) {
    if (is.na(rows[j])) next
    h <- histories[[as.character(p$person_id[rows[j]])]]
    if (!length(h)) next
    h <- h[order(tx$time[h])]
    prior <- if (left) tx$time[h] < times[j] else tx$time[h] <= times[j]
    out[j] <- if (any(prior)) tail(tx$destination[h[prior]], 1) else tx$origin[h[1]]
  }
  out
}
