#' Specify individual group transitions
#'
#' Construct or attach competing origin-to-destination hazards. Only the
#' individual changes group; ancestry, parentage and union histories are retained.
#' Duration is years since the latest recorded transition, or age without one.
#' Missing directed pairs, strata and gaps have zero risk. Restarts retain
#' transition histories and resume durations from them.
#'
#' @param weights NULL (default) for unweighted selection, or a
#'   [create_weight_model()] using person variables and `group_dest`, which
#'   chooses who in a cell moves without changing the cell's input rate.
#' @param rates A nonempty data frame requiring `group_origin`, `group_dest`,
#'   `rate`, and either `age_start`/`age_end` or `dur_start`/`dur_end` throughout.
#'   Half-open bands use years; rates are finite nonnegative hazards per
#'   person-year. Ends may be `Inf`. Bands cannot overlap within a stratum.
#'   Self-transitions are invalid. Optional `sex` is m/f. Optional `union_type`
#'   is unpartnered/partnered or unpartnered/cohab/marriage. Omitted columns
#'   apply universally. Destinations may initially be empty; other configured
#'   group-specific processes must cover all possible groups. Optional
#'   `time_start`/`time_end` specify half-open absolute calendar-year bands,
#'   including fractional years. Every active slice must be valid and the table
#'   must cover the requested run. Rate changes do not reset group duration.
#' @param model A `hammel_model`.
#' @return `group_transition()` returns a process with a rates tibble;
#'   `set_group_transition()` returns the updated model. Grouped simulation
#'   output always includes a `transitions` tibble with `person_id`, `time`,
#'   `group_origin`, and `group_dest`, including when empty.
#' @export
group_transition <- function(rates, weights = NULL) {
  process <- new_rate_process("group_transition", rates)
  validate_group_transition(process)
  validate_process_weights(weights, "transition")
  warn_cancelled_weight_terms(weights, weight_strata("transition", rates), "group transition")
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname group_transition
#' @export
set_group_transition <- function(model, rates, weights = NULL) {
  set_process(model, group_transition(rates, weights))
}

validate_group_transition <- function(process) {
  if (!identical(process$name, "group_transition") ||
      !identical(process$rate_unit, "person_year"))
    stop("Invalid group-transition specification or units.", call. = FALSE)
  r <- process$rates
  check_stratified_columns(r, c("group_origin", "group_dest", "sex", "union_type"), "Group-transition")
  if (validate_calendar_rates(r, function(x) {
    copy <- process; copy$rates <- x; validate_group_transition(copy)
  })) return(invisible(process))
  check_columns(r, c("group_origin", "group_dest", "rate"), "rates")
  age <- all(c("age_start", "age_end") %in% names(r))
  dur <- all(c("dur_start", "dur_end") %in% names(r))
  if (!nrow(r) || age == dur ||
      any(c("time_start", "time_end", "union_status", "group", "marital_status") %in% names(r)) ||
      (age && any(c("dur_start", "dur_end") %in% names(r))) ||
      (dur && any(c("age_start", "age_end") %in% names(r))))
    stop("Use one common age or duration band type for all group transitions.", call. = FALSE)
  for (column in c("group_origin", "group_dest")) {
    check_groups(data.frame(group = r[[column]]), column)
  }
  if (any(as.character(r$group_origin) == as.character(r$group_dest)))
    stop("Group transitions cannot have the same origin and destination.", call. = FALSE)
  start <- r[[if (age) "age_start" else "dur_start"]]
  end <- r[[if (age) "age_end" else "dur_end"]]
  check_numeric_column(start, "band start")
  check_numeric_column(r$rate, "rate")
  if (any(start < 0) || any(r$rate < 0) || !is.numeric(end) ||
      !is.null(dim(end)) || anyNA(end) || any(end <= start))
    stop("Transition rates and band starts must be nonnegative with valid band ends.", call. = FALSE)
  if ("sex" %in% names(r) && (!is.character(r$sex) || anyNA(r$sex) || any(!r$sex %in% c("m", "f"))))
    stop("Transition sex must be m or f.", call. = FALSE)
  if ("union_type" %in% names(r) && (!is.character(r$union_type) || anyNA(r$union_type) ||
      any(!r$union_type %in% c("unpartnered", "partnered", "cohab", "marriage")) ||
      ("partnered" %in% r$union_type && any(r$union_type %in% c("cohab", "marriage")))))
    stop("Transition union_type must use unpartnered/partnered or unpartnered/cohab/marriage.", call. = FALSE)
  keys <- intersect(c("group_origin", "group_dest", "sex", "union_type"), names(r))
  strata <- do.call(interaction, c(lapply(r[keys], as.factor), list(drop = TRUE, lex.order = TRUE)))
  for (rows in split(seq_len(nrow(r)), strata)) {
    rows <- rows[order(start[rows])]
    if (length(rows) > 1L && any(start[rows[-1L]] < head(cummax(end[rows]), -1L)))
      stop("Group-transition bands must not overlap within a stratum.", call. = FALSE)
  }
  invisible(process)
}

empty_transitions <- function(groups = character()) {
  tibble::tibble(person_id = numeric(), time = numeric(),
    group_origin = groups[integer()], group_dest = groups[integer()])
}

normalize_population_transitions <- function(population) {
  if ("group" %in% names(population$persons)) {
    if (is.null(population$transitions))
      population$transitions <- empty_transitions(population$persons$group)
    else population$transitions <- tibble::as_tibble(population$transitions)
  }
  population
}

validate_transition_history <- function(population) {
  tx <- population$transitions
  if (is.null(tx)) return(invisible(population))
  p <- population$persons
  check_columns(tx, c("person_id", "time", "group_origin", "group_dest"), "transitions")
  if (nrow(tx) && !"group" %in% names(p))
    stop("Transition histories require person groups.", call. = FALSE)
  check_numeric_column(tx$person_id, "transitions$person_id", id = TRUE)
  check_numeric_column(tx$time, "transitions$time")
  for (col in c("group_origin", "group_dest")) check_groups(data.frame(group = tx[[col]]), col)
  i <- match(tx$person_id, p$person_id)
  if (anyNA(i) || any(tx$time < p$birth_time[i]) ||
      any(tx$time >= p$death_time[i], na.rm = TRUE) ||
      (!is.null(population$end_time) && any(tx$time > population$end_time)) ||
      any(as.character(tx$group_origin) == as.character(tx$group_dest)))
    stop("Transition history has invalid people, dates or self-transitions.", call. = FALSE)
  for (rows in split(seq_len(nrow(tx)), tx$person_id)) {
    rows <- rows[order(tx$time[rows])]
    if (anyDuplicated(tx$time[rows]) ||
        (length(rows) > 1L && any(as.character(tx$group_dest[head(rows, -1L)]) !=
                                as.character(tx$group_origin[tail(rows, -1L)]))) ||
        as.character(tail(tx$group_dest[rows], 1)) != as.character(p$group[i[rows[1]]]))
      stop("Transition history must chain chronologically to each person's current group.", call. = FALSE)
  }
  invisible(population)
}

# Private observation copies carry indexed histories to all native counters.
prepare_group_observations <- function(persons, population, by_group) {
  tx <- if (inherits(population, "hammel_population")) population$transitions else NULL
  if (!isTRUE(by_group) || is.null(tx) || !nrow(tx)) return(persons)
  groups <- extend_identity_groups(unique(persons$group), c(as.character(tx$group_origin), as.character(tx$group_dest)))
  attr(persons, ".hammel_group_values") <- groups
  attr(persons, ".hammel_group_history") <- list(
    row = match(tx$person_id, persons$person_id) - 1L, time = tx$time,
    from = match(as.character(tx$group_origin), as.character(groups)) - 1L,
    to = match(as.character(tx$group_dest), as.character(groups)) - 1L)
  persons
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
    out[j] <- if (any(prior)) tail(tx$group_dest[h[prior]], 1) else tx$group_origin[h[1]]
  }
  out
}

prepare_union_group_observations <- function(unions, population, groups,
                                              pair = FALSE, typed = FALSE) {
  tx <- population$transitions
  if (is.null(tx) || !nrow(tx) || !nrow(unions)) return(unions)
  persons <- prepare_group_observations(population$persons, population, TRUE)
  identities <- observation_groups(persons, TRUE)
  indices <- match(persons$group, identities) - 1L
  # Map unordered pair labels once; native history traversal uses integer keys.
  pair_indices <- if (pair) {
    grid <- expand.grid(w = seq_along(identities), h = seq_along(identities))
    matrix(match(canonical_union_pairs(identities[grid$w], identities[grid$h]), groups) - 1L,
      nrow = length(identities))
  } else matrix(integer(), 0L, 0L)
  history <- union_group_history_cpp(persons, unions, indices, pair_indices,
    pair, typed, length(groups))
  if (length(history$row)) attr(unions, ".hammel_group_history") <- history
  unions
}
