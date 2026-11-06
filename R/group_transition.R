#' Specify individual group transitions
#'
#' Construct or attach competing origin-to-destination hazards for one group
#' dimension. Only the individual changes category; ancestry, parentage, union
#' histories and other dimensions are retained. Duration is years since the
#' person's latest recorded transition in this dimension, or age without one.
#' Missing directed pairs, strata and gaps have zero risk. Restarts retain
#' transition histories and resume durations from them.
#'
#' Each dimension has its own transition process: calling
#' `set_group_transition()` again for the same dimension replaces it, and other
#' dimensions' processes are kept. Transitions at the same instant in
#' different dimensions follow the declaration order of [set_groups()].
#'
#' @param rates A nonempty data frame requiring `origin`, `destination`,
#'   `rate`, and either `age_start`/`age_end` or `dur_start`/`dur_end` throughout.
#'   `origin` and `destination` are categories of the transition's dimension.
#'   Half-open bands use years; rates are finite nonnegative hazards per
#'   person-year. Ends may be `Inf`. Bands cannot overlap within a stratum.
#'   Self-transitions are invalid. Optional `sex` is m/f. Optional `union_type`
#'   is unpartnered/partnered or unpartnered/cohab/marriage. Columns named by
#'   other declared group dimensions give separate schedules by those
#'   categories; the transition's own dimension cannot be such a column.
#'   Omitted columns apply universally. Destinations may initially be empty;
#'   other configured group-specific processes must cover all possible
#'   categories. Optional `time_start`/`time_end` specify half-open absolute
#'   calendar-year bands, including fractional years. Every active slice must
#'   be valid and the table must cover the requested run. Rate changes do not
#'   reset duration.
#' @param group The name of the declared group dimension (see [set_groups()])
#'   whose categories change. `NULL` (default) means the model's only
#'   dimension, and is an error when several are declared.
#' @param weights NULL (default) for unweighted selection, or a
#'   [create_weight_model()] using person variables and `destination`, which
#'   chooses who in a cell moves without changing the cell's input rate.
#' @param model A `hammel_model`.
#' @return `group_transition()` returns a process with a rates tibble;
#'   `set_group_transition()` returns the updated model. Grouped simulation
#'   output always includes a `transitions` tibble with `person_id`, `time`,
#'   `group` (the dimension), `origin`, and `destination`, including when empty.
#' @export
group_transition <- function(rates, group = NULL, weights = NULL) {
  process <- new_rate_process("group_transition", rates)
  process$group <- check_transition_group(group)
  validate_group_transition(process)
  validate_process_weights(weights, "transition")
  warn_cancelled_weight_terms(weights, weight_strata("transition", rates, group), "group transition")
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname group_transition
#' @export
set_group_transition <- function(model, rates, group = NULL, weights = NULL) {
  set_process(model, group_transition(rates, group = group, weights = weights))
}

check_transition_group <- function(group) {
  if (is.null(group)) return(NULL)
  if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group))
    stop("`group` must be NULL or the name of one declared group dimension.", call. = FALSE)
  group
}

# A model's transition processes, one per dimension, keyed by dimension name
# (".default" for a process whose dimension is the model's only one).
transition_processes <- function(model) {
  x <- model$processes$group_transition
  if (is.null(x)) list() else x
}

# Store a transition process with the model's others, replacing any for the
# same dimension. A process without `group` takes the only declared dimension.
store_transition_process <- function(model, process) {
  declaration <- model_groups(model)
  if (is.null(process$group) && length(declaration) == 1L) process$group <- names(declaration)
  if (is.null(process$group) && length(declaration) > 1L)
    stop(sprintf("Name the dimension whose categories change with `group`: one of %s.",
                 paste0("\"", names(declaration), "\"", collapse = ", ")), call. = FALSE)
  processes <- transition_processes(model)
  processes[[if (is.null(process$group)) ".default" else process$group]] <- process
  model$processes$group_transition <- processes
  model
}

# Resolve each transition process to a declared dimension, in declaration order.
resolve_transition_processes <- function(processes, declaration) {
  if (!length(processes)) return(list())
  if (!length(declaration))
    stop("Group transitions require a declared group dimension; use set_groups().", call. = FALSE)
  dims <- names(declaration)
  resolved <- list()
  for (key in names(processes)) {
    x <- processes[[key]]
    group <- x$group
    if (is.null(group)) {
      if (length(dims) != 1L)
        stop("A group transition without `group` needs exactly one declared group dimension.", call. = FALSE)
      group <- dims
    }
    if (!group %in% dims)
      stop(sprintf("Group transitions use `%s`, which is not a declared group dimension.", group), call. = FALSE)
    if (!is.null(resolved[[group]]))
      stop(sprintf("Group `%s` has more than one transition process.", group), call. = FALSE)
    x$group <- group
    if (group %in% names(x$rates))
      stop(sprintf("Transition rates for `%s` cannot be stratified by `%s` itself; use `origin`.",
                   group, group), call. = FALSE)
    resolved[[group]] <- x
  }
  resolved[intersect(dims, names(resolved))]
}

validate_group_transition <- function(process) {
  if (!identical(process$name, "group_transition") ||
      !identical(process$rate_unit, "person_year"))
    stop("Invalid group-transition specification or units.", call. = FALSE)
  r <- process$rates
  old <- intersect(c("group_origin", "group_dest"), names(r))
  if (length(old))
    stop("Transition rates name categories with `origin` and `destination`, not ",
         paste0("`", old, "`", collapse = " and "), ".", call. = FALSE)
  check_stratified_columns(r, c("origin", "destination", "sex", "union_type"), "Group-transition")
  if (validate_calendar_rates(r, function(x) {
    copy <- process; copy$rates <- x; validate_group_transition(copy)
  })) return(invisible(process))
  check_columns(r, c("origin", "destination", "rate"), "rates")
  age <- all(c("age_start", "age_end") %in% names(r))
  dur <- all(c("dur_start", "dur_end") %in% names(r))
  if (!nrow(r) || age == dur ||
      any(c("time_start", "time_end", "union_status", "marital_status") %in% names(r)) ||
      (age && any(c("dur_start", "dur_end") %in% names(r))) ||
      (dur && any(c("age_start", "age_end") %in% names(r))))
    stop("Use one common age or duration band type for all group transitions.", call. = FALSE)
  for (column in c("origin", "destination")) {
    check_groups(data.frame(group = r[[column]]), column)
  }
  if (any(as.character(r$origin) == as.character(r$destination)))
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
  keys <- c(intersect(c("origin", "destination", "sex", "union_type"), names(r)), candidate_strata(r))
  strata <- do.call(interaction, c(lapply(r[keys], as.factor), list(drop = TRUE, lex.order = TRUE)))
  for (rows in split(seq_len(nrow(r)), strata)) {
    rows <- rows[order(start[rows])]
    if (length(rows) > 1L && any(start[rows[-1L]] < head(cummax(end[rows]), -1L)))
      stop("Group-transition bands must not overlap within a stratum.", call. = FALSE)
  }
  invisible(process)
}

empty_transitions <- function() {
  tibble::tibble(person_id = numeric(), time = numeric(), group = character(),
                 origin = character(), destination = character())
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

# Private observation copies carry indexed histories to all native counters.
prepare_group_observations <- function(persons, population, by_group) {
  if (!isTRUE(by_group)) return(persons)
  # A combined view of several dimensions carries indexed codes directly.
  combined <- if (inherits(population, "hammel_population")) population[[".hammel_combined"]]
  if (!is.null(combined)) {
    attr(persons, ".hammel_group_values") <- combined$values
    h <- combined$history
    if (nrow(h)) attr(persons, ".hammel_group_history") <- list(
      row = match(h$person_id, persons$person_id) - 1L, time = h$time, from = h$from, to = h$to)
    return(persons)
  }
  tx <- if (inherits(population, "hammel_population")) population$transitions else NULL
  if (is.null(tx) || !nrow(tx)) return(persons)
  groups <- extend_identity_groups(unique(persons$group), c(as.character(tx$origin), as.character(tx$destination)))
  attr(persons, ".hammel_group_values") <- groups
  attr(persons, ".hammel_group_history") <- list(
    row = match(tx$person_id, persons$person_id) - 1L, time = tx$time,
    from = match(as.character(tx$origin), as.character(groups)) - 1L,
    to = match(as.character(tx$destination), as.character(groups)) - 1L)
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
    out[j] <- if (any(prior)) tail(tx$destination[h[prior]], 1) else tx$origin[h[1]]
  }
  out
}

# Attach union-level group histories for `mode` 0 (wife), 1 (husband) or 2
# (both: wife + identities * husband), offset by type when `typed`.
prepare_union_group_observations <- function(unions, population, mode, typed,
                                              identity_count, group_count) {
  tx <- population$transitions
  combined <- population[[".hammel_combined"]]
  moves <- if (is.null(combined)) nrow(tx) else nrow(combined$history)
  if (!length(moves) || !moves || !nrow(unions)) return(unions)
  persons <- prepare_group_observations(population$persons, population, TRUE)
  identities <- observation_groups(persons, TRUE)
  indices <- match(persons$group, identities) - 1L
  history <- union_group_history_cpp(persons, unions, indices, as.integer(mode),
    typed, as.integer(identity_count), as.integer(group_count))
  if (length(history$row)) attr(unions, ".hammel_group_history") <- history
  unions
}
