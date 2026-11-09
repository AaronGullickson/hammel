# Construct and validate group-transition processes and household move rules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

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
#' @param household `NULL` (default): each person moves individually. A
#'   [household()]: a move carries the mover's household (D140). One clock
#'   per household: a woman in a union carries it, while her partner and
#'   dependents have none; everyone else has their own. When a carrier moves,
#'   their partner and the dependents attached to either of them move to the
#'   same destination, from wherever they are (members already there stay).
#'   Requires age bands.
#' @param model A `hammel_model`.
#' @return `group_transition()` returns a process with a rates tibble;
#'   `set_group_transition()` returns the updated model. Grouped simulation
#'   output always includes a `transitions` tibble with `person_id`, `time`,
#'   `group` (the dimension), `origin`, `destination`, `move_id` (shared by the
#'   rows of one move) and `initiator_id` (the person whose clock fired),
#'   including when empty.
#' @export
group_transition <- function(rates, group = NULL, weights = NULL, household = NULL) {
  process <- new_rate_process("group_transition", rates)
  process$group <- check_transition_group(group)
  if (!is.null(household) && !inherits(household, "hammel_household"))
    stop("`household` must be NULL or created with household().", call. = FALSE)
  process["household"] <- list(household)
  validate_group_transition(process)
  validate_process_weights(weights, "transition")
  warn_cancelled_weight_terms(weights, weight_strata("transition", rates, group), "group transition")
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname group_transition
#' @export
set_group_transition <- function(model, rates, group = NULL, weights = NULL, household = NULL) {
  set_process(model, group_transition(rates, group = group, weights = weights, household = household))
}

#' Define who moves together in household moves
#'
#' A household rule describes a household for [set_group_transition()] (D140).
#' There is no household identifier or residence record: when someone moves,
#' the moving party is assembled from current union and parent links. A
#' person's household is their current partner (married or cohabiting) and the
#' dependents attached to either of them. A dependent is alive, unpartnered,
#' younger than `age_independence`, and attached to their mother while she is
#' alive, otherwise to their father while he is alive. People with neither
#' parent alive are independent; grandparents, stepparents and other relatives
#' do not take dependents. A person under `age_independence` who leaves a union is
#' a dependent again if a parent is alive.
#'
#' @param age_independence Age in years at which a child stops being a
#'   dependent: below it, an unpartnered child with a living parent is a
#'   dependent. Defaults to 18.
#' @return A `hammel_household`, to pass as `household` to
#'   [set_group_transition()]. Each transition process may use its own rule.
#' @export
#' @examples
#' household(age_independence = 18)
household <- function(age_independence = 18) {
  if (!is.numeric(age_independence) || length(age_independence) != 1L || !is.finite(age_independence) ||
      age_independence < 0)
    stop("`age_independence` must be a finite nonnegative number of years.", call. = FALSE)
  structure(list(age_independence = as.numeric(age_independence)), class = "hammel_household")
}

# HELPER FUNCTIONS ---------------------------------------------------------

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
  if (!is.null(process$household) && dur)
    stop("Household moves need age bands: dependents reach their own clock at a set age.", call. = FALSE)
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
