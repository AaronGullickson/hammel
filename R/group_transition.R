# Construct and validate group-transition processes and household move rules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify individual group transitions
#'
#' @description
#'
#' Constructs a group-specific transition process specification for a simulation
#' model. The convenience function `set_group_tansition` is equivalent to
#' `set_process(model, group_transition(...))` and allows for direct piping from
#' a model.
#'
#' @details
#'
#' The primary and only required input to the function is a transition rate
#' schedule for the `rates` argument. This data frame must contain either
#' `age_start` and `age_end` columns or `dur_start` and `dur_end` columns to
#' define the half-open age/duration intervals (inclusive at the start and
#' exclusive at the end), an `origin` and `destination` column to define the
#' origin and destination categories for the transition and a `rate` column for
#' the transition rates. Transition rate are defined as the number of
#' transitions from the origin to destination category per person-year spent
#' in the origin category.
#'
#' For a set of K categories, there will be K(K-1) possible transition rates at
#' every age/duration band. Any cases left unspecified will be assumed to have a
#' zero rate. The rates are not assumed to be symmetric (e.g.A > B is not the
#' same as B < A), but users can provide schedules that make these rates
#'
#' Each group dimension has its own transition process: and users must call
#' `set_group_transition()` separately for each group dimension for which they
#' desire a process. Columns named by declared group dimensions besides the
#' focal group supply separate schedules for each combination of their
#' categories; dimensions without a column share the same schedule.
#'
#' Optional `time_start`/`time_end` columns select half-open calendar intervals
#' in absolute years, including fractional years. Both columns are required
#' together; starts are finite and ends exceed starts (possibly Inf). Active
#' slices must satisfy the strata and coverage rules and cover the requested
#' simulation run.
#'
#' Transition rates can also be stratified by `union_type` containing either
#' `"partnered"` and `"unpartnered"`, or all three of `"unpartnered"`,
#' `"cohab"`, `"marriage"`, or by `sex` ("m" and "f").
#'
#' @param rates A nonempty data frame defining the mortality rate schedules.
#'   Must contain at least `origin`, `destination`, `rate`, and either
#'   `age_start`/`age_end` or `dur_start`/`dur_end` throughout. See details for
#'   format and description of stratification variables that can be applied.
#' @param group The name of the declared group dimension (see [set_groups()])
#'   whose categories are affected by this transition rate schedule. `NULL`
#'   (default) means the model's only dimension, and is an error when several
#'   are declared.
#' @param weights NULL (default) for unweighted selection, or a [weight_model()]
#'   using person variables and `destination`, which chooses who in a cell moves
#'   without changing the cell's input rate.
#' @param household `NULL` (default): each person moves individually. A
#'   [household()]: a move carries the mover's household with them on a
#'   transition. If the household is set, then the rates apply to the household
#'   collectively or individuals if they are not in households. Can only be used
#'   with age bands.
#' @param model A `hammel_model` object, used by `set_group_transition()`.
#' @return `group_transition()` returns a process with a rates tibble;
#'   `set_group_transition()` returns the updated model.
#' @export
#' @examples
#' # create grouped population
#' pop_start <- generate_uniform_population(
#'   1000,
#'   groups = list(ethnicity = c(A = 1, B = 1))
#' )
#'
#' # create transition matrix with more A > B than B > A
#' ethnic_trans <- tibble::tribble(
#'  ~age_start, ~age_end, ~origin, ~destination, ~rate,
#'  0,          Inf,      "A",     "B",          0.05,
#'  0,          Inf,      "B",     "A",          0.02
#' )
#'
#' # create simulation model and run it
#' model <- hammel_model() |>
#'   set_population(pop_start) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_union_formation("on_demand") |>
#'   set_groups(ethnicity = group(ancestry = TRUE)) |>
#'   set_group_transition(ethnic_trans)
#' pop_end <- run_simulation(model, 100)
#'
#' # check population sizes - group A will grow relative to group B
#' calculate_population_size(pop_end, seq(0, 100, 20), by = "ethnicity")
group_transition <- function(
  rates,
  group = NULL,
  weights = NULL,
  household = NULL
) {
  new_group_transition(rates, group, weights, household)
}

#' @rdname group_transition
#' @export
set_group_transition <- function(
  model,
  rates,
  group = NULL,
  weights = NULL,
  household = NULL
) {
  attach_process(model, new_group_transition(rates, group, weights, household))
}

#' Define who moves together in household transitions
#'
#' @description
#'
#' Define a set of rules for determining household membership for transitions.
#'
#' @details
#'
#' A household rule describes a household for [set_group_transition()]. The
#' simulation does not track households independently, but determines who
#' belongs to a household when a household would move and assembles the moving
#' party.
#'
#' A person's household is their current partner (married or cohabiting) and the
#' dependents attached to either of them. A dependent is alive, unpartnered,
#' younger than `age_independence`, and attached to their mother while she is
#' alive, otherwise to their father while he is alive. People with neither
#' parent alive are independent. A person under `age_independence` who leaves a
#' union is a dependent again if a parent is alive.
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
  if (!is.numeric(age_independence) || length(age_independence) != 1L ||
        !is.finite(age_independence) ||
        age_independence < 0) {
    abort_hammel(
      c(
        "{.arg age_independence} must be a single nonnegative number of years.",
        "x" = "It is {describe_value(age_independence)}."
      ),
      "argument",
      argument = "age_independence"
    )
  }
  structure(
    list(age_independence = as.numeric(age_independence)),
    class = "hammel_household"
  )
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a group transition process
#'
#' @description
#'
#' Summary method for class `"hammel_group_transition"`
#'
#' @param object A `hammel_group_transition` object returned by
#'   [group_transition()].
#' @param categories Optional character vector of the dimension's
#'   categories, used to count possible transitions. `NULL` (default) uses the
#'   categories that appear in the `origin` or `destination` columns.
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_group_transition`
#' containing:
#' \describe{
#'   \item{group}{The name of the group dimension whose categories change, or
#'     `NULL` when the process applies to the model's only dimension.}
#'   \item{categories}{The categories used to define possible transitions.}
#'   \item{transitions}{A tibble of the `origin` and `destination` pairs that
#'     have a positive rate in at least one band or stratum.}
#'   \item{missing_transitions}{A tibble of the possible `origin` and
#'     `destination` pairs between distinct categories that have no positive
#'     rate in any band or stratum.}
#'   \item{rate_strata}{A list describing the transition rate schedule apart
#'     from its `origin` and `destination` columns, containing
#'     `stratified_by` (stratum column names), `strata` (a tibble of observed
#'     stratum combinations), `bands` (a tibble of exact band combinations),
#'     and `bands_by_stratum` (whether band combinations differ by stratum).
#'     Shared bands are stored once; otherwise `bands` includes the stratum
#'     columns. Band endpoints define half-open intervals `[start, end)`.}
#'   \item{weights}{A summary of the weight model used for selection within
#'     rate cells: a `summary.hammel_weight_model`, a list of `periods` and
#'     `models` for a calendar schedule, or `NULL` when selection is
#'     unweighted.}
#'   \item{age_independence}{The household age of independence, in years,
#'     or `NULL` when transitions move individuals rather than households.}
#' }
#'
#' @seealso [group_transition()], [household()]
#' @examples
#' rates <- tibble::tribble(
#'   ~age_start, ~age_end, ~origin, ~destination, ~rate,
#'   0,          Inf,      "A",     "B",          0.05,
#'   0,          Inf,      "B",     "A",          0.02
#' )
#' process <- group_transition(rates, household = household())
#' summary(process)
#' @export
summary.hammel_group_transition <- function(object, categories = NULL, ...) {
  rates <- object$rates
  strata <- intersect(
    names(rates), c(builtin_strata, candidate_strata(rates))
  )
  if (is.null(categories)) {
    categories <- unique(c(rates$origin, rates$destination))
  }
  transitions <- unique(rates[rates$rate > 0, c("origin", "destination")])
  possible <- expand.grid(
    destination = categories,
    origin = categories,
    stringsAsFactors = FALSE
  )
  possible <- tibble::as_tibble(possible[c("origin", "destination")])
  possible <- possible[possible$origin != possible$destination, ]
  covered <- paste(possible$origin, possible$destination) %in%
    paste(transitions$origin, transitions$destination)

  structure(
    list(
      group = object$group,
      categories = categories,
      transitions = transitions,
      missing_transitions = possible[!covered, ],
      rate_strata = summarize_rate_strata(
        rates,
        strata = setdiff(strata, c("origin", "destination"))
      ),
      weights = summarize_model_argument(object$weights),
      age_independence = object$household$age_independence
    ),
    class = "summary.hammel_group_transition"
  )
}

#' @rdname summary.hammel_group_transition
#' @param x A `summary.hammel_group_transition` object returned by
#'   `summary.hammel_group_transition()`.
#' @export
print.summary.hammel_group_transition <- function(x, ...) {

  cli::cli_h2("Group transition process")

  list <- open_summary_list()
  print_stat(
    "Group dimension",
    if (is.null(x$group)) "the model's only dimension" else x$group
  )
  print_transition_details(x)
  cli::cli_end(list)
  invisible(x)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Print the bullets of a group transition summary that do not depend on how
# its dimension is identified. Shared by the process and model summaries.
print_transition_details <- function(x) {
  n <- length(x$categories)
  print_stat(
    paste0(
      "Transitions (", nrow(x$transitions), " of ", n * (n - 1L),
      " possible)"
    ),
    format_transition_pairs(x$transitions)
  )
  print_stat(
    "Missing transitions",
    format_transition_pairs(x$missing_transitions)
  )
  print_rate_strata(x$rate_strata)
  household <- if (is.null(x$age_independence)) {
    "no"
  } else {
    paste0("yes, age of independence ", x$age_independence, " years")
  }
  print_stat("Household moves", household)
  print_model_argument("Weight model", x$weights)
  invisible(x)
}

# Format origin and destination pairs as arrows, showing at most `max`
# pairs followed by a count of the rest.
format_transition_pairs <- function(pairs, max = 10L) {
  if (!nrow(pairs)) {
    return("none")
  }
  labels <- paste(pairs$origin, cli::symbol$arrow_right, pairs$destination)
  labels <- gsub(" ", "\u00a0", labels, fixed = TRUE)
  if (length(labels) <= max) {
    return(paste(labels, collapse = ", "))
  }
  paste0(
    paste(labels[seq_len(max)], collapse = ", "),
    ", and ", length(labels) - max, " more"
  )
}

# Build a validated transition process for group_transition() and
# set_group_transition(); errors name `call`.
new_group_transition <- function(
  rates,
  group,
  weights,
  household,
  call = rlang::caller_env()
) {
  # Check that `rates` is a table with unique column names before
  # new_rate_process() converts it; validation checks its contents.
  check_columns(rates, character(), "rates", call = call)
  process <- new_rate_process("group_transition", rates)
  process$group <- check_transition_group(group, call = call)
  if (!is.null(household) && !inherits(household, "hammel_household")) {
    abort_hammel(
      c(
        "{.arg household} must be {.code NULL} or a household rule.",
        "x" = "It is {describe_value(household)}.",
        "i" = "Create one with {.fn household}."
      ),
      "specification",
      argument = "household",
      call = call
    )
  }
  process["household"] <- list(household)
  validate_group_transition(process, call = call)
  validate_process_weights(weights, "transition", call = call)
  warn_cancelled_weight_terms(
    weights,
    weight_strata("transition", rates, group),
    "group transition"
  )
  process$weights <- normalize_model_argument(weights)
  process
}

check_transition_group <- function(group, call = rlang::caller_env()) {
  if (is.null(group)) {
    return(NULL)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) ||
        !nzchar(group)) {
    abort_hammel(
      c(
        "{.arg group} must be {.code NULL} or the name of one group dimension.",
        "x" = "It is {describe_value(group)}."
      ),
      "argument",
      argument = "group",
      call = call
    )
  }
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
store_transition_process <- function(
  model,
  process,
  call = rlang::caller_env()
) {
  declaration <- model_groups(model)
  if (is.null(process$group) && length(declaration) == 1L) {
    process$group <- names(declaration)
  }
  if (is.null(process$group) && length(declaration) > 1L) {
    dimensions <- names(declaration)
    abort_hammel(
      c(
        "{.arg group} must name the dimension whose categories change.",
        "x" = "The model declares several dimensions: {.val {dimensions}}."
      ),
      "groups",
      argument = "group",
      call = call
    )
  }
  processes <- transition_processes(model)
  processes[[if (is.null(process$group)) ".default" else process$group]] <-
    process
  model$processes$group_transition <- processes
  model
}

# Resolve each transition process to a declared dimension, in declaration order.
resolve_transition_processes <- function(
  processes,
  declaration,
  call = rlang::caller_env()
) {
  if (!length(processes)) {
    return(list())
  }
  if (!length(declaration)) {
    abort_hammel(
      c(
        "Group transitions need a declared group dimension.",
        "i" = "Declare dimensions with {.fn set_groups}."
      ),
      "groups",
      call = call
    )
  }
  dims <- names(declaration)
  resolved <- list()
  for (key in names(processes)) {
    x <- processes[[key]]
    group <- x$group
    if (is.null(group)) {
      if (length(dims) != 1L) {
        abort_hammel(
          c(
            paste(
              "A group transition without {.arg group} needs exactly one",
              "declared group dimension."
            ),
            "x" = "The model declares {.val {dims}}.",
            "i" = "Name the dimension with {.arg group}."
          ),
          "groups",
          call = call
        )
      }
      group <- dims
    }
    if (!group %in% dims) {
      abort_hammel(
        c(
          "Group transitions use {.field {group}}, which isn't declared.",
          "i" = paste(
            "Declared dimensions: {.val {dims}}. Declare it with",
            "{.fn set_groups}."
          )
        ),
        "groups",
        group = group,
        call = call
      )
    }
    if (!is.null(resolved[[group]])) {
      abort_hammel(
        "Group {.field {group}} has more than one transition process.",
        "groups",
        group = group,
        call = call
      )
    }
    x$group <- group
    if (group %in% names(x$rates)) {
      abort_hammel(
        c(
          paste(
            "Transition rates for {.field {group}} can't be stratified by",
            "{.field {group}}."
          ),
          "i" = "Use {.field origin} for the category people move from."
        ),
        "groups",
        group = group,
        call = call
      )
    }
    resolved[[group]] <- x
  }
  resolved[intersect(dims, names(resolved))]
}

# Validate in two passes: row checks on the whole table, then checks within
# each calendar period.
validate_group_transition <- function(process, call = rlang::caller_env()) {
  if (!identical(process$name, "group_transition")) {
    abort_hammel(
      "{.arg process} must be a group transition process.",
      "specification",
      call = call
    )
  }
  r <- process$rates
  old <- intersect(c("group_origin", "group_dest"), names(r))
  if (length(old)) {
    abort_hammel(
      c(
        paste(
          "Transition rates name categories with {.field origin} and",
          "{.field destination}."
        ),
        "x" = "Found {.field {old}}.",
        "i" = "Rename these columns."
      ),
      "rates",
      argument = "rates",
      columns = old,
      call = call
    )
  }
  check_stratified_columns(
    r,
    c("origin", "destination", "sex", "union_type"),
    "Group transition",
    call = call
  )
  axis <- check_transition_rows(r, process$household, call = call)
  validate_transition_schedules(r, axis, call = call)
  invisible(process)
}

# Row checks for a transition table. Returns its band axis, "age" or "dur".
check_transition_rows <- function(r, household, call = rlang::caller_env()) {
  check_columns(r, c("origin", "destination", "rate"), "rates", call = call)
  if (!nrow(r)) {
    abort_hammel(
      c(
        "{.arg rates} must contain at least one transition.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  calendar_rate_segments(r, call = call)
  status <- intersect(c("union_status", "marital_status"), names(r))
  if (length(status)) {
    abort_hammel(
      c(
        "Transition rates can't have {.field {status}} columns.",
        "i" = "Stratify by union status with {.field union_type}."
      ),
      "rates",
      argument = "rates",
      columns = status,
      call = call
    )
  }
  band_columns <- intersect(
    c("age_start", "age_end", "dur_start", "dur_end"),
    names(r)
  )
  age <- all(c("age_start", "age_end") %in% band_columns)
  dur <- all(c("dur_start", "dur_end") %in% band_columns)
  if (!is.null(household) && dur) {
    abort_hammel(
      c(
        "Household moves need age bands.",
        "i" = paste(
          "Dependents leave the household at an age, so give",
          "{.field age_start} and {.field age_end} instead of duration bands."
        )
      ),
      "specification",
      argument = "household",
      call = call
    )
  }
  if (age == dur || length(band_columns) != 2L) {
    abort_hammel(
      c(
        "{.arg rates} must use either age bands or duration bands.",
        "x" = if (length(band_columns)) {
          "Band columns: {.field {band_columns}}."
        } else {
          "It has no band columns."
        },
        "i" = paste(
          "Give {.field age_start} and {.field age_end}, or {.field dur_start}",
          "and {.field dur_end}, for every transition."
        )
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  for (column in c("origin", "destination")) {
    check_group_labels(r[[column]], paste0("rates$", column), call = call)
  }
  rows <- which(as.character(r$origin) == as.character(r$destination))
  if (length(rows)) {
    abort_hammel(
      c(
        "Transitions must move between different categories.",
        "x" = paste(
          "{.field origin} equals {.field destination} at",
          "{format_rows(rows)}."
        )
      ),
      "rates",
      argument = "rates",
      rows = rows,
      call = call
    )
  }
  axis <- if (age) "age" else "dur"
  check_bands(r, axis, call = call)
  check_numeric_column(r$rate, "rates$rate", call = call)
  check_nonnegative(r$rate, "rates$rate", call = call)
  if ("sex" %in% names(r)) {
    check_labels(r$sex, "rates$sex", c("m", "f"), call = call)
  }
  if ("union_type" %in% names(r)) {
    check_labels(
      r$union_type, "rates$union_type",
      c("unpartnered", "partnered", "cohab", "marriage"),
      call = call
    )
  }
  axis
}

# Within each calendar period: one union-type scheme and bands that do not
# overlap within a stratum.
validate_transition_schedules <- function(r, axis, call = rlang::caller_env()) {
  if (validate_calendar_rates(
        r,
        function(x) {
          validate_transition_schedules(x, axis, call = call)
        },
        call = call
      )) {
    return(invisible(r))
  }
  if ("union_type" %in% names(r) && "partnered" %in% r$union_type &&
        any(r$union_type %in% c("cohab", "marriage"))) {
    types <- unique(r$union_type)
    abort_hammel(
      c(
        paste(
          "Column {.field union_type} of {.arg rates} must use",
          "{.val unpartnered} and {.val partnered}, or {.val unpartnered},",
          "{.val cohab} and {.val marriage}."
        ),
        "x" = "It mixes them: {.val {types}}."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  keys <- c(
    intersect(c("origin", "destination", "sex", "union_type"), names(r)),
    candidate_strata(r)
  )
  strata <- do.call(
    interaction, c(
      lapply(r[keys], as.factor), list(drop = TRUE, lex.order = TRUE)
    )
  )
  for (rows in split(seq_len(nrow(r)), strata)) {
    values <- vapply(
      keys,
      function(column) {
        as.character(r[[column]][rows[1]])
      },
      character(1)
    )
    schedule <- paste0(" for ", paste(keys, "=", values, collapse = ", "))
    check_band_overlap(r, rows, axis, schedule = schedule, call = call)
  }
  invisible(r)
}
