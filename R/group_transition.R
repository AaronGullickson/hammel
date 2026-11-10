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
  process <- new_rate_process("group_transition", rates)
  process$group <- check_transition_group(group)
  if (!is.null(household) && !inherits(household, "hammel_household")) {
    stop("`household` must be NULL or created with household().", call. = FALSE)
  }
  process["household"] <- list(household)
  validate_group_transition(process)
  validate_process_weights(weights, "transition")
  warn_cancelled_weight_terms(
    weights,
    weight_strata("transition", rates, group),
    "group transition"
  )
  process$weights <- normalize_model_argument(weights)
  process
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
  set_process(
    model,
    group_transition(
      rates,
      group = group,
      weights = weights,
      household = household
    )
  )
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
    stop(
      "`age_independence` must be a finite nonnegative number of years.",
      call. = FALSE
    )
  }
  structure(
    list(age_independence = as.numeric(age_independence)),
    class = "hammel_household"
  )
}

# HELPER FUNCTIONS ---------------------------------------------------------

check_transition_group <- function(group) {
  if (is.null(group)) {
    return(NULL)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) ||
        !nzchar(group)) {
    stop(
      "`group` must be NULL or the name of one declared group dimension.",
      call. = FALSE
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
store_transition_process <- function(model, process) {
  declaration <- model_groups(model)
  if (is.null(process$group) && length(declaration) == 1L) {
    process$group <- names(declaration)
  }
  if (is.null(process$group) && length(declaration) > 1L) {
    stop(
      sprintf(
        "Name the dimension whose categories change with `group`: one of %s.",
        paste0("\"", names(declaration), "\"", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  processes <- transition_processes(model)
  processes[[if (is.null(process$group)) ".default" else process$group]] <-
    process
  model$processes$group_transition <- processes
  model
}

# Resolve each transition process to a declared dimension, in declaration order.
resolve_transition_processes <- function(processes, declaration) {
  if (!length(processes)) {
    return(list())
  }
  if (!length(declaration)) {
    stop(
      paste0(
        "Group transitions require a declared group ",
        "dimension; use set_groups()."
      ),
      call. = FALSE
    )
  }
  dims <- names(declaration)
  resolved <- list()
  for (key in names(processes)) {
    x <- processes[[key]]
    group <- x$group
    if (is.null(group)) {
      if (length(dims) != 1L) {
        stop(
          paste0(
            "A group transition without `group` needs exactly ",
            "one declared group dimension."
          ),
          call. = FALSE
        )
      }
      group <- dims
    }
    if (!group %in% dims) {
      stop(
        sprintf(
          paste0(
            "Group transitions use `%s`, which is not a declared ",
            "group dimension."
          ),
          group
        ),
        call. = FALSE
      )
    }
    if (!is.null(resolved[[group]])) {
      stop(
        sprintf("Group `%s` has more than one transition process.", group),
        call. = FALSE
      )
    }
    x$group <- group
    if (group %in% names(x$rates)) {
      stop(
        sprintf(
          paste0(
            "Transition rates for `%s` cannot be stratified by ",
            "`%s` itself; use `origin`."
          ),
          group, group
        ),
        call. = FALSE
      )
    }
    resolved[[group]] <- x
  }
  resolved[intersect(dims, names(resolved))]
}

validate_group_transition <- function(process) {
  if (!identical(process$name, "group_transition")) {
    stop("Invalid group-transition specification.", call. = FALSE)
  }
  r <- process$rates
  old <- intersect(c("group_origin", "group_dest"), names(r))
  if (length(old)) {
    stop(
      paste0(
        "Transition rates name categories with `origin` and ",
        "`destination`, not "
      ),
      paste0("`", old, "`", collapse = " and "), ".",
      call. = FALSE
    )
  }
  check_stratified_columns(
    r,
    c("origin", "destination", "sex", "union_type"),
    "Group-transition"
  )
  if (validate_calendar_rates(
        r,
        function(x) {
       copy <- process
       copy$rates <- x
       validate_group_transition(copy)
     }
      )) {
    return(invisible(process))
  }
  check_columns(r, c("origin", "destination", "rate"), "rates")
  age <- all(c("age_start", "age_end") %in% names(r))
  dur <- all(c("dur_start", "dur_end") %in% names(r))
  if (!is.null(process$household) && dur) {
    stop(
      paste0(
        "Household moves need age bands: dependents reach ",
        "their own clock at a set age."
      ),
      call. = FALSE
    )
  }
  if (!nrow(r) || age == dur ||
        any(
          c("time_start", "time_end", "union_status", "marital_status") %in%
            names(r)
        ) ||
        (age && any(c("dur_start", "dur_end") %in% names(r))) ||
        (dur && any(c("age_start", "age_end") %in% names(r)))) {
    stop(
      "Use one common age or duration band type for all group transitions.",
      call. = FALSE
    )
  }
  for (column in c("origin", "destination")) {
    check_groups(data.frame(group = r[[column]]), column)
  }
  if (any(as.character(r$origin) == as.character(r$destination))) {
    stop(
      "Group transitions cannot have the same origin and destination.",
      call. = FALSE
    )
  }
  start <- r[[if (age) "age_start" else "dur_start"]]
  end <- r[[if (age) "age_end" else "dur_end"]]
  check_numeric_column(start, "band start")
  check_numeric_column(r$rate, "rate")
  if (any(start < 0) || any(r$rate < 0) || !is.numeric(end) ||
        !is.null(dim(end)) ||
        anyNA(end) ||
        any(end <= start)) {
    stop(
      paste0(
        "Transition rates and band starts must be ",
        "nonnegative with valid band ends."
      ),
      call. = FALSE
    )
  }
  if ("sex" %in% names(r) &&
        (!is.character(r$sex) || anyNA(r$sex) ||
           any(!r$sex %in% c("m", "f")))) {
    stop("Transition sex must be m or f.", call. = FALSE)
  }
  if ("union_type" %in% names(r) &&
        (!is.character(r$union_type) || anyNA(r$union_type) ||
           any(
             !r$union_type %in%
                c("unpartnered", "partnered", "cohab","marriage")
           ) ||
           ("partnered" %in% r$union_type &&
              any(r$union_type %in% c("cohab", "marriage"))))) {
    stop(
      paste0(
        "Transition union_type must use ",
        "unpartnered/partnered or unpartnered/cohab/marriage."
      ),
      call. = FALSE
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
    rows <- rows[order(start[rows])]
    if (length(rows) > 1L &&
          any(start[rows[-1L]] < head(cummax(end[rows]), -1L))) {
      stop(
        "Group-transition bands must not overlap within a stratum.",
        call. = FALSE
      )
    }
  }
  invisible(process)
}
