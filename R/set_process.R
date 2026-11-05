#' Set a demographic process specification
#'
#' Store a process specification under its name, replacing an existing process
#' with the same name. Other processes are retained. Built-in rate schedules
#' are validated again when attached. No simulation is executed.
#'
#' `set_mortality()` and `set_fertility()` construct the corresponding process
#' and delegate to `set_process()`.
#'
#' @param model A `hammel_model` object.
#' @param process A `hammel_process` specification.
#' @param rates An age-specific rate table passed to [mortality()] or [fertility()].
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   rate multipliers by group. Missing groups use 1.0; values are not normalized.
#'   Scales shared or group-specific fertility input rates before cell
#'   accounting. Without population groups, simulation warns and ignores it.
#'   Fertility also accepts a named list of vectors for `partnered` and/or
#'   `unpartnered`, or cohab/marriage when rates include those `union_type`
#'   strata. List names must occur in the schedules; omitted entries use 1.0.
#' @param hetfert NULL (default) for uniform mother selection, or a
#'   [het_fert()] specification for persistent lognormal selection weights.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [create_weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate. In
#'   fertility they combine with `hetfert` on the log scale.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
#' @param inheritance A group identity model from [create_group_inheritance()],
#'   or a tibble with numeric `time_start`, `time_end` and list-column `model`
#'   containing complete inheritance models. Intervals are half-open absolute
#'   years, must not overlap and must cover the requested run. Select the active
#'   model at each birth; older identities and ancestry remain unchanged.
#'   Base data frames with list columns are accepted. All declared possible
#'   identities across scheduled models join the simulation-wide registry.
#'   Static inputs may also be a method name requiring no additional parameters: `"mother"`, `"father"`,
#'   `"random_parent"`, `"same_sex"`, `"opposite_sex"`, or `"ancestry"`.
#'   Applies to fertility and defaults to selecting either parent with equal
#'   probability, with maternal fallback. Assignment is independent of recorded
#'   ancestry proportions.
#'
#' @return An updated `hammel_model`. The input model is not modified.
#' @export
#' @examples
#' mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008))
#' fertility_rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
#' model <- hammel_model() |>
#'   set_mortality(mortality_rates) |>
#'   set_fertility(fertility_rates)
#' model$processes
set_process <- function(model, process) {
  model <- attach_process(model, process)
  model
}

# Validate and store a process without user-facing advisory warnings; used by
# set_process() and by simulate()'s re-validation of edited specifications.
attach_process <- function(model, process) {
  check_model(model)
  if (!inherits(process, "hammel_process") || !is.list(process)) {
    stop("`process` must be a hammel_process object.", call. = FALSE)
  }
  name <- process$name
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name))) {
    stop("`name` must be a single nonempty string.", call. = FALSE)
  }
  if (name %in% c("mortality", "fertility")) {
    if (!identical(process$rate_unit, "person_year")) {
      stop("Built-in process rates must use person-year units.", call. = FALSE)
    }
    validate_rates(process$rates, mortality = name == "mortality")
    validate_process_weights(process$weights, "person")
    if (name == "fertility") {
      validate_fertility_group_prop_hazard(process$group_prop_hazard, process$rates)
      process$inheritance <- normalize_group_inheritance(process$inheritance)
    } else {
      validate_group_prop_hazard(process$group_prop_hazard)
    }
    if (name == "fertility") {
      validate_birth_interval(process$birth_interval)
      if (!"ffab" %in% names(process)) process$ffab <- 0.4886
      validate_ffab(process$ffab)
      validate_het_fert(process$hetfert)
      # Specifications saved while midpoint compensation existed (before D101
      # became the only method) may carry a spacing_method field.
      if (identical(process$spacing_method, "midpoint"))
        warning(paste('`spacing_method = "midpoint"` is no longer available;',
                      "cell-accounted fertility is used."), call. = FALSE)
      process$spacing_method <- NULL
    }
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_formation") {
    validate_union_formation(process)
    validate_process_weights(process$weights_conversion, "union")
    if (!is.null(process$rates)) process$rates <- tibble::as_tibble(process$rates)
    for (field in c("matching", "father_matching", "cohab_matching"))
      process[field] <- list(normalize_model_argument(process[[field]]))
  }
  if (name == "group_transition") {
    validate_group_transition(process)
    validate_process_weights(process$weights, "transition")
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_dissolution") {
    validate_union_dissolution(process)
    validate_process_weights(process$weights, "union")
    process$rates <- tibble::as_tibble(process$rates)
  }
  model$processes[name] <- list(process)
  model
}

#' @rdname set_process
#' @export
set_mortality <- function(model, rates, group_prop_hazard = NULL, weights = NULL) {
  set_process(model, mortality(rates, group_prop_hazard = group_prop_hazard, weights = weights))
}

#' @rdname set_process
#' @export
set_fertility <- function(model, rates, birth_interval = 0.75, group_prop_hazard = NULL,
                          inheritance = create_group_inheritance(), ffab = 0.4886,
                          hetfert = NULL, weights = NULL) {
  set_process(model, fertility(rates, birth_interval = birth_interval,
                              group_prop_hazard = group_prop_hazard,
                              inheritance = inheritance, ffab = ffab, hetfert = hetfert,
                              weights = weights))
}

#' Specify age-specific demographic rates
#'
#' Construct a process specification without running a simulation. Rates are
#' piecewise constant over age intervals `[age_start, age_end)`, in years,
#' and expressed per person-year. Overlapping intervals are rejected; adjacent
#' intervals are allowed. Input order and custom columns are preserved; rate
#' schedules are returned as tibbles. A `source` attribute on the input table
#' is copied to the process's `source` field.
#'
#' An optional `group` column supplies separate schedules by nonmissing character
#' label or finite numeric code. Without it, rates apply to all groups. Each
#' mortality group needs both sexes, and each conditional fertility group needs
#' both statuses. Overlap is checked within groups.
#' Simulation requires schedules for every population group. Without population
#' groups, simulation warns and uses the first group alphabetically (radix order
#' of character representations, including numeric codes).
#'
#' Mortality requires separate schedules for both sexes. Gaps and uncovered
#' ages are not filled or interpreted as zero rates. Coverage must be specified
#' before simulation. Optional `time_start`/`time_end` columns select half-open
#' calendar intervals in absolute years, including fractional years, on mortality
#' and fertility tables. Both columns are required together; starts are finite
#' and ends exceed starts (possibly Inf). Active slices must satisfy the existing
#' strata/coverage rules and cover the requested run. Rates are constant within
#' intervals, with no interpolation. Calendar changes precede same-time events,
#' preserving histories, birth spacing, parity and duration origins. Inheritance
#' also accepts whole-model calendar schedules; other fertility arguments remain
#' static.
#'
#' Fertility optionally accepts `union_type` containing both `"partnered"`
#' and `"unpartnered"`, or all three `"unpartnered"`, `"cohab"`, `"marriage"`.
#' Shared partnered rates cover either partnership type. Rates use woman-years
#' in the indicated type; the old fertility input `union_status` is rejected.
#' Without that column, one overall schedule applies to all women. A partial
#' status table or any other label is rejected; use explicit zero rates to
#' suppress births in a status. Overlap is checked separately within each status.
#' On-demand union formation requires overall rates; rate-based formation accepts
#' either form and warns when an overall schedule is used.
#'
#' Optional `parity_start`/`parity_end` use half-open current maternal parity
#' bands, covering zero through Inf without gaps or overlap in every group/type
#' stratum. Parity is counted from recorded maternal births (including deceased
#' children) and incremented only on successful births; no custom parity column
#' overrides history. A pooled parity band applies its rate to each contained
#' actual parity.
#'
#' Birth spacing uses cell-accounted fertility (D101). Women are grouped into
#' cells: one-year (or shorter) age sub-bands of each group/union-type schedule
#' after group multipliers, split at whole ages and band edges, and one actual
#' parity when rates are parity-specific. Each cell produces births at its input
#' rate times all its living members, including those inside the birth
#' interval, and each birth goes to an eligible member, uniformly by default
#' or proportional to persistent weights when `hetfert` is enabled. Observed
#' rates then equal inputs in expectation for any population composition,
#' including after restarts and union-type changes, while a cell has an eligible
#' member. Intensity in cells with members but none eligible is dropped and
#' reported in `simulation$unrealized_births`. An individual woman's hazard
#' depends on recent births by others in her cell; this is a calibration
#' device. Upper-open age bands stay one cell when `birth_interval * rate` is
#' below 0.001 and are otherwise split yearly. Zero-rate bands form no cells.
#'
#' For parity-specific schedules with a positive `birth_interval`, `fertility()`
#' warns when positive rates fall in one-year age cells where no woman can be
#' eligible: births are singletons at least `birth_interval` apart, starting
#' from the schedule's first positive rate, so parity k is first eligible at the
#' age found by following each parity's first positive rate. Union types are
#' pooled and groups checked separately. Partly reachable cells and higher
#' actual parities within an open parity band (always unreachable at the
#' youngest ages) are not reported.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`. Starts must be finite and nonnegative; ends must exceed starts
#'   and may be `Inf` for upper-open bands. Rates must
#'   be finite and nonnegative. Mortality additionally
#'   requires character `sex` values `"m"` and `"f"`, with both represented.
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative. Applies to fertility
#'   only; see cell-accounted fertility below.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   rate multipliers, with unique nonempty group names. Missing groups use 1.0;
#'   zero suppresses the risk. Values are not normalized. At execution, multiply
#'   the selected input schedule by the person's group factor, including all
#'   sex/status strata. Fertility multiplication precedes cell accounting.
#'   Without population groups, warn and ignore multipliers. Stored rate tables
#'   and their source metadata are unchanged; rate-table summaries use those
#'   stored rates, without applying population-dependent multipliers.
#'   Fertility with `union_type` rates also accepts a named list such as
#'   `list(partnered = c(A = 1.1, B = 0.9), unpartnered = c(A = 1.05))`.
#'   An omitted status or group uses 1.0. A single vector applies to all
#'   types (or overall fertility). Type lists require `union_type` rates and
#'   their names must match the supplied schedule types. Three-type schedules
#'   accept lists named unpartnered/cohab/marriage; shared partnered schedules
#'   accept unpartnered/partnered. Omitted types and groups have multiplier 1.
#' @param hetfert NULL (default) for uniform mother selection, or a
#'   [het_fert()] specification for persistent lognormal selection weights.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [create_weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate. In
#'   fertility they combine with `hetfert` on the log scale.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
#' @param inheritance A group identity model from [create_group_inheritance()],
#'   or a tibble with numeric `time_start`, `time_end` and list-column `model`
#'   containing complete inheritance models. Intervals are half-open absolute
#'   years, must not overlap and must cover the requested run. Select the active
#'   model at each birth; older identities and ancestry remain unchanged.
#'   Base data frames with list columns are accepted. All declared possible
#'   identities across scheduled models join the simulation-wide registry.
#'   Static inputs may also be a method name requiring no additional parameters: `"mother"`, `"father"`,
#'   `"random_parent"`, `"same_sex"`, `"opposite_sex"`, or `"ancestry"`.
#'   Fertility defaults to selecting either parent with equal probability,
#'   with maternal fallback. Other parental, ancestry,
#'   hierarchy, and custom rules assign identity after child sex and ancestry
#'   are calculated. Possible identities are checked against grouped schedules
#'   and categorical matching levels before execution. Applies to fertility only.
#' @return A named list containing `name`, `rates`, and `rate_unit`, plus
#'   `birth_interval` and `ffab` for fertility, inheriting
#'   from `hammel_mortality` or `hammel_fertility` and `hammel_process`.
#' @export
#' @examples
#' deaths <- mortality(data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008)))
#' births <- fertility(data.frame(age_start = 15, age_end = 50, rate = 0.1))
#' model <- hammel_model() |>
#'   set_process(deaths) |>
#'   set_process(births)
mortality <- function(rates, group_prop_hazard = NULL, weights = NULL) {
  validate_rates(rates, mortality = TRUE)
  validate_group_prop_hazard(group_prop_hazard)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(weights, weight_strata("mortality", rates, group_prop_hazard), "mortality")
  process <- new_rate_process("mortality", rates)
  process$group_prop_hazard <- group_prop_hazard
  process$weights <- weights
  process
}

#' @rdname mortality
#' @export
fertility <- function(rates, birth_interval = 0.75, group_prop_hazard = NULL,
                      inheritance = create_group_inheritance(), ffab = 0.4886,
                          hetfert = NULL, weights = NULL) {
  validate_rates(rates, mortality = FALSE)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(weights, weight_strata("fertility", rates, group_prop_hazard), "fertility")
  validate_fertility_group_prop_hazard(group_prop_hazard, rates)
  inheritance <- normalize_group_inheritance(inheritance)
  validate_birth_interval(birth_interval)
  validate_ffab(ffab)
  validate_het_fert(hetfert)
  process <- new_rate_process("fertility", rates)
  process$inheritance <- if (is.null(inheritance)) create_group_inheritance() else inheritance
  process$birth_interval <- birth_interval
  process$ffab <- ffab
  process$hetfert <- hetfert
  process$weights <- weights
  process$group_prop_hazard <- group_prop_hazard
  warn_unreachable_parity_cells(process$rates, birth_interval)
  process
}

# Positive parity-specific rates in one-year age cells (the fertility cells of
# D101) that no woman can occupy while eligible, given the
# birth interval: singleton births, each at least `interval` apart, starting at
# the schedule's first positive rate. Union types are pooled, because women move
# between them; groups are checked separately. Pooled parity bands are checked
# at their first parity only: higher actual parities in a pooled band reuse one
# schedule and always have unreachable young cells. Cells reachable for part of
# the year are not reported. Returns one row per unreachable age range.
unreachable_parity_cells <- function(rates, interval) {
  empty <- data.frame(group = character(), parity = numeric(), age_start = numeric(),
                      age_end = numeric(), earliest_eligible = numeric())
  if (!all(c("parity_start", "parity_end") %in% names(rates)) || !(interval > 0)) return(empty)
  out <- list()
  # Earliest next birth from eligibility age `e` given positive rows `r`.
  next_birth <- function(r, e) {
    ok <- r$age_end > e
    if (any(ok)) min(pmax(r$age_start[ok], e)) else Inf
  }
  for (rows in split_rate_groups(rates)) {
    r <- rates[rows, , drop = FALSE]
    r <- r[r$rate > 0, , drop = FALSE]
    if (!nrow(r)) next
    group <- if ("group" %in% names(r)) as.character(r$group[1]) else NA_character_
    starts <- sort(unique(rates$parity_start[rows]))
    ends <- vapply(starts, function(k) rates$parity_end[rows][match(k, rates$parity_start[rows])], 0)
    eligible <- 0 # earliest age at which a woman at the current parity is eligible
    for (b in seq_along(starts)) {
      k <- starts[b]
      at_k <- r[r$parity_start == k, , drop = FALSE]
      if (k > 0) {
        # Report the band's first parity only (pooled bands apply one schedule).
        for (i in seq_len(nrow(at_k))) {
          lo <- at_k$age_start[i]
          if (!is.finite(eligible)) {
            out[[length(out) + 1L]] <- data.frame(group = group, parity = k,
              age_start = lo, age_end = at_k$age_end[i], earliest_eligible = Inf)
            next
          }
          # Cells split at whole ages; report the run of cells ending no later
          # than the earliest eligible age as one range.
          end <- if (at_k$age_end[i] <= eligible) at_k$age_end[i] else floor(eligible)
          if (end > lo) out[[length(out) + 1L]] <- data.frame(group = group,
            parity = k, age_start = lo, age_end = end, earliest_eligible = eligible)
        }
      }
      # Carry eligibility through every parity of this band to the next band.
      if (b < length(starts)) {
        steps <- ends[b] - k
        while (steps > 0 && is.finite(eligible)) {
          birth <- next_birth(at_k, eligible)
          if (!is.finite(birth)) { eligible <- Inf; break }
          inside <- at_k$age_start <= birth & birth < at_k$age_end
          # Within an upper-open positive row every later step adds the interval.
          if (any(inside & is.infinite(at_k$age_end))) {
            eligible <- birth + steps * interval
            break
          }
          eligible <- birth + interval
          steps <- steps - 1
        }
      }
    }
  }
  if (length(out)) do.call(rbind, out) else empty
}

warn_unreachable_parity_cells <- function(rates, interval) {
  # The static reachability heuristic assumes the same schedules throughout a
  # life. Applying it to calendar slices could falsely reject reachable histories.
  if (has_calendar_rates(rates)) return(invisible(NULL))
  cells <- unreachable_parity_cells(rates, interval)
  if (!nrow(cells)) return(invisible(cells))
  keys <- unique(cells[c("group", "parity")])
  parts <- vapply(seq_len(nrow(keys)), function(i) {
    d <- cells[(is.na(cells$group) & is.na(keys$group[i]) | cells$group %in% keys$group[i]) &
               cells$parity == keys$parity[i], ]
    sprintf("%sparity %g at ages [%g, %g)%s",
            if (is.na(keys$group[i])) "" else paste0("group ", keys$group[i], ", "),
            keys$parity[i], min(d$age_start), max(d$age_end),
            if (is.finite(d$earliest_eligible[1]))
              sprintf(" (earliest eligible age %g)", d$earliest_eligible[1])
            else " (unreachable at any age)")
  }, character(1))
  shown <- head(parts, 5)
  if (length(parts) > 5) shown <- c(shown, sprintf("and %d more", length(parts) - 5))
  warning(sprintf(paste0("With `birth_interval` = %g, some positive parity-specific ",
    "fertility rates fall in one-year age cells where no woman can be eligible, so ",
    "those births cannot occur: %s."), interval, paste(shown, collapse = "; ")),
    call. = FALSE)
  invisible(cells)
}

new_rate_process <- function(name, rates) {
  fields <- list(name = name, rates = tibble::as_tibble(rates), rate_unit = "person_year")
  if (!is.null(attr(rates, "source"))) fields$source <- attr(rates, "source")
  structure(fields,
            class = c(paste0("hammel_", name), "hammel_process"))
}

validate_rates <- function(rates, mortality) {
  if (validate_calendar_rates(rates, function(x) validate_rates(x, mortality)))
    return(invisible(rates))
  columns <- c("age_start", "age_end", "rate")
  if (mortality) columns <- c("sex", columns)
  check_columns(rates, columns, "rates")
  if ("group" %in% names(rates)) {
    for (rows in split_rate_groups(rates)) {
      schedule <- rates[rows, setdiff(names(rates), "group"), drop = FALSE]
      validate_rates(schedule, mortality)
    }
    if (!nrow(rates)) stop("`rates` must contain at least one age interval.", call. = FALSE)
    return(invisible(rates))
  }
  if (nrow(rates) == 0L) {
    stop("`rates` must contain at least one age interval.", call. = FALSE)
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  check_numeric_column(rates$rate, "rates$rate")
  ends <- rates$age_end
  if (!is.numeric(ends) || !is.null(dim(ends)) || anyNA(ends) ||
      any(ends <= rates$age_start) || any(rates$age_start < 0)) {
    stop("Age intervals must have nonnegative starts and greater ends (possibly Inf).",
         call. = FALSE)
  }
  if (any(rates$rate < 0)) {
    stop("`rates$rate` must be nonnegative.", call. = FALSE)
  }
  if (mortality) {
    if (!is.character(rates$sex) || !is.null(dim(rates$sex)) ||
        anyNA(rates$sex) || any(!rates$sex %in% c("m", "f")) ||
        !all(c("m", "f") %in% rates$sex)) {
      stop('Mortality rates must include both "m" and "f" sex schedules.',
           call. = FALSE)
    }
    schedules <- split(seq_len(nrow(rates)), rates$sex)
  } else {
    if ("sex" %in% names(rates)) {
      stop("Fertility schedules do not take a sex column in this initial interface.",
           call. = FALSE)
    }
    if ("union_status" %in% names(rates))
      stop("Fertility input `union_status` has been replaced by `union_type`.", call. = FALSE)
    if ("union_type" %in% names(rates)) {
      status <- rates$union_type
      if (!is.character(status) || !is.null(dim(status)) || anyNA(status) ||
          !(setequal(status, c("partnered", "unpartnered")) ||
            setequal(status, c("unpartnered", "cohab", "marriage")))) {
        stop('Fertility `union_type` must contain unpartnered/partnered or unpartnered/cohab/marriage, with no other or missing values.',
             call. = FALSE)
      }
      schedules <- split(seq_len(nrow(rates)), status)
    } else {
      schedules <- list(seq_len(nrow(rates)))
    }
  }
  if (!mortality && any(c("parity_start", "parity_end") %in% names(rates))) {
    check_columns(rates, c("parity_start", "parity_end"), "rates")
    check_numeric_column(rates$parity_start, "rates$parity_start")
    upper <- rates$parity_end
    if (!is.numeric(upper) || !is.null(dim(upper)) || anyNA(upper) ||
        any(rates$parity_start < 0 | rates$parity_start != floor(rates$parity_start)) ||
        any(upper <= rates$parity_start) || any(is.finite(upper) & upper != floor(upper)) ||
        any(rates$parity_start > .Machine$integer.max) ||
        any(is.finite(upper) & upper > .Machine$integer.max))
      stop("Parity bands need nonnegative integer starts and greater integer ends or Inf within engine count capacity.", call. = FALSE)
    schedules <- unlist(lapply(schedules, function(rows) {
      bands <- unique(rates[rows, c("parity_start", "parity_end"), drop = FALSE])
      bands <- bands[order(bands$parity_start), , drop = FALSE]
      if (bands$parity_start[1] != 0 || tail(bands$parity_end, 1) != Inf ||
          (nrow(bands) > 1L && any(head(bands$parity_end, -1) != tail(bands$parity_start, -1))))
        stop("Parity bands must cover [0, Inf) without gaps or overlaps in every group/union-type schedule.", call. = FALSE)
      lapply(seq_len(nrow(bands)), function(i) rows[
        rates$parity_start[rows] == bands$parity_start[i] &
        rates$parity_end[rows] == bands$parity_end[i]])
    }), recursive = FALSE)
  }
  for (rows in schedules) {
    rows <- rows[order(rates$age_start[rows])]
    if (length(rows) > 1L &&
        any(rates$age_start[rows[-1L]] <
            cummax(rates$age_end[rows])[-length(rows)])) {
      stop("Age intervals must not overlap within a rate schedule.", call. = FALSE)
    }
  }
  invisible(rates)
}

validate_ffab <- function(probability) {
  if (!is.numeric(probability) || !is.null(dim(probability)) ||
      length(probability) != 1L || !is.finite(probability) ||
      probability < 0 || probability > 1) {
    stop("`ffab` must be a finite numeric probability in [0, 1].", call. = FALSE)
  }
  invisible(probability)
}
