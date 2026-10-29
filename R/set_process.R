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
#'   Scales shared or group-specific input rates before fertility birth-spacing
#'   compensation. Without population groups, simulation warns and ignores it.
#'   Fertility also accepts a named list of vectors for `partnered` and/or
#'   `unpartnered`, or cohab/marriage when rates include those `union_type`
#'   strata. List names must occur in the schedules; omitted entries use 1.0.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
#' @param spacing_method How eligible fertility is adjusted for the birth
#'   interval: `"midpoint"` (default) uses precomputed midpoint compensation;
#'   `"cell"` uses cell-accounted fertility (D101), described in [fertility()].
#' @param inheritance A group identity model from [create_group_inheritance()],
#'   or a method name requiring no additional parameters: `"mother"`, `"father"`,
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
      if (!"spacing_method" %in% names(process)) process$spacing_method <- "midpoint"
      validate_spacing_method(process$spacing_method)
    }
    process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_formation") {
    validate_union_formation(process)
    if (!is.null(process$rates)) process$rates <- tibble::as_tibble(process$rates)
  }
  if (name == "union_dissolution") {
    validate_union_dissolution(process)
    process$rates <- tibble::as_tibble(process$rates)
  }
  model$processes[name] <- list(process)
  model
}

#' @rdname set_process
#' @export
set_mortality <- function(model, rates, group_prop_hazard = NULL) {
  set_process(model, mortality(rates, group_prop_hazard = group_prop_hazard))
}

#' @rdname set_process
#' @export
set_fertility <- function(model, rates, birth_interval = 0.75, group_prop_hazard = NULL,
                          inheritance = create_group_inheritance(), ffab = 0.4886,
                          spacing_method = c("midpoint", "cell")) {
  set_process(model, fertility(rates, birth_interval = birth_interval,
                              group_prop_hazard = group_prop_hazard,
                              inheritance = inheritance, ffab = ffab,
                              spacing_method = spacing_method))
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
#' both statuses. Overlap and fertility compensation are checked within groups.
#' Simulation requires schedules for every population group. Without population
#' groups, simulation warns and uses the first group alphabetically (radix order
#' of character representations, including numeric codes).
#'
#' Mortality requires separate schedules for both sexes. Gaps and uncovered
#' ages are not filled or interpreted as zero rates. Coverage and fertility
#' eligibility must be specified before simulation. Calendar-time-varying
#' schedules are not yet supported; time interval columns are rejected.
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
#' overrides history. First-birth rates are unadjusted. Higher-parity spacing
#' compensation uses the preceding actual parity's unadjusted schedule in the
#' same group/type. In a final open parity band, entry and subsequent births
#' can therefore have different adjusted hazards. This SOCSIM-style approximation
#' does not guarantee exact recovery of arbitrary conditional rates.
#'
#' For finite fertility bands, one adjusted hazard is calculated at each original
#' band's midpoint: the input rate divided by one minus the input fertility integral
#' over the preceding `birth_interval`. The integral accounts for overlaps with
#' other input bands within the same status schedule. Execution rejects
#' nonpositive midpoint eligibility. Birth spacing persists across status changes;
#' compensation for the prior status history is not implemented.
#' Upper-open age bands use the interior-limit adjustment `f / (1 - b * f)`
#' for input rate f and birth interval b, rejecting nonpositive eligibility.
#' This approximates the first b years after entry, where preceding bands may
#' matter. This approximates the age-varying adjustment; initialization is not calibrated.
#'
#' @param rates A nonempty data frame containing numeric `age_start`, `age_end`,
#'   and `rate`. Starts must be finite and nonnegative; ends must exceed starts
#'   and may be `Inf` for upper-open bands. Rates must
#'   be finite and nonnegative. Mortality additionally
#'   requires character `sex` values `"m"` and `"f"`, with both represented.
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative. Applies to fertility
#'   only; eligible hazards use midpoint birth-spacing compensation.
#' @param group_prop_hazard Optional named numeric vector of finite nonnegative
#'   rate multipliers, with unique nonempty group names. Missing groups use 1.0;
#'   zero suppresses the risk. Values are not normalized. At execution, multiply
#'   the selected input schedule by the person's group factor, including all
#'   sex/status strata. Fertility multiplication precedes spacing compensation.
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
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
#' @param spacing_method `"midpoint"` (default) or `"cell"`. Midpoint uses the
#'   compensation described above. `"cell"` (D101) instead groups women into
#'   cells: one-year (or shorter) age sub-bands of each group/union-type
#'   schedule after group multipliers, split at whole ages and band edges, and
#'   one actual parity when rates are parity-specific (a pooled parity band
#'   applies its rate to each contained parity). Each cell produces births at its input rate times
#'   all its living members, including those inside the birth interval, and
#'   each birth goes to a uniformly chosen eligible member. Observed rates then
#'   equal inputs in expectation for any population composition, while a cell
#'   has an eligible member. Intensity in cells with members but none eligible
#'   is dropped and reported in `simulation$unrealized_births`. Midpoint
#'   feasibility checks do not apply. An individual woman's hazard depends on
#'   recent births by others in her cell; this is a calibration device.
#'
#'   For parity-specific schedules with a positive `birth_interval`, `fertility()`
#'   warns when positive rates fall in one-year age cells where no woman can be
#'   eligible under either method: births are singletons at least
#'   `birth_interval` apart, starting from the schedule's first positive rate,
#'   so parity k is first eligible at the age found by following each parity's
#'   first positive rate. Union types are pooled and groups checked separately.
#'   Partly reachable cells and higher actual parities within an open parity
#'   band (always unreachable at the youngest ages) are not reported.
#' @param inheritance A group identity model from [create_group_inheritance()],
#'   or a method name requiring no additional parameters: `"mother"`, `"father"`,
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
mortality <- function(rates, group_prop_hazard = NULL) {
  validate_rates(rates, mortality = TRUE)
  validate_group_prop_hazard(group_prop_hazard)
  process <- new_rate_process("mortality", rates)
  process$group_prop_hazard <- group_prop_hazard
  process
}

#' @rdname mortality
#' @export
fertility <- function(rates, birth_interval = 0.75, group_prop_hazard = NULL,
                      inheritance = create_group_inheritance(), ffab = 0.4886,
                      spacing_method = c("midpoint", "cell")) {
  if (identical(spacing_method, c("midpoint", "cell"))) spacing_method <- "midpoint"
  validate_spacing_method(spacing_method)
  validate_rates(rates, mortality = FALSE)
  validate_fertility_group_prop_hazard(group_prop_hazard, rates)
  inheritance <- normalize_group_inheritance(inheritance)
  validate_birth_interval(birth_interval)
  validate_ffab(ffab)
  process <- new_rate_process("fertility", rates)
  process$inheritance <- if (is.null(inheritance)) create_group_inheritance() else inheritance
  process$birth_interval <- birth_interval
  process$spacing_method <- spacing_method
  process$ffab <- ffab
  process$group_prop_hazard <- group_prop_hazard
  warn_unreachable_parity_cells(process$rates, birth_interval)
  process
}

# Positive parity-specific rates in one-year age cells (the cells used by
# spacing_method = "cell") that no woman can occupy while eligible, given the
# birth interval: singleton births, each at least `interval` apart, starting at
# the schedule's first positive rate. Union types are pooled, because women move
# between them; groups are checked separately. Open parity bands are checked at
# their first parity only: higher actual parities in an open band always have
# unreachable young cells. Cells reachable for part of the year are not
# reported. Returns a data frame of unreachable cells (possibly empty).
unreachable_parity_cells <- function(rates, interval) {
  empty <- data.frame(group = character(), parity = numeric(), age_start = numeric(),
                      age_end = numeric(), earliest_eligible = numeric())
  if (!all(c("parity_start", "parity_end") %in% names(rates)) || !(interval > 0)) return(empty)
  out <- list()
  for (rows in split_rate_groups(rates)) {
    r <- rates[rows, , drop = FALSE]
    r <- r[r$rate > 0, , drop = FALSE]
    if (!nrow(r)) next
    group <- if ("group" %in% names(r)) as.character(r$group[1]) else NA_character_
    finite_ends <- r$parity_end[is.finite(r$parity_end)]
    top <- max(c(r$parity_start, finite_ends - 1))
    eligible <- 0 # earliest age at which a woman at the current parity is eligible
    for (k in 0:top) {
      at_k <- r[r$parity_start <= k & k < r$parity_end, , drop = FALSE]
      if (k > 0) {
        for (i in seq_len(nrow(at_k))) {
          lo <- at_k$age_start[i]
          if (!is.finite(eligible)) {
            # Unreachable at every age (a lower parity has no positive rates):
            # report the whole band, which may be upper-open.
            out[[length(out) + 1L]] <- data.frame(group = group, parity = k,
              age_start = lo, age_end = at_k$age_end[i], earliest_eligible = Inf)
            next
          }
          stop_at <- min(at_k$age_end[i], eligible)
          while (lo < stop_at) {
            hi <- min(at_k$age_end[i], floor(lo) + 1)
            if (hi <= eligible) out[[length(out) + 1L]] <- data.frame(group = group,
              parity = k, age_start = lo, age_end = hi, earliest_eligible = eligible)
            lo <- hi
          }
        }
      }
      # Earliest next birth: the first positive-rate age at parity k from the
      # earliest eligible age onward; no such age makes later parities unreachable.
      reachable <- at_k$age_end > eligible
      eligible <- if (any(reachable))
        min(pmax(at_k$age_start[reachable], eligible)) + interval else Inf
    }
  }
  if (length(out)) do.call(rbind, out) else empty
}

warn_unreachable_parity_cells <- function(rates, interval) {
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

validate_spacing_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("midpoint", "cell")) {
    stop('`spacing_method` must be "midpoint" or "cell".', call. = FALSE)
  }
  invisible(method)
}

new_rate_process <- function(name, rates) {
  fields <- list(name = name, rates = tibble::as_tibble(rates), rate_unit = "person_year")
  if (!is.null(attr(rates, "source"))) fields$source <- attr(rates, "source")
  structure(fields,
            class = c(paste0("hammel_", name), "hammel_process"))
}

validate_rates <- function(rates, mortality) {
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
  if (any(c("time_start", "time_end") %in% names(rates))) {
    stop("Calendar-time-varying rate schedules are not yet supported.", call. = FALSE)
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
