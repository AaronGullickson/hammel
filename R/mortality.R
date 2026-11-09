#' Specify age-specific demographic rates
#'
#' Construct a process specification without running a simulation. Rates are
#' piecewise constant over age intervals `[age_start, age_end)`, in years,
#' and expressed per person-year. Overlapping intervals are rejected; adjacent
#' intervals are allowed. Input order and custom columns are preserved; rate
#' schedules are returned as tibbles. A `source` attribute on the input table
#' is copied to the process's `source` field.
#'
#' Columns named by declared group dimensions (see [set_groups()]), such as
#' `race` or `region`, supply separate schedules for each combination of their
#' categories; dimensions without a column share the schedule. Each mortality
#' stratum needs both sexes, and each conditional fertility stratum needs both
#' statuses. Overlap is checked within strata. Simulation requires schedules for
#' every combination of the categories of the dimensions a table uses; a
#' table using an undeclared persons column is an error.
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
#' member. A cell with members but none eligible keeps its intensity; when it
#' would produce a birth, no child is born and the event is counted in
#' `simulation$unrealized_births`. An individual woman's hazard
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
#'   Use [stratify_rates()] to add proportional group differences.
#' @param birth_interval Minimum maternal birth-to-birth interval in years,
#'   defaulting to `0.75`. Must be finite and nonnegative. Applies to fertility
#'   only; see cell-accounted fertility below.
#' @param hetfert NULL (default) for uniform mother selection, or a
#'   [het_fert()] specification for persistent lognormal selection weights.
#' @param weights NULL (default) for unweighted selection within rate cells, or
#'   a [create_weight_model()] specification. Weights choose which member of a
#'   cell experiences an event without changing the cell's input rate. In
#'   fertility they combine with `hetfert` on the log scale.
#' @param ffab Probability that a newborn is female, default `0.4886`.
#'   Must be a finite numeric scalar in `[0, 1]`. Applies to fertility only.
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
#' @order 1
mortality <- function(rates, weights = NULL) {
  validate_rates(rates, mortality = TRUE)
  validate_process_weights(weights, "person")
  warn_cancelled_weight_terms(weights, weight_strata("mortality", rates), "mortality")
  process <- new_rate_process("mortality", rates)
  process$weights <- normalize_model_argument(weights)
  process
}

#' @rdname set_process
#' @export
#' @order 1
set_mortality <- function(model, rates, weights = NULL) {
  set_process(model, mortality(rates, weights = weights))
}
