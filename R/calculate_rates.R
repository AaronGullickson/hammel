# Calculate observed demographic event rates from population histories and
# exposure.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Calculate observed demographic event rates
#'
#' Calculate mortality, fertility, union formation, marriage, dissolution,
#' cohabitation conversion, and group-transition rates from recorded population
#' histories. Each function counts its events and divides by its own risk
#' exposure; the denominators are not interchangeable. These functions analyze
#' histories without running a simulation or changing the population.
#'
#' @param population One independent history. All functions accept a
#'   `hammel_population` with the required histories and a known endpoint;
#'   observation cannot extend beyond it. Mortality and fertility also accept
#'   bare persons data frames when union stratification is not requested.
#'   Dissolution also accepts a bare union table. For bare tables, the caller
#'   must ensure complete follow-up. See each calculator's section for its
#'   required histories and table columns.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#'   Cells are half-open: events at internal boundaries enter the cell starting
#'   there; the final calendar boundary is excluded.
#' @param age_breaks Strictly increasing nonnegative age boundaries in years. A
#'   finite final boundary gains Inf, creating an open final band starting
#'   there. Required for mortality, fertility, union formation, and marriage.
#'   Conversion and group transitions require exactly one of `age_breaks` and
#'   `dur_breaks`; conversion uses the woman's attained age.
#' @param dur_breaks Strictly increasing nonnegative duration boundaries in
#'   years. A finite final boundary gains Inf. Required for dissolution;
#'   optional for union formation; an alternative to age for conversion and
#'   group transitions. The duration origin and union-formation mixed
#'   age/duration output are described in the corresponding sections below.
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over. Accepted variables and their event-time interpretation differ
#'   by calculator; see the sections below. Declared group dimensions must be
#'   named explicitly, for example `"race"`, `"race_w"`, or `"race_h"`.
#' @param conf_level NULL (default) to omit confidence intervals, or a finite
#'   numeric scalar strictly between 0 and 1. Adds `conf_low` and `conf_high`
#'   from two-sided exact [stats::poisson.test()] intervals for observed event
#'   counts over exposure. Bounds are NA at zero exposure. For fertility, these
#'   describe the observed rate, not the unrealized-adjusted rate.
#' @param pool_dissolved Logical, default `FALSE`. For calculators accepting
#'   this argument, `TRUE` with `"previous_union"` in `by` combines separated
#'   and widowed as `"dissolved"`, pooling counts and exposure before
#'   calculating rates and confidence intervals. Union formation adds
#'   `"previous_union"` automatically when `dur_breaks` is supplied.
#' @param parity_breaks For fertility only: NULL to pool parities, or strictly
#'   increasing nonnegative integer boundaries starting at zero. Bands are
#'   half-open. A finite final boundary receives Inf, as for age breaks.
#' @param unrealized For fertility only: logical, default `FALSE`. With `TRUE`,
#'   the population must come from `simulate()`; add `unrealized_births`
#'   (counted births that could not occur because their fertility cell had
#'   members but no eligible woman, from `simulation$unrealized_births`) and
#'   `rate_fertility_adjusted = (births + unrealized_births) / person_years`.
#'   Values are `NA` where they cannot be attributed exactly: periods before the
#'   most recent simulation's start (restarts keep only the latest run's
#'   record), breaks that cut through a recorded one-year period or fertility
#'   cell, and `by` or parity strata the fertility schedule did not use
#'   (`"previous_union"` always).
#' @param group The declared group dimension whose transitions are counted.
#'   `NULL` (default) means the only declared dimension.
#' @return A tibble containing calendar bands, age or duration bands where
#'   applicable, requested strata, event counts, exposure, and the observed
#'   rate. Rates are NA at zero exposure; zero events with positive exposure
#'   give zero. Confidence limits are included only when requested. Person-based
#'   exposure uses `person_years` (woman-years for fertility); union-based
#'   exposure uses `union_years`. Each calculator's section specifies its
#'   count/rate columns, row ordering, and additional output.
#'
#' @section Choosing a calculator:
#' \tabular{lll}{
#' Calculator \tab Events \tab Exposure\cr
#' Mortality \tab Deaths \tab Person-years alive\cr
#' Fertility \tab Maternal births \tab Woman-years alive\cr
#' Union formation \tab Partnership starts \tab Unpartnered person-years\cr
#' Marriage \tab Marriage starts and conversions \tab Nonmarried person-years\cr
#' Dissolution \tab Dissolutions \tab Union-years\cr
#' Cohabitation conversion \tab Conversions \tab Cohabitation union-years\cr
#' Group transitions \tab Directed group moves \tab Person-years in origin
#' category
#' }
#'
#' @section Person-state stratification:
#' `by` lists categorical variables to split exposure and events over. Each is
#' the person's state at the time, changing only at recorded events:
#' * A declared group dimension's name (such as `"race"`; see [set_groups()]):
#'   the current category, following that dimension's transition history.
#'   Several dimensions may be listed; each gets its own column, with categories
#'   in order of first appearance.
#' * `"union_type"`: current union state, `"unpartnered"`, `"cohab"` or
#'   `"marriage"` for typed histories, or `"unpartnered"` and `"partnered"` for
#'   untyped ones.
#' * `"previous_union"`: how the most recent ended union ended, `"none"`,
#'   `"separated"` or `"widowed"`, including during later unions. Conversion is
#'   not an ending; cohabitation endings follow the population's cohabitation
#'   policy.
#'
#' Union variables require a `hammel_population` with union histories;
#' `"previous_union"` (and typed `"union_type"`) require nonoverlapping
#' histories for each person. Each event counts under the state just before it
#' in the simulation's event order: deaths see no same-time changes, so a death
#' that ends a union counts as partnered; births count in a union starting at
#' the same time and in one ending then by dissolution or conversion, but not
#' one ending by the partner's death, and use the mother's group before any
#' same-time transition.
#'
#' For mortality and fertility, rows vary by time fastest, then age (and parity
#' when requested), then sex, then the `by` variables in the order given,
#' earlier ones varying faster. Strata sum to pooled counts and exposure;
#' recompute rates after pooling.
#'
#' For mortality and fertility, `by` selects declared group dimensions,
#' `"union_type"`, and `"previous_union"`. Group transitions additionally allow
#' `"sex"`. Union entry and union-level rates have the distinct interpretations
#' described in their own sections. Do not sum shared exposure across
#' union-entry types or group-transition destinations.
#'
#' @section Mortality:
#' `calculate_mortality_rates()`: Report deaths per person-year, with exact
#' half-open period and age splitting.
#'
#' \strong{Output:} A tibble containing period, sex, age, the `by` columns,
#' person_years, deaths and rate_mortality.
#'
#' @section Fertility:
#' `calculate_fertility_rates()`: Report maternal births per woman-year,
#' including time excluded by minimum birth spacing. Only female rows are
#' returned. Current parity is reconstructed from all recorded maternal births,
#' including pre-window births and deceased children; custom person parity
#' columns are ignored. Each birth counts at the parity immediately before it
#' and advances subsequent exposure. Simultaneous births progress successively
#' in recorded person-row order.
#'
#' `by` works as in [calculate_mortality_rates()]: births count under the
#' mother's state just before them.
#'
#' \strong{Output:} A tibble with period, sex (female), age, optional
#' parity_start and parity_end, the `by` columns, person_years, births and
#' rate_fertility, and with `unrealized = TRUE`, unrealized_births and
#' rate_fertility_adjusted. Unknown-mother births are in the period-level
#' unassigned_births attribute.
#'
#' @section Union formation:
#' `calculate_union_rates()`: Count union starts and unpartnered person-years by
#' age, period, sex, and optionally `previous_union`: `"none"`, `"separated"`,
#' or `"widowed"`. The most recent episode's ending determines the category;
#' each start counts in the category immediately before formation. Partnered
#' time contributes no exposure. Men describe realized matching outcomes, not a
#' male search-rate input. Typed entries pool cohabitation and direct marriage;
#' conversion is excluded because partners remain partnered. Previous-category
#' classification follows the population's recorded cohabitation policy, default
#' FALSE without metadata. With `dur_breaks`, only none rows use age;
#' separated/widowed rows use duration since the most recent ending;
#' `pool_dissolved` can combine them.
#'
#' Histories must be complete for the intended interpretation. `"none"` means no
#' recorded earlier union; omitted history cannot establish never-partnered
#' status. Overlapping episodes are rejected. Touching episodes are allowed;
#' tied starts follow union ID order. Time and age cells are half-open, with
#' events on internal boundaries assigned to the cell beginning there. Every
#' episode counts, including zero-duration episodes. Rates are formations per
#' unpartnered person-year, NA at zero exposure. Search attempts that fail to
#' form unions are not counted, so these rates can differ from search hazards.
#'
#' \strong{\code{population}:} A `hammel_population` with person and union
#' histories and a known endpoint. Observation cannot extend beyond that
#' endpoint.
#'
#' \strong{\code{by}:} `NULL` (default) or a character vector of variables to
#' split results over: declared group dimensions (each person's own categories
#' before formation), `"previous_union"` (`"none"`, `"separated"` or
#' `"widowed"`, the category immediately before formation) and `"union_type"`
#' (the type of union entered, `"cohab"` or `"marriage"`; requires typed
#' histories). Both `union_type` rows share the full unpartnered exposure, which
#' must not be summed across them. Rows vary by time fastest, then age (or
#' duration), then sex, then the `by` variables in the order given.
#'
#' \strong{\code{dur_breaks}:} Optional nonnegative duration boundaries in
#' years. When supplied, separated/widowed rows use years since the latest union
#' ending; none rows still use `age_breaks`. Adds `"previous_union"` to `by`
#' (last, unless listed). A finite final boundary gains Inf. Adds
#' `dur_start`/`dur_end`, with NA in the unused age/duration pair. Duration rows
#' include all ages; this is not an age-by-duration surface.
#'
#' \strong{Output:} A tibble with `time_start`, `time_end`, `sex`, `age_start`,
#' `age_end`, duration columns when `dur_breaks` is supplied, the `by` columns,
#' `person_years`, `unions`, and `rate_union`. Empty cells are retained. With
#' non-NULL `conf_level`, also includes `conf_low` and `conf_high`.
#'
#' @section Marriage:
#' `calculate_marriage_rates()`: Count direct marriages and conversions to
#' marriage over nonmarried person-years: unpartnered and cohabiting time both
#' contribute exposure. Married time contributes none. Each marriage episode
#' counts once per spouse, including a marriage beginning at conversion. Untyped
#' unions are interpreted as marriages. Men describe realized outcomes, not male
#' search-rate inputs.
#'
#' Complete, nonoverlapping histories are required. Cells are half-open; starts
#' on boundaries enter the cell beginning there. Exposure is censored at death
#' and the observation endpoint. Rates are NA when exposure is zero. This
#' age-specific measure does not provide duration strata. With `by =
#' "previous_union"`, cohabiting exposure and conversion retain the category
#' held immediately before that cohabitation began. Ending categories follow the
#' recorded cohabitation policy, as in [calculate_union_rates()]. For entry into
#' any partnership over unpartnered exposure, use [calculate_union_rates()].
#'
#' \strong{\code{by}:} `NULL` (default) or a character vector of variables to
#' split results over: declared group dimensions and `"previous_union"`, as in
#' [calculate_union_rates()].
#'
#' \strong{Output:} A tibble with `time_start`, `time_end`, `sex`, `age_start`,
#' `age_end`, the `by` columns, `person_years`, `marriages`, and
#' `rate_marriage`. Empty cells are retained. Person-years include cohabiting
#' time. With non-NULL `conf_level`, also includes `conf_low` and `conf_high`.
#'
#' @section Dissolution:
#' `calculate_dissolution_rates()`: Split each union episode across calendar
#' periods and elapsed-duration bands, counting union-years once per episode,
#' not once per spouse. Death and dissolution both stop exposure, but only
#' `end_reason = "dissolution"` counts toward dissolutions. Missing ends are
#' censored at the observation endpoint. Repeated episodes are separate
#' observations; distinct overlapping records are not merged. Analyze one
#' independent history per call, then combine counts.
#'
#' Exposure and event cells are half-open. Dissolutions at internal boundaries
#' belong to the cell beginning there, and events at the last time boundary are
#' excluded. Zero-duration episodes may contribute events with no exposure.
#' Rates are dissolutions per union-year: zero with positive exposure and no
#' events, and NA when exposure is zero. This is a cause-specific rate, not the
#' probability of dissolution in the presence of competing partner deaths.
#'
#' \strong{\code{population}:} A `hammel_population`, or a union data frame
#' containing unique positive whole-number `union_id`, `start_time`, `end_time`,
#' and `end_reason`. Ends and reasons are NA for ongoing episodes; ended reasons
#' are `"dissolution"`, `"death_husband"`, `"death_wife"`, or `"conversion"`.
#' Population objects need a known endpoint; observation cannot extend beyond
#' it. For bare tables, ensure complete follow-up through the requested
#' endpoint.
#'
#' \strong{\code{dur_breaks}:} Strictly increasing nonnegative duration
#' boundaries in years. If the last boundary is finite, append Inf for an open
#' final band.
#'
#' \strong{\code{by}:} `NULL` (default) or a character vector of variables to
#' split results over: `"<group>_w"` and `"<group>_h"` for declared group
#' dimensions, such as `"race_w"` and `"region_h"` (the wife's and husband's
#' current categories, following transition histories; directional, so sum
#' A-wife/B-husband and B-wife/A-husband rows for an unordered pair; any
#' combination of dimensions and partners) and `"union_type"` (cohab or
#' marriage, each with its own exposure; conversion stops cohabitation exposure
#' without counting as dissolution, and duration resets at the start of the
#' linked marriage). Population inputs require declared groups for the group
#' variables; bare union tables require `group_w`, `group_h` or `union_type`
#' columns. A dissolution counts under the groups before any same-time
#' transition.
#'
#' \strong{Output:} A tibble with `time_start`, `time_end`, `dur_start`,
#' `dur_end`, the `by` columns, `union_years`, `dissolutions`, and
#' `rate_dissolution`. Rows vary by time fastest, then duration, then the `by`
#' variables in the order given. With non-NULL `conf_level`, also includes
#' `conf_low` and `conf_high`.
#'
#' @section Cohabitation conversion:
#' `calculate_cohab_convert_rate()`: Count conversions over cohabitation
#' union-years, once per episode. Marriage episodes contribute neither exposure
#' nor events. Conversion, dissolution and partner death end cohabitation
#' exposure; only `end_reason = "conversion"` counts as an event. Ongoing
#' episodes are censored at the observation endpoint. Cells are half-open;
#' boundary events enter the cell beginning there, and events at the final time
#' boundary are excluded. Zero-duration episodes can contribute events without
#' exposure. Rates are NA when exposure is zero.
#'
#' \strong{\code{population}:} A `hammel_population` with typed union histories
#' and a known endpoint. Observation cannot extend beyond it.
#'
#' \strong{\code{dur_breaks}:} Optional strictly increasing nonnegative duration
#' boundaries, in years since cohabitation began. Supply exactly one of
#' `dur_breaks` and `age_breaks`. A finite final boundary gains an open-ended
#' band through Inf.
#'
#' \strong{\code{by}:} `NULL` (default) or a character vector of `"<group>_w"`
#' (the woman's current category in a declared group dimension, such as
#' `"race_w"`, matching conversion input schedules) and `"<group>_h"` (the
#' man's), in any combination of dimensions. Rows vary by time fastest, then the
#' axis band, then the `by` variables in the order given.
#'
#' \strong{Output:} A tibble with `time_start`, `time_end`, `dur_start`,
#' `dur_end` (or `age_start`, `age_end`), the `by` columns, `union_years`,
#' `conversions` and `rate_cohab_convert`. Empty cells are retained. Rates are
#' conversions per cohabitation union-year, not conversion probabilities with
#' competing risks. With non-NULL `conf_level`, also includes `conf_low` and
#' `conf_high`.
#'
#' @section Group transitions:
#' `calculate_transition_rates()`: Count individual moves and person-years in
#' their origin category of one group dimension in half-open Lexis cells. Each
#' outgoing destination shares its origin's exposure: do not sum denominators
#' across destinations. Duration is years since the latest entry into a category
#' of this dimension, or age before its first recorded transition; return moves
#' reset duration and union changes do not. Events use the exiting group's
#' duration. Death and competing moves end exposure. Rates are cause-specific
#' hazards, not transition probabilities. Internal-boundary events enter the
#' cell starting there; the final calendar boundary is excluded.
#'
#' \strong{\code{population}:} A grouped `hammel_population` with a known
#' endpoint and complete person, transition and (when requested) union
#' histories.
#'
#' \strong{\code{dur_breaks}:} Optional strictly increasing nonnegative
#' group-duration boundaries. Supply exactly one of age_breaks and dur_breaks. A
#' finite final axis boundary gains Inf.
#'
#' \strong{\code{by}:} `NULL` (default) or a character vector of variables to
#' split results over: `"sex"`, `"union_type"` (current union state; typed
#' histories use unpartnered/cohab/marriage, untyped unpartnered/partnered) and
#' `"previous_union"`, as in [calculate_mortality_rates()], and other declared
#' group dimensions (the person's category there at the time; only moves in
#' `group` count, and duration is unaffected by moves in other dimensions).
#' Union variables require exclusive histories. A transition counts under the
#' union state after same-time union changes, which precede transitions.
#'
#' \strong{Output:} A tibble with `time_start`, `time_end`, age or duration
#' bounds, `origin`, `destination`, the `by` columns, `person_years`,
#' `transitions`, and `rate_transition` (NA at zero exposure). All distinct
#' directed pairs among current and historical identities are retained,
#' including zero-event cells. Fewer than two known groups returns a typed empty
#' tibble. Rows vary by time fastest, then the axis band, then destination, then
#' origin, then the `by` variables in the order given. With non-NULL
#' `conf_level`, also includes `conf_low` and `conf_high`.
#'
#' @seealso [calculate_population_size()], [calculate_tfr()],
#'   [calculate_patfr()], [calculate_life_table()]
#' @name calculate_rates
#' @rdname calculate_rates
#' @order 1
#' @export
#' @examples
#'
#' # Dissolution
#' unions <- data.frame(
#'   union_id = 1:3, start_time = c(0, 1, 2),
#'   end_time = c(3, 4, NA),
#'   end_reason = c("dissolution", "death_husband", NA)
#' )
#' calculate_dissolution_rates(unions, c(0, 2, 5), c(0, 2, Inf))
#'
#' # Group transitions
#' persons <- data.frame(
#'   person_id = 1, sex = "f", birth_time = 0,
#'   death_time = NA_real_, mother_id = NA_real_,
#'   father_id = NA_real_, group = "B"
#' )
#' population <- hammel_population(persons,
#'   end_time = 5,
#'   transitions = data.frame(
#'     person_id = 1, time = 2,
#'     origin = "A", destination = "B"
#'   ), groups = "group"
#' )
#' calculate_transition_rates(population, c(0, 5), age_breaks = c(0, Inf))
#'
#' # Person and partnership rates from a supplied history:
#' persons <- data.frame(
#'   person_id = 1:3, sex = c("f", "m", "f"),
#'   birth_time = c(-25, -27, 2), death_time = c(NA_real_, 4, NA_real_),
#'   mother_id = c(NA_real_, NA_real_, 1),
#'   father_id = c(NA_real_, NA_real_, 2)
#' )
#' unions <- data.frame(
#'   union_id = 1:2, wife_id = 1, husband_id = 2,
#'   start_time = c(1, 3), end_time = c(3, 4),
#'   end_reason = c("conversion", "death_husband"),
#'   union_type = c("cohab", "marriage"), previous_cohab_id = c(NA_real_, 1)
#' )
#' population <- hammel_population(persons, unions, end_time = 5)
#' calculate_mortality_rates(population, c(0, 5), c(0, Inf))
#' calculate_fertility_rates(population, c(0, 5), c(0, Inf))
#' calculate_union_rates(population, c(0, 5), c(0, Inf))
#' calculate_marriage_rates(population, c(0, 5), c(0, Inf))
#' calculate_dissolution_rates(population, c(0, 5), c(0, Inf))
#' calculate_cohab_convert_rate(population, c(0, 5), dur_breaks = c(0, Inf))
calculate_mortality_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  conf_level = NULL,
  pool_dissolved = FALSE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "union_type", "previous_union"),
    "calculate_mortality_rates", alias$dims
  )
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(obs, by, time_breaks, age_breaks,
    sex = TRUE,
    pool_dissolved = pool_dissolved
  )
  out <- person_rate_rows(
    obs, by, counts, time_breaks, age_breaks, c("m", "f"),
    pool_dissolved
  )
  out$deaths <- counts$deaths[attr(out, "index")]
  out$rate_mortality <- as.numeric(ifelse(out$person_years > 0,
    out$deaths / out$person_years, NA_real_
  ))
  attr(out, "index") <- NULL
  finish_strata(
    add_rate_confidence(out, "deaths", "person_years", conf_level),
    alias
  )
}

#' @rdname calculate_rates
#' @order 2
#' @export
calculate_fertility_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  parity_breaks = NULL,
  unrealized = FALSE,
  conf_level = NULL,
  pool_dissolved = FALSE
) {
  validate_rate_conf_level(conf_level)
  if (!is.logical(unrealized) || length(unrealized) != 1L || is.na(unrealized))
    stop("`unrealized` must be TRUE or FALSE.", call. = FALSE)
  if (unrealized) {
    record <- if (inherits(
      population,
      "hammel_population"
    )) population$simulation else NULL
    if (is.null(record$unrealized_births))
      stop(paste(
        "`unrealized = TRUE` requires a population returned by simulate()",
        "with an unrealized-births record."
      ), call. = FALSE)
  }
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "union_type", "previous_union"),
    "calculate_fertility_rates", alias$dims
  )
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  if (!is.null(parity_breaks)) {
    check_demographic_breaks(parity_breaks, "parity_breaks", ages = TRUE)
    finite <- parity_breaks[is.finite(parity_breaks)]
    if (parity_breaks[1] != 0 || any(finite != floor(finite)))
      stop(paste0(
        "`parity_breaks` must start at zero and contain ",
        "integer boundaries or Inf."
      ), call. = FALSE)
    if (is.finite(tail(parity_breaks, 1))) parity_breaks <- c(
      parity_breaks,
      Inf
    )
  }
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(obs, by, time_breaks, age_breaks,
    sex = TRUE,
    parity_breaks = if (is.null(parity_breaks)) numeric() else parity_breaks,
    pool_dissolved = pool_dissolved
  )
  out <- person_rate_rows(
    obs, by, counts, time_breaks, age_breaks, "f",
    pool_dissolved,
    parity_breaks
  )
  out$births <- counts$births[attr(out, "index")]
  out$rate_fertility <- as.numeric(ifelse(out$person_years > 0,
    out$births / out$person_years, NA_real_
  ))
  attr(out, "index") <- NULL
  if (unrealized) out <- add_unrealized_births(
    out, population$simulation,
    alias
  )
  attr(out, "unassigned_births") <- tibble::tibble(
    time_start = head(
      time_breaks, -1
    ),
    time_end = tail(time_breaks, -1), births = counts$unassigned_births
  )
  out <- finish_strata(add_rate_confidence(
    out, "births", "person_years",
    conf_level
  ), alias)
  out
}

#' @rdname calculate_rates
#' @order 6
#' @export
calculate_cohab_convert_rate <- function(
  population,
  time_breaks,
  dur_breaks = NULL,
  age_breaks = NULL,
  by = NULL,
  conf_level = NULL
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by, c("_w", "_h"))
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group_w", "group_h"), "calculate_cohab_convert_rate",
    alias$dims
  )
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a hammel_population with typed union histories.",
      call. = FALSE
    )
  if (is.null(dur_breaks) == is.null(age_breaks))
    stop("Supply exactly one of `dur_breaks` and `age_breaks`.", call. = FALSE)
  validate_population(population)
  if (!"union_type" %in% names(population$unions))
    stop(paste0(
      "Conversion observation requires typed union ",
      "histories with `union_type`."
    ), call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(breaks, if (duration) "dur_breaks" else "age_breaks",
    ages = TRUE
  )
  if (is.finite(tail(breaks, 1))) breaks <- c(breaks, Inf)
  if (is.null(population$end_time) || tail(time_breaks, 1) >
    population$end_time)
    stop("Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  unions <- population$unions[population$unions$union_type == "cohab", ]
  wives <- match(unions$wife_id, population$persons$person_id)
  origins <- if (duration) unions$start_time else population$persons$birth_time[
    wives
  ]
  nt <- length(time_breaks) - 1L
  na <- length(breaks) - 1L
  strata <- union_strata(population, unions, by,
    stride = nt * na,
    typed_split = FALSE
  )
  axis <- if (duration) "dur" else "age"
  rows <- build_rows(c(
    list(band_dim(time_breaks, "time", 1), band_dim(
      breaks,
      axis, nt
    )),
    strata$dims
  ))
  counts <- dissolution_counts_cpp(strata$unions, time_breaks, breaks,
    strata$index, strata$count,
    origins = origins, event_reason = "conversion"
  )
  result <- tibble::as_tibble(rows[c("time_start", "time_end", paste0(axis, c(
    "_start", "_end"
  )), by)])
  result$union_years <- counts$union_years[rows$.index]
  result$conversions <- counts$dissolutions[rows$.index]
  result$rate_cohab_convert <- as.numeric(ifelse(result$union_years > 0,
    result$conversions / result$union_years, NA_real_
  ))
  result <- finish_strata(
    result, alias, "conversions", "union_years",
    "rate_cohab_convert"
  )
  add_rate_confidence(result, "conversions", "union_years", conf_level)
}

#' @rdname calculate_rates
#' @order 5
#' @export
calculate_dissolution_rates <- function(
  population,
  time_breaks,
  dur_breaks,
  by = NULL,
  conf_level = NULL
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by, c("_w", "_h"))
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group_w", "group_h", "union_type"),
    "calculate_dissolution_rates", alias$dims
  )
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
  if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
  if (inherits(population, "hammel_population")) {
    validate_population(population)
    if (is.null(population$end_time) || tail(time_breaks, 1) >
      population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.",
        call. = FALSE
      )
    unions <- population$unions
  } else {
    unions <- population
    check_columns(
      unions, c("union_id", "start_time", "end_time", "end_reason"),
      "unions"
    )
    check_numeric_column(unions$union_id, "unions$union_id", id = TRUE)
    if (anyDuplicated(unions$union_id))
      stop("`unions$union_id` must be unique; analyze replications separately.",
        call. = FALSE
      )
    check_numeric_column(unions$start_time, "unions$start_time")
    check_numeric_column(unions$end_time, "unions$end_time", missing = TRUE)
    if (any(unions$end_time < unions$start_time, na.rm = TRUE))
      stop("Union end times must not precede start times.", call. = FALSE)
    reason <- unions$end_reason
    if (!is.character(reason) || !is.null(dim(reason)) ||
      any(!is.na(reason) & !reason %in% c(
        "dissolution", "death_husband",
        "death_wife", "conversion"
      )) ||
      any(is.na(reason) != is.na(unions$end_time)))
      stop(paste0(
        "Ended unions require a valid end reason; ongoing ",
        "unions require NA end and reason."
      ), call. = FALSE)
  }
  nt <- length(time_breaks) - 1L
  nd <- length(dur_breaks) - 1L
  strata <- union_strata(population, unions, by, stride = nt * nd)
  rows <- build_rows(c(
    list(band_dim(time_breaks, "time", 1), band_dim(
      dur_breaks, "dur", nt
    )),
    strata$dims
  ))
  counts <- dissolution_counts_cpp(
    strata$unions, time_breaks, dur_breaks,
    strata$index, strata$count
  )
  result <- tibble::as_tibble(rows[c(
    "time_start", "time_end", "dur_start",
    "dur_end", by
  )])
  result$union_years <- counts$union_years[rows$.index]
  result$dissolutions <- counts$dissolutions[rows$.index]
  result$rate_dissolution <- as.numeric(ifelse(result$union_years > 0,
    result$dissolutions / result$union_years, NA_real_
  ))
  result <- finish_strata(
    result, alias, "dissolutions", "union_years",
    "rate_dissolution"
  )
  add_rate_confidence(result, "dissolutions", "union_years", conf_level)
}

#' @rdname calculate_rates
#' @order 4
#' @export
calculate_marriage_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  conf_level = NULL,
  pool_dissolved = FALSE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "previous_union"), "calculate_marriage_rates",
    alias$dims
  )
  finish_strata(calculate_entry_rates(population, time_breaks, age_breaks, by,
    dur_breaks = NULL,
    marriage_measure = TRUE, conf_level = conf_level,
    pool_dissolved = pool_dissolved
  ), alias)
}

#' @rdname calculate_rates
#' @order 7
#' @export
calculate_transition_rates <- function(
  population,
  time_breaks,
  age_breaks = NULL,
  dur_breaks = NULL,
  by = NULL,
  conf_level = NULL,
  pool_dissolved = FALSE,
  group = NULL
) {
  validate_rate_conf_level(conf_level)
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a grouped hammel_population.", call. = FALSE)
  validate_population(population)
  if (!declared_group(population))
    stop("`calculate_transition_rates()` requires a declared group dimension.",
      call. = FALSE
    )
  dims <- names(population$groups)
  if (is.null(group)) {
    if (length(dims) != 1L)
      stop(sprintf(
        paste0(
          "Name the dimension whose transitions are counted ",
          "with `group`: one of %s."
        ),
        paste0("\"", dims, "\"", collapse = ", ")
      ), call. = FALSE)
    group <- dims
  }
  if (!is.character(group) || length(group) != 1L || !group %in% dims)
    stop("`group` must name one declared group dimension.", call. = FALSE)
  if (group %in% by)
    stop(sprintf(
      paste0(
        "Transitions in `%s` are split by origin and ",
        "destination, not by `%s` in `by`."
      ),
      group, group
    ), call. = FALSE)
  # Other dimensions in `by` combine with this one: the internal group is
  # the combination, and only moves in this dimension count.
  external <- check_by(
    by, c("sex", "union_type", "previous_union", setdiff(
      dims, group
    )),
    "calculate_transition_rates"
  )
  others <- intersect(external, dims)
  alias <- NULL
  if (length(others)) {
    alias <- observation_alias(population, by, extra = group)
    population <- alias$population
    by <- setdiff(alias$by, "group")
  } else population <- dimension_view(population, group)
  by <- check_by(
    by, c("sex", "union_type", "previous_union"),
    "calculate_transition_rates", alias$dims
  )
  check_pool_dissolved(pool_dissolved, by)
  if (is.null(age_breaks) == is.null(dur_breaks))
    stop("Supply exactly one of `age_breaks` and `dur_breaks`.", call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(breaks, if (duration) "dur_breaks" else "age_breaks",
    ages = TRUE
  )
  if (is.finite(tail(breaks, 1))) breaks <- c(breaks, Inf)
  obs <- person_observation(population, by, tail(time_breaks, 1),
    grouped =
      TRUE,
    validate = FALSE
  )
  union_vars <- any(c("union_type", "previous_union") %in% by)
  size <- if (is.null(alias)) 0L else length(alias$combined$categories[[1]])
  counts <- run_person_lexis(obs, setdiff(by, "sex"), time_breaks, breaks,
    axis = if (duration) 1L else 0L, sex = "sex" %in% by, group = TRUE,
    pool_dissolved = pool_dissolved, exclusive = union_vars,
    group_target = if (size) obs$groups %% size
  )
  s <- counts$strides
  groups <- obs$groups
  dims <- c(
    list(
      band_dim(time_breaks, "time", 1),
      band_dim(breaks, if (duration) "dur" else "age", s[["a"]]),
      label_dim("destination", groups, counts$cells),
      label_dim("origin", groups, s[["group"]])
    ),
    person_by_dims(by, obs, counts, pool_dissolved)
  )
  rows <- build_rows(dims)
  if (size) {
    # Combined codes: keep moves in this dimension only (other dimensions
    # unchanged), then label each part.
    o <- rows$origin; d <- rows$destination
    rows <- rows[o %% size != d %% size & o %/% size == d %/% size, ,
      drop =
        FALSE
    ]
  } else rows <- rows[rows$origin != rows$destination, , drop = FALSE]
  dest <- rows$.key3 - 1L
  axis <- if (duration) "dur" else "age"
  out <- tibble::as_tibble(rows[c(
    "time_start", "time_end", paste0(axis, c(
      "_start", "_end"
    )),
    "origin", "destination", by
  )])
  out$person_years <- counts$person_years[rows$.index - dest * counts$cells]
  out$transitions <- counts$transitions[rows$.index]
  out$rate_transition <- as.numeric(ifelse(out$person_years > 0,
    out$transitions / out$person_years, NA_real_
  ))
  if (size) {
    decode <- alias$combined$decode
    code <- out$origin
    for (d in others) out[[d]] <- decode[[d]][code + 1L]
    out$origin <- decode[[group]][code + 1L]
    out$destination <- decode[[group]][out$destination + 1L]
    out <- out[c(
      "time_start", "time_end", paste0(axis, c("_start", "_end")),
      "origin", "destination",
      external, "person_years", "transitions", "rate_transition"
    )]
    # The first `by` variable varies fastest; ties keep the calculator's order.
    keys <- lapply(rev(external), function(b) {
      x <- as.character(out[[b]])
      match(x, if (b %in% others) as.character(alias$combined$categories[[
        b
      ]]) else unique(x))
    })
    out <- out[do.call(order, c(keys, list(seq_len(nrow(out))))), ,
      drop =
        FALSE
    ]
  }
  add_rate_confidence(out, "transitions", "person_years", conf_level)
}

#' @rdname calculate_rates
#' @order 3
#' @export
calculate_union_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  dur_breaks = NULL,
  conf_level = NULL,
  pool_dissolved = FALSE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "previous_union", "union_type"),
    "calculate_union_rates", alias$dims
  )
  finish_strata(calculate_entry_rates(population, time_breaks, age_breaks, by,
    dur_breaks,
    marriage_measure = FALSE, conf_level = conf_level,
    pool_dissolved = pool_dissolved
  ), alias)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Observed event rates and their rate-specific calculation helpers.

# Shared validation of period and age breaks; a finite final age gains Inf.
person_breaks <- function(time_breaks, age_breaks) {
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(age_breaks, "age_breaks", ages = TRUE)
  if (is.finite(tail(age_breaks, 1))) age_breaks <- c(age_breaks, Inf)
  age_breaks
}

# Rows for person rates: time, sex, age, optional parity, then `by`.
person_rate_rows <- function(
  obs,
  by,
  counts,
  time_breaks,
  age_breaks,
  sexes,
  pool_dissolved,
  parity_breaks = NULL
) {
  s <- counts$strides
  sex_dim <- label_dim("sex", sexes, s[["sex"]], code = match(sexes, c(
    "m",
    "f"
  )) - 1L)
  dims <- c(
    list(band_dim(time_breaks, "time", 1), band_dim(
      age_breaks, "age",
      s[["a"]]
    ), sex_dim),
    if (!is.null(parity_breaks)) list(band_dim(parity_breaks, "parity", s[[
      "parity"
    ]])),
    person_by_dims(by, obs, counts, pool_dissolved)
  )
  rows <- build_rows(dims)
  columns <- c(
    "time_start", "time_end", "sex", "age_start", "age_end",
    if (!is.null(parity_breaks)) c("parity_start", "parity_end"), by
  )
  out <- tibble::as_tibble(rows[columns])
  out$person_years <- counts$person_years[rows$.index]
  attr(out, "index") <- rows$.index
  out
}

# Attribute recorded unrealized births (D101) to observation rows. A record is
# one run-anchored annual period (or a fragment split at a fertility-rate
# change) and one fertility cell; it is added to
# the row whose period, age, parity, group and union type contain it. Rows that
# a record only partly overlaps, or that a stratum the schedule did not use
# would have to split, become NA, as do periods before the run start.
add_unrealized_births <- function(out, simulation, alias = list()) {
  u <- simulation$unrealized_births
  # Records name each group dimension; the table uses the internal `group`
  # (a combined code when several dimensions are requested).
  if (!is.null(alias$combined)) u$group <- encode_combined(alias, u)
  else for (from in names(alias$rename)) names(u)[names(u) == alias$rename[[
    from
  ]]] <- from
  n <- nrow(out)
  value <- numeric(n)
  covered <- out$time_start >= simulation$start_time
  unknown <- !covered
  if (nrow(u) && n) {
    # Interval indices of each record's bounds within each observation axis.
    within <- function(lo, hi, starts, ends) {
      breaks <- sort(unique(c(starts, ends)))
      i <- findInterval(lo, breaks)
      j <- findInterval(hi, breaks, left.open = TRUE)
      list(lo = breaks[pmax(i, 1L)], ok = i >= 1L & i == j & i < length(breaks))
    }
    t <- within(u$period_start, u$period_end, out$time_start, out$time_end)
    a <- within(u$age_start, u$age_end, out$age_start, out$age_end)
    keys <- list(t$lo, a$lo)
    row_keys <- list(out$time_start, out$age_start)
    ok <- t$ok & a$ok
    # Records outside every observation row contribute nothing anywhere.
    inside <- u$period_start < max(out$time_end) & u$period_end > min(
      out$time_start
    ) &
      u$age_start < max(out$age_end) & u$age_end > min(out$age_start)
    if ("parity_start" %in% names(out)) {
      if ("parity_start" %in% names(u)) {
        p <- within(
          u$parity_start, u$parity_end, out$parity_start,
          out$parity_end
        )
        keys <- c(keys, list(p$lo)); row_keys <- c(row_keys, list(
          out$parity_start
        ))
        ok <- ok & p$ok
      } else ok[] <- FALSE # pooled-parity schedule cannot be split by parity
    }
    if ("group" %in% names(out)) {
      keys <- c(keys, list(as.character(u$group)))
      row_keys <- c(row_keys, list(as.character(out$group)))
    }
    # The fertility schedule never splits by previous union.
    if ("previous_union" %in% names(out)) ok[] <- FALSE
    if ("union_type" %in% names(out)) {
      if ("union_type" %in% names(u)) {
        ok <- ok & u$union_type %in% out$union_type
        keys <- c(keys, list(u$union_type)); row_keys <- c(row_keys, list(
          out$union_type
        ))
      } else ok[] <- FALSE # overall schedule cannot be split by union type
    }
    key <- do.call(paste, c(keys, sep = "\r"))
    row_key <- do.call(paste, c(row_keys, sep = "\r"))
    row <- match(key, row_key)
    exact <- inside & ok & !is.na(row)
    sums <- tapply(u$unrealized_births[exact], row[exact], sum)
    value[as.integer(names(sums))] <- as.numeric(sums)
    # Records that overlap rows without fitting one row make those rows unknown.
    for (r in which(inside & !exact)) {
      hit <- out$time_start < u$period_end[r] & out$time_end > u$period_start[
        r
      ] &
        out$age_start < u$age_end[r] & out$age_end > u$age_start[r]
      if ("parity_start" %in% names(out) && "parity_start" %in% names(u))
        hit <- hit & out$parity_start < u$parity_end[r] & out$parity_end >
          u$parity_start[r]
      if ("group" %in% names(out)) hit <- hit & as.character(out$group) ==
        as.character(u$group[r])
      if ("union_type" %in% names(out) && "union_type" %in% names(u) &&
        u$union_type[r] %in% out$union_type)
        hit <- hit & out$union_type == u$union_type[r]
      unknown <- unknown | hit
    }
  }
  value[unknown] <- NA_real_
  out$unrealized_births <- value
  out$rate_fertility_adjusted <- as.numeric(ifelse(out$person_years > 0,
    (out$births + out$unrealized_births) / out$person_years, NA_real_
  ))
  out
}

check_demographic_breaks <- function(breaks, name, ages = FALSE) {
  if (!is.numeric(breaks) || !is.null(dim(breaks)) || length(breaks) < 2L ||
    anyNA(breaks) ||
    any(!is.finite(head(breaks, -1))) || any(diff(breaks) <= 0) ||
    (!ages && !is.finite(tail(breaks, 1))) || (ages && breaks[1] < 0)) {
    stop(sprintf(
      "`%s` must be strictly increasing %sboundaries.", name,
      if (ages) {
        if (name == "dur_breaks") "nonnegative duration " else
          "nonnegative age "
      } else "finite time "
    ), call. = FALSE)
  }
}

# Union entry (and marriage entry) rates from the shared engine. Formation
# kinds: 0 cohabitation, 1 direct marriage, 2 marriage by conversion.
calculate_entry_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by,
  dur_breaks,
  marriage_measure,
  conf_level = NULL,
  pool_dissolved = FALSE
) {
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a hammel_population with union histories.",
      call. = FALSE
    )
  if (!is.null(dur_breaks)) {
    check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
    if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
    if (!"previous_union" %in% by) by <- c(by, "previous_union")
  }
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  entry_type <- "union_type" %in% by
  obs <- person_observation(population, by, tail(time_breaks, 1),
    needs_unions =
      TRUE,
    requires_typed = entry_type
  )
  # union_type here splits formations only; exposure is not split by it.
  state_by <- setdiff(by, "union_type")
  run <- function(axis, breaks, mask) run_person_lexis(obs, state_by,
    time_breaks, breaks,
    axis = axis, sex = TRUE, utype_filter = if (marriage_measure) 2L else 1L,
    previous_mask = mask, pool_dissolved = pool_dissolved, exclusive = TRUE
  )
  kinds <- if (marriage_measure) list(c(1L, 2L)) else if (entry_type) list(
    0L,
    1L
  ) else list(c(0L, 1L))
  block <- function(counts, breaks, axis_name, previous_codes, source) {
    s <- counts$strides
    dims <- c(
      list(
        band_dim(time_breaks, "time", 1), band_dim(
          breaks, axis_name,
          s[["a"]]
        ),
        label_dim("sex", c("m", "f"), s[["sex"]])
      ),
      lapply(by, function(name) if (name == "union_type")
        label_dim("union_type", c("cohab", "marriage"), counts$cells) else
        person_by_dims(name, obs, counts, pool_dissolved, previous_codes)[[1]])
    )
    rows <- build_rows(dims)
    rows$.source <- rep(source, nrow(rows))
    rows
  }
  if (is.null(dur_breaks)) {
    counts <- list(run(0L, age_breaks, c(TRUE, TRUE, TRUE)))
    rows <- block(counts[[1]], age_breaks, "age", NULL, 1L)
  } else {
    counts <- list(
      run(0L, age_breaks, c(TRUE, FALSE, FALSE)),
      run(2L, dur_breaks, c(FALSE, TRUE, TRUE))
    )
    later <- if (pool_dissolved) 1L else 1:2
    a <- block(counts[[1]], age_breaks, "age", 0L, 1L)
    b <- block(counts[[2]], dur_breaks, "dur", later, 2L)
    a$dur_start <- a$dur_end <- rep(NA_real_, nrow(a))
    b$age_start <- b$age_end <- rep(NA_real_, nrow(b))
    rows <- order_rows_keep(rbind(a, b[names(a)]))
  }
  exposure <- numeric(nrow(rows))
  events <- numeric(nrow(rows))
  for (k in seq_along(counts)) {
    here <- rows$.source == k
    x <- counts[[k]]
    cells <- x$cells
    # Formation kind offsets: the union_type dimension codes cohab 0, marriage
    # 1.
    type <- if (entry_type) (rows$.index[here] - 1) %/% cells else rep(0, sum(
      here
    ))
    base <- rows$.index[here] - type * cells
    exposure[here] <- x$person_years[base]
    for (j in seq_along(kinds)) {
      use <- if (entry_type) type == j - 1L else rep(TRUE, length(type))
      for (kind in kinds[[j]])
        events[here][use] <- events[here][use] + x$formations[base[use] +
          cells * kind]
    }
  }
  columns <- c(
    "time_start", "time_end", "sex", "age_start", "age_end",
    if (!is.null(dur_breaks)) c("dur_start", "dur_end"), by
  )
  result <- tibble::as_tibble(rows[columns])
  result$person_years <- exposure
  result$unions <- events
  result$rate_union <- as.numeric(ifelse(result$person_years > 0,
    result$unions / result$person_years, NA_real_
  ))
  if (marriage_measure) {
    names(result)[names(result) == "unions"] <- "marriages"
    names(result)[names(result) == "rate_union"] <- "rate_marriage"
  }
  add_rate_confidence(
    result, if (marriage_measure) "marriages" else "unions",
    "person_years", conf_level
  )
}

# Attach union-level group histories for `mode` 0 (wife), 1 (husband) or 2
# (both: wife + identities * husband), offset by type when `typed`.
prepare_union_group_observations <- function(
  unions,
  population,
  mode,
  typed,
  identity_count,
  group_count
) {
  tx <- population$transitions
  combined <- population[[".hammel_combined"]]
  moves <- if (is.null(combined)) nrow(tx) else nrow(combined$history)
  if (!length(moves) || !moves || !nrow(unions)) return(unions)
  persons <- prepare_group_observations(population$persons, population, TRUE)
  identities <- observation_groups(persons, TRUE)
  indices <- match(persons$group, identities) - 1L
  history <- union_group_history_cpp(
    persons, unions, indices, as.integer(mode),
    typed, as.integer(identity_count), as.integer(group_count)
  )
  if (length(history$row)) attr(unions, ".hammel_group_history") <- history
  unions
}

# Run the person-level engine. Returns the native counts and the strides of
# each dimension in its arrays, for mapping to output rows.
run_person_lexis <- function(
  obs,
  by,
  time_breaks,
  axis_breaks,
  axis = 0L,
  parity_breaks = numeric(),
  sex = FALSE,
  group = "group" %in% by,
  utype_filter = 0L,
  previous_mask = c(TRUE, TRUE, TRUE),
  pool_dissolved = FALSE,
  exclusive = NULL,
  group_target = NULL
) {
  utype <- "union_type" %in% by
  previous <- "previous_union" %in% by
  if (is.null(exclusive)) exclusive <- obs$typed || previous || utype_filter >
    0L
  map <- if (pool_dissolved) c(0L, 1L, 1L) else 0:2
  options <- list(
    typed = obs$typed, exclusive = exclusive,
    counts_cohab = obs$counts_cohab,
    axis = as.integer(axis), utype_filter = as.integer(utype_filter),
    previous_mask = previous_mask, sex = sex, group = group, union_type = utype,
    previous_union = previous, previous_map = map
  )
  if (!is.null(group_target)) options$group_target <- as.integer(group_target)
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  counts <- person_lexis_cpp(
    obs$persons, obs$unions, obs$group_indices, ng,
    time_breaks, axis_breaks, if (length(parity_breaks)) parity_breaks else
      numeric(1), options
  )
  nt <- length(time_breaks) - 1L
  na <- length(axis_breaks) - 1L
  sizes <- c(
    t = nt, a = na, sex = if (sex) 2L else 1L,
    parity = if (length(parity_breaks) > 1L) length(parity_breaks) - 1L else 1L,
    group = if (group) ng else 1L, union_type = if (utype) length(
      union_type_labels(obs$typed)
    ) else 1L,
    previous_union = if (previous) length(unique(map)) else 1L
  )
  strides <- cumprod(c(1, head(sizes, -1)))
  names(strides) <- names(sizes)
  counts$strides <- strides
  counts$cells <- prod(sizes)
  counts
}

# Order rows combined from several blocks by their keys, the last dimension
# varying slowest.
order_rows_keep <- function(rows) {
  keys <- rev(grep("^\\.key", names(rows), value = TRUE))
  rows[do.call(order, unname(as.list(rows[keys]))), , drop = FALSE]
}

band_dim <- function(breaks, prefix, stride) {
  n <- length(breaks) - 1L
  columns <- data.frame(breaks[seq_len(n)], breaks[seq_len(n) + 1L])
  names(columns) <- paste0(prefix, c("_start", "_end"))
  dim_of(columns, stride = stride)
}

# Group strata for union-level rates (dissolution, conversion): the wife's and
# husband's groups and, optionally, union type. Returns native indices, the
# native group count and the output dimensions for the `by` variables.
union_strata <- function(
  population,
  unions,
  by,
  stride,
  typed_split = "union_type" %in% by
) {
  wives <- "group_w" %in% by
  husbands <- "group_h" %in% by
  mode <- if (wives && husbands) 2L else if (husbands) 1L else 0L
  identities <- NULL
  w <- h <- rep(0L, nrow(unions))
  if (wives || husbands) {
    if (inherits(population, "hammel_population")) {
      if (!declared_group(population))
        stop(paste0(
          "Splitting by \"group_w\" or \"group_h\" requires a ",
          "declared group dimension."
        ), call. = FALSE)
      persons <- prepare_group_observations(
        population$persons, population,
        TRUE
      )
      identities <- observation_groups(persons, TRUE)
      w <- match(
        persons$group[match(unions$wife_id, persons$person_id)],
        identities
      ) - 1L
      h <- match(
        persons$group[match(unions$husband_id, persons$person_id)],
        identities
      ) - 1L
    } else {
      needed <- c(if (wives) "group_w", if (husbands) "group_h")
      missing <- setdiff(needed, names(unions))
      if (length(missing))
        stop(sprintf(
          "Splitting a union table by %s requires union column(s) %s.",
          paste0("\"", needed, "\"", collapse = " and "),
          paste0("`", missing, "`", collapse = ", ")
        ), call. = FALSE)
      values <- c(if (wives) as.character(unions$group_w), if (
        husbands) as.character(unions$group_h))
      if (anyNA(values) || any(!nzchar(trimws(values))))
        stop(paste0(
          "Union `group_w` and `group_h` must contain ",
          "nonmissing, nonempty labels."
        ), call. = FALSE)
      identities <- sort(unique(values), method = "radix")
      if (wives) w <- match(as.character(unions$group_w), identities) - 1L
      if (husbands) h <- match(as.character(unions$group_h), identities) - 1L
    }
  }
  ng <- if (is.null(identities)) 1L else length(identities)
  combos <- if (mode == 2L) ng * ng else if (wives || husbands) ng else 1L
  index <- if (mode == 2L) w + ng * h else if (husbands) h else w
  types <- c("cohab", "marriage")
  if (typed_split) {
    if (!"union_type" %in% names(unions))
      stop("Splitting by \"union_type\" requires a union `union_type` column.",
        call. = FALSE
      )
    if (!is.character(unions$union_type) || !is.null(dim(unions$union_type)) ||
      anyNA(unions$union_type) || any(!unions$union_type %in% types))
      stop(paste0(
        "`union_type` must contain only cohab or marriage ",
        "without missing values."
      ), call. = FALSE)
    index <- index + combos * (match(unions$union_type, types) - 1L)
  }
  if ((wives || husbands) && inherits(population, "hammel_population"))
    unions <- prepare_union_group_observations(
      unions, population, mode,
      typed_split, ng, combos
    )
  dims <- lapply(by, function(name) switch(name,
    group_w = label_dim("group_w", identities, stride),
    group_h = label_dim("group_h", identities, stride * if (mode == 2L) ng else
      1L),
    union_type = label_dim("union_type", types, stride * combos)
  ))
  list(
    unions = unions, index = as.integer(index),
    count = combos * if (typed_split) 2L else 1L, dims = dims
  )
}

# Encode per-dimension labels (columns named `<dim><suffix>`) as combined codes.
encode_combined <- function(alias, data, suffix = "") {
  combined <- alias$combined
  sizes <- vapply(combined$categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  code <- 0
  for (k in seq_along(combined$dims)) {
    column <- paste0(combined$dims[k], suffix)
    if (!column %in% names(data)) return(rep(NA_integer_, nrow(data)))
    code <- code + (match(as.character(data[[column]]), as.character(
      combined$categories[[k]]
    )) - 1L) *
      strides[k]
  }
  as.integer(code)
}

# Shared by public event-rate summaries, after their final aggregation.
validate_rate_conf_level <- function(conf_level) {
  if (!is.null(conf_level) && (!is.numeric(conf_level) ||
    length(conf_level) != 1L || !is.null(dim(conf_level)) ||
    !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1))
    stop(
      paste0(
        "`conf_level` must be NULL or a finite numeric ",
        "scalar strictly between 0 and 1."
      ),
      call. = FALSE
    )
  invisible(NULL)
}

add_rate_confidence <- function(table, events, exposure, conf_level) {
  validate_rate_conf_level(conf_level)
  if (is.null(conf_level)) return(table)
  table$conf_low <- rep(NA_real_, nrow(table))
  table$conf_high <- rep(NA_real_, nrow(table))
  for (i in which(table[[exposure]] > 0)) {
    limits <- stats::poisson.test(
      x = table[[events]][i], T = table[[exposure]][
        i
      ],
      conf.level = conf_level
    )$conf.int
    table$conf_low[i] <- limits[1L]
    table$conf_high[i] <- limits[2L]
  }
  table
}
