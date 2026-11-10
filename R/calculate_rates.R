# Calculate observed demographic event rates from population histories and
# exposure.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Calculate observed demographic event rates
#'
#' @description
#'
#' Calculate age-specific mortality, fertility, union formation, dissolution,
#' cohabitation conversion, and group transition rates from simulated data by
#' time intervals. Rates can be stratified across a variety of dimensions.
#'
#' @details
#'
#' Rates are calculated by counting person-years and events within a given Lexis
#' triangle from a Lexis surface. Some rates allow for either `age_breaks` or
#' `dur_breaks` but one must always be provided. `time_breaks` are also always
#' required but users can specify a single time interval for the entire
#' simulation.
#'
#' All breaks defined by the function use half-open intervals that are inclusive
#' at the lower bound and exclusive at the upper bound.
#'
#' For more details on each rate, see below.
#'
#' @param population A `hammel_population` object returned from
#'   [run_simulation()].
#' @param time_breaks Strictly increasing finite calendar boundaries in years or
#'   fractions of years that define the time intervals of the Lexis surface.
#' @param age_breaks Strictly increasing nonnegative age boundaries in years or
#'   fractions of years that define the age intervals of the Lexis surface. If
#'   the last boundary is finite, an open-ended age interval extended to Inf is
#'   added. Required for mortality, fertility, and marriage. Cohabitation
#'   conversion and group transitions require exactly one of `age_breaks` and
#'   `dur_breaks`. Union formation can take both (see details below).
#' @param dur_breaks Strictly increasing nonnegative duration boundaries in
#'   years or fractions of years that define the duration intervals of the Lexis
#'   surface. If the last boundary is finite, an open-ended duration interval
#'   extended to Inf is added. Required for union dissolution; optional for
#'   union formation; an alternative to age for cohabitation conversion and
#'   group transitions.
#' @param parity_breaks For fertility only: NULL to pool parities, or strictly
#'   increasing nonnegative integer boundaries starting at zero that define
#'   categories for stratifying by parity.
#' @param by Character vector of variables to stratify rates across. Defaults to
#'   `NULL` for no stratification. Accepted stratification variables differ by
#'   calculator; see the sections below.
#' @param conf_level NULL (default) to omit confidence intervals, or a finite
#'   numeric scalar strictly between 0 and 1 indicating the level of confidence
#'   to use for rate confidence intervals. Adds `conf_low` and `conf_high` from
#'   two-sided exact [stats::poisson.test()] intervals for observed event counts
#'   over exposure. Bounds are NA at zero exposure. For fertility, these
#'   describe the observed rate, not the unrealized-adjusted rate.
#' @param use_detailed_previous_union Logical value that indicates whether
#'   stratification by `"previous_union"` will use `"separated"` and `"widowed"`
#'   or a combined `"dissolved"` category.
#' @param unrealized For fertility only: logical, default `FALSE`. When `TRUE`,
#'   adds a column for `rate_fertility_adjusted` that gives the fertility rate
#'   when unrealized births (counted births that could not be realized because
#'   their fertility cell had members but no eligible women).
#' @param group For group transitions only. The declared group dimension whose
#'   transitions are counted. `NULL` (default) means the only declared dimension
#'   and will only work when one group is defined.
#' @return A tibble containing calendar bands, age or duration bands where
#'   applicable, requested strata, event counts, exposure, and the observed
#'   rate. Rates are NA at zero exposure; zero events with positive exposure
#'   give zero. Confidence limits are included only when requested. Person-based
#'   exposure uses `person_years` (woman-years for fertility); union-based
#'   exposure uses `union_years`. Each calculator's section specifies its
#'   count/rate columns, and additional output.
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
#' @section Stratifying with `by`: The values `by` accepts in each calculator.
#'   Group means the name of a declared group dimension, such as `"race"` or
#'   `"region"`.
#'
#'   Person-level rates count events and exposure for individuals:
#'
#'   | Function | Group | Sex | `union_type` | `previous_union` |
#'   |:---|:---:|:---:|:---:|:---:|
#'   | `calculate_mortality_rates()` | Yes | Always | Yes | Yes |
#'   | `calculate_fertility_rates()` | Yes | Women | Yes | Yes |
#'   | `calculate_union_rates()` | Yes | Always | Yes | Yes (a) |
#'   | `calculate_marriage_rates()` | Yes | Always | No | Yes |
#'   | `calculate_transition_rates()` | Yes (b) | Yes | Yes | Yes |
#'
#'   Union-level rates count events and exposure for unions:
#'
#'   | Function | Group | `union_type` | `previous_union` |
#'   |:---|:---:|:---:|:---:|
#'   | `calculate_dissolution_rates()` | Yes (c) | Yes | No |
#'   | `calculate_cohab_convert_rate()` | Yes (c) | No | No |
#'
#'   Always: results are always split by sex, so `sex` is not a `by` value.
#'   Women: fertility is reported for women only.
#'
#'   (a) Supplying `dur_breaks` adds `previous_union` automatically.
#'
#'   (b) Other declared dimensions only, not the one whose transitions are
#'   counted (the `group` argument).
#'
#'   (c) The partners' groups: `"<dim>_w"` for the wife's current category,
#'   `"<dim>_h"` for the husband's, or both for the couple's combination, such
#'   as `c("race_w", "race_h")`. A bare dimension name is ambiguous for a union.
#'
#'   (d)
#'
#' @section Mortality: Use `calculate_mortality_rates()`: Reports deaths per
#'   person-year lived.

#' @section Fertility: Use `calculate_fertility_rates()`: Reports maternal
#'   births per woman-year, including time excluded by minimum birth spacing.
#'   Only female rows are returned. Current parity is reconstructed from all
#'   recorded maternal births, including pre-window births and deceased
#'   children. Setting `unrealized = TRUE` will calculate a separate
#'   `rate_fertility_adjusted` column that calculates fertility rates as if
#'   unrealized births were realized.
#'
#' @section Union formation: Use `calculate_union_rates()`: Reports union starts
#'   (both marriage and cohabitation) per unpartnered person-year lived.
#'   Marriages that occur because of the conversion of a cohabitation are not
#'   counted as events in this function. To calculate cohabitation and direct
#'   marriage separately, use `by = "union_type`.
#'
#'   The `age_breaks` argument is always required for this function, but
#'   `dur_breaks` can also be used when stratified by `previous_union` to get
#'   duration-based rates for union formation since the previous marriage for
#'   individuals who have been previously married.
#'
#' @section Marriage: Use `calculate_marriage_rates()`: Reports direct marriages
#'   and conversions to marriage from cohabitation over nonmarried person-years.
#'   Unpartnered and cohabiting time both contribute exposure. Married time
#'   contributes none.
#'
#' @section Dissolution: Use `calculate_dissolution_rates()`: Reports union
#'   dissolutions over union-years lived. To calculate dissolution separately by
#'   marriage and cohabitation, use `by = "union_type"`.
#'
#'   Dissolution can only be calculated using `dur_breaks` where duration is
#'   measured since the formation of the current union.
#'
#' @section Cohabitation conversion: Use `calculate_cohab_convert_rate()`:
#'   Reports conversions from cohabitation to marriage over cohabitation
#'   union-years.
#'
#' @section Group transitions: use `calculate_transition_rates()`: Reports state
#'   transitions from origin to destination category per person-years lived in
#'   origin category.
#'
#'   When more than one group is in use in a simulation, the `group` argument
#'   must be used to specify which group transition rates are desired.
#'   Transition rates can be stratified by other groups besides the focal group.
#'
#'   Transition rates can be defined by `age_breaks` or `dur_breaks`. If
#'   `dur_breaks` is used, The duration for the first transition is equivalent
#'   to age.
#'
#' @seealso [calculate_population_size()], [calculate_tfr()],
#'   [calculate_patfr()], [calculate_life_table()]
#' @name calculate_rates
#' @rdname calculate_rates
#' @order 1
#' @export
#' @examples
#'
#' # run a simulation
#' pop_start <- generate_uniform_population(1000, end_time = 0)
#' model <- hammel_model() |>
#'   set_population(pop_start) |>
#'   set_fertility(usa_2000_fertility_union) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_union_formation("rate", rate = usa_2000_marriage)
#' pop_end <- run_simulation(model, 100)
#'
#' calculate_mortality_rates(
#'   pop_end,
#'   time_breaks = seq(0, 100, 20),
#'   age_breaks = c(0, 1, seq(5, 100, 5))
#' )
#'
#' calculate_fertility_rates(
#'   pop_end,
#'   time_breaks = seq(0, 100, 20),
#'   age_breaks = seq(15, 40, 5),
#'   by = "union_type"
#' )
calculate_mortality_rates <- function(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "union_type", "previous_union"),
    "calculate_mortality_rates", alias$dims
  )
  check_detailed_previous_union(use_detailed_previous_union, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(
    obs, by, time_breaks, age_breaks,
    sex = TRUE,
    use_detailed_previous_union = use_detailed_previous_union
  )
  out <- person_rate_rows(
    obs, by, counts, time_breaks, age_breaks, c("m", "f"),
    use_detailed_previous_union
  )
  out$deaths <- counts$deaths[attr(out, "index")]
  out$rate_mortality <- as.numeric(
    ifelse(out$person_years > 0, out$deaths / out$person_years, NA_real_)
  )
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
  use_detailed_previous_union = TRUE
) {
  validate_rate_conf_level(conf_level)
  if (!is.logical(unrealized) || length(unrealized) != 1L ||
        is.na(unrealized)) {
    stop("`unrealized` must be TRUE or FALSE.", call. = FALSE)
  }
  if (unrealized && is.null(population$unrealized_births)) {
    stop(
      paste(
        "`unrealized = TRUE` requires a population returned by",
        "run_simulation(), which records unrealized births."
      ),
      call. = FALSE
    )
  }
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "union_type", "previous_union"),
    "calculate_fertility_rates", alias$dims
  )
  check_detailed_previous_union(use_detailed_previous_union, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  if (!is.null(parity_breaks)) {
    check_demographic_breaks(parity_breaks, "parity_breaks", ages = TRUE)
    finite <- parity_breaks[is.finite(parity_breaks)]
    if (parity_breaks[1] != 0 || any(finite != floor(finite))) {
      stop(
        paste0(
          "`parity_breaks` must start at zero and contain ",
          "integer boundaries or Inf."
        ),
        call. = FALSE
      )
    }
    if (is.finite(tail(parity_breaks, 1))) {
      parity_breaks <- c(parity_breaks, Inf)
    }
  }
  obs <- person_observation(population, by, tail(time_breaks, 1))
  counts <- run_person_lexis(
    obs, by, time_breaks, age_breaks,
    sex = TRUE,
    parity_breaks = if (is.null(parity_breaks)) numeric() else parity_breaks,
    use_detailed_previous_union = use_detailed_previous_union
  )
  out <- person_rate_rows(
    obs, by, counts, time_breaks, age_breaks, "f",
    use_detailed_previous_union, parity_breaks
  )
  out$births <- counts$births[attr(out, "index")]
  out$rate_fertility <- as.numeric(
    ifelse(out$person_years > 0, out$births / out$person_years, NA_real_)
  )
  attr(out, "index") <- NULL
  if (unrealized) {
    out <- add_unrealized_births(out, population, alias)
  }
  attr(out, "unassigned_births") <- tibble::tibble(
    time_start = head(time_breaks, -1),
    time_end = tail(time_breaks, -1),
    births = counts$unassigned_births
  )
  out <- finish_strata(
    add_rate_confidence(out, "births", "person_years", conf_level), alias
  )
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
  if (!inherits(population, "hammel_population")) {
    stop(
      "`population` must be a hammel_population with typed union histories.",
      call. = FALSE
    )
  }
  if (is.null(dur_breaks) == is.null(age_breaks)) {
    stop("Supply exactly one of `dur_breaks` and `age_breaks`.", call. = FALSE)
  }
  validate_population(population)
  if (!"union_type" %in% names(population$unions)) {
    stop(
      paste0(
        "Conversion observation requires typed union ",
        "histories with `union_type`."
      ),
      call. = FALSE
    )
  }
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(
    breaks, if (duration) "dur_breaks" else "age_breaks",
    ages = TRUE
  )
  if (is.finite(tail(breaks, 1))) {
    breaks <- c(breaks, Inf)
  }
  if (is.null(population$end_time) ||
        tail(time_breaks, 1) > population$end_time) {
    stop(
      "Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  }
  unions <- population$unions[population$unions$union_type == "cohab", ]
  wives <- match(unions$wife_id, population$persons$person_id)
  origins <- if (duration) {
    unions$start_time
  } else {
    population$persons$birth_time[wives]
  }
  nt <- length(time_breaks) - 1L
  na <- length(breaks) - 1L
  strata <- union_strata(
    population, unions, by,
    stride = nt * na,
    typed_split = FALSE
  )
  axis <- if (duration) "dur" else "age"
  rows <- build_rows(
    c(
      list(band_dim(time_breaks, "time", 1), band_dim(breaks, axis, nt)),
      strata$dims
    )
  )
  counts <- dissolution_counts_cpp(
    strata$unions, time_breaks, breaks,
    strata$index, strata$count,
    origins = origins,
    event_reason = "conversion"
  )
  result <- tibble::as_tibble(
    rows[c("time_start", "time_end", paste0(axis, c("_start", "_end")), by)]
  )
  result$union_years <- counts$union_years[rows$.index]
  result$conversions <- counts$dissolutions[rows$.index]
  result$rate_cohab_convert <- as.numeric(
    ifelse(
      result$union_years > 0,
      result$conversions / result$union_years, NA_real_
    )
  )
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
  if (is.finite(tail(dur_breaks, 1))) {
    dur_breaks <- c(dur_breaks, Inf)
  }
  if (!inherits(population, "hammel_population")) {
    stop("`population` must be a hammel_population.", call. = FALSE)
  }
  validate_population(population)
  if (is.null(population$end_time) ||
        tail(time_breaks, 1) > population$end_time) {
    stop(
      "Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  }
  unions <- population$unions
  nt <- length(time_breaks) - 1L
  nd <- length(dur_breaks) - 1L
  strata <- union_strata(population, unions, by, stride = nt * nd)
  rows <- build_rows(
    c(
      list(band_dim(time_breaks, "time", 1), band_dim(dur_breaks, "dur", nt)),
      strata$dims
    )
  )
  counts <- dissolution_counts_cpp(
    strata$unions, time_breaks, dur_breaks,
    strata$index, strata$count
  )
  result <- tibble::as_tibble(
    rows[c("time_start", "time_end", "dur_start", "dur_end", by)]
  )
  result$union_years <- counts$union_years[rows$.index]
  result$dissolutions <- counts$dissolutions[rows$.index]
  result$rate_dissolution <- as.numeric(
    ifelse(
      result$union_years > 0,
      result$dissolutions / result$union_years, NA_real_
    )
  )
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
  use_detailed_previous_union = TRUE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "previous_union"), "calculate_marriage_rates",
    alias$dims
  )
  finish_strata(
    calculate_entry_rates(
      population, time_breaks, age_breaks, by,
      dur_breaks = NULL,
      marriage_measure = TRUE,
      conf_level = conf_level,
      use_detailed_previous_union = use_detailed_previous_union
    ), alias
  )
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
  use_detailed_previous_union = TRUE,
  group = NULL
) {
  validate_rate_conf_level(conf_level)
  if (!inherits(population, "hammel_population")) {
    stop("`population` must be a grouped hammel_population.", call. = FALSE)
  }
  validate_population(population)
  if (!declared_group(population)) {
    stop(
      "`calculate_transition_rates()` requires a declared group dimension.",
      call. = FALSE
    )
  }
  dims <- names(population$groups)
  if (is.null(group)) {
    if (length(dims) != 1L) {
      stop(
        sprintf(
          paste0(
            "Name the dimension whose transitions are counted ",
            "with `group`: one of %s."
          ),
          paste0("\"", dims, "\"", collapse = ", ")
        ),
        call. = FALSE
      )
    }
    group <- dims
  }
  if (!is.character(group) || length(group) != 1L || !group %in% dims) {
    stop("`group` must name one declared group dimension.", call. = FALSE)
  }
  if (group %in% by) {
    stop(
      sprintf(
        paste0(
          "Transitions in `%s` are split by origin and ",
          "destination, not by `%s` in `by`."
        ),
        group, group
      ),
      call. = FALSE
    )
  }
  # Other dimensions in `by` combine with this one: the internal group is
  # the combination, and only moves in this dimension count.
  external <- check_by(
    by, c("sex", "union_type", "previous_union", setdiff(dims, group)),
    "calculate_transition_rates"
  )
  others <- intersect(external, dims)
  alias <- NULL
  if (length(others)) {
    alias <- observation_alias(population, by, extra = group)
    population <- alias$population
    by <- setdiff(alias$by, "group")
  } else {
    population <- dimension_view(population, group)
  }
  by <- check_by(
    by, c("sex", "union_type", "previous_union"),
    "calculate_transition_rates", alias$dims
  )
  check_detailed_previous_union(use_detailed_previous_union, by)
  if (is.null(age_breaks) == is.null(dur_breaks)) {
    stop("Supply exactly one of `age_breaks` and `dur_breaks`.", call. = FALSE)
  }
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(
    breaks, if (duration) "dur_breaks" else "age_breaks",
    ages = TRUE
  )
  if (is.finite(tail(breaks, 1))) {
    breaks <- c(breaks, Inf)
  }
  obs <- person_observation(
    population, by, tail(time_breaks, 1),
    grouped = TRUE,
    validate = FALSE
  )
  union_vars <- any(c("union_type", "previous_union") %in% by)
  size <- if (is.null(alias)) 0L else length(alias$combined$categories[[1]])
  counts <- run_person_lexis(
    obs, setdiff(by, "sex"), time_breaks, breaks,
    axis = if (duration) 1L else 0L,
    sex = "sex" %in% by,
    group = TRUE,
    use_detailed_previous_union = use_detailed_previous_union,
    exclusive = union_vars,
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
    person_by_dims(by, obs, counts, use_detailed_previous_union)
  )
  rows <- build_rows(dims)
  if (size) {
    # Combined codes: keep moves in this dimension only (other dimensions
    # unchanged), then label each part.
    o <- rows$origin
    d <- rows$destination
    rows <- rows[o %% size != d %% size & o %/% size == d %/% size, ,
      drop = FALSE
    ]
  } else {
    rows <- rows[rows$origin != rows$destination, , drop = FALSE]
  }
  dest <- rows$.key3 - 1L
  axis <- if (duration) "dur" else "age"
  out <- tibble::as_tibble(
    rows[c(
      "time_start", "time_end", paste0(axis, c("_start", "_end")),
      "origin", "destination", by
    )]
  )
  out$person_years <- counts$person_years[rows$.index - dest * counts$cells]
  out$transitions <- counts$transitions[rows$.index]
  out$rate_transition <- as.numeric(
    ifelse(out$person_years > 0, out$transitions / out$person_years, NA_real_)
  )
  if (size) {
    decode <- alias$combined$decode
    code <- out$origin
    for (d in others) {
      out[[d]] <- decode[[d]][code + 1L]
    }
    out$origin <- decode[[group]][code + 1L]
    out$destination <- decode[[group]][out$destination + 1L]
    out <- out[c(
      "time_start", "time_end", paste0(axis, c("_start", "_end")),
      "origin", "destination",
      external, "person_years", "transitions", "rate_transition"
    )]
    # The first `by` variable varies fastest; ties keep the calculator's order.
    keys <- lapply(
      rev(external), function(b) {
        x <- as.character(out[[b]])
        match(
          x, if (b %in% others) {
            as.character(alias$combined$categories[[b]])
          } else {
            unique(x)
          }
        )
      }
    )
    out <- out[do.call(order, c(keys, list(seq_len(nrow(out))))), ,
      drop = FALSE
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
  dur_breaks = NULL,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
) {
  validate_rate_conf_level(conf_level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(
    by, c("group", "previous_union", "union_type"),
    "calculate_union_rates", alias$dims
  )
  finish_strata(
    calculate_entry_rates(
      population, time_breaks, age_breaks, by,
      dur_breaks,
      marriage_measure = FALSE,
      conf_level = conf_level,
      use_detailed_previous_union = use_detailed_previous_union
    ), alias
  )
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Observed event rates and their rate-specific calculation helpers.

# Shared validation of period and age breaks; a finite final age gains Inf.
person_breaks <- function(time_breaks, age_breaks) {
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(age_breaks, "age_breaks", ages = TRUE)
  if (is.finite(tail(age_breaks, 1))) {
    age_breaks <- c(age_breaks, Inf)
  }
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
  use_detailed_previous_union,
  parity_breaks = NULL
) {
  s <- counts$strides
  sex_dim <- label_dim(
    "sex", sexes, s[["sex"]],
    code = match(sexes, c("m", "f")) - 1L
  )
  dims <- c(
    list(
      band_dim(time_breaks, "time", 1), band_dim(
        age_breaks, "age",
        s[["a"]]
      ), sex_dim
    ),
    if (!is.null(parity_breaks)) {
      list(band_dim(parity_breaks, "parity", s[["parity"]]))
    },
    person_by_dims(by, obs, counts, use_detailed_previous_union)
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
# would have to split, become NA, as do periods before the record starts. NA
# parity or union type in a record means that run's schedule pooled it.
add_unrealized_births <- function(out, population, alias = list()) {
  u <- population$unrealized_births
  # Records name each group dimension; the table uses the internal `group`
  # (a combined code when several dimensions are requested).
  if (!is.null(alias$combined)) {
    u$group <- encode_combined(alias, u)
  } else {
    for (from in names(alias$rename)) {
      names(u)[names(u) == alias$rename[[from]]] <- from
    }
  }
  n <- nrow(out)
  value <- numeric(n)
  covered <- out$time_start >= population$unrealized_since
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
        keys <- c(keys, list(p$lo))
        row_keys <- c(row_keys, list(out$parity_start))
        ok <- ok & !is.na(u$parity_start) & p$ok %in% TRUE
      } else {
        ok[] <- FALSE
      } # pooled-parity schedule cannot be split by parity
    }
    if ("group" %in% names(out)) {
      keys <- c(keys, list(as.character(u$group)))
      row_keys <- c(row_keys, list(as.character(out$group)))
    }
    # The fertility schedule never splits by previous union.
    if ("previous_union" %in% names(out)) {
      ok[] <- FALSE
    }
    if ("union_type" %in% names(out)) {
      if ("union_type" %in% names(u)) {
        ok <- ok & u$union_type %in% out$union_type
        keys <- c(keys, list(u$union_type))
        row_keys <- c(row_keys, list(out$union_type))
      } else {
        ok[] <- FALSE
      } # overall schedule cannot be split by union type
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
      if ("parity_start" %in% names(out) && "parity_start" %in% names(u) &&
            !is.na(u$parity_start[r])) {
        hit <- hit & out$parity_start < u$parity_end[r] & out$parity_end >
          u$parity_start[r]
      }
      if ("group" %in% names(out)) {
        hit <- hit & as.character(out$group) == as.character(u$group[r])
      }
      if ("union_type" %in% names(out) && "union_type" %in% names(u) &&
            u$union_type[r] %in% out$union_type) {
        hit <- hit & out$union_type == u$union_type[r]
      }
      unknown <- unknown | hit
    }
  }
  value[unknown] <- NA_real_
  out$unrealized_births <- value
  out$rate_fertility_adjusted <- as.numeric(
    ifelse(
      out$person_years > 0,
      (out$births + out$unrealized_births) / out$person_years, NA_real_
    )
  )
  out
}

check_demographic_breaks <- function(breaks, name, ages = FALSE) {
  if (!is.numeric(breaks) || !is.null(dim(breaks)) || length(breaks) < 2L ||
        anyNA(breaks) ||
        any(!is.finite(head(breaks, -1))) ||
        any(diff(breaks) <= 0) ||
        (!ages && !is.finite(tail(breaks, 1))) ||
        (ages && breaks[1] < 0)) {
    stop(
      sprintf(
        "`%s` must be strictly increasing %sboundaries.", name,
        if (ages) {
          if (name == "dur_breaks") {
            "nonnegative duration "
          } else {
            "nonnegative age "
          }
        } else {
          "finite time "
        }
      ),
      call. = FALSE
    )
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
  use_detailed_previous_union = TRUE
) {
  if (!inherits(population, "hammel_population")) {
    stop(
      "`population` must be a hammel_population with union histories.",
      call. = FALSE
    )
  }
  if (!is.null(dur_breaks)) {
    check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
    if (is.finite(tail(dur_breaks, 1))) {
      dur_breaks <- c(dur_breaks, Inf)
    }
    if (!"previous_union" %in% by) {
      by <- c(by, "previous_union")
    }
  }
  check_detailed_previous_union(use_detailed_previous_union, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  entry_type <- "union_type" %in% by
  obs <- person_observation(
    population, by, tail(time_breaks, 1),
    requires_typed = entry_type
  )
  # union_type here splits formations only; exposure is not split by it.
  state_by <- setdiff(by, "union_type")
  run <- function(axis, breaks, mask) {
    run_person_lexis(
      obs, state_by,
      time_breaks, breaks,
      axis = axis,
      sex = TRUE,
      utype_filter = if (marriage_measure) 2L else 1L,
      previous_mask = mask,
      use_detailed_previous_union = use_detailed_previous_union,
      exclusive = TRUE
    )
  }
  kinds <- if (marriage_measure) {
    list(c(1L, 2L))
  } else if (entry_type) {
    list(0L, 1L)
  } else {
    list(c(0L, 1L))
  }
  block <- function(counts, breaks, axis_name, previous_codes, source) {
    s <- counts$strides
    dims <- c(
      list(
        band_dim(time_breaks, "time", 1), band_dim(breaks, axis_name, s[["a"]]),
        label_dim("sex", c("m", "f"), s[["sex"]])
      ),
      lapply(
        by, function(name) {
          if (name == "union_type") {
            label_dim("union_type", c("cohab", "marriage"), counts$cells)
          } else {
            person_by_dims(
              name, obs, counts, use_detailed_previous_union, previous_codes
            )[[1]]
          }
        }
      )
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
    later <- if (use_detailed_previous_union) 1:2 else 1L
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
    type <- if (entry_type) {
      (rows$.index[here] - 1) %/% cells
    } else {
      rep(0, sum(here))
    }
    base <- rows$.index[here] - type * cells
    exposure[here] <- x$person_years[base]
    for (j in seq_along(kinds)) {
      use <- if (entry_type) type == j - 1L else rep(TRUE, length(type))
      for (kind in kinds[[j]]) {
        events[here][use] <- events[here][use] + x$formations[base[use] +
          cells * kind]
      }
    }
  }
  columns <- c(
    "time_start", "time_end", "sex", "age_start", "age_end",
    if (!is.null(dur_breaks)) c("dur_start", "dur_end"), by
  )
  result <- tibble::as_tibble(rows[columns])
  result$person_years <- exposure
  result$unions <- events
  result$rate_union <- as.numeric(
    ifelse(
      result$person_years > 0,
      result$unions / result$person_years, NA_real_
    )
  )
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
  if (!length(moves) || !moves || !nrow(unions)) {
    return(unions)
  }
  persons <- prepare_group_observations(population$persons, population, TRUE)
  identities <- observation_groups(persons, TRUE)
  indices <- match(persons$group, identities) - 1L
  history <- union_group_history_cpp(
    persons, unions, indices, as.integer(mode),
    typed, as.integer(identity_count), as.integer(group_count)
  )
  if (length(history$row)) {
    attr(unions, ".hammel_group_history") <- history
  }
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
  use_detailed_previous_union = TRUE,
  exclusive = NULL,
  group_target = NULL
) {
  utype <- "union_type" %in% by
  previous <- "previous_union" %in% by
  if (is.null(exclusive)) {
    exclusive <- obs$typed || previous || utype_filter > 0L
  }
  map <- if (use_detailed_previous_union) 0:2 else c(0L, 1L, 1L)
  options <- list(
    typed = obs$typed,
    exclusive = exclusive,
    counts_cohab = obs$counts_cohab,
    axis = as.integer(axis),
    utype_filter = as.integer(utype_filter),
    previous_mask = previous_mask,
    sex = sex,
    group = group,
    union_type = utype,
    previous_union = previous,
    previous_map = map
  )
  if (!is.null(group_target)) {
    options$group_target <- as.integer(group_target)
  }
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  counts <- person_lexis_cpp(
    obs$persons, obs$unions, obs$group_indices, ng,
    time_breaks, axis_breaks, if (length(parity_breaks)) {
      parity_breaks
    } else {
      numeric(1)
    }, options
  )
  nt <- length(time_breaks) - 1L
  na <- length(axis_breaks) - 1L
  sizes <- c(
    t = nt,
    a = na,
    sex = if (sex) 2L else 1L,
    parity = if (length(parity_breaks) > 1L) length(parity_breaks) - 1L else 1L,
    group = if (group) ng else 1L,
    union_type = if (utype) {
      length(union_type_labels(obs$typed))
    } else {
      1L
    },
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
    if (!declared_group(population)) {
      stop(
        paste0(
          "Splitting by \"group_w\" or \"group_h\" requires a ",
          "declared group dimension."
        ),
        call. = FALSE
      )
    }
    persons <- prepare_group_observations(population$persons, population, TRUE)
    identities <- observation_groups(persons, TRUE)
    w <- match(
      persons$group[match(unions$wife_id, persons$person_id)],
      identities
    ) - 1L
    h <- match(
      persons$group[match(unions$husband_id, persons$person_id)],
      identities
    ) - 1L
  }
  ng <- if (is.null(identities)) 1L else length(identities)
  combos <- if (mode == 2L) ng * ng else if (wives || husbands) ng else 1L
  index <- if (mode == 2L) w + ng * h else if (husbands) h else w
  types <- c("cohab", "marriage")
  if (typed_split) {
    if (!"union_type" %in% names(unions)) {
      stop(
        "Splitting by \"union_type\" requires a union `union_type` column.",
        call. = FALSE
      )
    }
    if (!is.character(unions$union_type) || !is.null(dim(unions$union_type)) ||
          anyNA(unions$union_type) ||
          any(!unions$union_type %in% types)) {
      stop(
        paste0(
          "`union_type` must contain only cohab or marriage ",
          "without missing values."
        ),
        call. = FALSE
      )
    }
    index <- index + combos * (match(unions$union_type, types) - 1L)
  }
  if (wives || husbands) {
    unions <- prepare_union_group_observations(
      unions, population, mode,
      typed_split, ng, combos
    )
  }
  dims <- lapply(
    by, function(name) {
      switch(
        name,
        group_w = label_dim("group_w", identities, stride),
        group_h = label_dim(
          "group_h", identities, stride * if (mode == 2L) {
            ng
          } else {
            1L
          }
        ),
        union_type = label_dim("union_type", types, stride * combos)
      )
    }
  )
  list(
    unions = unions,
    index = as.integer(index),
    count = combos * if (typed_split) 2L else 1L,
    dims = dims
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
    if (!column %in% names(data)) {
      return(rep(NA_integer_, nrow(data)))
    }
    code <- code + (
      match(
        as.character(data[[column]]),
        as.character(combined$categories[[k]])
      ) - 1L
    ) * strides[k]
  }
  as.integer(code)
}

# Shared by public event-rate summaries, after their final aggregation.
validate_rate_conf_level <- function(conf_level) {
  if (!is.null(conf_level) &&
      (!is.numeric(conf_level) || length(conf_level) != 1L ||
       !is.null(dim(conf_level)) ||
       !is.finite(conf_level) ||
       conf_level <= 0 ||
       conf_level >= 1)) {
    stop(
      paste0(
        "`conf_level` must be NULL or a finite numeric ",
        "scalar strictly between 0 and 1."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

add_rate_confidence <- function(table, events, exposure, conf_level) {
  validate_rate_conf_level(conf_level)
  if (is.null(conf_level)) {
    return(table)
  }
  table$conf_low <- rep(NA_real_, nrow(table))
  table$conf_high <- rep(NA_real_, nrow(table))
  for (i in which(table[[exposure]] > 0)) {
    limits <- stats::poisson.test(
      x = table[[events]][i],
      T = table[[exposure]][i],
      conf.level = conf_level
    )$conf.int
    table$conf_low[i] <- limits[1L]
    table$conf_high[i] <- limits[2L]
  }
  table
}
