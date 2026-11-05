#' Calculate union exposure and formation rates
#'
#' Count union starts and unpartnered person-years by age, period, sex, and
#' optionally `previous_union`: `"none"`, `"separated"`, or `"widowed"`. The most recent
#' episode's ending determines the category; each start counts in the category
#' immediately before formation. Partnered time contributes no exposure.
#' Men describe realized matching outcomes, not a male search-rate input.
#' Typed entries pool cohabitation and direct marriage; conversion is excluded
#' because partners remain partnered. Previous-category classification follows
#' the population's recorded cohabitation policy, default FALSE without metadata.
#' With `dur_breaks`, only none rows use age; separated/widowed rows use
#' duration since the most recent ending; `pool_dissolved` can combine them.
#'
#' Histories must be complete for the intended interpretation. `"none"` means
#' no recorded earlier union; omitted history cannot establish never-partnered
#' status. Overlapping episodes are rejected. Touching episodes are allowed;
#' tied starts follow union ID order. Time and age cells are half-open, with
#' events on internal boundaries assigned to the cell beginning there. Every
#' episode counts, including zero-duration episodes. Rates are formations per
#' unpartnered person-year, NA at zero exposure. Search attempts that fail to
#' form unions are not counted, so these rates can differ from search hazards.
#'
#' @param population A `hammel_population` with person and union histories and
#'   a known endpoint. Observation cannot extend beyond that endpoint.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param age_breaks Strictly increasing nonnegative age boundaries in years.
#'   A finite final boundary gains an open-ended band through Inf.
#' @param by `NULL` (default) or a character vector of variables to split
#'   results over: `"group"` (each person's own group before formation),
#'   `"previous_union"` (`"none"`, `"separated"` or `"widowed"`, the category
#'   immediately before formation) and `"union_type"` (the type of union
#'   entered, `"cohab"` or `"marriage"`; requires typed histories). Both
#'   `union_type` rows share the full unpartnered exposure, which must not be
#'   summed across them. Rows vary by time fastest, then age (or duration),
#'   then sex, then the `by` variables in the order given.
#' @param pool_dissolved Logical, default FALSE. With TRUE, return `"none"`
#'   and `"dissolved"`, pooling separated and widowed events and exposure before
#'   calculating rates and confidence intervals. Requires `"previous_union"` in
#'   `by` (added automatically with `dur_breaks`).
#' @param dur_breaks Optional nonnegative duration boundaries in years. When
#'   supplied, separated/widowed rows use years since the latest union ending;
#'   none rows still use `age_breaks`. Adds `"previous_union"` to `by` (last,
#'   unless listed). A finite final boundary gains Inf.
#'   Adds `dur_start`/`dur_end`, with NA in the unused age/duration pair.
#'   Duration rows include all ages; this is not an age-by-duration surface.
#' @inheritParams calculate_mortality_rates
#' @return A tibble with `time_start`, `time_end`, `sex`, `age_start`, `age_end`,
#'   duration columns when `dur_breaks` is supplied, the `by` columns,
#'   `person_years`, `unions`, and `rate_union`. Empty cells are retained.
#'   With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.
#' @export
calculate_union_rates <- function(population, time_breaks, age_breaks, by = NULL,
                                  dur_breaks = NULL, conf.level = NULL, pool_dissolved = FALSE) {
  validate_rate_conf_level(conf.level)
  alias <- observation_alias(population, by)
  population <- alias$population
  by <- alias$by
  by <- check_by(by, c("group", "previous_union", "union_type"), "calculate_union_rates")
  rename_group_columns(calculate_entry_rates(population, time_breaks, age_breaks, by, dur_breaks,
    marriage_measure = FALSE, conf.level = conf.level, pool_dissolved = pool_dissolved), alias$rename)
}

# Union entry (and marriage entry) rates from the shared engine. Formation
# kinds: 0 cohabitation, 1 direct marriage, 2 marriage by conversion.
calculate_entry_rates <- function(population, time_breaks, age_breaks, by, dur_breaks,
                                  marriage_measure, conf.level = NULL, pool_dissolved = FALSE) {
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a hammel_population with union histories.", call. = FALSE)
  if (!is.null(dur_breaks)) {
    check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
    if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
    if (!"previous_union" %in% by) by <- c(by, "previous_union")
  }
  check_pool_dissolved(pool_dissolved, by)
  age_breaks <- person_breaks(time_breaks, age_breaks)
  entry_type <- "union_type" %in% by
  obs <- person_observation(population, by, tail(time_breaks, 1), needs_unions = TRUE,
                            requires_typed = entry_type)
  # union_type here splits formations only; exposure is not split by it.
  state_by <- setdiff(by, "union_type")
  run <- function(axis, breaks, mask) run_person_lexis(obs, state_by, time_breaks, breaks,
    axis = axis, sex = TRUE, utype_filter = if (marriage_measure) 2L else 1L,
    previous_mask = mask, pool_dissolved = pool_dissolved, exclusive = TRUE)
  kinds <- if (marriage_measure) list(c(1L, 2L)) else if (entry_type) list(0L, 1L) else list(c(0L, 1L))
  block <- function(counts, breaks, axis_name, previous_codes, source) {
    s <- counts$strides
    dims <- c(list(band_dim(time_breaks, "time", 1), band_dim(breaks, axis_name, s[["a"]]),
      label_dim("sex", c("m", "f"), s[["sex"]])),
      lapply(by, function(name) if (name == "union_type")
        label_dim("union_type", c("cohab", "marriage"), counts$cells) else
        person_by_dims(name, obs, counts, pool_dissolved, previous_codes)[[1]]))
    rows <- build_rows(dims)
    rows$.source <- rep(source, nrow(rows))
    rows
  }
  if (is.null(dur_breaks)) {
    counts <- list(run(0L, age_breaks, c(TRUE, TRUE, TRUE)))
    rows <- block(counts[[1]], age_breaks, "age", NULL, 1L)
  } else {
    counts <- list(run(0L, age_breaks, c(TRUE, FALSE, FALSE)),
                   run(2L, dur_breaks, c(FALSE, TRUE, TRUE)))
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
    # Formation kind offsets: the union_type dimension codes cohab 0, marriage 1.
    type <- if (entry_type) (rows$.index[here] - 1) %/% cells else rep(0, sum(here))
    base <- rows$.index[here] - type * cells
    exposure[here] <- x$person_years[base]
    for (j in seq_along(kinds)) {
      use <- if (entry_type) type == j - 1L else rep(TRUE, length(type))
      for (kind in kinds[[j]])
        events[here][use] <- events[here][use] + x$formations[base[use] + cells * kind]
    }
  }
  columns <- c("time_start", "time_end", "sex", "age_start", "age_end",
               if (!is.null(dur_breaks)) c("dur_start", "dur_end"), by)
  result <- tibble::as_tibble(rows[columns])
  result$person_years <- exposure
  result$unions <- events
  result$rate_union <- as.numeric(ifelse(result$person_years > 0,
    result$unions / result$person_years, NA_real_))
  if (marriage_measure) {
    names(result)[names(result) == "unions"] <- "marriages"
    names(result)[names(result) == "rate_union"] <- "rate_marriage"
  }
  add_rate_confidence(result, if (marriage_measure) "marriages" else "unions",
    "person_years", conf.level)
}
