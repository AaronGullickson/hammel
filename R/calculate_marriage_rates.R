#' Calculate marriage exposure and formation rates by previous union
#'
#' Count union starts and unpartnered person-years by age, period, sex, and
#' optionally `previous_union`: `"none"`, `"dissolved"`, or `"widowed"`. The most recent
#' episode's ending determines the category; each start counts in the category
#' immediately before formation. Partnered time contributes no exposure.
#' Men describe realized matching outcomes, not a male search-rate input.
#' With `dur_breaks`, only none rows use age; dissolved/widowed rows use
#' duration since the most recent ending and always remain separate.
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
#' @param by_group Logical, default FALSE. With TRUE, require person groups and
#'   return separate rows by each individual's group. No pooled row is added.
#' @param by_previous_union Logical, default FALSE. With TRUE, separate exposure
#'   and starts by previous union category and include `previous_union`. Otherwise
#'   pool counts and exposure before calculating rates, and omit that column.
#'   Supplying `dur_breaks` automatically enables this option.
#' @param dur_breaks Optional nonnegative duration boundaries in years. When
#'   supplied, dissolved/widowed rows use years since the latest union ending;
#'   none rows still use `age_breaks`. Automatically enables `by_previous_union`
#'   regardless of its supplied value. A finite final boundary gains Inf.
#'   Adds `dur_start`/`dur_end`, with NA in the unused age/duration pair.
#'   Duration rows include all ages; this is not an age-by-duration surface.
#' @return A tibble with `time_start`, `time_end`, `sex`, `age_start`, `age_end`,
#'   `person_years`, `unions`, and `rate_marriage`, plus `group` and
#'   `previous_union` when requested, and duration columns when `dur_breaks` is
#'   supplied. Empty cells are retained.
#' @export
calculate_marriage_rates <- function(population, time_breaks, age_breaks,
                                     by_group = FALSE, by_previous_union = FALSE,
                                     dur_breaks = NULL) {
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a hammel_population with union histories.", call. = FALSE)
  if (!is.logical(by_previous_union) || length(by_previous_union) != 1L ||
      !is.null(dim(by_previous_union)) || is.na(by_previous_union))
    stop("`by_previous_union` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(dur_breaks)) {
    check_demographic_breaks(dur_breaks, "dur_breaks", ages = TRUE)
    if (is.finite(tail(dur_breaks, 1))) dur_breaks <- c(dur_breaks, Inf)
    by_previous_union <- TRUE
  }
  validate_population(population)
  check_demographic_breaks(time_breaks, "time_breaks")
  check_demographic_breaks(age_breaks, "age_breaks", ages = TRUE)
  if (is.finite(tail(age_breaks, 1))) age_breaks <- c(age_breaks, Inf)
  if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
    stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
  persons <- population$persons
  groups <- observation_groups(persons, by_group)
  group_count <- if (by_group) length(groups) else 1L
  indices <- if (by_group) match(persons$group, groups) - 1L else rep(0L, nrow(persons))
  if (is.null(dur_breaks)) {
    grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
      a = seq_len(length(age_breaks) - 1L), sex = c("m", "f"),
      previous = if (by_previous_union) c("none", "dissolved", "widowed") else "pooled",
      g = seq_len(group_count), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    result <- tibble::tibble(time_start = time_breaks[grid$t],
      time_end = time_breaks[grid$t + 1L], sex = grid$sex,
      age_start = age_breaks[grid$a], age_end = age_breaks[grid$a + 1L])
  } else {
    axis_grid <- function(breaks, categories, duration) {
      grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
        a = seq_len(length(breaks) - 1L), sex = c("m", "f"), previous = categories,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      tibble::tibble(time_start = time_breaks[grid$t], time_end = time_breaks[grid$t + 1L],
        sex = grid$sex, age_start = if (duration) NA_real_ else breaks[grid$a],
        age_end = if (duration) NA_real_ else breaks[grid$a + 1L],
        dur_start = if (duration) breaks[grid$a] else NA_real_,
        dur_end = if (duration) breaks[grid$a + 1L] else NA_real_,
        previous_union = grid$previous)
    }
    one_group <- rbind(axis_grid(age_breaks, "none", FALSE),
      axis_grid(dur_breaks, c("dissolved", "widowed"), TRUE))
    result <- one_group[rep(seq_len(nrow(one_group)), times = group_count), ]
    grid <- list(previous = result$previous_union,
      g = rep(seq_len(group_count), each = nrow(one_group)))
  }
  if (by_previous_union) result$previous_union <- grid$previous
  if (by_group) result$group <- groups[grid$g]
  counts <- marriage_counts_cpp(persons, population$unions, time_breaks, age_breaks,
                               indices, group_count, dur_breaks)
  if (!by_previous_union) {
    cells <- (length(time_breaks) - 1L) * (length(age_breaks) - 1L) * 2L
    pool <- function(values) {
      values <- matrix(values, nrow = cells * 3L, ncol = group_count)
      rows <- seq_len(cells)
      as.numeric(values[rows, , drop = FALSE] +
        values[rows + cells, , drop = FALSE] + values[rows + 2L * cells, , drop = FALSE])
    }
    counts <- lapply(counts, pool)
  }
  result$person_years <- counts$person_years
  result$unions <- counts$unions
  result$rate_marriage <- as.numeric(ifelse(result$person_years > 0,
    result$unions / result$person_years, NA_real_))
  result
}
