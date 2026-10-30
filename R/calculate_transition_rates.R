#' Calculate age- or duration-specific group-transition rates
#'
#' Count individual moves and person-years in their origin group in half-open
#' Lexis cells. Each outgoing destination shares its origin's exposure: do not
#' sum denominators across destinations. Duration is years since the latest
#' group entry, or age before the first recorded transition; return moves reset
#' duration and union changes do not. Events use the exiting group's duration.
#' Death and competing moves end exposure. Rates are cause-specific hazards,
#' not transition probabilities. Internal-boundary events enter the cell
#' starting there; the final calendar boundary is excluded.
#'
#' @param population A grouped `hammel_population` with a known endpoint and
#'   complete person, transition and (when requested) union histories.
#' @param time_breaks Strictly increasing finite calendar boundaries in years.
#' @param age_breaks Optional strictly increasing nonnegative age boundaries.
#' @param dur_breaks Optional strictly increasing nonnegative group-duration
#'   boundaries. Supply exactly one of age_breaks and dur_breaks. A finite final
#'   axis boundary gains Inf.
#' @param by_sex Logical, default FALSE. Disaggregate by `sex` (m/f).
#' @param by_union_type Logical, default FALSE. Disaggregate exposure and events
#'   by current union type. Typed histories use unpartnered/cohab/marriage;
#'   untyped histories use unpartnered/partnered. Requires exclusive histories.
#' @return A tibble with `time_start`, `time_end`, age or duration bounds,
#'   `group_origin`, `group_dest`, optional `sex` and `union_type`, `person_years`,
#'   `transitions`, and `rate_transition` (NA at zero exposure). All distinct
#'   directed pairs among current and historical identities are retained,
#'   including zero-event cells. Fewer than two known groups returns a typed
#'   empty tibble. Period varies fastest within axis bands and pair blocks.
#' @export
#' @examples
#' persons <- data.frame(person_id = 1, sex = "f", birth_time = 0,
#'   death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_, group = "B")
#' population <- hammel_population(persons, end_time = 5,
#'   transitions = data.frame(person_id = 1, time = 2,
#'     group_origin = "A", group_dest = "B"))
#' calculate_transition_rates(population, c(0, 5), age_breaks = c(0, Inf))
calculate_transition_rates <- function(population, time_breaks, age_breaks = NULL,
                                      dur_breaks = NULL, by_sex = FALSE,
                                      by_union_type = FALSE) {
  if (!inherits(population, "hammel_population"))
    stop("`population` must be a grouped hammel_population.", call. = FALSE)
  validate_population(population)
  if (!"group" %in% names(population$persons))
    stop("`calculate_transition_rates()` requires a persons `group` column.", call. = FALSE)
  for (name in c("by_sex", "by_union_type")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || !is.null(dim(value)) || is.na(value))
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  if (is.null(age_breaks) == is.null(dur_breaks))
    stop("Supply exactly one of `age_breaks` and `dur_breaks`.", call. = FALSE)
  check_demographic_breaks(time_breaks, "time_breaks")
  duration <- !is.null(dur_breaks)
  breaks <- if (duration) dur_breaks else age_breaks
  check_demographic_breaks(breaks, if (duration) "dur_breaks" else "age_breaks", ages = TRUE)
  if (is.finite(tail(breaks, 1))) breaks <- c(breaks, Inf)
  if (is.null(population$end_time) || tail(time_breaks, 1) > population$end_time)
    stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
  persons <- prepare_group_observations(population$persons, population, TRUE)
  groups <- observation_groups(persons, TRUE)
  tx <- population$transitions
  if (is.null(tx)) tx <- empty_transitions(persons$group)
  indexed <- data.frame(row = match(tx$person_id, persons$person_id) - 1L,
    time = tx$time, from = match(as.character(tx$group_origin), as.character(groups)) - 1L,
    to = match(as.character(tx$group_dest), as.character(groups)) - 1L)
  typed <- "union_type" %in% names(population$unions)
  types <- if (!by_union_type) "pooled" else if (typed)
    c("unpartnered", "cohab", "marriage") else c("unpartnered", "partnered")
  pairs <- expand.grid(dest = seq_along(groups), origin = seq_along(groups), KEEP.OUT.ATTRS = FALSE)
  pairs <- pairs[pairs$dest != pairs$origin, , drop = FALSE]
  grid <- expand.grid(t = seq_len(length(time_breaks) - 1L),
    a = seq_len(length(breaks) - 1L), pair = seq_len(nrow(pairs)),
    s = seq_len(if (by_sex) 2L else 1L), k = seq_along(types), KEEP.OUT.ATTRS = FALSE)
  counts <- transition_counts_cpp(persons, population$unions, indexed,
    match(persons$group, groups) - 1L, length(groups), time_breaks, breaks,
    duration, by_sex, by_union_type)
  # Native arrays store origin exposure once, while event arrays retain direction.
  nt <- length(time_breaks) - 1L
  na <- length(breaks) - 1L
  ns <- if (by_sex) 2L else 1L
  origin <- pairs$origin[grid$pair] - 1L
  dest <- pairs$dest[grid$pair] - 1L
  cell <- grid$t - 1L + nt * (grid$a - 1L + na * (origin + length(groups) *
    (grid$s - 1L + ns * (grid$k - 1L))))
  event_cell <- grid$t - 1L + nt * (grid$a - 1L + na * (dest + length(groups) *
    (origin + length(groups) * (grid$s - 1L + ns * (grid$k - 1L)))))
  out <- tibble::tibble(time_start = time_breaks[grid$t], time_end = time_breaks[grid$t + 1L])
  axis <- if (duration) "dur" else "age"
  out[[paste0(axis, "_start")]] <- breaks[grid$a]
  out[[paste0(axis, "_end")]] <- breaks[grid$a + 1L]
  out$group_origin <- groups[origin + 1L]
  out$group_dest <- groups[dest + 1L]
  if (by_sex) out$sex <- c("m", "f")[grid$s]
  if (by_union_type) out$union_type <- types[grid$k]
  out$person_years <- counts$person_years[cell + 1]
  out$transitions <- counts$transitions[event_cell + 1]
  out$rate_transition <- as.numeric(ifelse(out$person_years > 0, out$transitions / out$person_years, NA_real_))
  out
}
