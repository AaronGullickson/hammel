# Reconstruct wife-side matching choice sets with recorded weight offsets.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Extract matching choice sets from simulation results
#'
#' @description
#'
#' DRAFT (agent-written; needs the user's review). Build wife-side
#' conditional-logit choice sets for unions formed in a simulation run, with
#' each candidate man's log matching weight as an offset.
#'
#' @details
#'
#' The run must keep matching weights, with
#' `run_simulation(..., union_weight_sample = )` above 0. Each choice set is
#' one recorded union: its husband (`chosen = TRUE`) and a uniform sample of
#' the other men eligible when it formed, who were alive, at least the
#' union process's `min_male_age` and in no live union. Known-kin exclusions
#' and hard group restrictions (`NA` in `group_exogamy`) are not applied;
#' drop restricted pairs with the category columns if needed.
#'
#' Each candidate's `log_weight` is the calibrated weight of his matching
#' type (categories and rounded ancestry profile) at that search, 0 when his
#' model has no calibrated weights. With it as an offset, a conditional logit
#' of `chosen` on the matching model's terms estimates its coefficients:
#' for example, `survival::clogit(chosen ~ age_diff + I(age_diff^2) +
#' offset(log_weight) + strata(choice_id), data = sets)`.
#'
#' Sampling uses R's random-number generator; set a seed to reproduce it.
#'
#' @param population A `hammel_population` returned by [run_simulation()] with
#'   `union_weight_sample` above 0.
#' @param unions The number of recorded unions to sample, or all of them if
#'   fewer.
#' @param alternatives The number of other eligible men sampled for each
#'   union, or all of them if fewer.
#' @param union_types The union types to sample from: `"marriage"`,
#'   `"cohabitation"` or both.
#' @return A tibble with one row per candidate: `choice_id` (the union ID),
#'   `union_type`, `time`, `chosen`, `wife_id`, `man_id`, `birth_time_w`,
#'   `birth_time_h`, each matching dimension's category at `time` with `_w`
#'   and `_h` suffixes, rounded ancestry shares named
#'   `<dimension>_<category>_w`/`_h` when matching used ancestry, and
#'   `log_weight`.
#' @export
extract_choice_sets <- function(
  population,
  unions = 1000,
  alternatives = 50,
  union_types = c("marriage", "cohabitation")
) {
  check_analysis_population(population)
  record <- population$matching_weights
  if (is.null(record)) {
    abort_hammel(
      c(
        "{.arg population} has no recorded matching weights.",
        "i" = paste(
          "Run the simulation with {.arg union_weight_sample} above 0, for",
          "example {.code run_simulation(model, years, union_weight_sample =",
          "0.05)}."
        )
      ),
      "population",
      argument = "population"
    )
  }
  check_count(unions, "unions")
  check_count(alternatives, "alternatives")
  check_labels(union_types, "union_types", c("marriage", "cohabitation"))

  kept <- record$unions[record$unions$union_type %in% union_types, ]
  if (nrow(kept) > unions) {
    kept <- kept[sample.int(nrow(kept), unions), ]
  }
  persons <- population$persons
  all_unions <- population$unions
  row <- match(kept$union_id, all_unions$union_id)
  wives <- match(all_unions$wife_id[row], persons$person_id)
  husbands <- match(all_unions$husband_id[row], persons$person_id)

  men <- which(persons$sex == "m")
  sets <- lapply(
    seq_len(nrow(kept)), function(k) {
      t <- kept$time[k]
      live <- all_unions$start_time < t &
        (is.na(all_unions$end_time) | all_unions$end_time > t)
      partnered <- persons$person_id[men] %in% all_unions$husband_id[live]
      eligible <- men[
        t - persons$birth_time[men] >= record$min_male_age &
          (is.na(persons$death_time[men]) | persons$death_time[men] > t) &
          !partnered & men != husbands[k]
      ]
      if (length(eligible) > alternatives) {
        eligible <- eligible[sample.int(length(eligible), alternatives)]
      }
      list(man = c(husbands[k], eligible), k = k)
    }
  )
  size <- vapply(sets, function(x) length(x$man), integer(1))
  k <- rep(seq_len(nrow(kept)), size)
  man <- unlist(lapply(sets, `[[`, "man"))
  wife <- wives[k]
  time <- kept$time[k]
  out <- tibble::tibble(
    choice_id = kept$union_id[k],
    union_type = kept$union_type[k],
    time = time,
    chosen = sequence(size) == 1L,
    wife_id = persons$person_id[wife],
    man_id = persons$person_id[man],
    birth_time_w = persons$birth_time[wife],
    birth_time_h = persons$birth_time[man]
  )

  # Matching categories just before each search: transitions at the same
  # time follow searches.
  dims <- record$dimensions
  man_key <- character(length(man))
  for (d in dims) {
    view <- dimension_view(population, d)
    histories <- if (!is.null(view$transitions)) {
      split(seq_len(nrow(view$transitions)), view$transitions$person_id)
    }
    w <- groups_at(view, wife, time, left = TRUE, histories = histories)
    h <- groups_at(view, man, time, left = TRUE, histories = histories)
    out[[paste0(d, "_w")]] <- as.character(w)
    out[[paste0(d, "_h")]] <- as.character(h)
    man_key <- paste(man_key, h, sep = "\r")
  }
  # Rounded ancestry shares, as matching saw them.
  for (j in seq_along(record$ancestry$dims)) {
    d <- record$ancestry$dims[j]
    components <- record$ancestry$components[[j]]
    table <- population$ancestry[[d]]
    shares <- as.matrix(table[, components, drop = FALSE])
    units <- inspect_round_profile_cpp(shares, record$ancestry$units)
    position <- match(persons$person_id, table$person_id)
    for (g in seq_along(components)) {
      name <- paste0(d, "_", components[g])
      out[[paste0(name, "_w")]] <- units[position[wife], g] /
        record$ancestry$units
      out[[paste0(name, "_h")]] <- units[position[man], g] /
        record$ancestry$units
      man_key <- paste(man_key, units[position[man], g], sep = "\r")
    }
  }

  # Offsets: the man's type entry for the union's calibration, else 0.
  type_id <- record$types$type_id[match(man_key, record$type_keys)]
  calibration <- kept$calibration_id[k]
  weight <- record$weights$log_weight[
    match(
      paste(calibration, type_id),
      paste(record$weights$calibration_id, record$weights$type_id)
    )
  ]
  weight[is.na(weight)] <- 0
  out$log_weight <- weight
  out
}

# INTERNAL HELPERS ---------------------------------------------------------

# The engine's sampled weight record (D200) as tibbles: each type's matching
# categories and rounded ancestry shares, the kept unions and their
# calibrations, and the stored log weights by calibration and type.
format_matching_weights <- function(
  record,
  dimensions,
  registries,
  layout,
  min_male_age
) {
  dims <- dimensions[record$type_dims + 1L]
  codes <- record$type_codes
  units <- record$type_units
  types <- tibble::tibble(type_id = seq_len(nrow(codes)))
  key <- character(nrow(codes))
  for (j in seq_along(dims)) {
    labels <- as.character(registries[[dims[j]]])[codes[, j] + 1L]
    types[[dims[j]]] <- labels
    key <- paste(key, labels, sep = "\r")
  }
  ancestry <- NULL
  if (!is.null(layout)) {
    ancestry <- layout[c("dims", "components", "units")]
    names <- unlist(
      lapply(
        seq_along(layout$dims), function(j) {
          paste0(layout$dims[j], "_", layout$components[[j]])
        }
      )
    )
    for (g in seq_along(names)) {
      types[[names[g]]] <- units[, g] / layout$units
      key <- paste(key, units[, g], sep = "\r")
    }
  }
  list(
    types = types,
    unions = tibble::tibble(
      union_id = record$union_id,
      union_type = c("marriage", "cohabitation")[record$cohab + 1L],
      time = record$time,
      calibration_id = record$calibration_id + 1L
    ),
    weights = tibble::tibble(
      calibration_id = record$weight_calibration + 1L,
      type_id = record$weight_type + 1L,
      log_weight = record$log_weight
    ),
    min_male_age = min_male_age,
    dimensions = dims,
    ancestry = ancestry,
    type_keys = key
  )
}
