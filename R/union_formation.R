# Construct and validate union-formation processes and search rate schedules.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify a union-formation process
#'
#' @description
#'
#' Constructs a union formation process specification for a simulation model.
#' The convenience function `set_union_formation` is equivalent to
#' `set_process(model, union_formation(...))` and allows for direct piping from
#' a model.
#'
#' @details
#'
#' The union formation process is female-centric. Women initiate searches and
#' select male partners. Two different approaches to union search initiation are
#' provided for the required `initiation` argument:
#'
#' - The "on_demand" approach will initiate a union search when an unpartnered
#' woman gives birth. This approach should not be used with fertility rates that
#' are stratified by `union_type` because unpartnered fertility is what creates
#' the search. This approach is useful if the user has no particular interest in
#' union formation processes but wants to assign parentage correctly for kinship
#' determination.
#' - The "rate" approach uses female marriage rates defined by
#' the `rates` argument. The input should be a marriage rate schedule that may
#' be stratified by group dimensions and/or `previous_union`. When stratified by
#' `previous_union`, remarriages may be duration-specific rather than
#' age-specific. This approach should generally be used with fertility rates
#' that are stratified by `union_type`, although it is not required. Failure to
#' stratify by `union_type` will lead to independence between fertility and
#' union formation.
#'
#' Once a search is initiated, a woman will select a male partner from among all
#' eligible male partners. Men related to the women are considered ineligible
#' out to a kinship depth defined by `kinship_generations` and men below
#' `min_male_age` are also ineligible. the decision itself is determined by the
#' provided `matching_model` which uses a conditional logit formulation to score
#' potential suitors by age differences and group exogamy. See
#' [matching_model()] for more details.
#'
#' It is possible for this search to fail if no eligible men are available. When
#' using the "on_demand" method, the `cancel_birth_if_no_match` argument
#' controls whether the birth continues if no match is found. If set to FALSE,
#' the birth continues but with `NA` values for the father. If set to TRUE,
#' observed fertility rates may deviate from input rates.
#'
#' When unpartnered fertility is allowed, the `matching_model` will also be used
#' to identify the father of a birth to an unpartnered woman. Users can
#' optionally specify a separate `father_matching` model for this case.
#'
#' Cohabitation can be added to the union formation process by specifying
#' `cohab_prob` and `union_convert_rate`. The `cohab_prob` argument expects a
#' data frame that is similar to a rate process object. It contains age-specific
#' bands defined by `age_start` and `age_end`, but also contains a `prob` column
#' rather than a `rate` column. This `prob` column defines the probability that
#' a union formed by a woman of the given age is a cohabitation instead of a
#' marriage. Like rate processes, `cohab_prob` can be stratified by time period
#' using `time_start` and `time_end`, group dimensions, or by `previous_union`.
#'
#' The `union_convert_rate` is not required, but is necessary for cohabitations
#' to be able to convert into marriages. It follows the standard rate process
#' format and can be either duration or age dependent. It can also be stratified
#' by time period and group dimension.
#'
#' By default, partner selection for cohabitation will use the `matching_model`
#' argument. Users can specify a matching model specific to cohabitations with
#' the `cohab_matching` argument.
#'
#' @param initiation Either `"on_demand"` for birth-triggered searches or
#'   `"rate"` for scheduled searches. No default is selected.
#' @param matching A specification created by [matching_model()], or a tibble
#'   with numeric `time_start`, `time_end`, and list-column `model` holding
#'   complete matching objects. Defaults to `age_gap_matching_model(0, 5)`:
#'   preferred age gap zero years, with selection weight halving five years on
#'   either side.
#' @param cancel_birth_if_no_match Logical, default `TRUE`. Cancel a birth when
#'   its on-demand partner search or rate-mode nonmarital father selection
#'   fails. Set `FALSE` to allow births with unknown fathers. Historical missing
#'   father IDs are retained. Partner scarcity can reduce realized fertility
#'   below inputs.
#' @param father_matching Optional matching specification for fathers of
#'   unpartnered births in rate mode, created by [matching_model()], or a
#'   `time_start`/`time_end`/`model` schedule as for `matching`. Defaults to
#'   `NULL`, using the general `matching` model active at that birth.
#' @param min_male_age Minimum male age in continuous years for partner
#'   eligibility, default `15`. Finite nonnegative scalar; zero disables this
#'   age exclusion.
#' @param kinship_generations Generation bound on both sides of a known common
#'   ancestor defining partner ineligiblity, default `2`. Excludes parents,
#'   grandparents, siblings, aunts/uncles, nieces/nephews and first cousins,
#'   including half relationships, when known. Zero imposes no cross-person kin
#'   exclusion. See [extract_kin()] for the shared distance definitions.
#' @param rates For rate mode, a nonempty data frame defining the marriage rate
#'   schedules. Must contain at least numeric `age_start`, `age_end`, and
#'   `rate`. See details for description of stratification variables that can be
#'   applied.
#' @param diagnostics Logical, default `FALSE`. When `TRUE`, adds additional
#'   technical variables from matching model to each union record.
#' @param cohab_prob Optional age-specific probability table defining the
#'   probability that a union is a cohabitation. See details for specification.
#' @param cohab_matching Optional matching model or whole-model schedule for
#'   cohabitation; NULL shares the general `matching` model.
#' @param union_convert_rate Optional age- or duration-specific conversion
#'   table, defining conversion rates of cohabitations into marriages. See
#'   details for specification.
#' @param weights NULL (default) for unweighted selection of searching women, or
#'   a [weight_model()] using person variables, which chooses which unpartnered
#'   woman in a cell searches without changing the cell's input search rate.
#'   Rate mode only.
#' @param weights_conversion NULL (default) for unweighted selection, or a
#'   [weight_model()] using union variables, which chooses which cohabitation in
#'   a cell converts to marriage without changing the cell's input rate.
#' @param cohab_as_previous_union Logical, default FALSE. FALSE
#'   only marriage defines the `previous_union` value; TRUE cohabitation or
#'   marriage define the `previous_union` value.
#' @param model A `hammel_model` object.
#' @return `union_formation()` returns a `hammel_union_formation` process;
#'   `set_union_formation()` returns an updated model.
#' @export
#' @examples
#' union_process <- union_formation(
#'   initiation = "rate",
#'   rates = usa_2000_marriage,
#'   matching = age_gap_matching_model(3, 5)
#' )
#'
#' model <- hammel_model() |>
#'   set_fertility(usa_2000_fertility_union) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_process(union_process)
union_formation <- function(
  initiation,
  rates = NULL,
  matching = age_gap_matching_model(0, 5),
  father_matching = NULL,
  cohab_matching = NULL,
  cohab_prob = NULL,
  union_convert_rate = NULL,
  weights = NULL,
  weights_conversion = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  cancel_birth_if_no_match = TRUE,
  cohab_as_previous_union = FALSE,
  diagnostics = FALSE
) {
  process <- structure(
    list(
      name = "union_formation",
      initiation = initiation,
      matching = matching,
      min_male_age = min_male_age,
      kinship_generations = kinship_generations,
      father_matching = father_matching,
      cancel_birth_if_no_match = cancel_birth_if_no_match,
      diagnostics = diagnostics,
      rates = rates,
      cohab_prob = cohab_prob,
      cohab_matching = cohab_matching,
      union_convert_rate = union_convert_rate,
      cohab_as_previous_union = cohab_as_previous_union,
      weights_conversion = normalize_model_argument(weights_conversion),
      weights = normalize_model_argument(weights)
    ),
    class = c("hammel_union_formation", "hammel_process")
  )
  validate_union_formation(process)
  validate_process_weights(weights_conversion, "union")
  validate_process_weights(weights, "person")
  if (!is.null(weights) && !identical(initiation, "rate")) {
    stop("Search `weights` require rate-based union formation.", call. = FALSE)
  }
  if (!is.null(weights) && !is.null(rates)) {
    warn_cancelled_weight_terms(
      weights, weight_strata("search", rates),
      "search"
    )
  }
  if (!is.null(weights_conversion) && !is.null(union_convert_rate)) {
    warn_cancelled_weight_terms(
      weights_conversion, weight_strata("conversion", union_convert_rate),
      "conversion"
    )
  }
  for (field in c("matching", "father_matching", "cohab_matching")) {
    process[field] <- list(normalize_model_argument(process[[field]]))
  }
  for (field in c("cohab_prob", "union_convert_rate")) {
    if (!is.null(process[[field]])) {
      process[[field]] <- tibble::as_tibble(process[[field]])
    }
  }
  if (!is.null(rates)) {
    process$rates <- tibble::as_tibble(rates)
    if (!is.null(attr(rates, "source"))) {
      process$source <- attr(rates, "source")
    }
  }
  process
}

#' @rdname union_formation
#' @export
set_union_formation <- function(
  model,
  initiation,
  rates = NULL,
  matching = age_gap_matching_model(0, 5),
  father_matching = NULL,
  cohab_matching = NULL,
  cohab_prob = NULL,
  union_convert_rate = NULL,
  weights = NULL,
  weights_conversion = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  cancel_birth_if_no_match = TRUE,
  cohab_as_previous_union = FALSE,
  diagnostics = FALSE
) {
  set_process(
    model, union_formation(
      initiation = initiation,
      rates = rates,
      matching = matching,
      father_matching = father_matching,
      cohab_matching = cohab_matching,
      cohab_prob = cohab_prob,
      union_convert_rate = union_convert_rate,
      weights = weights,
      weights_conversion = weights_conversion,
      min_male_age = min_male_age,
      kinship_generations = kinship_generations,
      cancel_birth_if_no_match = cancel_birth_if_no_match,
      cohab_as_previous_union = cohab_as_previous_union,
      diagnostics = diagnostics
    )
  )
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a union formation process
#'
#' @description
#' 
#' Summary method for class `"hammel_union_formation"`
#'
#' @param object A `hammel_union_formation` object returned by [union_formation()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_union_formation`
#' containing:
#' \describe{
#'   \item{method}{The union initiation method: `"rate"` or `"on_demand"`.}
#'   \item{min_male_age}{The minimum eligible male age, in years.}
#'   \item{kinship_generations}{The configured number of generations used
#'     for kinship exclusion when selecting partners.}
#'   \item{cohab_as_previous_union}{Whether cohabitations count as previous
#'     unions when determining previous-union status.}
#'   \item{rate_strata}{A list describing the union-search rate schedule, or
#'     `NULL` when no rate schedule is configured. The list contains
#'     `stratified_by` (stratum column names), `strata` (a tibble of observed
#'     stratum combinations), `bands` (a tibble of exact band combinations),
#'     and `bands_by_stratum` (whether band combinations differ by stratum).
#'     Shared bands are stored once; otherwise `bands` includes the stratum
#'     columns. Band endpoints define half-open intervals `[start, end)`.}
#'   \item{cohabitation}{`NULL` when cohabitation is disabled; otherwise a
#'     list containing `cohab_prob_range` (the minimum and maximum supplied
#'     cohabitation probabilities), `cohab_rate_strata` (strata and bands of
#'     the cohabitation probability schedule), and `convert_rate_strata`
#'     (strata and bands of the conversion rate schedule, or `NULL` when no
#'     conversion schedule is configured). Both schedule summaries use the
#'     same structure as `rate_strata`.}
#' }
#'
#' @seealso [union_formation()]
#' @examples
#' process <- union_formation("rate", usa_2000_marriage)
#' summary(process)
#' @export
summary.hammel_union_formation <- function(object, ...) {

  structure(
    list(
      method = object$initiation,
      min_male_age = object$min_male_age,
      kinship_generations = object$kinship_generations,
      cohab_as_previous_union = object$cohab_as_previous_union,
      rate_strata = summarize_rate_strata(object$rates),
      cohabitation = if (is.null(object$cohab_prob)) {
        NULL
      } else {
        list(
          cohab_prob_range = range(object$cohab_prob$prob),
          cohab_rate_strata = summarize_rate_strata(object$cohab_prob),
          convert_rate_strata = summarize_rate_strata(object$union_convert_rate)
        )
      }
    ),
    class = "summary.hammel_union_formation"
  )
}

#' @rdname summary.hammel_union_formation
#' @param x A `summary.hammel_union_formation` object returned by
#'   `summary.hammel_union_formation()`.
#' @export
print.summary.hammel_union_formation <- function(x, ...) {

  cat(cli::style_bold("Union Formation Process Summary", "\n\n"))

  print_stat("Search Method:", x$method)

  if(!is.null(x$rate_strata)) {
    print_rate_strata(x$rate_strata)
  }

  cat("\n")

  cat(cli::style_bold("Male eligibility restrictions", "\n"))
  print_stat("Minimum male age:", x$min_male_age)
  print_stat("Generational proximity:", x$kinship_generations)
  
  cat("\n")

  cat(cli::style_bold("Cohabitation", "\n"))
  if(is.null(x$cohabitation)) {
    cat(cli::col_br_red("Disabled"), "\n")
  } else {
    print_stat(
      "Probability range:",
      paste(x$cohabitation$cohab_prob_range, collapse = "\u002d")
    )
    print_stat(
      "Counts as previous union",
      ifelse(x$cohab_as_previous_union, "yes", "no")
    )
    print_rate_strata(x$cohabitation$cohab_rate_strata)
    cat("\n")
    if(!is.null(x$cohabitation$convert_rate_strata)) {
      cat(cli::style_italic("Cohabitation conversion rates"))
      print_rate_strata(x$cohabitation$convert_rate_strata)
    } else {
      cat(cli::style_italic("No cohabitation conversion to marriage allowed"))
    }
  }


}

# HELPER FUNCTIONS ---------------------------------------------------------

validate_union_formation <- function(process) {
  validate_cohabitation(process)
  if (!identical(process$name, "union_formation")) {
    stop("Invalid union-formation specification.", call. = FALSE)
  }
  if (!is.character(process$initiation) || length(process$initiation) != 1L ||
        is.na(process$initiation) ||
        !process$initiation %in% c("rate", "on_demand")) {
    stop('`initiation` must be "rate" or "on_demand".', call. = FALSE)
  }
  age <- process$min_male_age
  if (!is.numeric(age) || length(age) != 1L || !is.null(dim(age)) ||
        !is.finite(age) ||
        age < 0) {
    stop(
      "`min_male_age` must be a finite nonnegative number of years.",
      call. = FALSE
    )
  }
  diagnostic <- process$diagnostics
  if (!is.logical(diagnostic) || length(diagnostic) != 1L ||
        !is.null(dim(diagnostic)) ||
        is.na(diagnostic)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  cancel <- process$cancel_birth_if_no_match
  if (!is.logical(cancel) || length(cancel) != 1L || !is.null(dim(cancel)) ||
        is.na(cancel)) {
    stop("`cancel_birth_if_no_match` must be TRUE or FALSE.", call. = FALSE)
  }
  validate_kinship_generations(process$kinship_generations)
  validate_model_argument(process$matching, validate_matching_model, "matching")
  if (!is.null(process$father_matching)) {
    validate_model_argument(
      process$father_matching, validate_matching_model,
      "father_matching"
    )
  }
  if (process$initiation == "on_demand") {
    if (!is.null(process$rates)) {
      stop("On-demand union formation does not take rates.", call. = FALSE)
    }
  } else {
    validate_union_search_rates(process$rates)
  }
  invisible(process)
}

expand_previous_union_rates <- function(rates) {
  if (!"previous_union" %in% names(rates)) {
    return(rates)
  }
  changed <- FALSE
  expand_slice <- function(x) {
    copies <- lapply(
      split_rate_groups(x), function(rows) {
        copy <- x[rows[x$previous_union[rows] == "dissolved"], , drop = FALSE]
        if (nrow(copy)) {
          changed <<- TRUE
        }
        do.call(
          rbind, lapply(
            c("separated", "widowed"), function(label) {
              out <- copy
              out$previous_union <- rep(label, nrow(copy))
              out
            }
          )
        )
      }
    )
    do.call(
      rbind, c(list(x[x$previous_union != "dissolved", , drop = FALSE]), copies)
    )
  }
  if (!has_calendar_rates(rates)) {
    return(expand_slice(rates))
  }
  pieces <- lapply(
    calendar_rate_segments(rates), function(segment) {
      x <- segment$rates
      if (!nrow(x)) {
        return(NULL)
      }
      x$time_start <- segment$start
      x$time_end <- segment$end
      x <- x[, names(rates), drop = FALSE]
      expand_slice(x)
    }
  )
  if (changed) do.call(rbind, pieces) else rates
}

validate_union_search_rates <- function(rates) {
  check_stratified_columns(rates, "previous_union", "Union-search")
  if (has_calendar_rates(rates)) {
    validate_calendar_rates(rates, validate_union_search_rates)
    effective <- expand_previous_union_rates(rates)
    if (all(c("previous_union", "dur_start") %in% names(effective)) &&
          any(
            vapply(
              split(!is.na(effective$dur_start), effective$previous_union),
              function(x) length(unique(x)) != 1L,
              logical(1)
            )
          )) {
      stop(
        paste0(
          "Each previous_union category must use the same ",
          "age or duration clock across calendar periods."
        ),
        call. = FALSE
      )
    }
    return(invisible(rates))
  }
  check_columns(rates, c("age_start", "age_end", "rate"), "rates")
  check_numeric_column(rates$rate, "rates$rate")
  if (!nrow(rates) || any(rates$rate < 0) ||
        any(c("sex", "time_start", "time_end") %in% names(rates))) {
    stop(
      paste0(
        "Use nonnegative female search rates in valid age ",
        "or duration bands, without sex or calendar columns."
      ),
      call. = FALSE
    )
  }
  duration_columns <- c("dur_start", "dur_end") %in% names(rates)
  if (any(duration_columns) && !all(duration_columns)) {
    stop("Supply both `dur_start` and `dur_end`, or neither.", call. = FALSE)
  }
  boundaries <- c(
    "age_start", "age_end", if (all(duration_columns)) {
      c("dur_start", "dur_end")
    }
  )
  for (column in boundaries) {
    x <- rates[[column]]
    # Bare NA is logical in R; permit all-NA placeholders for an unused axis.
    if (!(is.numeric(x) || (is.logical(x) && all(is.na(x)))) ||
          !is.null(dim(x)) ||
          any(is.nan(x))) {
      stop(
        sprintf(
          paste0(
            "`rates$%s` must contain numeric boundaries or NA ",
            "for an unused clock."
          ), column
        ),
        call. = FALSE
      )
    }
  }
  age <- !is.na(rates$age_start) & !is.na(rates$age_end)
  duration <- if (all(duration_columns)) {
    !is.na(rates$dur_start) & !is.na(rates$dur_end)
  } else {
    rep(FALSE, nrow(rates))
  }
  age_empty <- is.na(rates$age_start) & is.na(rates$age_end)
  duration_empty <- if (all(duration_columns)) {
    is.na(rates$dur_start) & is.na(rates$dur_end)
  } else {
    rep(TRUE, nrow(rates))
  }
  if (any(!((age & duration_empty) | (duration & age_empty)))) {
    stop(
      paste0(
        "Each row requires exactly one complete age or ",
        "duration band; unused boundaries must both be NA."
      ),
      call. = FALSE
    )
  }
  starts <- rates$age_start
  ends <- rates$age_end
  if (any(duration)) {
    if (!"previous_union" %in% names(rates)) {
      stop(
        paste0(
          "Duration search bands require `previous_union` ",
          "('dissolved', 'separated', or 'widowed')."
        ),
        call. = FALSE
      )
    }
    starts[duration] <- rates$dur_start[duration]
    ends[duration] <- rates$dur_end[duration]
  }
  if (any(!is.finite(starts)) || any(starts < 0) || any(ends <= starts)) {
    stop(
      paste0(
        "Use valid nonnegative age or duration bands with ",
        "finite starts and ends greater than starts."
      ),
      call. = FALSE
    )
  }
  strata <- split_rate_groups(rates, candidate_strata(rates))
  if ("previous_union" %in% names(rates)) {
    previous <- rates$previous_union
    categories <- c("none", "separated", "widowed", "dissolved")
    if (!is.character(previous) || !is.null(dim(previous)) || anyNA(previous) ||
          any(!previous %in% categories)) {
      stop(
        paste0(
          "`previous_union` must contain 'none', ",
          "'separated', 'widowed', or pooled 'dissolved'."
        ),
        call. = FALSE
      )
    }
    if (any(
      vapply(
        strata,
        function(rows) {
          !(setequal(previous[rows], c("none", "dissolved")) ||
            setequal(previous[rows], c("none", "separated", "widowed")))
        },
        logical(1)
      )
    )) {
      stop(
        paste0(
          "Supply either 'none' and 'dissolved', or all ",
          "three detailed `previous_union` categories: ",
          "'none', 'separated', 'widowed', within each group ",
          "and calendar period. Do not mix pooled and ",
          "detailed categories."
        ),
        call. = FALSE
      )
    }
    if (any(duration & previous == "none")) {
      stop(
        "`previous_union = 'none'` requires age-specific bands.",
        call. = FALSE
      )
    }
    effective <- expand_previous_union_rates(rates)
    effective_duration <- if (all(duration_columns)) {
      !is.na(effective$dur_start)
    } else {
      rep(FALSE, nrow(effective))
    }
    if (any(
          vapply(
            split(effective_duration, effective$previous_union),
            function(x) length(unique(x)) != 1L,
            logical(1)
          )
        )) {
      stop(
        paste0(
          "Each `previous_union` category must use the same ",
          "age or duration clock for all bands and groups."
        ),
        call. = FALSE
      )
    }
    strata <- unlist(
      lapply(strata, function(rows) split(rows, previous[ rows ])),
      recursive = FALSE
    )
  }
  for (rows in strata) {
    rows <- rows[order(starts[rows])]
    if (length(rows) > 1L &&
          any(tail(starts[rows], -1) < head(cummax(ends[rows]), -1))) {
      stop(
        "Union-search age or duration bands must not overlap.",
        call. = FALSE
      )
    }
  }
  invisible(rates)
}
