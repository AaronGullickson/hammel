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
  new_union_formation(
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
  attach_process(
    model, new_union_formation(
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
#'   \item{cancel_birth_if_no_match}{Whether a birth is cancelled when its
#'     on-demand partner search or rate-mode nonmarital father selection
#'     fails, rather than occurring with an unknown father.}
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
#'     the cohabitation probability schedule), `convert_rate_strata`
#'     (strata and bands of the conversion rate schedule, or `NULL` when no
#'     conversion schedule is configured), `matching` (a summary of the
#'     cohabitation matching model, or `NULL` when the general matching
#'     model is used) and `conversion_weights` (a summary of the conversion
#'     weight model, or `NULL` when unweighted). Both schedule summaries use the
#'     same structure as `rate_strata`.}
#'   \item{matching}{A summary of the matching model: a
#'     `summary.hammel_matching_model`, or a list of `periods` and `models`
#'     for a calendar schedule.}
#'   \item{father_matching}{A summary of the matching model for fathers of
#'     unpartnered births, as for `matching`, or `NULL` when the general
#'     matching model is used.}
#'   \item{search_weights}{A summary of the weight model choosing which
#'     woman searches: a `summary.hammel_weight_model`, a list of `periods`
#'     and `models` for a calendar schedule, or `NULL` when searching women
#'     are selected without weights.}
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
      cancel_birth_if_no_match = object$cancel_birth_if_no_match,
      cohab_as_previous_union = object$cohab_as_previous_union,
      rate_strata = summarize_rate_strata(object$rates),
      search_weights = summarize_model_argument(object$weights),
      matching = summarize_model_argument(object$matching),
      father_matching = summarize_model_argument(object$father_matching),
      cohabitation = if (is.null(object$cohab_prob)) {
        NULL
      } else {
        list(
          cohab_prob_range = range(object$cohab_prob$prob),
          cohab_rate_strata = summarize_rate_strata(object$cohab_prob),
          convert_rate_strata = summarize_rate_strata(
            object$union_convert_rate
          ),
          matching = summarize_model_argument(object$cohab_matching),
          conversion_weights = summarize_model_argument(
            object$weights_conversion
          )
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

  cli::cli_h2("Union formation process")

  list <- open_summary_list()
  print_stat("Search method", x$method)
  print_rate_strata(x$rate_strata)
  rate_mode <- identical(x$method, "rate")
  if (rate_mode) {
    print_model_argument("Search weight model", x$search_weights)
  }
  print_model_argument("Matching model", x$matching)
  if (rate_mode) {
    print_model_argument(
      "Father matching model",
      x$father_matching,
      none = "same as matching model"
    )
  }
  print_stat(
    "Births without a matched father",
    if (x$cancel_birth_if_no_match) "cancelled" else "allowed"
  )

  eligibility <- open_nested_stat("Male eligibility restrictions")
  print_stat("Minimum male age", x$min_male_age)
  print_stat("Generational proximity", x$kinship_generations)
  cli::cli_end(eligibility)

  if (is.null(x$cohabitation)) {
    print_stat("Cohabitation", "disabled")
  } else {
    cohabitation <- open_nested_stat("Cohabitation", "enabled")
    print_stat(
      "Probability range",
      paste(x$cohabitation$cohab_prob_range, collapse = "\u002d")
    )
    print_stat(
      "Counts as previous union",
      ifelse(x$cohab_as_previous_union, "yes", "no")
    )
    print_rate_strata(x$cohabitation$cohab_rate_strata)
    print_model_argument(
      "Matching model",
      x$cohabitation$matching,
      none = "same as general matching model"
    )
    if (is.null(x$cohabitation$convert_rate_strata)) {
      print_stat("Conversion to marriage", "not allowed")
    } else {
      conversion <- open_nested_stat("Conversion to marriage rates")
      print_rate_strata(x$cohabitation$convert_rate_strata)
      print_model_argument(
        "Weight model",
        x$cohabitation$conversion_weights
      )
      cli::cli_end(conversion)
    }
    cli::cli_end(cohabitation)
  }
  cli::cli_end(list)
  invisible(x)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Build a validated union-formation process for union_formation() and
# set_union_formation(); errors name `call`.
new_union_formation <- function(
  initiation,
  rates,
  matching,
  father_matching,
  cohab_matching,
  cohab_prob,
  union_convert_rate,
  weights,
  weights_conversion,
  min_male_age,
  kinship_generations,
  cancel_birth_if_no_match,
  cohab_as_previous_union,
  diagnostics,
  call = rlang::caller_env()
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
  validate_union_formation(process, call = call)
  validate_process_weights(weights_conversion, "union", call = call)
  validate_process_weights(weights, "person", call = call)
  if (!is.null(weights) && !identical(initiation, "rate")) {
    abort_hammel(
      c(
        "Search {.arg weights} need rate-based union formation.",
        "x" = "{.arg initiation} is {.val {initiation}}.",
        "i" = paste(
          "Use {.code initiation = \"rate\"} with search {.arg rates}, or",
          "remove {.arg weights}."
        )
      ),
      "specification",
      argument = "weights",
      call = call
    )
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


validate_union_formation <- function(process, call = rlang::caller_env()) {
  if (!identical(process$name, "union_formation")) {
    abort_hammel(
      "{.arg process} must be a union formation process.",
      "specification",
      call = call
    )
  }
  initiation <- process$initiation
  if (!is.character(initiation) || length(initiation) != 1L ||
        is.na(initiation) ||
        !initiation %in% c("rate", "on_demand")) {
    abort_hammel(
      c(
        "{.arg initiation} must be {.val rate} or {.val on_demand}.",
        "x" = "It is {describe_value(initiation)}."
      ),
      "argument",
      argument = "initiation",
      call = call
    )
  }
  validate_cohabitation(process, call = call)
  age <- process$min_male_age
  if (!is.numeric(age) || length(age) != 1L || !is.null(dim(age)) ||
        !is.finite(age) ||
        age < 0) {
    abort_hammel(
      c(
        "{.arg min_male_age} must be a single nonnegative number of years.",
        "x" = "It is {describe_value(age)}."
      ),
      "argument",
      argument = "min_male_age",
      call = call
    )
  }
  check_flag(process$diagnostics, "diagnostics", call = call)
  check_flag(
    process$cancel_birth_if_no_match, "cancel_birth_if_no_match",
    call = call
  )
  validate_kinship_generations(process$kinship_generations, call = call)
  validate_model_argument(
    process$matching, validate_matching_model, "matching",
    call = call
  )
  if (!is.null(process$father_matching)) {
    validate_model_argument(
      process$father_matching, validate_matching_model, "father_matching",
      call = call
    )
  }
  if (initiation == "on_demand") {
    if (!is.null(process$rates)) {
      abort_hammel(
        c(
          "On-demand union formation doesn't use search {.arg rates}.",
          "i" = paste(
            "Remove {.arg rates}, or use {.code initiation = \"rate\"} for",
            "rate-based searches."
          )
        ),
        "specification",
        argument = "rates",
        call = call
      )
    }
  } else if (is.null(process$rates)) {
    abort_hammel(
      c(
        "Rate-based union formation needs search {.arg rates}.",
        "i" = paste(
          "Supply female search rates by age, or use",
          "{.code initiation = \"on_demand\"}."
        )
      ),
      "specification",
      argument = "rates",
      call = call
    )
  } else {
    validate_union_search_rates(process$rates, call = call)
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

# Validate a union-search table, or a cohabitation-probability table whose
# probabilities are copied into `rate`. `name` and `value` name the user's
# table and value column in messages. Row checks run on the whole table, then
# schedule checks within each calendar period.
validate_union_search_rates <- function(
  rates,
  name = "rates",
  value = "rate",
  call = rlang::caller_env()
) {
  label <- if (name == "rates") "Union search" else "Cohabitation probability"
  check_stratified_columns(rates, "previous_union", label, call = call)
  duration <- check_search_rows(rates, name, value, call = call)
  validate_search_schedules(rates, name, call = call)
  if ("previous_union" %in% names(rates)) {
    effective <- expand_previous_union_rates(rates)
    uses_duration <- if ("dur_start" %in% names(effective)) {
      !is.na(effective$dur_start)
    } else {
      rep(FALSE, nrow(effective))
    }
    clocks <- vapply(
      split(uses_duration, effective$previous_union),
      function(x) length(unique(x)) != 1L,
      logical(1)
    )
    if (any(clocks)) {
      mixed <- names(clocks)[clocks]
      abort_hammel(
        c(
          paste(
            "Each {.field previous_union} category in {.arg {name}} must use",
            "one clock: age bands or duration bands."
          ),
          "x" = "{.val {mixed}} use{?s/} both.",
          "i" = "The clock must be the same in every group and period."
        ),
        "rates",
        argument = name,
        call = call
      )
    }
  }
  invisible(rates)
}

# Row checks for search and cohabitation-probability tables. Returns which
# rows use duration bands.
check_search_rows <- function(rates, name, value, call = rlang::caller_env()) {
  check_columns(rates, c("age_start", "age_end", "rate"), name, call = call)
  if (!nrow(rates)) {
    abort_hammel(
      c(
        "{.arg {name}} must contain at least one band.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = name,
      call = call
    )
  }
  calendar_rate_segments(rates, name = name, call = call)
  for (column in candidate_strata(rates)) {
    check_group_labels(
      rates[[column]], paste0(name, "$", column),
      call = call
    )
  }
  values <- paste0(name, "$", value)
  check_numeric_column(rates$rate, values, call = call)
  check_nonnegative(rates$rate, values, call = call)
  duration <- check_clock_rows(rates, name, call = call)
  previous <- if ("previous_union" %in% names(rates)) {
    rates$previous_union
  }
  if (any(duration) && is.null(previous)) {
    abort_hammel(
      c(
        paste(
          "Duration bands in {.arg {name}} need a {.field previous_union}",
          "column."
        ),
        "x" = "Duration bands at {format_rows(which(duration))}.",
        "i" = paste(
          "Durations count from the end of the previous union, so give",
          "{.val separated}, {.val widowed} or {.val dissolved}."
        )
      ),
      "rates",
      argument = name,
      rows = which(duration),
      call = call
    )
  }
  if (!is.null(previous)) {
    check_labels(
      previous, paste0(name, "$previous_union"),
      c("none", "separated", "widowed", "dissolved"),
      call = call
    )
    rows <- which(duration & previous == "none")
    if (length(rows)) {
      abort_hammel(
        c(
          "Rows with {.code previous_union = \"none\"} must use age bands.",
          "x" = "Duration bands at {format_rows(rows)}.",
          "i" = "Women with no previous union have no union duration."
        ),
        "rates",
        argument = name,
        rows = rows,
        call = call
      )
    }
  }
  invisible(duration)
}

# Each row of a search, cohabitation or conversion table uses exactly one
# clock: a complete age band or a complete duration band, with the other
# pair NA. Missing band columns count as NA. Returns which rows use
# duration bands.
check_clock_rows <- function(rates, name, call = rlang::caller_env()) {
  has_duration <- c("dur_start", "dur_end") %in% names(rates)
  if (any(has_duration) && !all(has_duration)) {
    absent <- c("dur_start", "dur_end")[!has_duration]
    abort_hammel(
      c(
        paste(
          "{.arg {name}} must have both {.field dur_start} and",
          "{.field dur_end}, or neither."
        ),
        "x" = "Missing: {.field {absent}}."
      ),
      "table",
      argument = name,
      columns = absent,
      call = call
    )
  }
  columns <- c("age_start", "age_end", "dur_start", "dur_end")
  bands <- lapply(
    stats::setNames(columns, columns),
    function(column) {
      if (column %in% names(rates)) {
        rates[[column]]
      } else {
        rep(NA_real_, nrow(rates))
      }
    }
  )
  for (column in intersect(columns, names(rates))) {
    x <- bands[[column]]
    # Bare NA is logical in R; permit all-NA placeholders for an unused axis.
    if (!(is.numeric(x) || (is.logical(x) && all(is.na(x)))) ||
          !is.null(dim(x)) ||
          any(is.nan(x))) {
      subject <- describe_column(paste0(name, "$", column))
      abort_hammel(
        c(
          "{subject} must contain numbers, or {.code NA} for an unused clock.",
          "x" = "It is {describe_value(x)}."
        ),
        "table",
        argument = paste0(name, "$", column),
        call = call
      )
    }
  }
  complete <- function(axis) {
    !is.na(bands[[paste0(axis, "_start")]]) &
      !is.na(bands[[paste0(axis, "_end")]])
  }
  empty <- function(axis) {
    is.na(bands[[paste0(axis, "_start")]]) &
      is.na(bands[[paste0(axis, "_end")]])
  }
  age <- complete("age")
  duration <- complete("dur")
  rows <- which(!((age & empty("dur")) | (duration & empty("age"))))
  if (length(rows)) {
    abort_hammel(
      c(
        paste(
          "Each row of {.arg {name}} must have one complete age or duration",
          "band."
        ),
        "x" = "Rows without exactly one: {format_rows(rows)}.",
        "i" = paste(
          "Give {.field age_start} and {.field age_end}, or {.field dur_start}",
          "and {.field dur_end}, and set the other pair to {.code NA}."
        )
      ),
      "rates",
      argument = name,
      rows = rows,
      call = call
    )
  }
  starts <- ifelse(duration, bands$dur_start, bands$age_start)
  ends <- ifelse(duration, bands$dur_end, bands$age_end)
  rows <- which(!is.finite(starts) | starts < 0 | ends <= starts)
  if (length(rows)) {
    abort_hammel(
      c(
        paste(
          "Bands in {.arg {name}} must have finite, nonnegative starts and",
          "greater ends."
        ),
        "x" = "Invalid bands at {format_rows(rows)}.",
        "i" = "Use {.code Inf} for a band with no upper limit."
      ),
      "rates",
      argument = name,
      rows = rows,
      call = call
    )
  }
  duration
}

# Schedule checks for search and cohabitation-probability tables within each
# calendar period: complete previous-union categories in each group stratum
# and bands that do not overlap.
validate_search_schedules <- function(rates, name, call = rlang::caller_env()) {
  if (validate_calendar_rates(
        rates,
        function(x) {
          validate_search_schedules(x, name, call = call)
        },
        call = call
      )) {
    return(invisible(rates))
  }
  strata_columns <- candidate_strata(rates)
  strata <- split_rate_groups(rates, strata_columns, call = call)
  stratum_text <- function(rows) {
    if (!length(strata_columns)) {
      return("")
    }
    values <- vapply(
      strata_columns,
      function(column) {
        as.character(rates[[column]][rows[1]])
      },
      character(1)
    )
    paste0(" for ", paste(strata_columns, "=", values, collapse = ", "))
  }
  previous <- if ("previous_union" %in% names(rates)) {
    rates$previous_union
  }
  if (!is.null(previous)) {
    for (rows in strata) {
      present <- unique(previous[rows])
      if (!(setequal(present, c("none", "dissolved")) ||
              setequal(present, c("none", "separated", "widowed")))) {
        stratum <- stratum_text(rows)
        subject <- describe_column(paste0(name, "$previous_union"))
        abort_hammel(
          c(
            "{subject} must use pooled or detailed categories.",
            "x" = "Categories{stratum}: {.val {present}}.",
            "i" = paste(
              "Use {.val none} and {.val dissolved}, or {.val none},",
              "{.val separated} and {.val widowed}, in every group and",
              "period, without mixing the two."
            )
          ),
          "rates",
          argument = name,
          call = call
        )
      }
    }
    strata <- unlist(
      lapply(strata, function(rows) split(rows, previous[rows])),
      recursive = FALSE
    )
  }
  for (rows in strata) {
    axis <- if ("dur_start" %in% names(rates) &&
                  !is.na(rates$dur_start[rows[1]])) {
      "dur"
    } else {
      "age"
    }
    schedule <- stratum_text(rows)
    if (!is.null(previous)) {
      schedule <- paste0(
        if (nzchar(schedule)) paste0(schedule, ",") else " for",
        " previous_union = ", previous[rows[1]]
      )
    }
    check_band_overlap(rates, rows, axis, name, schedule, call = call)
  }
  invisible(rates)
}
