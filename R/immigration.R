# Specify external immigration streams, validate them, and prepare their
# migrants and arrival records for simulation.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify external immigration process
#'
#' @description
#'
#' `immigration()` constructs an immigration process from one or more named
#' streams of migrants arriving from outside the simulated population. Each
#' stream is described by [immigration_stream()]. The convenience function
#' `set_immigration()` is equivalent to
#' `set_process(model, immigration(...))` and allows for direct piping from a
#' model.
#'
#' @details
#'
#' A model has at most one immigration process. Each call to
#' `set_immigration()` defines the complete set of streams, replacing any
#' streams set earlier. Calling `set_immigration(model)` without streams
#' removes immigration from the model.
#'
#' Each named argument provides one [immigration_stream()]. Every stream must
#' be named and streams cannot be named `"native"` or `"birthplace_ancestry"`.
#'
#' Simulations with an immigration process will automatically create a
#' `birthplace` group whose categories match the stream names, plus a
#' "native" category. This group can be used in the same way as any
#' other group (e.g. to stratify rate schedules, apply weights in weight
#' models). To track ancestry from this group, use the
#' `birthplace_ancestry = TRUE` argument of this function.
#'
#' Simulated populations of models with immigration also gain an `arrival_time`
#' column in `persons`, which is `NA` for people who did not immigrate, and a
#' `pre_migration` column in `unions` and `transitions` that flags unions and
#' transitions that occurred to migrants before entry. Additionally, simulations
#' with immigration will contain an `external_relatives` tibble that includes
#' relevant relatives of the migrating person for complete history. These
#' relatives are partners from former unions and children that did not
#' migrate with the person. These external relatives are not updated within
#' the simulation and a `censor_time` is provided for their last date of
#' observation.
#'
#' This additional returned information is provided
#' in order to be able to reconstruct full histories for migrants. Events
#' occurring prior to migration never show up in calculated rates, nor do
#' person years of exposure for migrants prior to arrival.
#'
#' @param ... Named [immigration_stream()] objects, one per stream. The
#'   argument names are the stream names.
#' @param birthplace_ancestry Whether to track ancestry in the automatic
#'   `birthplace` group. Defaults to `FALSE`.
#' @param model A `hammel_model` object, used by `set_immigration()`.
#' @return `immigration()` returns an inspectable process inheriting from
#'   `hammel_immigration` and `hammel_process`, with `name`, a named list
#'   of `streams` and `birthplace_ancestry`. `set_immigration()` returns an
#'   updated `hammel_model`.
#' @seealso [immigration_stream()], [migrant_assignment()], [set_process()]
#' @export
#' @examples
#' template <- generate_uniform_population(500)
#' model <- hammel_model() |>
#'   set_population(generate_uniform_population(1000)) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_union_formation("on_demand") |>
#'   set_immigration(
#'     other_place = immigration_stream(data.frame(count = 20), template)
#'   )
#' result <- run_simulation(model, years = 10, seed = 1)
#' table(result$persons$birthplace)
#' @order 1
immigration <- function(..., birthplace_ancestry = FALSE) {
  new_immigration(list(...), birthplace_ancestry)
}

#' @rdname immigration
#' @export
#' @order 2
set_immigration <- function(model, ..., birthplace_ancestry = FALSE) {
  attach_process(model, new_immigration(list(...), birthplace_ancestry))
}

#' Describe one stream of immigrants
#'
#' @description
#'
#' Describe a stream of migrants for [immigration()]: the template population
#' migrants are drawn from, how many arrive, who is selected and how groups
#' the template lacks are assigned.
#'
#' @details
#'
#' Migrant arrivals follow a Poisson process. Users can enter the `intensity`
#' rate for this Poisson process in one of two ways. With a `count` intensity,
#' the expected number of arrivals per year is the count and this does not
#' depend on the simulated population's size. With a `rate` intensity, it
#' is the rate times the number of people living in the simulated population
#' at each moment, including earlier immigrants, so immigration grows and
#' shrinks with the population.
#'
#' To use either approach, users provide a data.frame with either a `count`
#' or `rate` column. This data.frame can be banded by `time_start` and
#' `time_end` to allow for changing immigration intensity over time. It
#' cannot be stratified in any other way. For constant intensity, users
#' provide a tibble with one row and one column (either `count` or `rate`).
#'
#' Migrants are sampled from the provided `template` population with
#' replacement. Users can apply a [weight_model()] with the `weights`
#' argument to target individuals more likely to be migrants (e.g. working
#' age adults).
#'
#' By default, only the living members of the template's population at its end
#' time are sampled. If the provided template population has a range of time
#' that covers the current simulation, users can instead sample from living
#' persons at different cut points that correspond to the time in the
#' simulated population. To do so, users set `step` which determines how
#' frequently the template population is updated to match the simulation
#' time. For example, with `step = 5`, every five years of the simulation
#' the living population of the template available for sampling will be updated
#' to the simulation time. More frequent steps will slow the simulation.
#'
#' By default, the template's living people are divided into households as
#' for household moves in [set_group_transition()]: a woman in a union and
#' her partner, or an unpartnered adult, together with the dependents of
#' either. Each arrival draws one household, with replacement, represented
#' by its carrier: the woman in a union, or the unpartnered adult. `intensity`
#' rates are adjusted by household size to still produce correct rates for
#' individual people.
#'
#' With `household = NULL`, the stream clones individuals instead: anyone
#' drawn keeps their age, sex and groups but none of their history. Clones
#' arrive never partnered, with no children, unknown parents and group
#' durations starting at arrival. This suits templates built for the purpose,
#' such as young adults sampled as labor migrants; for templates with histories,
#' household streams are more realistic.
#'
#' If the template shares a group dimension with the simulated population, 
#' this group dimension will be carried along with migrants along with any
#' group transition history. Categories must be identically defined in shared
#' group dimensions. Template group dimensions that the simulation does
#' not declare are ignored. Each simulation group dimension that the template
#' lacks needs a [migrant_assignment()] in `assignment`.
#'
#' @param intensity A data frame with exactly one nonnegative intensity
#'   column: `count`, the expected arrivals per year, or `rate`, the expected
#'   arrivals per person-year lived in the simulated population. Without
#'   calendar columns it has one row, applying throughout. Optional
#'   `time_start` and `time_end` columns give piecewise-constant intensities
#'   over half-open calendar intervals, which must not overlap and must cover
#'   the simulation run; use a zero to pause a stream.
#' @param weights `NULL` (default) for unweighted selection of migrants, or
#'   a [weight_model()] specification. The model may use `age` (at the
#'   template's `end_time`, or at the cut with `step`), `sex` and the
#'   template's groups shared with the simulation. For household streams,
#'   weights apply to each household through its carrier, and the model may
#'   also use `partnered` (0/1), `household_size`, `n_dependents` and
#'   `partner_age` (`NA` for unpartnered carriers; in a formula the term then
#'   contributes nothing). Weights change who is selected, not how many
#'   arrive. To vary weights over time, provide a data frame with
#'   `time_start`, `time_end` and a list column `model` of weight models
#'   covering the run, as for other processes.
#' @param template A `hammel_population` object from which migrants are drawn.
#' @param assignment `NULL` (default) or a list of [migrant_assignment()]
#'   objects named by the simulation's group dimensions that the template
#'   lacks.
#' @param household A [household()] rule (default) to draw whole households
#'   with their histories, or `NULL` to clone individuals without histories.
#' @param step `NULL` (default) to draw migrants from the template's living
#'   population at its `end_time`, or a positive number of years to cut the
#'   template over time: arrivals in each `step`-year interval from the run's
#'   start are drawn from the template as it stood at the interval's start.
#'   The template must have been simulated on the same calendar over every
#'   cut.
#' @return A `hammel_immigration_stream` object.
#' @seealso [immigration()], [migrant_assignment()], [weight_model()]
#' @export
#' @examples
#' template <- generate_uniform_population(500)
#' immigration_stream(data.frame(count = 20), template)
#'
#' # Immigration at 5 per 1,000 person-years lived in the population.
#' immigration_stream(data.frame(rate = 0.005), template)
#'
#' # Fewer arrivals after 2010, and younger migrants are more likely.
#' immigration_stream(
#'   data.frame(
#'     time_start = c(2000, 2010),
#'     time_end = c(2010, Inf),
#'     count = c(40, 20)
#'   ),
#'   template,
#'   weights = weight_model(~age, c(age = log(0.97)))
#' )
immigration_stream <- function(
  intensity,
  template,
  weights = NULL,
  assignment = NULL,
  household = hammel::household(),
  step = NULL
) {
  new_immigration_stream(
    intensity, template, weights, assignment, household, step
  )
}

#' Assign migrants to a group category
#'
#' @description
#'
#' Describe how migrants in an [immigration_stream()] are assigned a category
#' of a simulation group dimension that their template population lacks.
#'
#' @details
#'
#' Three methods are available:
#'
#' - `"fixed"` assigns every migrant the category `fixed_value`.
#' - `"probability"` draws each migrant's category from `probs`, a numeric
#'   vector of probabilities named by category. Unlisted categories have
#'   probability zero, and the probabilities must sum to one.
#' - `"custom"` calls `fun(members)` at each arrival. `members` is a tibble
#'   with one row per arriving person, holding `person_id` (the template
#'   ID), `age` (at arrival), `sex`, the template's groups shared with the
#'   simulation and `arrival_time`. For household streams it also holds
#'   `role` (`"carrier"`, `"partner"` or `"dependent"`), `household` (a
#'   number shared by the household's members) and `partnered`
#'   (0/1), `household_size`, `n_dependents` and `partner_age` (`NA` for
#'   unpartnered carriers). Rows come in household order, carrier first.
#'   `fun` returns one category per row, so members
#'   of a household can receive different categories. Random choices are
#'   made inside `fun`, for example with [sample()]; to keep a household
#'   together, draw once and repeat the result.
#'
#' The `"fixed"` and `"probability"` methods assign every member of a household
#' the same category, drawn once per household. A `"custom"` function is called
#' for every arriving household, so it is slower than the other methods, and is
#' not needed for groups the template already declares, which members keep.
#'
#' @param method One of `"fixed"`, `"probability"` or `"custom"`.
#' @param fixed_value For `"fixed"`, one category label.
#' @param probs For `"probability"`, probabilities named by category. Must sum 
#'   to one.
#' @param fun For `"custom"`, a function of `members` returning one category
#'   per member, as described in details.
#' @return A `hammel_migrant_assignment` object.
#' @seealso [immigration_stream()], [group_inheritance()]
#' @export
#' @examples
#' migrant_assignment("fixed", fixed_value = "north")
#' migrant_assignment("probability", probs = c(north = 0.6, south = 0.4))
#' # Generation by age at arrival, person by person.
#' migrant_assignment(
#'   "custom",
#'   fun = function(members) {
#'     ifelse(members$age < 13, "1.5", "1")
#'   }
#' )
#' # One random region per household, likelier north for young carriers.
#' migrant_assignment(
#'   "custom",
#'   fun = function(members) {
#'     young <- members$age[1] < 30
#'     region <- sample(
#'       c("north", "south"), 1,
#'       prob = if (young) c(0.8, 0.2) else c(0.5, 0.5)
#'     )
#'     rep(region, nrow(members))
#'   }
#' )
migrant_assignment <- function(
  method,
  fixed_value = NULL,
  probs = NULL,
  fun = NULL
) {
  model <- structure(
    list(
      method = method,
      fixed_value = fixed_value,
      probs = probs,
      fun = fun
    ),
    class = "hammel_migrant_assignment"
  )
  validate_migrant_assignment(model)
  model
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize an immigration process
#'
#' @description
#'
#' Summary method for class `"hammel_immigration"`
#'
#' @param object A `hammel_immigration` object returned by [immigration()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_immigration`
#' containing `birthplace_ancestry` (whether ancestry is tracked in the
#' `birthplace` group) and `streams`, a named list with one element per
#' stream holding
#' `template_size` (persons in the template), `eligible` (households, or for
#' cloning streams template people, that can be drawn), `household` (the
#' household rule's `age_independence`, or `NULL` for cloning streams),
#' `step` (the years between template cuts, or `NULL`),
#' `intensity` (a tibble with a `count` or `rate` column),
#' `weights` (a weight model summary, or `NULL` for uniform draws) and
#' `assignment` (the method assigning each group dimension, named by
#' dimension).
#'
#' @seealso [immigration()]
#' @examples
#' process <- immigration(
#'   origin = immigration_stream(
#'     data.frame(count = 10), generate_uniform_population(100)
#'   )
#' )
#' summary(process)
#' @export
summary.hammel_immigration <- function(object, ...) {
  streams <- lapply(
    object$streams,
    function(stream) {
      intensity <- tibble::as_tibble(stream$intensity)
      parties <- migrant_parties(stream$template, stream$household)
      list(
        template_size = nrow(stream$template$persons),
        eligible = length(parties$party_start) - 1L,
        household = stream$household$age_independence,
        step = stream$step,
        intensity = intensity,
        weights = summarize_model_argument(stream$weights),
        assignment = vapply(stream$assignment, `[[`, "", "method")
      )
    }
  )
  structure(
    list(
      streams = streams,
      birthplace_ancestry = isTRUE(object$birthplace_ancestry)
    ),
    class = "summary.hammel_immigration"
  )
}

#' @rdname summary.hammel_immigration
#' @param x A `summary.hammel_immigration` object returned by
#'   `summary.hammel_immigration()`.
#' @export
print.summary.hammel_immigration <- function(x, ...) {
  cli::cli_h2("Immigration process")
  list <- open_summary_list()
  print_stat("Streams", length(x$streams))
  print_stat(
    "Birthplace ancestry tracked",
    if (isTRUE(x$birthplace_ancestry)) "yes" else "no"
  )
  for (name in names(x$streams)) {
    stream <- x$streams[[name]]
    item <- open_nested_stat("Stream", name)
    print_stat(
      "Template",
      paste0(
        stream$template_size, " persons, ", stream$eligible,
        if (is.null(stream$household)) " eligible" else " households"
      )
    )
    if (!is.null(stream$household)) {
      print_stat(
        "Households",
        paste("dependents under", stream$household, "years")
      )
    }
    if (!is.null(stream$step)) {
      print_stat("Template cuts", paste("every", stream$step, "years"))
    }
    intensity <- intensity_table(stream$intensity)
    unit <- if (intensity$kind[1] == "rate") {
      " per person-year"
    } else {
      " arrivals per year"
    }
    print_stat(
      "Intensity",
      if ("time_start" %in% names(stream$intensity)) {
        paste0(
          intensity$intensity, unit, " in ",
          format_band(intensity$time_start, intensity$time_end),
          collapse = "; "
        )
      } else {
        paste0(intensity$intensity, unit)
      }
    )
    print_model_argument("Weight model", stream$weights, none = "uniform")
    print_stat(
      "Assigned groups",
      if (length(stream$assignment)) {
        paste0(
          names(stream$assignment), " (", stream$assignment, ")",
          collapse = ", "
        )
      } else {
        "none"
      }
    )
    cli::cli_end(item)
  }
  cli::cli_end(list)
  invisible(x)
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Columns that record arrivals in `persons` (D185); the stream is the
# arrival's `birthplace` (D190).
arrival_columns <- "arrival_time"

# The automatic group dimension of models with immigration (D190), its
# category for people born in the population, and names streams can't take.
birthplace_dimension <- "birthplace"
native_birthplace <- "native"
reserved_stream_names <- c(native_birthplace, "birthplace_ancestry")

assignment_methods <- c("fixed", "probability", "custom")

# Value columns of intensity tables (D194).
intensity_kinds <- c("count", "rate")

# Carrier variables of household streams, for weights and assignment (D189).
household_variables <- c(
  "partnered", "household_size", "n_dependents", "partner_age"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Build a validated immigration process from a list of streams; errors name
# `call`.
new_immigration <- function(
  streams,
  birthplace_ancestry = FALSE,
  call = rlang::caller_env()
) {
  process <- structure(
    list(
      name = "immigration",
      streams = streams,
      birthplace_ancestry = birthplace_ancestry
    ),
    class = c("hammel_immigration", "hammel_process")
  )
  validate_immigration(process, call = call)
  process
}

validate_immigration <- function(process, call = rlang::caller_env()) {
  check_flag(
    process$birthplace_ancestry, "birthplace_ancestry",
    call = call
  )
  streams <- process$streams
  if (!is.list(streams) || inherits(streams, "hammel_immigration_stream")) {
    abort_hammel(
      c(
        "Immigration streams must be a list of immigration streams.",
        "x" = "It is {describe_value(streams)}."
      ),
      "specification",
      argument = "...",
      call = call
    )
  }
  if (!length(streams)) {
    return(invisible(process))
  }
  stream_names <- names(streams)
  if (is.null(stream_names)) {
    stream_names <- rep("", length(streams))
  }
  unnamed <- which(is.na(stream_names) | !nzchar(trimws(stream_names)))
  if (length(unnamed)) {
    abort_hammel(
      c(
        "Every immigration stream must be named.",
        "x" = "Unnamed {cli::qty(length(unnamed))}stream{?s}: {unnamed}.",
        "i" = "Name each stream, as in {.code origin = immigration_stream()}."
      ),
      "specification",
      argument = "...",
      call = call
    )
  }
  reserved <- intersect(stream_names, reserved_stream_names)
  if (length(reserved)) {
    abort_hammel(
      c(
        "Immigration streams can't be named {.val {reserved}}.",
        "i" = paste(
          "Stream names are categories of the {.field birthplace} group,",
          "where {.val native} means born in the population."
        )
      ),
      "specification",
      argument = "...",
      call = call
    )
  }
  repeated <- unique(stream_names[duplicated(stream_names)])
  if (length(repeated)) {
    abort_hammel(
      c(
        "Immigration stream names must be unique.",
        "x" = "Repeated: {.val {repeated}}."
      ),
      "specification",
      argument = "...",
      call = call
    )
  }
  for (name in stream_names) {
    if (!inherits(streams[[name]], "hammel_immigration_stream")) {
      abort_hammel(
        c(
          "Stream {.val {name}} must be an immigration stream.",
          "x" = "It is {describe_value(streams[[name]])}.",
          "i" = "Create streams with {.fn immigration_stream}."
        ),
        "specification",
        argument = name,
        call = call
      )
    }
    with_error_context(
      validate_immigration_stream(streams[[name]], call = call),
      paste0("In immigration stream \"", name, "\".")
    )
  }
  invisible(process)
}

# Build a validated stream for immigration_stream(); errors name `call`.
new_immigration_stream <- function(
  intensity,
  template,
  weights,
  assignment,
  household,
  step = NULL,
  call = rlang::caller_env()
) {
  if (missing(intensity)) {
    abort_hammel(
      c(
        "{.arg intensity} is required.",
        "i" = "Give arrivals per year or per person-year in a table."
      ),
      "argument",
      argument = "intensity",
      call = call
    )
  }
  if (missing(template)) {
    abort_hammel(
      c(
        "{.arg template} is required.",
        "i" = "Give the population migrants are drawn from."
      ),
      "argument",
      argument = "template",
      call = call
    )
  }
  stream <- structure(
    list(
      template = template,
      intensity = intensity,
      weights = normalize_model_argument(weights),
      assignment = assignment,
      household = household,
      step = step
    ),
    class = "hammel_immigration_stream"
  )
  validate_immigration_stream(stream, call = call)
  histories <- c(
    unions = NROW(template$unions), transitions = NROW(template$transitions)
  )
  if (is.null(household) && any(histories > 0)) {
    held <- names(histories)[histories > 0]
    warn_hammel(
      c(
        "Migrants won't bring the template's {.field {held}}.",
        "i" = paste(
          "A stream without a household rule clones people without their",
          "histories; use {.code household = household()} to import them."
        )
      ),
      "no_effect",
      argument = "household",
      call = call
    )
  }
  stream
}

validate_immigration_stream <- function(stream, call = rlang::caller_env()) {
  check_population(stream$template, "template", call = call)
  check_intensity(stream$intensity, call = call)
  weights <- stream$weights
  if (!is.null(weights) && !is_model_schedule(weights) &&
        !inherits(weights, "hammel_weight_model")) {
    abort_hammel(
      c(
        "{.arg weights} must be {.code NULL} or a weight model.",
        "x" = "It is {describe_value(weights)}.",
        "i" = "Create one with {.fn weight_model}."
      ),
      "specification",
      argument = "weights",
      call = call
    )
  }
  validate_process_weights(weights, "migrant", call = call)
  step <- stream$step
  if (!is.null(step) &&
        (!is.numeric(step) || length(step) != 1L || !is.null(dim(step)) ||
           !is.finite(step) || step <= 0)) {
    abort_hammel(
      c(
        paste(
          "{.arg step} must be {.code NULL} or a single positive number of",
          "years."
        ),
        "x" = "It is {describe_value(step)}."
      ),
      "argument",
      argument = "step",
      call = call
    )
  }
  household <- stream$household
  used <- unlist(lapply(model_entries(weights), function(model) {
    if (!is.null(model$formula)) all.vars(model$formula)
  }))
  clone_only <- intersect(used, household_variables)
  if (is.null(household) && length(clone_only)) {
    abort_hammel(
      c(
        "Weights for a cloning stream can't use {.field {clone_only}}.",
        "i" = paste(
          "Household variables describe a household's carrier; use",
          "{.code household = household()} to draw households."
        )
      ),
      "specification",
      argument = "weights",
      variables = clone_only,
      call = call
    )
  }
  if (!is.null(household) && !inherits(household, "hammel_household")) {
    abort_hammel(
      c(
        "{.arg household} must be {.code NULL} or a household rule.",
        "x" = "It is {describe_value(household)}.",
        "i" = "Create one with {.fn household}."
      ),
      "specification",
      argument = "household",
      call = call
    )
  }
  assignment <- stream$assignment
  if (is.null(assignment)) {
    return(invisible(stream))
  }
  dims <- names(assignment)
  if (!is.list(assignment) ||
        inherits(assignment, "hammel_migrant_assignment") ||
        !length(assignment) || is.null(dims) || anyNA(dims) ||
        any(!nzchar(dims)) || anyDuplicated(dims) ||
        !all(
          vapply(
            assignment, inherits, logical(1), "hammel_migrant_assignment"
          )
        )) {
    abort_hammel(
      c(
        paste(
          "{.arg assignment} must be a list of migrant assignments named by",
          "group dimension."
        ),
        "x" = "It is {describe_value(assignment)}.",
        "i" = paste(
          "For example, {.code list(region = migrant_assignment(\"fixed\",",
          "fixed_value = \"north\"))}."
        )
      ),
      "specification",
      argument = "assignment",
      call = call
    )
  }
  for (rule in assignment) {
    validate_migrant_assignment(rule, call = call)
  }
  invisible(stream)
}

# A stream's intensity (D194): a table with exactly one of a `count` column
# (arrivals per year) or a `rate` column (arrivals per destination
# person-year), and optional calendar columns that do not overlap.
check_intensity <- function(intensity, call = rlang::caller_env()) {
  if (!is.data.frame(intensity)) {
    abort_hammel(
      c(
        paste(
          "{.arg intensity} must be a table with a {.field count} or a",
          "{.field rate} column."
        ),
        "x" = "It is {describe_value(intensity)}.",
        "i" = paste(
          "For a constant flow of 200 arrivals a year, use",
          "{.code tibble(count = 200)}."
        )
      ),
      "argument",
      argument = "intensity",
      call = call
    )
  }
  kinds <- intersect(intensity_kinds, names(intensity))
  if (length(kinds) != 1L) {
    abort_hammel(
      c(
        paste(
          "{.arg intensity} must have exactly one of {.field count} or",
          "{.field rate}."
        ),
        "x" = if (length(kinds)) "It has both." else "It has neither.",
        "i" = paste(
          "Use {.field count} for arrivals per year, or {.field rate} for",
          "arrivals per person-year lived in the population."
        )
      ),
      "rates",
      argument = "intensity",
      call = call
    )
  }
  calendar <- intersect(c("time_start", "time_end"), names(intensity))
  extra <- setdiff(names(intensity), c(calendar, kinds))
  if (length(extra)) {
    abort_hammel(
      c(
        "{.arg intensity} can only have calendar and intensity columns.",
        "x" = "Unexpected: {.field {extra}}.",
        "i" = paste(
          "Arrivals can't be stratified. Use weights and the template to",
          "decide who arrives."
        )
      ),
      "rates",
      argument = "intensity",
      columns = extra,
      call = call
    )
  }
  if (!nrow(intensity)) {
    abort_hammel(
      c(
        "{.arg intensity} must contain at least one row.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = "intensity",
      call = call
    )
  }
  column <- paste0("intensity$", kinds)
  check_numeric_column(intensity[[kinds]], column, call = call)
  check_nonnegative(intensity[[kinds]], column, call = call)
  if (!length(calendar)) {
    if (nrow(intensity) != 1L) {
      abort_hammel(
        c(
          "An {.arg intensity} without calendar columns must have one row.",
          "x" = "It has {nrow(intensity)} rows.",
          "i" = "Add {.field time_start} and {.field time_end} to vary it."
        ),
        "rates",
        argument = "intensity",
        call = call
      )
    }
    return(invisible(intensity))
  }
  check_bands(intensity, "time", "intensity", minimum = -Inf, call = call)
  rows <- order(intensity$time_start)
  starts <- intensity$time_start[rows]
  ends <- intensity$time_end[rows]
  overlap <- which(utils::tail(starts, -1) < utils::head(ends, -1))
  if (length(overlap)) {
    bad <- sort(rows[c(overlap, overlap + 1L)])
    abort_hammel(
      c(
        "Calendar intervals in {.arg intensity} must not overlap.",
        "x" = "Overlapping intervals at {format_rows(unique(bad))}."
      ),
      "rates",
      argument = "intensity",
      rows = unique(bad),
      call = call
    )
  }
  invisible(intensity)
}

validate_migrant_assignment <- function(model, call = rlang::caller_env()) {
  method <- model$method
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% assignment_methods) {
    abort_hammel(
      c(
        "{.arg method} must be one of {.or {.val {assignment_methods}}}.",
        "x" = "It is {describe_value(method)}."
      ),
      "argument",
      argument = "method",
      call = call
    )
  }
  arguments <- c(fixed = "fixed_value", probability = "probs", custom = "fun")
  unused <- setdiff(arguments, arguments[[method]])
  unused <- unused[!vapply(model[unused], is.null, logical(1))]
  if (length(unused)) {
    abort_hammel(
      c(
        "The {.val {method}} method doesn't use {.arg {unused}}.",
        "i" = "Supply only {.arg {arguments[[method]]}}."
      ),
      "argument",
      argument = unused[1],
      call = call
    )
  }
  value <- model[[arguments[[method]]]]
  if (is.null(value)) {
    abort_hammel(
      "The {.val {method}} method needs {.arg {arguments[[method]]}}.",
      "argument",
      argument = arguments[[method]],
      call = call
    )
  }
  if (method == "fixed" &&
        ((!is.character(value) && !is.numeric(value)) ||
           length(value) != 1L || is.na(value))) {
    abort_hammel(
      c(
        "{.arg fixed_value} must be one category label.",
        "x" = "It is {describe_value(value)}."
      ),
      "argument",
      argument = "fixed_value",
      call = call
    )
  }
  if (method == "probability") {
    labels <- names(value)
    if (!is.numeric(value) || !length(value) || !is.null(dim(value)) ||
          any(!is.finite(value)) || any(value < 0) || is.null(labels) ||
          anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
      abort_hammel(
        c(
          paste(
            "{.arg probs} must be nonnegative probabilities named by",
            "category."
          ),
          "x" = "It is {describe_value(value)}."
        ),
        "argument",
        argument = "probs",
        call = call
      )
    }
    if (abs(sum(value) - 1) > 1e-8) {
      total <- sum(value)
      abort_hammel(
        c(
          "{.arg probs} must sum to one.",
          "x" = "They sum to {total}."
        ),
        "argument",
        argument = "probs",
        call = call
      )
    }
  }
  if (method == "custom" && !is.function(value)) {
    abort_hammel(
      c(
        "{.arg fun} must be a function.",
        "x" = "It is {describe_value(value)}."
      ),
      "argument",
      argument = "fun",
      call = call
    )
  }
  invisible(model)
}

# A stream's parties (D187). Without a household rule, each living person
# is a one-person party to clone (D188). With one, the template's living
# persons form households under D140's rule at its end_time: a carrier (the
# woman in a current union, or an unpartnered person who is not a
# dependent), their partner and the dependents of either, recursively.
# Returns the members' template rows, party by party (carrier, partner, then
# dependents by birth), the zero-based start of each party and its end, and
# links within parties as member offsets (-1 none): parents, and each
# woman's current partner with the union's template row.
migrant_parties <- function(template, household) {
  persons <- template$persons
  n <- nrow(persons)
  if (is.null(household)) {
    rows <- which(is.na(persons$death_time))
    none <- rep(-1L, length(rows))
    return(
      list(
        rows = rows,
        party_start = c(0L, seq_along(rows)),
        mother = none,
        father = none,
        partner = none,
        union_row = rep(NA_integer_, length(rows))
      )
    )
  }
  alive <- is.na(persons$death_time)
  unions <- template$unions
  current <- which(is.na(unions$end_time))
  wife <- match(unions$wife_id[current], persons$person_id)
  husband <- match(unions$husband_id[current], persons$person_id)
  partner <- rep(NA_integer_, n)
  partner[wife] <- husband
  partner[husband] <- wife
  union_row <- rep(NA_integer_, n)
  union_row[wife] <- current
  # A dependent attaches to their living mother, otherwise their living
  # father; parents who are external relatives are absent.
  mother <- match(persons$mother_id, persons$person_id)
  father <- match(persons$father_id, persons$person_id)
  living_mother <- !is.na(mother) & alive[mother]
  living_father <- !is.na(father) & alive[father]
  attached <- ifelse(
    living_mother, mother, ifelse(living_father, father, NA_integer_)
  )
  dependent <- alive & is.na(partner) & !is.na(attached) &
    persons$birth_time + household$age_independence > template$end_time
  head <- seq_len(n)
  repeat {
    up <- dependent[head]
    if (!any(up)) {
      break
    }
    head[up] <- attached[head[up]]
  }
  male_head <- !is.na(partner[head]) & persons$sex[head] == "m"
  carrier <- ifelse(male_head, partner[head], head)
  members <- which(alive)
  # Carrier first, then the partner, then dependents by birth, so parents
  # precede their children.
  rank <- ifelse(
    members == carrier[members], 0L, ifelse(dependent[members], 2L, 1L)
  )
  members <- members[order(
    carrier[members], rank, persons$birth_time[members], members
  )]
  party <- match(carrier[members], unique(carrier[members]))
  first <- match(seq_len(max(c(0L, party))), party)
  offset <- function(target) {
    at <- match(target, members)
    same <- !is.na(at) & party[at] == party[seq_along(at)]
    ifelse(same, at - first[party], -1L)
  }
  list(
    rows = members,
    party_start = as.integer(c(first, length(members) + 1L) - 1L),
    mother = as.integer(offset(mother[members])),
    father = as.integer(offset(father[members])),
    partner = as.integer(
      ifelse(persons$sex[members] == "f", offset(partner[members]), -1L)
    ),
    union_row = as.integer(union_row[members])
  )
}

# The declaration and population of a run with the automatic birthplace
# dimension (D190), present when the model has immigration or the population
# already has it. Categories are "native", the streams, and any earlier
# categories; newborns are always native. People without a recorded
# birthplace are native. `user_groups` is the model's set_groups()
# declaration, which can't declare birthplace alongside immigration.
with_birthplace <- function(
  population,
  declaration,
  process,
  user_groups = NULL,
  call = rlang::caller_env()
) {
  name <- birthplace_dimension
  existing <- population$groups[[name]]
  if (!isTRUE(existing$automatic)) {
    existing <- NULL
  }
  if (is.null(process) && is.null(existing)) {
    return(list(population = population, declaration = declaration))
  }
  if (!is.null(process) && name %in% names(user_groups)) {
    abort_hammel(
      c(
        "Models with immigration can't declare a {.field birthplace} group.",
        "i" = paste(
          "{.field birthplace} is created automatically: {.val native} for",
          "people born in the population, and the stream name for immigrants."
        )
      ),
      "groups",
      group = name,
      call = call
    )
  }
  persons <- population$persons
  if (!name %in% names(persons)) {
    persons[[name]] <- rep(native_birthplace, nrow(persons))
  }
  values <- as.character(persons[[name]])
  categories <- unique(
    c(
      native_birthplace, as.character(existing$categories),
      names(process$streams), values[!is.na(values)]
    )
  )
  persons[[name]] <- values
  population$persons <- persons
  ancestry <- if (is.null(process)) {
    isTRUE(existing$ancestry)
  } else {
    isTRUE(process$birthplace_ancestry)
  }
  entry <- group(
    categories = categories,
    ancestry = ancestry,
    inheritance = group_inheritance("fixed", fixed_value = native_birthplace)
  )
  entry$automatic <- TRUE
  declaration[[name]] <- entry
  list(population = population, declaration = declaration)
}

# A stream's intensity as calendar intervals and its kind ("count" or
# "rate"); a table without calendar columns applies throughout.
intensity_table <- function(intensity) {
  kind <- intersect(intensity_kinds, names(intensity))
  if (!"time_start" %in% names(intensity)) {
    return(
      data.frame(
        time_start = -Inf, time_end = Inf,
        intensity = as.numeric(intensity[[kind]]),
        kind = kind
      )
    )
  }
  rows <- order(intensity$time_start)
  data.frame(
    time_start = as.numeric(intensity$time_start[rows]),
    time_end = as.numeric(intensity$time_end[rows]),
    intensity = as.numeric(intensity[[kind]][rows]),
    kind = kind
  )
}

# Give each ancestry table of a destination every category of its dimension,
# so migrants' shares in categories nobody yet has can be recorded.
complete_ancestry_columns <- function(population, registries) {
  for (d in names(population$ancestry)) {
    table <- population$ancestry[[d]]
    categories <- as.character(registries[[d]])
    known <- setdiff(names(table), c("person_id", "unknown"))
    for (category in setdiff(categories, known)) {
      table[[category]] <- 0
    }
    known <- c(intersect(categories, known), setdiff(known, categories))
    known <- c(known, setdiff(categories, known))
    population$ancestry[[d]] <- table[c("person_id", known, "unknown")]
  }
  population
}

# Engine input for the immigration process (D185): for each stream, its
# intensity intervals and, for every eligible template person, age, sex,
# destination categories (-1 where assigned at arrival), ancestry shares
# (NA in blocks filled at arrival), cumulative selection probabilities and
# assignment probabilities. `components` are each destination dimension's
# ancestry columns and `offsets` their first engine column.
prepare_immigration <- function(
  process,
  declaration,
  registries,
  components,
  offsets,
  start,
  end,
  counts_cohab = FALSE,
  cohab_enabled = FALSE,
  call = rlang::caller_env()
) {
  dimensions <- names(declaration)
  streams <- lapply(
    names(process$streams),
    function(name) {
      with_error_context(
        prepare_stream(
          process$streams[[name]], name, declaration, registries, components,
          start, end, counts_cohab, cohab_enabled,
          call = call
        ),
        paste0("In immigration stream \"", name, "\".")
      )
    }
  )
  # The engine column of each category's ancestry share, by dimension.
  category_component <- lapply(
    seq_along(dimensions),
    function(k) {
      if (!length(components[[k]])) {
        return(integer())
      }
      as.integer(
        offsets[k] +
          match(as.character(registries[[k]]), components[[k]]) - 1L
      )
    }
  )
  oldest <- max(
    0, unlist(lapply(streams, function(x) {
      vapply(x$slices, function(y) max(c(0, y$age)), numeric(1))
    }))
  )
  list(
    streams = streams,
    category_component = category_component,
    min_birth = start - oldest
  )
}

prepare_stream <- function(
  stream,
  name,
  declaration,
  registries,
  components,
  start,
  end,
  counts_cohab,
  cohab_enabled,
  call = rlang::caller_env()
) {
  intensity <- intensity_table(stream$intensity)
  if ("time_start" %in% names(stream$intensity)) {
    check_time_coverage(
      intensity$time_start, intensity$time_end, start, end,
      "The stream's intensity",
      call = call
    )
  }
  # Template cuts (D196): one truncated template per cut, used from its cut
  # time; without cuts, the template at its end_time.
  cuts <- template_cuts(stream, start, end, call = call)
  templates <- if (is.null(stream$step)) {
    list(stream$template)
  } else {
    lapply(cuts, function(cut) population_at(stream$template, cut))
  }
  # Slices start at every cut and weight-schedule boundary in the run
  # (D197); each uses its latest cut's template and its active weights.
  weights <- stream$weights
  check_model_coverage(weights, start, end, "stream's weight", call = call)
  edges <- if (is_model_schedule(weights)) {
    c(weights$time_start, weights$time_end)
  }
  starts <- sort(unique(c(cuts, edges[edges > start & edges < end])))
  template_of <- findInterval(starts, cuts)
  slices <- lapply(
    seq_along(starts),
    function(k) {
      period <- NULL
      if (is_model_schedule(weights)) {
        row <- which(
          weights$time_start <= starts[k] & starts[k] < weights$time_end
        )
        period <- c(weights$time_start[row], weights$time_end[row])
      }
      prepare_slice(
        stream, templates[[template_of[k]]], name, declaration, registries,
        components, counts_cohab, cohab_enabled,
        weights = model_at(weights, starts[k]),
        cut = if (is.null(stream$step)) NULL else cuts[template_of[k]],
        period = period,
        call = call
      )
    }
  )
  list(
    time_start = intensity$time_start,
    time_end = intensity$time_end,
    intensity = intensity$intensity,
    per_person = identical(intensity$kind[1], "rate"),
    cuts = as.numeric(starts),
    slices = slices,
    templates = templates[template_of]
  )
}

# The cut times of a stream (D196): the run's start and every `step` years
# after it, each within the template's simulated history; without `step`,
# the run's start alone.
template_cuts <- function(stream, start, end, call = rlang::caller_env()) {
  step <- stream$step
  if (is.null(step)) {
    return(start)
  }
  cuts <- start + step * (seq_len(max(1, ceiling((end - start) / step))) - 1)
  template <- stream$template
  first <- template$unrealized_since
  if (is.null(first)) {
    abort_hammel(
      c(
        "Template cuts need a template with a simulated history.",
        "i" = paste(
          "Use a template from {.fn run_simulation}, or {.code step = NULL}",
          "to draw from the template's end."
        )
      ),
      "population",
      argument = "step",
      call = call
    )
  }
  outside <- cuts[cuts < first | cuts > template$end_time]
  if (length(outside)) {
    abort_hammel(
      c(
        "The template's history must cover every cut.",
        "x" = paste(
          "Cuts at {outside} fall outside its simulated period, {first} to",
          "{template$end_time}."
        ),
        "i" = paste(
          "Cuts run from the run's start every {.arg step} years, on the",
          "template's calendar."
        )
      ),
      "coverage",
      argument = "step",
      call = call
    )
  }
  cuts
}

# One template's engine input for a stream: its parties, selection and
# assignment probabilities, expected party size and imported state.
prepare_slice <- function(
  stream,
  template,
  name,
  declaration,
  registries,
  components,
  counts_cohab,
  cohab_enabled,
  weights = stream$weights,
  cut = NULL,
  period = NULL,
  call = rlang::caller_env()
) {
  rules <- stream$assignment
  if (birthplace_dimension %in% names(rules)) {
    abort_hammel(
      c(
        "Assignments can't set {.field birthplace}.",
        "i" = "Each migrant's birthplace is their stream's name."
      ),
      "specification",
      argument = "assignment",
      group = birthplace_dimension,
      call = call
    )
  }
  # A migrant's birthplace is their stream (D190); the template's is ignored.
  if (birthplace_dimension %in% names(declaration)) {
    rules[[birthplace_dimension]] <- migrant_assignment(
      "fixed",
      fixed_value = name
    )
  }
  parties <- migrant_parties(template, stream$household)
  rows <- parties$rows
  if (!length(rows)) {
    abort_hammel(
      c(
        "The stream's template has nobody who can migrate.",
        "i" = "Migrants are drawn from people alive at the template's end."

      ),
      "population",
      argument = "template",
      call = call
    )
  }
  persons <- template$persons[rows, , drop = FALSE]
  dimensions <- names(declaration)
  template_dims <- setdiff(names(template$groups), birthplace_dimension)
  shared <- intersect(dimensions, template_dims)
  assigned <- setdiff(dimensions, template_dims)
  template_registries <- if (length(template_dims)) {
    group_registries(template)
  } else {
    list()
  }
  for (d in shared) {
    ours <- as.character(registries[[d]])
    theirs <- as.character(template_registries[[d]])
    if (!setequal(ours, theirs)) {
      abort_hammel(
        c(
          paste(
            "Group {.field {d}} must have the same categories in the",
            "template and the simulation."
          ),
          "x" = paste(
            "The simulation has {.val {ours}}; the template has",
            "{.val {theirs}}."
          ),
          "i" = paste(
            "Declare the same categories for both with",
            "{.code group(categories = ...)}."
          )
        ),
        "groups",
        group = d,
        call = call
      )
    }
  }
  overriding <- intersect(names(rules), template_dims)
  if (length(overriding)) {
    abort_hammel(
      c(
        "The template already has {.field {overriding}}.",
        "i" = paste(
          "Assignments are only for simulation groups the template lacks;",
          "migrants keep their template categories."
        )
      ),
      "specification",
      argument = "assignment",
      group = overriding,
      call = call
    )
  }
  unknown <- setdiff(names(rules), dimensions)
  if (length(unknown)) {
    abort_hammel(
      c(
        "The simulation doesn't declare {.field {unknown}}.",
        "i" = "Assign only declared group dimensions."
      ),
      "groups",
      argument = "assignment",
      group = unknown,
      call = call
    )
  }
  missing_rules <- setdiff(assigned, names(rules))
  if (length(missing_rules)) {
    abort_hammel(
      c(
        "The template lacks {.field {missing_rules}}.",
        "i" = paste(
          "Give each such group a {.fn migrant_assignment} in",
          "{.arg assignment}."
        )
      ),
      "specification",
      argument = "assignment",
      group = missing_rules,
      call = call
    )
  }
  # Every party member, and the carriers who represent parties in weights
  # and assignment (D187).
  people <- tibble::tibble(
    person_id = persons$person_id,
    age = template$end_time - persons$birth_time,
    sex = as.character(persons$sex)
  )
  for (d in shared) {
    people[[d]] <- persons[[d]]
  }
  starts <- parties$party_start
  carriers <- utils::head(starts, -1) + 1L
  members <- people[carriers, , drop = FALSE]
  # Household streams describe each household by its carrier (D189).
  if (!is.null(stream$household)) {
    partner <- parties$partner[carriers]
    size <- diff(starts)
    members$partnered <- as.integer(partner >= 0)
    members$household_size <- as.integer(size)
    members$n_dependents <- as.integer(size - 1L - members$partnered)
    members$partner_age <- ifelse(
      partner >= 0, people$age[carriers + pmax(partner, 0L)], NA_real_
    )
  }
  probabilities <- migrant_probabilities(
    weights, members, registries[shared], name,
    cut = cut,
    period = period,
    call = call
  )
  # Intensity is in persons per year; the clock counts parties (D187).
  expected_size <- if (is.null(stream$household)) {
    1
  } else {
    sum(probabilities * diff(starts))
  }
  n <- nrow(people)
  categories <- matrix(-1L, nrow = n, ncol = length(dimensions))
  assignment <- vector("list", length(dimensions))
  arriving <- NULL
  for (k in seq_along(dimensions)) {
    d <- dimensions[k]
    if (d %in% shared) {
      categories[, k] <- match(
        as.character(people[[d]]), as.character(registries[[d]])
      ) - 1L
    } else if (rules[[d]]$method == "custom") {
      if (is.null(arriving)) {
        arriving <- arriving_members(people, parties, members, stream)
      }
      assignment[[k]] <- custom_assignment(
        rules[[d]], arriving, as.character(registries[[d]]), d,
        call = call
      )
    } else {
      assignment[[k]] <- assignment_probabilities(
        rules[[d]], nrow(members), as.character(registries[[d]]), d,
        call = call
      )
    }
  }
  live <- template$unions[parties$union_row, , drop = FALSE]
  list(
    expected_size = expected_size,
    age = as.numeric(people$age),
    female = people$sex == "f",
    categories = categories,
    ancestry = migrant_ancestry(
      template, rows, people, dimensions, shared, registries, components
    ),
    cumulative = c(utils::head(cumsum(probabilities), -1), 1),
    assignment = assignment,
    rows = rows,
    party_start = as.integer(starts),
    mother = parties$mother,
    father = parties$father,
    partner = parties$partner,
    union_since = template$end_time - as.numeric(live$start_time),
    union_cohab = if ("union_type" %in% names(live)) {
      !is.na(live$union_type) & live$union_type == "cohab"
    } else {
      rep(FALSE, n)
    },
    union_row = parties$union_row
  ) |>
    c(
      if (is.null(stream$household)) {
        clone_state(n, length(dimensions))
      } else {
        migrant_history_state(
          template, rows, dimensions, shared, counts_cohab, cohab_enabled,
          call = call
        )
      }
    )
}

# Clones' state at arrival (D188): never partnered, parity zero, and every
# group entered at arrival.
clone_state <- function(n, dimensions) {
  list(
    parity = integer(n),
    since_birth = rep(NA_real_, n),
    previous_union = integer(n),
    legal_previous = integer(n),
    since_union_end = rep(NA_real_, n),
    since_entry = matrix(NA_real_, nrow = n, ncol = dimensions)
  )
}

# A population as it stood at `time` (D196), for template cuts: people born
# by then who had arrived, unions started by then (live if they ended
# later), transitions by then with groups at that time, external relatives
# censored at it, and deaths after it undone. Internal for now.
population_at <- function(population, time) {
  if (time >= population$end_time) {
    return(population)
  }
  persons <- population$persons
  keep <- persons$birth_time <= time
  if ("arrival_time" %in% names(persons)) {
    keep <- keep &
      (is.na(persons$arrival_time) | persons$arrival_time <= time)
  }
  external <- population[["external_relatives"]]
  has_external <- !is.null(external) && nrow(external) > 0
  if (has_external) {
    external <- external[external$birth_time <= time, , drop = FALSE]
  }
  # Groups at the time: the origin of a person's first later move.
  tx <- population$transitions
  for (d in names(population$groups)) {
    if (is.null(tx) || !nrow(tx)) {
      break
    }
    later <- tx[tx$group == d & tx$time > time, , drop = FALSE]
    later <- later[order(later$person_id, later$time), , drop = FALSE]
    first <- later[!duplicated(later$person_id), , drop = FALSE]
    at <- match(first$person_id, persons$person_id)
    persons[[d]][at[!is.na(at)]] <- first$origin[!is.na(at)]
  }
  persons <- persons[keep, , drop = FALSE]
  persons$death_time[!is.na(persons$death_time) &
                       persons$death_time > time] <- NA
  ids <- c(persons$person_id, if (has_external) external$person_id)
  unions <- population$unions
  unions <- unions[unions$start_time <= time & unions$wife_id %in% ids &
                     unions$husband_id %in% ids, , drop = FALSE]
  ongoing <- !is.na(unions$end_time) & unions$end_time > time
  unions$end_time[ongoing] <- NA
  unions$end_reason[ongoing] <- NA
  if (has_external) {
    # Keep external relatives still linked to a remaining person.
    member <- unions$wife_id %in% persons$person_id |
      unions$husband_id %in% persons$person_id
    linked <- external$person_id %in% c(
      unions$wife_id[member], unions$husband_id[member],
      persons$mother_id, persons$father_id
    ) |
      external$mother_id %in% persons$person_id |
      external$father_id %in% persons$person_id
    external <- external[linked, , drop = FALSE]
    external$death_time[!is.na(external$death_time) &
                          external$death_time > time] <- NA
    external$censor_time <- ifelse(
      is.na(external$death_time),
      pmin(external$censor_time, time),
      external$death_time
    )
    ids <- c(persons$person_id, external$person_id)
    unions <- unions[unions$wife_id %in% ids & unions$husband_id %in% ids, ,
                     drop = FALSE]
    for (field in c("mother_id", "father_id")) {
      external[[field]][!external[[field]] %in% ids] <- NA
    }
    population$external_relatives <- external
  }
  for (field in c("mother_id", "father_id")) {
    persons[[field]][!persons[[field]] %in% ids] <- NA
  }
  if (!is.null(tx)) {
    population$transitions <- tx[tx$time <= time &
                                   tx$person_id %in% persons$person_id, ,
                                 drop = FALSE]
  }
  if (!is.null(population$ancestry)) {
    population$ancestry <- lapply(
      population$ancestry,
      function(shares) {
        shares[shares$person_id %in% persons$person_id, , drop = FALSE]
      }
    )
  }
  population$persons <- persons
  population$unions <- unions
  population$end_time <- time
  population
}

# Template people: persons and any external relatives, in one table with the
# columns histories need (D186). Template arrival records are ignored.
template_people <- function(template) {
  columns <- c(
    "person_id", "sex", "birth_time", "death_time", "mother_id",
    "father_id"
  )
  people <- as.data.frame(template$persons)[columns]
  external <- template$external_relatives
  if (!is.null(external) && nrow(external)) {
    people <- rbind(people, as.data.frame(external)[columns])
  }
  people
}

# Each eligible template person's state at the template endpoint, as elapsed
# times for the engine (D186): parity and time since the last birth (women),
# previous-union category, legal category and time since the last union
# ended, and time since entering each destination group (NA for groups the
# template lacks, whose entry is the arrival).
migrant_history_state <- function(
  template,
  rows,
  dimensions,
  shared,
  counts_cohab,
  cohab_enabled,
  call = rlang::caller_env()
) {
  persons <- template$persons[rows, , drop = FALSE]
  ids <- persons$person_id
  end <- template$end_time
  people <- template_people(template)
  unions <- template$unions
  mine <- unions$wife_id %in% ids | unions$husband_id %in% ids
  if (!cohab_enabled && "union_type" %in% names(unions) &&
        any(unions$union_type[mine] == "cohab")) {
    abort_hammel(
      c(
        paste(
          "The stream's template has cohabitations, but cohabitation isn't",
          "enabled."
        ),
        "i" = paste(
          "Migrants bring their union histories. Enable cohabitation with",
          "{.arg cohab_prob} in {.fn set_union_formation}."
        )
      ),
      "specification",
      argument = "template",
      call = call
    )
  }
  female <- persons$sex == "f"
  # Parity and the latest birth of each person, by numeric matching.
  mother <- match(people$mother_id, ids)
  births <- which(!is.na(mother))
  parity <- tabulate(mother[births], nbins = length(ids))
  last <- rep(NA_real_, length(ids))
  if (length(births)) {
    o <- births[order(mother[births], -people$birth_time[births])]
    first <- o[!duplicated(mother[o])]
    last[mother[first]] <- people$birth_time[first]
  }
  parity[!female] <- 0L
  last[!female] <- NA_real_
  previous <- call_native(
    previous_union_info_cpp(people, unions, counts_cohab),
    call = call
  )
  at <- match(ids, people$person_id)
  since_entry <- matrix(
    NA_real_,
    nrow = length(ids), ncol = length(dimensions)
  )
  tx <- template$transitions
  for (k in seq_along(dimensions)) {
    d <- dimensions[k]
    if (!d %in% shared) {
      next
    }
    entered <- persons$birth_time
    if (!is.null(tx) && nrow(tx)) {
      moves <- tx[tx$group == d & tx$person_id %in% ids, , drop = FALSE]
      if (nrow(moves)) {
        moves <- moves[order(moves$person_id, -moves$time), , drop = FALSE]
        latest <- moves[!duplicated(moves$person_id), , drop = FALSE]
        entered[match(latest$person_id, ids)] <- latest$time
      }
    }
    since_entry[, k] <- end - entered
  }
  list(
    parity = as.integer(parity),
    since_birth = end - last,
    previous_union = as.integer(previous$previous_union[at]),
    legal_previous = as.integer(previous$legal_previous[at]),
    since_union_end = end - previous$end_time[at],
    since_entry = since_entry
  )
}

# Selection probabilities of eligible template people: uniform, or
# proportional to exp(eta) under the stream's weight model.
migrant_probabilities <- function(
  weights,
  members,
  categories,
  name,
  cut = NULL,
  period = NULL,
  call = rlang::caller_env()
) {
  n <- nrow(members)
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  if (!is.null(weights$fun)) {
    bin <- list(process = "immigration", stream = name)
    if (!is.null(cut)) {
      bin$cut <- cut
    }
    if (!is.null(period)) {
      bin$time_start <- period[1]
      bin$time_end <- period[2]
    }
    eta <- weights$fun(members, bin)
    if (!is.numeric(eta) || !is.null(dim(eta)) || length(eta) != n ||
          any(!is.finite(eta))) {
      abort_hammel(
        c(
          "A custom immigration weight function must return log weights.",
          "x" = "It returned {describe_value(eta)} for {n} member{?s}.",
          "i" = paste(
            "Return one finite log weight per row of {.code members}."
          )
        ),
        "runtime",
        process = "immigration",
        call = call
      )
    }
  } else {
    columns <- weight_design(
      weights, "migrant",
      list(categories = categories, ancestry = character()),
      call = call
    )
    eta <- numeric(n)
    for (column in columns) {
      value <- rep(1, n)
      for (f in column$factors) {
        x <- members[[f$var]]
        value <- value * if (is.na(f$level)) {
          as.numeric(x)
        } else {
          as.numeric(as.character(x) == f$level)
        }
      }
      # A term with a missing value, such as the partner's age of an
      # unpartnered carrier, contributes nothing.
      term <- column$coefficient * value
      term[is.na(term)] <- 0
      eta <- eta + term
    }
  }
  w <- exp(eta - max(eta))
  w / sum(w)
}

# A fixed or probability assignment's probabilities over a dimension's
# categories, the same for each of a stream's `n` parties (D185).
assignment_probabilities <- function(
  rule,
  n,
  categories,
  dimension,
  call = rlang::caller_env()
) {
  probs <- if (rule$method == "fixed") {
    stats::setNames(1, as.character(rule$fixed_value))
  } else {
    rule$probs
  }
  unknown <- setdiff(names(probs), categories)
  if (length(unknown)) {
    abort_hammel(
      c(
        "Assigned categories must belong to {.field {dimension}}.",
        "x" = "Unknown: {.val {unknown}}.",
        "i" = "Its categories are {.val {categories}}."
      ),
      "groups",
      group = dimension,
      call = call
    )
  }
  row <- numeric(length(categories))
  row[match(names(probs), categories)] <- probs
  matrix(row, nrow = n, ncol = length(categories), byrow = TRUE)
}

# Every member a stream can bring, for custom assignment (D193): their own
# template ID, age, sex and shared groups, and for household streams their
# role and their household's variables (D189).
arriving_members <- function(people, parties, members, stream) {
  table <- as.data.frame(people)
  if (is.null(stream$household)) {
    return(table)
  }
  starts <- parties$party_start
  size <- diff(starts)
  party <- rep(seq_along(size), size)
  offset <- seq_len(nrow(table)) - 1L - starts[party]
  partner <- parties$partner[starts[party] + 1L]
  table$role <- ifelse(
    offset == 0L, "carrier", ifelse(offset == partner, "partner", "dependent")
  )
  for (v in household_variables) {
    table[[v]] <- members[[v]][party]
  }
  table
}

# A custom assignment, called by the engine at each arrival (D193) with the
# arriving members' rows, the arrival time and a household number. The
# user's function returns one category per member; the result is their
# zero-based category codes.
custom_assignment <- function(
  rule,
  arriving,
  categories,
  dimension,
  call = rlang::caller_env()
) {
  # Evaluate now: the caller's loop variables change before arrivals.
  force(rule)
  force(categories)
  force(dimension)
  force(call)
  household <- "role" %in% names(arriving)
  function(rows, time, party) {
    members <- arriving[rows, , drop = FALSE]
    rownames(members) <- NULL
    if (household) {
      members$household <- rep(party, nrow(members))
    }
    members$arrival_time <- rep(time, nrow(members))
    members <- tibble::as_tibble(members)
    result <- rule$fun(members)
    if (is.factor(result)) {
      result <- as.character(result)
    }
    if (!(is.character(result) || is.numeric(result)) ||
          !is.null(dim(result)) || length(result) != nrow(members) ||
          anyNA(result)) {
      abort_hammel(
        c(
          paste(
            "A custom assignment function must return one category per",
            "member."
          ),
          "x" = paste(
            "It returned {describe_value(result)} for {nrow(members)}",
            "member{?s}."
          )
        ),
        "runtime",
        group = dimension,
        call = call
      )
    }
    codes <- match(as.character(result), categories) - 1L
    if (anyNA(codes)) {
      unknown <- unique(as.character(result)[is.na(codes)])
      abort_hammel(
        c(
          "A custom assignment function returned an undeclared category.",
          "x" = "It returned {.val {unknown}}.",
          "i" = paste(
            "The categories of {.field {dimension}} are",
            "{.val {categories}}."
          )
        ),
        "runtime",
        group = dimension,
        call = call,
        categories = unknown
      )
    }
    as.integer(codes)
  }
}

# Migrants' ancestry in the engine's columns: template shares where the
# template tracks the dimension, all in the template category otherwise, and
# NA for assigned dimensions (filled at arrival).
migrant_ancestry <- function(
  template,
  rows,
  members,
  dimensions,
  shared,
  registries,
  components
) {
  n <- nrow(members)
  blocks <- lapply(
    seq_along(dimensions),
    function(k) {
      columns <- components[[k]]
      d <- dimensions[k]
      block <- matrix(0, nrow = n, ncol = length(columns))
      if (!length(columns)) {
        return(block)
      }
      if (!d %in% shared) {
        block[] <- NA_real_
        return(block)
      }
      shares <- template$ancestry[[d]]
      if (d %in% tracked_ancestry(template$groups) && !is.null(shares)) {
        shares <- shares[
          match(members$person_id, shares$person_id), ,
          drop = FALSE
        ]
        for (j in seq_along(columns)) {
          if (columns[j] %in% names(shares)) {
            block[, j] <- shares[[columns[j]]]
          }
        }
        return(block)
      }
      block[cbind(seq_len(n), match(as.character(members[[d]]), columns))] <- 1
      block
    }
  )
  do.call(cbind, c(list(matrix(numeric(), nrow = n, ncol = 0L)), blocks))
}

# Arrival columns of a population's persons: both or neither, consistent,
# and within each migrant's life and the population's endpoint.
check_arrivals <- function(persons, end_time, call = rlang::caller_env()) {
  if (!"arrival_time" %in% names(persons)) {
    return(invisible(persons))
  }
  arrival <- persons$arrival_time
  check_numeric_column(
    arrival, "persons$arrival_time",
    missing = TRUE,
    call = call
  )
  # An arrival's birthplace is their stream, never "native" (D190).
  if (birthplace_dimension %in% names(persons)) {
    place <- as.character(persons[[birthplace_dimension]])
    abort_population_rows(
      which(!is.na(arrival) & (is.na(place) | place == native_birthplace)),
      "Immigrants' {.field birthplace} must be their stream.",
      "Immigrants with birthplace {.val native} or {.code NA}",
      call = call
    )
  }
  abort_population_rows(
    which(!is.na(arrival) & arrival < persons$birth_time),
    "Arrivals must not precede births.",
    "{.field arrival_time} is before {.field birth_time}",
    call = call
  )
  abort_population_rows(
    which(!is.na(arrival) & !is.na(persons$death_time) &
            persons$death_time < arrival),
    "Deaths must not precede arrivals.",
    "{.field death_time} is before {.field arrival_time}",
    call = call
  )
  abort_population_rows(
    which(!is.na(arrival) & arrival > end_time),
    "Arrivals must not occur after {.field end_time} ({end_time}).",
    "Later arrivals",
    call = call
  )
  invisible(persons)
}

# External relatives (D186): people referenced by immigrants' histories who
# are not population members. Checks the specification's rules 1-7.
check_external_relatives <- function(x, call = rlang::caller_env()) {
  external <- x$external_relatives
  persons <- x$persons
  required <- c(
    "person_id", "sex", "birth_time", "death_time", "mother_id",
    "father_id", "censor_time"
  )
  check_columns(external, required, "external_relatives", call = call)
  check_numeric_column(
    external$person_id, "external_relatives$person_id",
    id = TRUE,
    call = call
  )
  check_unique_ids(
    external$person_id, "external_relatives$person_id",
    call = call
  )
  shared <- which(external$person_id %in% persons$person_id)
  abort_population_rows(
    shared,
    "External relatives must not also be in {.field persons}.",
    "IDs also in {.field persons}",
    table = "external_relatives",
    call = call
  )
  check_labels(
    external$sex, "external_relatives$sex", c("m", "f"),
    call = call
  )
  for (field in c("birth_time", "death_time", "censor_time")) {
    check_numeric_column(
      external[[field]], paste0("external_relatives$", field),
      missing = field == "death_time",
      call = call
    )
  }
  ids <- c(persons$person_id, external$person_id)
  for (field in c("mother_id", "father_id")) {
    parents <- external[[field]]
    check_numeric_column(
      parents, paste0("external_relatives$", field),
      missing = TRUE,
      id = TRUE,
      call = call
    )
    abort_population_rows(
      which(
        !is.na(parents) &
          (!parents %in% ids | parents == external$person_id)
      ),
      "{.field {field}} must refer to another person or be {.code NA}.",
      "Unknown or self references in {.field {field}}",
      table = "external_relatives",
      call = call
    )
  }
  died <- !is.na(external$death_time)
  abort_population_rows(
    which(died & external$death_time < external$birth_time),
    "Deaths must not precede births.",
    "{.field death_time} is before {.field birth_time}",
    table = "external_relatives",
    call = call
  )
  abort_population_rows(
    which(
      (died & external$censor_time != external$death_time) |
        (!died & (external$censor_time < external$birth_time |
                    external$censor_time > x$end_time))
    ),
    paste(
      "{.field censor_time} must equal {.field death_time}, or for the",
      "living lie between {.field birth_time} and {.field end_time}."
    ),
    "Invalid {.field censor_time}",
    table = "external_relatives",
    call = call
  )
  unions <- x$unions
  member_unions <- unions$wife_id %in% persons$person_id |
    unions$husband_id %in% persons$person_id
  referenced <- external$person_id %in% c(
    unions$wife_id[member_unions], unions$husband_id[member_unions],
    persons$mother_id, persons$father_id
  ) |
    external$mother_id %in% persons$person_id |
    external$father_id %in% persons$person_id
  abort_population_rows(
    which(!referenced),
    paste(
      "Each external relative must be a partner, parent or child of someone",
      "in {.field persons}."
    ),
    "Unreferenced external relatives",
    table = "external_relatives",
    call = call
  )
  current <- unions[is.na(unions$end_time), , drop = FALSE]
  abort_population_rows(
    which(external$person_id %in% c(current$wife_id, current$husband_id)),
    "External relatives must not be in a union without an end time.",
    "Current partners",
    table = "external_relatives",
    call = call
  )
  for (d in names(x$groups)) {
    categories <- x$groups[[d]]$categories
    if (!d %in% names(external) || is.null(categories)) {
      next
    }
    values <- external[[d]]
    abort_population_rows(
      which(
        !is.na(values) & !as.character(values) %in% as.character(categories)
      ),
      "External relatives' {.field {d}} must be a declared category or NA.",
      "Undeclared categories",
      table = "external_relatives",
      call = call
    )
  }
  invisible(x)
}

# `pre_migration` on unions and transitions: present exactly when persons
# has arrival records, logical and complete (D186).
check_pre_migration <- function(x, call = rlang::caller_env()) {
  arrivals <- all(arrival_columns %in% names(x$persons))
  for (table in c("unions", "transitions")) {
    rows <- x[[table]]
    if (is.null(rows)) {
      next
    }
    present <- "pre_migration" %in% names(rows)
    if (present != arrivals) {
      abort_hammel(
        c(
          paste(
            "{.field {table}} must have {.field pre_migration} exactly when",
            "{.field persons} has arrival records."
          ),
          "x" = if (arrivals) {
            paste(
              "{.field persons} has arrival records, but {.field {table}}",
              "has no {.field pre_migration}."
            )
          } else {
            paste(
              "{.field {table}} has {.field pre_migration} without arrival",
              "records."
            )
          }
        ),
        "population",
        argument = table,
        call = call
      )
    }
    if (present &&
          (!is.logical(rows$pre_migration) || anyNA(rows$pre_migration))) {
      abort_hammel(
        c(
          "{.field {table}$pre_migration} must be TRUE or FALSE.",
          "x" = "It is {describe_value(rows$pre_migration)}."
        ),
        "table",
        argument = paste0(table, "$pre_migration"),
        call = call
      )
    }
  }
  invisible(x)
}

# Bind two tables whose columns may differ; a column one lacks is NA there.
bind_columns <- function(old, new) {
  if (is.null(old)) {
    return(tibble::as_tibble(new))
  }
  columns <- union(names(old), names(new))
  for (column in setdiff(columns, names(old))) {
    old[[column]] <- new[[column]][rep(NA_integer_, nrow(old))]
  }
  for (column in setdiff(columns, names(new))) {
    new[[column]] <- old[[column]][rep(NA_integer_, nrow(new))]
  }
  tibble::as_tibble(
    rbind(as.data.frame(old)[columns], as.data.frame(new)[columns])
  )
}

# Add the imported histories of a run's arrivals (D186) to its result: each
# migrant's template unions, their partners and children as external
# relatives, their template transitions in shared groups, and an entry row
# at arrival for each assigned group. Dates shift by the arrival time minus
# the template endpoint; every ID is fresh.
import_histories <- function(
  result,
  arrivals,
  process,
  input,
  declaration,
  registries
) {
  if (is.null(arrivals) || !length(arrivals$person_id)) {
    return(result)
  }
  dimensions <- names(declaration)
  # Each arrival's categories at arrival, by dimension (engine codes).
  arrived_in <- matrix(
    arrivals$categories,
    ncol = length(dimensions), byrow = TRUE
  )
  next_person <- max(
    c(result$persons$person_id, result$external_relatives$person_id)
  )
  next_union <- max(c(0, result$unions$union_id))
  next_move <- max(c(0, result$transitions$move_id), na.rm = TRUE)
  external <- list()
  unions <- list()
  moves <- list()
  live_unions <- list()
  # Arrivals by stream and template cut (D196); each cut's template supplies
  # its arrivals' histories.
  batches <- unique(data.frame(s = arrivals$stream, k = arrivals$slice))
  for (batch in seq_len(nrow(batches))) {
    s <- batches$s[batch] + 1L
    k <- batches$k[batch] + 1L
    a <- which(arrivals$stream == s - 1L & arrivals$slice == k - 1L)
    template <- input$streams[[s]]$templates[[k]]
    people <- template_people(template)
    t_unions <- as.data.frame(template$unions)
    t_moves <- template$transitions
    if (is.null(t_moves)) {
      t_moves <- empty_transitions()
    }
    template_dims <- intersect(
      dimensions, setdiff(names(template$groups), birthplace_dimension)
    )
    # Clones import no records, and enter every group at arrival (D188).
    if (is.null(process$streams[[s]]$household)) {
      t_unions <- t_unions[0, , drop = FALSE]
      t_moves <- t_moves[0, , drop = FALSE]
      people <- people[0, , drop = FALSE]
      template_dims <- character()
    }
    assigned <- setdiff(dimensions, template_dims)
    pid <- arrivals$person_id[a]
    time <- arrivals$time[a]
    offset <- time - template$end_time
    party <- arrivals$party[a]
    rows <- input$streams[[s]]$slices[[k]]$rows[arrivals$row[a] + 1L]
    tp <- template$persons$person_id[rows]
    key <- as.character(tp)
    member_key <- paste(party, tp)
    # Records whose `ids` are an arrival's template person, as (arrival,
    # record) pairs in arrival and record order. Only the batch's arrivals
    # are indexed, not the whole template.
    by_arrival <- split(seq_along(tp), tp)
    pairs <- function(ids, records = seq_along(ids)) {
      hit <- which(ids %in% tp)
      found <- by_arrival[as.character(ids[hit])]
      arrival <- as.integer(unlist(found, use.names = FALSE))
      record <- as.integer(rep(records[hit], lengths(found)))
      o <- order(arrival, record)
      list(arrival = arrival[o], record = as.integer(record[o]))
    }
    parents <- c(people$mother_id, people$father_id)
    known <- !is.na(parents)
    # Each party's ended unions, once each; current unions arrive live from
    # the engine (D187).
    u <- pairs(
      c(t_unions$wife_id, t_unions$husband_id),
      rep(seq_len(nrow(t_unions)), 2L)
    )
    u_party <- party[u$arrival]
    keep <- !duplicated(paste(u_party, u$record)) &
      !is.na(t_unions$end_time[u$record])
    u_party <- u_party[keep]
    u_record <- u$record[keep]
    u_first <- u$arrival[keep]
    k <- pairs(parents[known], rep(seq_len(nrow(people)), 2L)[known])
    # Relatives: partners and children outside the party, once per party.
    rel_party <- c(
      u_party, u_party, party[k$arrival]
    )
    rel_id <- c(
      t_unions$wife_id[u_record], t_unions$husband_id[u_record],
      people$person_id[k$record]
    )
    child <- rep(
      c(FALSE, TRUE), c(2L * length(u_record), length(k$record))
    )
    keep <- !paste(rel_party, rel_id) %in% member_key &
      !duplicated(paste(rel_party, rel_id))
    rel_party <- rel_party[keep]
    rel_id <- rel_id[keep]
    child <- child[keep]
    new_id <- next_person + seq_along(rel_id)
    next_person <- next_person + length(rel_id)
    lookup_key <- c(paste(rel_party, rel_id), member_key)
    lookup_id <- c(new_id, pid)
    map <- function(group, id) {
      lookup_id[match(paste(group, id), lookup_key)]
    }
    party_offset <- offset[match(rel_party, party)]
    party_time <- time[match(rel_party, party)]
    if (length(rel_id)) {
      r <- match(rel_id, people$person_id)
      table <- data.frame(
        person_id = new_id,
        sex = people$sex[r],
        birth_time = people$birth_time[r] + party_offset,
        death_time = people$death_time[r] + party_offset,
        mother_id = as.numeric(
          ifelse(child, map(rel_party, people$mother_id[r]), NA)
        ),
        father_id = as.numeric(
          ifelse(child, map(rel_party, people$father_id[r]), NA)
        ),
        stringsAsFactors = FALSE
      )
      for (d in dimensions) {
        table[[d]] <- if (d %in% template_dims) {
          c(
            template$persons[[d]],
            if (!is.null(template$external_relatives)) {
              template$external_relatives[[d]]
            }
          )[r]
        } else {
          result$persons[[d]][rep(NA_integer_, length(r))]
        }
      }
      censor <- rep(NA_real_, nrow(people))
      if (!is.null(template$external_relatives)) {
        censor[nrow(template$persons) +
                 seq_len(nrow(template$external_relatives))] <-
          template$external_relatives$censor_time
      }
      table$censor_time <- ifelse(
        is.na(table$death_time),
        ifelse(is.na(censor[r]), party_time, censor[r] + party_offset),
        table$death_time
      )
      external[[length(external) + 1L]] <- table
    }
    # Members' parents outside the party are external relatives when they
    # were a member's previous partner (D187).
    at <- match(pid, result$persons$person_id)
    for (field in c("mother_id", "father_id")) {
      outside <- map(party, template$persons[[field]][rows])
      fill <- is.na(result$persons[[field]][at]) & !is.na(outside)
      result$persons[[field]][at[fill]] <- outside[fill]
    }
    imported_key <- character()
    imported_id <- numeric()
    if (length(u_record)) {
      ids <- next_union + seq_along(u_record)
      next_union <- next_union + length(u_record)
      shift <- offset[u_first]
      x <- t_unions[u_record, , drop = FALSE]
      imported_key <- paste(u_party, x$union_id)
      imported_id <- ids
      table <- data.frame(
        union_id = ids,
        wife_id = map(u_party, x$wife_id),
        husband_id = map(u_party, x$husband_id),
        start_time = x$start_time + shift,
        end_time = x$end_time + shift,
        end_reason = as.character(x$end_reason),
        stringsAsFactors = FALSE
      )
      if ("union_type" %in% names(result$unions)) {
        table$union_type <- if ("union_type" %in% names(x)) {
          as.character(x$union_type)
        } else {
          rep("marriage", nrow(x))
        }
      }
      if ("previous_cohab_id" %in% names(result$unions)) {
        table$previous_cohab_id <- if ("previous_cohab_id" %in% names(x)) {
          ids[match(paste(u_party, x$previous_cohab_id), imported_key)]
        } else {
          rep(NA_real_, nrow(x))
        }
      }
      table$pre_migration <- rep(TRUE, nrow(table))
      unions[[length(unions) + 1L]] <- table
    }
    # Unions that arrived live are pre-migration; a converted cohabitation
    # links to its imported record.
    live <- which(arrivals$union_party %in% party)
    if (length(live)) {
      x <- t_unions[arrivals$union_row[live], , drop = FALSE]
      previous <- if ("previous_cohab_id" %in% names(x)) {
        imported_id[match(
          paste(arrivals$union_party[live], x$previous_cohab_id), imported_key
        )]
      } else {
        rep(NA_real_, length(live))
      }
      live_unions[[length(live_unions) + 1L]] <- data.frame(
        union_id = arrivals$union_id[live],
        previous_cohab_id = previous
      )
    }
    if (length(dimensions)) {
      imported <- list(arrival = integer(), record = integer())
      if (nrow(t_moves)) {
        shared <- which(t_moves$group %in% template_dims)
        imported <- pairs(t_moves$person_id[shared], shared)
      }
      m <- t_moves[imported$record, , drop = FALSE]
      entry_arrival <- rep(seq_along(a), each = length(assigned))
      entry_group <- rep(assigned, times = length(a))
      entry_code <- arrived_in[cbind(
        a[entry_arrival], match(entry_group, dimensions)
      )]
      entry_destination <- vapply(
        seq_along(entry_group),
        function(j) {
          as.character(registries[[entry_group[j]]])[entry_code[j] + 1L]
        },
        ""
      )
      n_moves <- length(imported$record) + length(entry_arrival)
      if (n_moves) {
        moves[[length(moves) + 1L]] <- data.frame(
          person_id = pid[c(imported$arrival, entry_arrival)],
          time = c(m$time + offset[imported$arrival], time[entry_arrival]),
          group = c(as.character(m$group), entry_group),
          origin = c(
            as.character(m$origin), rep(NA_character_, length(entry_group))
          ),
          destination = c(as.character(m$destination), entry_destination),
          move_id = next_move + seq_len(n_moves),
          initiator_id = pid[c(imported$arrival, entry_arrival)],
          pre_migration = rep(
            c(TRUE, FALSE), c(length(imported$record), length(entry_arrival))
          ),
          stringsAsFactors = FALSE
        )
        next_move <- next_move + n_moves
      }
    }
  }
  if (length(external)) {
    result$external_relatives <- bind_columns(
      result$external_relatives, do.call(rbind, external)
    )
  }
  if (length(unions)) {
    result$unions <- bind_columns(result$unions, do.call(rbind, unions))
  }
  if (length(live_unions)) {
    live <- do.call(rbind, live_unions)
    at <- match(live$union_id, result$unions$union_id)
    if (!"pre_migration" %in% names(result$unions)) {
      result$unions$pre_migration <- rep(NA, nrow(result$unions))
    }
    result$unions$pre_migration[at] <- TRUE
    if ("previous_cohab_id" %in% names(result$unions)) {
      result$unions$previous_cohab_id[at] <- live$previous_cohab_id
    }
  }
  if (length(moves)) {
    result$transitions <- bind_columns(
      result$transitions, do.call(rbind, moves)
    )
  }
  result
}

# Observation copies of persons with any external relatives appended and
# flagged (D186), so native histories resolve their links. Their group is a
# placeholder: they are never counted.
with_external_rows <- function(persons, external) {
  if (is.null(external) || !nrow(external)) {
    return(persons)
  }
  extra <- persons[rep(NA_integer_, nrow(external)), , drop = FALSE]
  columns <- c(
    "person_id", "sex", "birth_time", "death_time", "mother_id", "father_id"
  )
  for (column in columns) {
    extra[[column]] <- external[[column]]
  }
  if ("group" %in% names(persons) && nrow(persons)) {
    extra$group <- persons$group[1]
  }
  persons$.hammel_external <- rep(FALSE, nrow(persons))
  extra$.hammel_external <- rep(TRUE, nrow(extra))
  out <- rbind(as.data.frame(persons), as.data.frame(extra))
  # rbind drops the observation attributes; they are added afterwards.
  out
}

# Persons followed by any external relatives (D186), for resolving partner
# and parent references. Columns external relatives lack are NA.
history_people <- function(population) {
  persons <- population$persons
  external <- population[["external_relatives"]]
  if (is.null(external) || !nrow(external)) {
    return(persons)
  }
  external <- external[setdiff(names(external), "censor_time")]
  bind_columns(persons, external)
}

# People for kinship: persons, then external relatives when included
# (D186). Without them, references to external relatives are unknown.
kin_people <- function(population, include_external) {
  persons <- population$persons
  external <- population[["external_relatives"]]
  columns <- c(
    "person_id", "sex", "birth_time", "death_time", "mother_id", "father_id"
  )
  people <- as.data.frame(persons)[columns]
  if (is.null(external) || !nrow(external)) {
    return(people)
  }
  if (!include_external) {
    for (field in c("mother_id", "father_id")) {
      people[[field]][people[[field]] %in% external$person_id] <- NA
    }
    return(people)
  }
  extra <- as.data.frame(external)[c(columns, "censor_time")]
  people$censor_time <- NA_real_
  rbind(people, extra)
}

# Unions as observed in the population, with `.origin` holding each union's
# start for duration axes. Pre-migration unions are left-truncated at their
# partners' arrival (D187): those that ended by then are dropped, and the
# rest are observed from it, so a union that arrived live contributes from
# arrival at its imported duration.
observed_unions <- function(population) {
  unions <- population$unions
  unions$.origin <- unions$start_time
  if (!"pre_migration" %in% names(unions)) {
    return(unions)
  }
  persons <- population$persons
  arrival <- function(id) {
    persons$arrival_time[match(id, persons$person_id)]
  }
  from <- pmax(
    arrival(unions$wife_id), arrival(unions$husband_id),
    na.rm = TRUE
  )
  pre <- unions$pre_migration & !is.na(from)
  keep <- !unions$pre_migration |
    (pre & (is.na(unions$end_time) | unions$end_time > from))
  unions$start_time[pre] <- from[pre]
  unions[keep, , drop = FALSE]
}

# With arrival records, unions and transitions record `pre_migration`;
# events of the run, and earlier histories without the column, are not
# pre-migration (D186).
mark_pre_migration <- function(population) {
  if (!all(arrival_columns %in% names(population$persons))) {
    return(population)
  }
  for (table in c("unions", "transitions")) {
    rows <- population[[table]]
    if (is.null(rows)) {
      next
    }
    flag <- rows[["pre_migration"]]
    if (is.null(flag)) {
      flag <- rep(FALSE, nrow(rows))
    }
    flag[is.na(flag)] <- FALSE
    rows$pre_migration <- as.logical(flag)
    population[[table]] <- rows
  }
  population
}

# Fill the arrival times of simulated persons from the engine's records.
# Earlier arrivals and other people keep their values; new rows are NA.
record_arrivals <- function(persons, arrivals) {
  if (!"arrival_time" %in% names(persons)) {
    persons$arrival_time <- rep(NA_real_, nrow(persons))
  }
  if (!is.null(arrivals) && length(arrivals$person_id)) {
    rows <- match(arrivals$person_id, persons$person_id)
    persons$arrival_time[rows] <- arrivals$time
  }
  persons
}
