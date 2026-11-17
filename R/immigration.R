# Specify external immigration streams, validate them, and prepare their
# migrants and arrival records for simulation.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify external immigration
#'
#' @description
#'
#' DRAFT: agent-written documentation that needs user review and editing.
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
#' Each named argument is one stream, and its name identifies the stream in
#' the `stream` column of simulated persons. Every stream must be named, and
#' names must be unique.
#'
#' A model has at most one immigration process. Each call to
#' `set_immigration()` defines the complete set of streams, replacing any
#' streams set earlier, as [set_groups()] replaces the whole group
#' declaration. Calling `set_immigration(model)` without streams removes
#' immigration from the model.
#'
#' Migrants arrive as households, or as clones of individuals when their
#' stream has no household rule (see [immigration_stream()]). A household
#' member keeps the age, sex and history of the template person they were
#' drawn from, with every date shifted so that ages and durations at the
#' template's `end_time` are unchanged at arrival.
#' Their template unions arrive as pre-migration unions, and their previous
#' partners and children become `external_relatives` of the simulated
#' population: records that inform histories, such as parity and
#' previous-union status, but are never simulated or counted. Their own
#' parents are unknown unless they arrive in the same household. Group
#' duration continues from the template, or for groups assigned by
#' [migrant_assignment()] starts at arrival, recorded as a `transitions` row
#' with no `origin`. A household's current union arrives with it, at its
#' imported duration; it is a pre-migration union, so rate calculations count
#' its dissolution and conversion only from arrival and never count its
#' formation.
#'
#' Simulated populations of models with immigration gain `arrival_time` and
#' `stream` columns in `persons`, which are `NA` for people who did not
#' immigrate, and a `pre_migration` column in `unions` and `transitions`.
#' External relatives have a `censor_time`, after which nothing is known about
#' them. Observation functions such as [calculate_life_table()] count a
#' migrant's exposure and events only from their arrival.
#'
#' @param ... Named [immigration_stream()] objects, one per stream. The
#'   argument names are the stream names.
#' @param model A `hammel_model` object, used by `set_immigration()`.
#' @return `immigration()` returns an inspectable process inheriting from
#'   `hammel_immigration` and `hammel_process`, with `name` and a named list
#'   of `streams`. `set_immigration()` returns an updated `hammel_model`.
#' @seealso [immigration_stream()], [migrant_assignment()], [set_process()]
#' @export
#' @examples
#' template <- generate_uniform_population(500)
#' model <- hammel_model() |>
#'   set_population(generate_uniform_population(1000)) |>
#'   set_mortality(usa_2000_mortality) |>
#'   set_fertility(usa_2000_fertility) |>
#'   set_union_formation("on_demand") |>
#'   set_immigration(origin = immigration_stream(template, intensity = 20))
#' result <- run_simulation(model, years = 10, seed = 1)
#' table(result$persons$stream, useNA = "ifany")
#' @order 1
immigration <- function(...) {
  new_immigration(list(...))
}

#' @rdname immigration
#' @export
#' @order 2
set_immigration <- function(model, ...) {
  attach_process(model, new_immigration(list(...)))
}

#' Describe one stream of immigrants
#'
#' @description
#'
#' DRAFT: agent-written documentation that needs user review and editing.
#'
#' Describe a stream of migrants for [immigration()]: the template population
#' migrants are drawn from, how many arrive, who is selected and how groups
#' the template lacks are assigned.
#'
#' @details
#'
#' Arrivals follow a Poisson process with the given intensity in migrants per
#' year. The expected number of arrivals over a period is the integral of the
#' intensity; the actual number is random.
#'
#' By default, the template's living people are divided into households as
#' for household moves in [set_group_transition()]: a woman in a union and
#' her partner, or an unpartnered adult, together with the dependents of
#' either. Each arrival draws one household, with replacement, represented
#' by its carrier: the woman in a union, or the unpartnered adult. Every
#' member arrives at the same time, keeps their age at the template's
#' `end_time` (so their birth time is the arrival time minus that age) and
#' their history, and shares the categories assigned by `assignment`.
#' `intensity` still counts people, so households arrive at `intensity`
#' divided by the expected household size, and arrivals come in groups, with
#' more variation than individual arrivals. Drawing does not change the
#' template.
#'
#' With `household = NULL`, the stream clones individuals instead: anyone
#' alive at the template's `end_time` can be drawn, and keeps their age, sex
#' and groups but none of their history. Clones arrive never partnered, with
#' no children, unknown parents and group durations starting at arrival.
#' This suits templates built for the purpose, such as young adults sampled
#' as labor migrants; for templates with histories, household streams are
#' more realistic, and a cloning stream warns that the histories will not
#' be imported.
#'
#' Group dimensions that both the template and the simulation declare must
#' have the same categories, and migrants keep their template category.
#' Template dimensions that the simulation does not declare are ignored. Each
#' simulation dimension that the template lacks needs a [migrant_assignment()]
#' in `assignment`. When the simulation tracks ancestry for a dimension,
#' migrants take their template ancestry shares if the template tracks them,
#' and otherwise have all their ancestry in their category.
#'
#' @param template A `hammel_population` from which migrants are drawn.
#' @param intensity The expected number of arrivals per year: a single
#'   nonnegative number, or a data frame with numeric `time_start`, `time_end`
#'   and `intensity` columns giving piecewise-constant intensities over
#'   half-open calendar intervals. Intervals must not overlap and must cover
#'   the simulation run; use a zero intensity to pause a stream.
#' @param weights `NULL` (default) to draw eligible template people uniformly,
#'   or a [weight_model()] whose formula or function may use `age` (at the
#'   template's `end_time`), `sex` and the template's groups shared with the
#'   simulation. Each eligible person is drawn with probability proportional
#'   to their weight; weights do not change the stream's intensity. For
#'   household streams, weights describe each household by its carrier, and
#'   may also use `partnered` (0/1), `household_size`, `n_dependents` and
#'   `partner_age` (`NA` for unpartnered carriers; in a formula the term then
#'   contributes nothing).
#' @param assignment `NULL` (default) or a list of [migrant_assignment()]
#'   objects named by the simulation's group dimensions that the template
#'   lacks.
#' @param household A [household()] rule (default) to draw whole households
#'   with their histories, or `NULL` to clone individuals without histories.
#' @return An `hammel_immigration_stream` object.
#' @seealso [immigration()], [migrant_assignment()], [weight_model()]
#' @export
#' @examples
#' template <- generate_uniform_population(500)
#' immigration_stream(template, intensity = 20)
#'
#' # Fewer arrivals after 2010, and younger migrants are more likely.
#' immigration_stream(
#'   template,
#'   intensity = data.frame(
#'     time_start = c(2000, 2010),
#'     time_end = c(2010, Inf),
#'     intensity = c(40, 20)
#'   ),
#'   weights = weight_model(~age, c(age = log(0.97)))
#' )
immigration_stream <- function(
  template,
  intensity,
  weights = NULL,
  assignment = NULL,
  household = hammel::household()
) {
  new_immigration_stream(template, intensity, weights, assignment, household)
}

#' Assign migrants to a group category
#'
#' @description
#'
#' DRAFT: agent-written documentation that needs user review and editing.
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
#' - `"custom"` draws each migrant's category from probabilities returned by
#'   `fun(members)`. `members` is a tibble with one row per eligible template
#'   person, holding `person_id` (the template ID), `age` (at the template's
#'   `end_time`), `sex` and the template's groups shared with the simulation.
#'   For household streams, the rows are household carriers, with the
#'   household's `partnered` (0/1), `household_size`, `n_dependents` and
#'   `partner_age` (`NA` for unpartnered carriers), and every member of a
#'   household receives the carrier's category.
#'   `fun` returns a data frame or matrix with one row per member and one
#'   column per category, named by category; each row's probabilities must
#'   sum to one.
#'
#' The probabilities are computed before the simulation starts; the category
#' is drawn when each migrant arrives.
#'
#' @param method One of `"fixed"`, `"probability"` or `"custom"`.
#' @param fixed_value For `"fixed"`, one category label.
#' @param probs For `"probability"`, probabilities named by category.
#' @param fun For `"custom"`, a function of `members` returning
#'   probabilities, as described in details.
#' @return A `hammel_migrant_assignment` object.
#' @seealso [immigration_stream()], [group_inheritance()]
#' @export
#' @examples
#' migrant_assignment("fixed", fixed_value = "north")
#' migrant_assignment("probability", probs = c(north = 0.6, south = 0.4))
#' migrant_assignment(
#'   "custom",
#'   fun = function(members) {
#'     young <- members$age < 30
#'     data.frame(
#'       north = ifelse(young, 0.8, 0.5),
#'       south = ifelse(young, 0.2, 0.5)
#'     )
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
#' DRAFT: agent-written documentation that needs user review and editing.
#'
#' Summary method for class `"hammel_immigration"`
#'
#' @param object A `hammel_immigration` object returned by [immigration()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_immigration`
#' containing `streams`, a named list with one element per stream holding
#' `template_size` (persons in the template), `eligible` (households, or for
#' cloning streams template people, that can be drawn), `household` (the
#' household rule's `age_independence`, or `NULL` for cloning streams),
#' `intensity` (a number, or a tibble of calendar intervals),
#' `weights` (a weight model summary, or `NULL` for uniform draws) and
#' `assignment` (the method assigning each group dimension, named by
#' dimension).
#'
#' @seealso [immigration()]
#' @examples
#' process <- immigration(
#'   origin = immigration_stream(generate_uniform_population(100), 10)
#' )
#' summary(process)
#' @export
summary.hammel_immigration <- function(object, ...) {
  streams <- lapply(
    object$streams,
    function(stream) {
      intensity <- stream$intensity
      if (is.data.frame(intensity)) {
        intensity <- tibble::as_tibble(
          intensity[c("time_start", "time_end", "intensity")]
        )
      }
      parties <- migrant_parties(stream$template, stream$household)
      list(
        template_size = nrow(stream$template$persons),
        eligible = length(parties$party_start) - 1L,
        household = stream$household$age_independence,
        intensity = intensity,
        weights = summarize_model_argument(stream$weights),
        assignment = vapply(stream$assignment, `[[`, "", "method")
      )
    }
  )
  structure(list(streams = streams), class = "summary.hammel_immigration")
}

#' @rdname summary.hammel_immigration
#' @param x A `summary.hammel_immigration` object returned by
#'   `summary.hammel_immigration()`.
#' @export
print.summary.hammel_immigration <- function(x, ...) {
  cli::cli_h2("Immigration process")
  list <- open_summary_list()
  print_stat("Streams", length(x$streams))
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
    intensity <- stream$intensity
    print_stat(
      "Intensity",
      if (is.data.frame(intensity)) {
        paste0(
          intensity$intensity, " per year in ",
          format_band(intensity$time_start, intensity$time_end),
          collapse = "; "
        )
      } else {
        paste(intensity, "per year")
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

# Columns that record arrivals in `persons` (D185).
arrival_columns <- c("arrival_time", "stream")

assignment_methods <- c("fixed", "probability", "custom")

# Carrier variables of household streams, for weights and assignment (D189).
household_variables <- c(
  "partnered", "household_size", "n_dependents", "partner_age"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Build a validated immigration process from a list of streams; errors name
# `call`.
new_immigration <- function(streams, call = rlang::caller_env()) {
  process <- structure(
    list(name = "immigration", streams = streams),
    class = c("hammel_immigration", "hammel_process")
  )
  validate_immigration(process, call = call)
  process
}

validate_immigration <- function(process, call = rlang::caller_env()) {
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
  template,
  intensity,
  weights,
  assignment,
  household,
  call = rlang::caller_env()
) {
  if (missing(intensity)) {
    abort_hammel(
      c(
        "{.arg intensity} is required.",
        "i" = "Give the expected number of arrivals per year."
      ),
      "argument",
      argument = "intensity",
      call = call
    )
  }
  stream <- structure(
    list(
      template = template,
      intensity = intensity,
      weights = weights,
      assignment = assignment,
      household = household
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
  if (is_model_schedule(weights)) {
    abort_hammel(
      c(
        "Immigration weights can't follow a calendar schedule.",
        "i" = "Give a single {.fn weight_model}."
      ),
      "specification",
      argument = "weights",
      call = call
    )
  }
  if (!is.null(weights) && !inherits(weights, "hammel_weight_model")) {
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
  household <- stream$household
  used <- if (!is.null(weights$formula)) all.vars(weights$formula)
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

# A stream's intensity: one nonnegative number, or calendar intervals with
# nonnegative intensities that do not overlap.
check_intensity <- function(intensity, call = rlang::caller_env()) {
  if (!is.data.frame(intensity)) {
    if (!is.numeric(intensity) || length(intensity) != 1L ||
          !is.null(dim(intensity)) || !is.finite(intensity) ||
          intensity < 0) {
      abort_hammel(
        c(
          paste(
            "{.arg intensity} must be a single nonnegative number or a data",
            "frame of calendar intervals."
          ),
          "x" = "It is {describe_value(intensity)}."
        ),
        "argument",
        argument = "intensity",
        call = call
      )
    }
    return(invisible(intensity))
  }
  required <- c("time_start", "time_end", "intensity")
  check_columns(intensity, required, "intensity", call = call)
  extra <- setdiff(names(intensity), required)
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
        "{.arg intensity} must contain at least one calendar interval.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = "intensity",
      call = call
    )
  }
  check_bands(intensity, "time", "intensity", minimum = -Inf, call = call)
  check_numeric_column(
    intensity$intensity, "intensity$intensity",
    call = call
  )
  check_nonnegative(intensity$intensity, "intensity$intensity", call = call)
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

# A stream's intensity as calendar intervals; a constant applies throughout.
intensity_table <- function(intensity) {
  if (!is.data.frame(intensity)) {
    return(
      data.frame(time_start = -Inf, time_end = Inf, intensity = intensity)
    )
  }
  rows <- order(intensity$time_start)
  data.frame(
    time_start = as.numeric(intensity$time_start[rows]),
    time_end = as.numeric(intensity$time_end[rows]),
    intensity = as.numeric(intensity$intensity[rows])
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
    0, vapply(streams, function(x) max(c(0, x$age)), numeric(1))
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
  template <- stream$template
  intensity <- intensity_table(stream$intensity)
  if (is.data.frame(stream$intensity)) {
    check_time_coverage(
      intensity$time_start, intensity$time_end, start, end,
      "The stream's intensity",
      call = call
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
  template_dims <- names(template$groups)
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
  rules <- stream$assignment
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
    stream$weights, members, registries[shared], name,
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
  for (k in seq_along(dimensions)) {
    d <- dimensions[k]
    if (d %in% shared) {
      categories[, k] <- match(
        as.character(people[[d]]), as.character(registries[[d]])
      ) - 1L
    } else {
      assignment[[k]] <- assignment_probabilities(
        rules[[d]],
        members[intersect(
          c("person_id", "age", "sex", shared, household_variables),
          names(members)
        )],
        as.character(registries[[d]]), d,
        call = call
      )
    }
  }
  live <- template$unions[parties$union_row, , drop = FALSE]
  list(
    time_start = intensity$time_start,
    time_end = intensity$time_end,
    intensity = intensity$intensity / expected_size,
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
  children <- !is.na(people$mother_id)
  key <- as.character(ids)
  mothers <- as.character(people$mother_id[children])
  counts <- c(table(mothers))
  parity <- as.integer(counts[key])
  parity[is.na(parity)] <- 0L
  last <- rep(NA_real_, length(ids))
  if (any(children)) {
    latest <- c(tapply(people$birth_time[children], mothers, max))
    last <- as.numeric(latest[key])
  }
  parity[!female] <- 0L
  last[!female] <- NA_real_
  previous <- call_native(
    previous_union_info_cpp(people, unions, counts_cohab),
    call = call
  )
  at <- match(ids, people$person_id)
  since_entry <- matrix(NA_real_, nrow = length(ids), ncol = length(dimensions))
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
        latest <- tapply(moves$time, moves$person_id, max)
        i <- match(names(latest), as.character(ids))
        entered[i] <- as.numeric(latest)
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
  call = rlang::caller_env()
) {
  n <- nrow(members)
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  if (!is.null(weights$fun)) {
    eta <- weights$fun(members, list(process = "immigration", stream = name))
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

# Each eligible template person's probabilities over a dimension's
# categories, as an n x K matrix.
assignment_probabilities <- function(
  rule,
  members,
  categories,
  dimension,
  call = rlang::caller_env()
) {
  n <- nrow(members)
  k <- length(categories)
  if (rule$method %in% c("fixed", "probability")) {
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
    row <- numeric(k)
    row[match(names(probs), categories)] <- probs
    return(matrix(row, nrow = n, ncol = k, byrow = TRUE))
  }
  result <- rule$fun(members)
  labels <- colnames(result)
  if ((!is.data.frame(result) && !is.matrix(result)) ||
        NROW(result) != n || is.null(labels)) {
    abort_hammel(
      c(
        paste(
          "A custom assignment function must return probabilities with one",
          "row per member and one column per category."
        ),
        "x" = "It returned {describe_value(result)} for {n} member{?s}."
      ),
      "runtime",
      group = dimension,
      call = call
    )
  }
  unknown <- setdiff(labels, categories)
  if (length(unknown) || anyDuplicated(labels)) {
    abort_hammel(
      c(
        "Custom assignment columns must be distinct categories.",
        "x" = if (length(unknown)) "Unknown: {.val {unknown}}.",
        "i" = "The categories of {.field {dimension}} are {.val {categories}}."
      ),
      "runtime",
      group = dimension,
      call = call
    )
  }
  values <- as.matrix(as.data.frame(result))
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0) ||
        any(abs(rowSums(values) - 1) > 1e-8)) {
    abort_hammel(
      c(
        "Custom assignment probabilities must be valid.",
        "x" = paste(
          "Each row must hold finite, nonnegative probabilities summing to",
          "one."
        )
      ),
      "runtime",
      group = dimension,
      call = call
    )
  }
  out <- matrix(0, nrow = n, ncol = k)
  out[, match(labels, categories)] <- values
  out
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
  present <- intersect(arrival_columns, names(persons))
  if (!length(present)) {
    return(invisible(persons))
  }
  if (length(present) == 1L) {
    absent <- setdiff(arrival_columns, present)
    abort_hammel(
      c(
        "Arrival records need both {.field arrival_time} and {.field stream}.",
        "x" = paste(
          "{.field persons} has {.field {present}} but not",
          "{.field {absent}}."
        )
      ),
      "population",
      argument = "persons",
      call = call
    )
  }
  arrival <- persons$arrival_time
  check_numeric_column(
    arrival, "persons$arrival_time",
    missing = TRUE,
    call = call
  )
  stream <- persons$stream
  if (!is.character(stream)) {
    abort_hammel(
      c(
        "{.field persons$stream} must contain stream names.",
        "x" = "It is {describe_value(stream)}."
      ),
      "table",
      argument = "persons$stream",
      call = call
    )
  }
  abort_population_rows(
    which(is.na(arrival) != is.na(stream)),
    paste(
      "{.field arrival_time} and {.field stream} must both be missing or",
      "both be present."
    ),
    "Mismatched arrival records",
    call = call
  )
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
  for (s in seq_along(process$streams)) {
    a <- which(arrivals$stream == s - 1L)
    if (!length(a)) {
      next
    }
    template <- process$streams[[s]]$template
    people <- template_people(template)
    t_unions <- as.data.frame(template$unions)
    t_moves <- template$transitions
    if (is.null(t_moves)) {
      t_moves <- empty_transitions()
    }
    template_dims <- intersect(dimensions, names(template$groups))
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
    rows <- input$streams[[s]]$rows[arrivals$row[a] + 1L]
    tp <- template$persons$person_id[rows]
    key <- as.character(tp)
    member_key <- paste(party, tp)
    # Records of each arrival's template person, as (arrival, record) pairs.
    pairs <- function(index) {
      found <- index[key]
      list(
        arrival = rep(seq_along(a), lengths(found)),
        record = as.integer(unlist(found, use.names = FALSE))
      )
    }
    unions_of <- split(
      rep(seq_len(nrow(t_unions)), 2L),
      as.character(c(t_unions$wife_id, t_unions$husband_id))
    )
    parents <- c(people$mother_id, people$father_id)
    known <- !is.na(parents)
    children_of <- split(
      rep(seq_len(nrow(people)), 2L)[known], as.character(parents[known])
    )
    # Each party's ended unions, once each; current unions arrive live from
    # the engine (D187).
    u <- pairs(unions_of)
    u_party <- party[u$arrival]
    keep <- !duplicated(paste(u_party, u$record)) &
      !is.na(t_unions$end_time[u$record])
    u_party <- u_party[keep]
    u_record <- u$record[keep]
    u_first <- u$arrival[keep]
    k <- pairs(children_of)
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
        mother_id = ifelse(child, map(rel_party, people$mother_id[r]), NA),
        father_id = ifelse(child, map(rel_party, people$father_id[r]), NA),
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
        imported <- pairs(
          split(shared, as.character(t_moves$person_id[shared]))
        )
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

# Fill the arrival columns of simulated persons from the engine's records.
# Earlier arrivals and other people keep their values; new rows are NA.
record_arrivals <- function(persons, arrivals, stream_names) {
  for (column in arrival_columns) {
    if (!column %in% names(persons)) {
      persons[[column]] <- if (column == "stream") {
        rep(NA_character_, nrow(persons))
      } else {
        rep(NA_real_, nrow(persons))
      }
    }
  }
  if (!is.null(arrivals) && length(arrivals$person_id)) {
    rows <- match(arrivals$person_id, persons$person_id)
    persons$arrival_time[rows] <- arrivals$time
    persons$stream[rows] <- stream_names[arrivals$stream + 1L]
  }
  persons
}
