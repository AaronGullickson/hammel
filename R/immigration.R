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
#' Migrants arrive individually and without histories. A migrant keeps the age
#' and sex of the template person they were drawn from, has unknown parents,
#' no union history and parity zero, and enters every applicable rate process
#' and matching pool at arrival. Simulated populations of models with
#' immigration gain `arrival_time` and `stream` columns in `persons`, which
#' are `NA` for people who did not immigrate. Observation functions such as
#' [calculate_life_table()] count a migrant's exposure and events only from
#' their arrival.
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
#' At each arrival one eligible person is drawn, with replacement, from the
#' template. Eligible people are alive at the template's `end_time` and have
#' no current partner. The migrant keeps their age at the template's
#' `end_time`, so their birth time is the arrival time minus that age. Drawing
#' does not change the template.
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
#'   to their weight; weights do not change the stream's intensity.
#' @param assignment `NULL` (default) or a list of [migrant_assignment()]
#'   objects named by the simulation's group dimensions that the template
#'   lacks.
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
  assignment = NULL
) {
  new_immigration_stream(template, intensity, weights, assignment)
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
#' `template_size` (persons in the template), `eligible` (template people who
#' can be drawn), `intensity` (a number, or a tibble of calendar intervals),
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
      list(
        template_size = nrow(stream$template$persons),
        eligible = length(eligible_migrants(stream$template)),
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
        stream$template_size, " persons, ", stream$eligible, " eligible"
      )
    )
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
      assignment = assignment
    ),
    class = "hammel_immigration_stream"
  )
  validate_immigration_stream(stream, call = call)
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

# Rows of the template's persons that can be drawn: alive at its end_time and
# without a current (unended) union.
eligible_migrants <- function(template) {
  persons <- template$persons
  unions <- template$unions
  current <- unions[is.na(unions$end_time), , drop = FALSE]
  partnered <- c(current$wife_id, current$husband_id)
  which(is.na(persons$death_time) & !persons$person_id %in% partnered)
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
  call = rlang::caller_env()
) {
  dimensions <- names(declaration)
  streams <- lapply(
    names(process$streams),
    function(name) {
      with_error_context(
        prepare_stream(
          process$streams[[name]], name, declaration, registries, components,
          start, end,
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
  rows <- eligible_migrants(template)
  if (!length(rows)) {
    abort_hammel(
      c(
        "The stream's template has nobody who can migrate.",
        "i" = paste(
          "Migrants are drawn from people alive at the template's",
          "{.field end_time} without a current partner."
        )
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
  members <- tibble::tibble(
    person_id = persons$person_id,
    age = template$end_time - persons$birth_time,
    sex = as.character(persons$sex)
  )
  for (d in shared) {
    members[[d]] <- persons[[d]]
  }
  probabilities <- migrant_probabilities(
    stream$weights, members, registries[shared], name,
    call = call
  )
  n <- nrow(members)
  categories <- matrix(-1L, nrow = n, ncol = length(dimensions))
  assignment <- vector("list", length(dimensions))
  for (k in seq_along(dimensions)) {
    d <- dimensions[k]
    if (d %in% shared) {
      categories[, k] <- match(
        as.character(members[[d]]), as.character(registries[[d]])
      ) - 1L
    } else {
      assignment[[k]] <- assignment_probabilities(
        rules[[d]], members[c("person_id", "age", "sex", shared)],
        as.character(registries[[d]]), d,
        call = call
      )
    }
  }
  list(
    time_start = intensity$time_start,
    time_end = intensity$time_end,
    intensity = intensity$intensity,
    age = as.numeric(members$age),
    female = members$sex == "f",
    categories = categories,
    ancestry = migrant_ancestry(
      template, rows, members, dimensions, shared, registries, components
    ),
    cumulative = c(utils::head(cumsum(probabilities), -1), 1),
    assignment = assignment
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
      eta <- eta + column$coefficient * value
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
