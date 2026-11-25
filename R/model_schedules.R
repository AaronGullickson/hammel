# Validate, select, and prepare calendar schedules of component models.

# HELPER FUNCTIONS ---------------------------------------------------------

# Public arguments retain direct model objects or ordinary inspectable tibbles.
is_model_schedule <- function(x) is.data.frame(x)

validate_model_argument <- function(
  x,
  validate,
  name,
  call = rlang::caller_env()
) {
  if (!is_model_schedule(x)) {
    return(validate(x, call = call))
  }
  check_columns(x, c("time_start", "time_end", "model"), name, call = call)
  if (!nrow(x)) {
    abort_hammel(
      c(
        "The {.arg {name}} schedule must contain at least one model.",
        "x" = "It has no rows."
      ),
      "specification",
      argument = name,
      call = call
    )
  }
  check_bands(
    x, "time", name,
    minimum = -Inf,
    kind = "specification",
    call = call
  )
  check_band_overlap(
    x, seq_len(nrow(x)), "time", name,
    kind = "specification",
    call = call
  )
  if (!is.list(x$model) || !is.null(dim(x$model))) {
    abort_hammel(
      c(
        "Column {.field model} of {.arg {name}} must be a list of models.",
        "x" = "It is {describe_value(x$model)}."
      ),
      "specification",
      argument = name,
      call = call
    )
  }
  for (i in seq_len(nrow(x))) {
    model <- x$model[[i]]
    if (is.null(model) || is_model_schedule(model) || is.character(model)) {
      abort_hammel(
        c(
          "Each row of {.arg {name}} must contain a complete model object.",
          "x" = "Row {i} contains {describe_value(model)}."
        ),
        "specification",
        argument = name,
        rows = i,
        call = call
      )
    }
    tryCatch(
      with_error_context(
        validate(model, call = call),
        paste0("In row ", i, " of the ", name, " schedule.")
      ),
      error = function(e) {
        if (inherits(e, "hammel_error")) {
          stop(e)
        }
        abort_hammel(
          "Row {i} of {.arg {name}} contains an invalid model.",
          "specification",
          argument = name,
          rows = i,
          parent = e,
          call = call
        )
      }
    )
  }
  invisible(x)
}

normalize_model_argument <- function(x) {
  if (is_model_schedule(x)) tibble::as_tibble(x) else x
}

model_entries <- function(x) {
  if (is.null(x)) list() else if (is_model_schedule(x)) x$model else list(x)
}

model_at <- function(x, time) {
  if (!is_model_schedule(x)) {
    return(x)
  }
  i <- which(x$time_start <= time & time < x$time_end)
  if (length(i) != 1L) {
    abort_internal("No unique active model in a schedule at time {time}.")
  }
  x$model[[i]]
}

check_model_coverage <- function(
  x,
  start,
  end,
  name,
  call = rlang::caller_env()
) {
  if (!is_model_schedule(x)) {
    return(invisible(NULL))
  }
  check_time_coverage(
    x$time_start, x$time_end, start, end,
    paste("The", name, "model schedule"),
    call = call
  )
}

collect_model_bank <- function(arguments) {
  bank <- list()
  for (argument in arguments) {
    for (model in model_entries(argument)) {
      if (!any(vapply(bank, identical, logical(1), y = model))) {
        bank[[length(bank) + 1L]] <- model
      }
    }
  }
  bank
}

# One shared identity registry and model bank for the entire native invocation.
# `inheritance` is a list of each dimension's inheritance argument; all
# dimensions' rules share one bank, with a route of bank ids per dimension.
prepare_model_schedules <- function(
  arguments,
  values,
  inheritance,
  registries,
  components,
  start,
  end,
  ancestry = NULL,
  call = rlang::caller_env()
) {
  matching_bank <- collect_model_bank(
    arguments[c("matching", "father_matching", "cohab_matching")]
  )
  prepared_matching <- lapply(
    matching_bank, prepare_matching_model,
    population_groups = values,
    ancestry = ancestry,
    call = call
  )
  times <- start
  for (x in c(arguments, inheritance)) {
    if (is_model_schedule(x)) {
      times <- c(times, x$time_start, x$time_end)
    }
  }
  times <- sort(unique(times[times >= start & (times < end | times == start)]))
  id <- function(model, bank) {
    which(vapply(bank, identical, logical(1), y = model))[1L] - 1L
  }
  matching_ids <- vapply(
    times, function(t) {
      id(model_at(arguments$matching, t), matching_bank)
    }, integer(1)
  )
  role_ids <- function(role) {
    if (is.null(arguments[[role]])) {
      matching_ids
    } else {
      vapply(
        times, function(t) {
          id(model_at(arguments[[role]], t), matching_bank)
        }, integer(1)
      )
    }
  }
  prepared_inheritance <- list()
  routes <- list()
  for (d in seq_along(inheritance)) {
    bank <- collect_model_bank(inheritance[d])
    offset <- length(prepared_inheritance)
    prepared_inheritance <- c(
      prepared_inheritance, lapply(
        bank,
        prepare_group_inheritance,
        values = registries[[d]],
        components = components[[d]],
        call = call
      )
    )
    routes[[d]] <- vapply(
      times, function(t) id(model_at(inheritance[[d]], t), bank), integer(1)
    ) + offset
  }
  list(
    values = values,
    matching = prepared_matching[[matching_ids[1] + 1L]],
    father_matching = prepared_matching[[role_ids("father_matching")[1] + 1L]],
    cohab_matching = prepared_matching[[role_ids("cohab_matching")[1] + 1L]],
    inheritance = lapply(routes, function(r) prepared_inheritance[[r[1] + 1L]]),
    config = list(
      matching_models = prepared_matching,
      inheritance_models = prepared_inheritance,
      inheritance_routes = routes,
      epochs = tibble::tibble(
        time_start = times,
        matching = matching_ids,
        father_matching = role_ids("father_matching"),
        cohab_matching = role_ids("cohab_matching")
      )
    )
  )
}
