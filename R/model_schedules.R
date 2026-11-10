# Validate, select, and prepare calendar schedules of component models.

# HELPER FUNCTIONS ---------------------------------------------------------

# Public arguments retain direct model objects or ordinary inspectable tibbles.
is_model_schedule <- function(x) is.data.frame(x)

validate_model_argument <- function(x, validate, name) {
  if (!is_model_schedule(x)) {
    return(validate(x))
  }
  check_columns(x, c("time_start", "time_end", "model"), name)
  check_numeric_column(x$time_start, paste0(name, "$time_start"))
  ends <- x$time_end
  if (!nrow(x) || !is.numeric(ends) || !is.null(dim(ends)) || anyNA(ends) ||
        any(ends <= x$time_start)) {
    stop(
      sprintf(
        "%s schedule needs finite starts and greater ends (possibly Inf).", name
      ),
      call. = FALSE
    )
  }
  rows <- order(x$time_start)
  if (length(rows) > 1L &&
        any(tail(x$time_start[rows], -1L) < head(cummax(ends[ rows ]), -1L))) {
    stop(sprintf("%s model intervals must not overlap.", name), call. = FALSE)
  }
  if (!is.list(x$model) || !is.null(dim(x$model))) {
    stop(
      sprintf("%s$model must be a list column of complete models.", name),
      call. = FALSE
    )
  }
  for (i in seq_len(nrow(x))) {
    model <- x$model[[i]]
    if (is.null(model) || is_model_schedule(model) || is.character(model)) {
      stop(
        sprintf(
          "%s$model[[%d]] must contain a complete model object.", name,
          i
        ),
        call. = FALSE
      )
    }
    tryCatch(
      validate(model),
      error = function(e) {
        stop(
          sprintf("%s$model[[%d]]: %s", name, i, conditionMessage(e)),
          call. = FALSE
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
    stop("No unique active model in schedule.", call. = FALSE)
  }
  x$model[[i]]
}

check_model_coverage <- function(x, start, end, name) {
  if (!is_model_schedule(x)) {
    return(invisible(NULL))
  }
  check_time_coverage(
    x$time_start, x$time_end, start, end,
    sprintf("%s model schedule has a coverage gap", name)
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
  end
) {
  matching_bank <- collect_model_bank(
    arguments[c("matching", "father_matching", "cohab_matching")]
  )
  prepared_matching <- lapply(
    matching_bank, prepare_matching_model,
    population_groups = values
  )
  times <- start
  for (x in c(arguments, inheritance)) {
    if (is_model_schedule(x)) {
      times <- c(times, x$time_start, x$time_end)
    }
  }
  times <- sort(unique(times[times >= start & (times < end | times == start)]))
  id <- function(model, bank) which(
    vapply(bank, identical, logical(1), y = model)
  )[1L] - 1L
  matching_ids <- vapply(
    times, function(t) id(
      model_at(arguments$matching, t),
      matching_bank
    ), integer(1)
  )
  role_ids <- function(role) {
    if (is.null(arguments[[role]])) {
      matching_ids
    } else {
      vapply(
        times, function(t) id(
          model_at(arguments[[role]], t),
          matching_bank
        ), integer(1)
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
        components = components[[d]]
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
