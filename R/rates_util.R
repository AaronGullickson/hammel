# Validate input rate schedules and prepare calendar and group strata for
# execution.

# HELPER FUNCTIONS ---------------------------------------------------------

# Shared validation and engine preparation of input rate schedules.

# RATE SUMMARY METADATA ---------------------------------------------------

# Summarize an already-validated rate (or probability) table. NULL means no
# table is configured. Strata retain their input types and observed order;
# no unobserved combinations are added. Callers may supply explicit stratum
# names when the input also carries columns that are not schedule strata.
#
# Returns stratum names, a tibble of observed strata, and a bands tibble.
# Bands retain joint age/duration/parity/calendar combinations, including
# unused NA clock pairs in mixed age/duration schedules. Identical band sets
# across every stratum are returned once, without stratum columns; otherwise
# the bands include their stratum columns and bands_by_stratum is TRUE.
# Rates and duplicate rows do not affect the comparison. Boundary values and
# half-open interval semantics are preserved, without formatting or rounding.
summarize_rate_strata <- function(rates, strata = NULL) {
  if (is.null(rates)) {
    return(NULL)
  }
  if (is.null(strata)) {
    strata <- intersect(
      names(rates), c(builtin_strata, candidate_strata(rates))
    )
  }
  band_columns <- intersect(
    c(
      "age_start", "age_end", "dur_start", "dur_end",
      "parity_start", "parity_end", "time_start", "time_end"
    ),
    names(rates)
  )
  rates <- tibble::as_tibble(rates)
  stratum_values <- unique(rates[strata])
  bands <- unique(rates[c(strata, band_columns)])
  common <- unique(rates[band_columns])

  # Each distinct stratum/band row is one observed combination. The full
  # Cartesian count is reached exactly when all strata share all bands.
  shared <- !length(strata) ||
    nrow(bands) == as.double(nrow(stratum_values)) * nrow(common)
  if (shared) {
    bands <- common
  }

  list(
    stratified_by = strata,
    strata = stratum_values,
    bands = bands,
    bands_by_stratum = !shared
  )
}

# CALENDAR SCHEDULES AND COVERAGE ------------------------------------------

# Shared calendar validation for flat rate tables. Each active slice follows
# the existing demographic validator; time gaps are checked against the run.
has_calendar_rates <- function(rates) {
  any(c("time_start", "time_end") %in% names(rates))
}

calendar_rate_segments <- function(
  rates,
  start = -Inf,
  end = Inf,
  name = "rates",
  call = rlang::caller_env()
) {
  if (!has_calendar_rates(rates)) {
    return(list(list(start = -Inf, end = Inf, rates = rates)))
  }
  check_columns(rates, c("time_start", "time_end"), name, call = call)
  if (!nrow(rates)) {
    abort_hammel(
      c(
        "{.arg {name}} must contain at least one calendar interval.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = name,
      call = call
    )
  }
  check_bands(rates, "time", name, minimum = -Inf, call = call)
  ends <- rates$time_end
  times <- sort(unique(c(rates$time_start, ends)))
  indices <- seq_len(length(times) - 1L)
  keep <- if (start == end) {
    times[indices] <= start & times[indices + 1L] > start
  } else {
    times[indices] < end & times[indices + 1L] > start
  }
  lapply(
    indices[keep], function(i) {
      time <- times[i]
      rows <- rates$time_start <= time & time < rates$time_end
      list(
        start = time,
        end = times[i + 1L],
        rates = rates[rows, setdiff(names(rates), c("time_start", "time_end")),
          drop = FALSE
        ]
      )
    }
  )
}

calendar_rate_slices <- function(
  rates,
  start = -Inf,
  end = Inf,
  call = rlang::caller_env()
) {
  lapply(calendar_rate_segments(rates, start, end, call = call), `[[`, "rates")
}

# Validate each calendar period of a rate table with `validate`, which takes
# the period's rows (and its interval when `with_interval`). An error names
# the period. Returns FALSE when the table has no calendar columns.
validate_calendar_rates <- function(
  rates,
  validate,
  with_interval = FALSE,
  call = rlang::caller_env()
) {
  if (!has_calendar_rates(rates)) {
    return(FALSE)
  }
  for (segment in calendar_rate_segments(rates, call = call)) {
    if (nrow(segment$rates)) {
      with_error_context(
        if (with_interval) {
          validate(segment$rates, c(segment$start, segment$end))
        } else {
          validate(segment$rates)
        },
        paste0(
          "In calendar period ",
          format_band(segment$start, segment$end), "."
        )
      )
    }
  }
  TRUE
}

# Every time from `start` up to `end` must fall in some interval
# [starts, ends). `subject` names the schedule, as in "Mortality rates".
check_time_coverage <- function(
  starts,
  ends,
  start,
  end,
  subject,
  call = rlang::caller_env()
) {
  times <- sort(unique(c(start, starts, ends)))
  times <- times[times >= start & (times < end | times == start)]
  for (time in times) {
    if (!any(starts <= time & time < ends)) {
      abort_hammel(
        c(
          "{subject} must cover the whole simulation.",
          "x" = "Nothing applies at time {time}.",
          "i" = paste(
            "The simulation runs from {start} to {end}. Extend the",
            "{.field time_start} and {.field time_end} intervals to cover",
            "{time}, or shorten {.arg years}."
          )
        ),
        "coverage",
        time = time,
        call = call
      )
    }
  }
  invisible(NULL)
}

check_calendar_coverage <- function(
  rates,
  start,
  end,
  name,
  call = rlang::caller_env()
) {
  if (is.null(rates) || !has_calendar_rates(rates)) {
    return(invisible(NULL))
  }
  calendar_rate_slices(rates, call = call) # validate bounds of edited inputs
  labels <- c(
    mortality = "Mortality rates",
    fertility = "Fertility rates",
    search = "Union search rates",
    dissolution = "Union dissolution rates",
    conversion = "Cohabitation conversion rates",
    cohabitation_probability = "Cohabitation probabilities"
  )
  subject <- if (name %in% names(labels)) {
    labels[[name]]
  } else {
    paste(name, "rates")
  }
  check_time_coverage(
    rates$time_start, rates$time_end, start, end, subject,
    call = call
  )
}

# GROUP REGISTRIES AND ENGINE STRATA ---------------------------------------

# Group dimensions in rate tables and the engine (D134). Each declared
# dimension has a registry: its ordered categories. A rate table stratifies by
# the dimensions it has columns for: the dimension's name for person-level
# processes, or `<name>_w` and `<name>_h` (wife's and husband's categories)
# for union-level ones. Its strata combine those categories in declaration
# order, wife before husband.

# Ordered categories of every declared dimension: declared categories, or
# the population's values plus transition destinations and categories with
# ancestry shares, in first-occurrence order.
# `transition_rates` is a list of rate tables keyed by dimension.
group_registries <- function(population, transition_rates = list()) {
  declaration <- population$groups
  persons <- population$persons
  registries <- stats::setNames(
    lapply(
      names(declaration), function(name) {
        categories <- declaration[[name]]$categories
        if (!is.null(categories)) {
          return(categories)
        }
        values <- unique(persons[[name]])
        # Transition histories and this dimension's transition rates.
        tx <- population$transitions
        if (!is.null(tx)) {
          tx <- tx[as.character(tx$group) == name, , drop = FALSE]
        }
        rates <- transition_rates[[name]]
        values <- extend_identity_groups(
          values, c(
            as.character(tx$origin),
            as.character(tx$destination), as.character(rates$origin),
            as.character(rates$destination)
          )
        )
        shares <- population[["ancestry"]][[name]]
        if (!is.null(shares)) {
          known <- setdiff(names(shares), c("person_id", "unknown"))
          positive <- known[vapply(
            shares[known], function(x) any(x > 0), logical(1)
          )]
          values <- extend_identity_groups(values, positive)
        }
        values
      }
    ), names(declaration)
  )
  attr(registries, "declared") <- vapply(
    declaration, function(g) !is.null(g$categories), logical(1)
  )
  registries
}

# The group columns of a rate table: list of (dimension, role) pairs, with role
# 0 for a person or wife and 1 for a husband.
rate_strata_columns <- function(rates, registries, union_level = FALSE) {
  out <- list()
  for (d in seq_along(registries)) {
    name <- names(registries)[d]
    if (!union_level) {
      if (name %in% names(rates)) {
        out[[length(out) + 1L]] <- list(column = name, dim = d, role = 0L)
      }
    } else {
      for (role in 0:1) {
        column <- paste0(name, if (role) "_h" else "_w")
        if (column %in% names(rates)) {
          out[[length(out) + 1L]] <- list(column = column, dim = d, role = role)
        }
      }
    }
  }
  out
}

# Index each rate row by its stratum and check that every stratum is covered
# in each active calendar slice. Returns the rates (with `.hammel_group` when
# stratified) and the engine specification of the process's strata.
# Without `coverage`, uncovered strata have zero risk (group transitions).
prepare_strata_rates <- function(
  rates,
  process,
  registries,
  start = -Inf,
  end = Inf,
  union_level = FALSE,
  coverage = TRUE,
  call = rlang::caller_env()
) {
  rates$.hammel_group <- NULL
  # Person-level tables use the dimension's name; union-level tables the
  # wife's and husband's columns.
  dims <- names(registries)
  wrong <- if (union_level) {
    c(intersect(c(dims, "group_pair"), names(rates)))
  } else {
    intersect(c(paste0(dims, "_w"), paste0(dims, "_h")), names(rates))
  }
  if (length(wrong)) {
    subject <- upper_first(process)
    hint <- if (union_level) {
      paste(
        "Union rates are stratified by the wife's and husband's groups:",
        "use {.field <group>_w} and {.field <group>_h} columns."
      )
    } else {
      "These rates are stratified by a group dimension's own name."
    }
    abort_hammel(
      c(
        "{subject} rates can't use {.field {wrong}} columns.",
        "i" = hint
      ),
      "groups",
      columns = wrong,
      call = call
    )
  }
  columns <- rate_strata_columns(rates, registries, union_level)
  spec <- list(
    dims = integer(),
    roles = integer(),
    strides = integer(),
    count = 1L
  )
  if (!length(columns)) {
    check_strata_overlap(rates, process, call = call)
    return(list(rates = rates, spec = spec))
  }
  counts <- vapply(columns, function(x) length(registries[[x$dim]]), integer(1))
  strides <- as.integer(cumprod(c(1, utils::head(counts, -1))))
  index <- integer(nrow(rates))
  declared <- attr(registries, "declared")
  keep <- rep(TRUE, nrow(rates))
  for (k in seq_along(columns)) {
    x <- columns[[k]]
    code <- match(
      as.character(rates[[x$column]]), as.character(registries[[x$dim]])
    ) - 1L
    # Rows for categories absent from an inferred registry are unused; with
    # declared categories they are errors.
    if (anyNA(code) && isTRUE(declared[[x$dim]])) {
      subject <- upper_first(process)
      dimension <- names(registries)[x$dim]
      unknown <- unique(as.character(rates[[x$column]])[is.na(code)])
      categories <- as.character(registries[[x$dim]])
      abort_hammel(
        c(
          paste(
            "{subject} rates use categories that aren't declared for",
            "group {.field {dimension}}."
          ),
          "x" = "Column {.field {x$column}} contains {.val {unknown}}.",
          "i" = "Declared categories: {.val {categories}}."
        ),
        "groups",
        group = dimension,
        categories = unknown,
        call = call
      )
    }
    keep <- keep & !is.na(code)
    index <- index + code * strides[k]
  }
  count <- as.integer(prod(counts))
  label <- function(k) {
    paste(
      vapply(
        seq_along(columns), function(j) {
          value <- as.character(registries[[columns[[j]]$dim]])[
            (k %/% strides[j]) %% counts[j] + 1L
          ]
          if (length(columns) == 1L) {
            value
          } else {
            paste0(columns[[j]]$column, "=", value)
          }
        }, ""
      ),
      collapse = ", "
    )
  }
  # Coverage is checked before rows for absent categories are dropped, so a
  # slice whose rows all name absent categories still reports every stratum
  # as missing (rather than looking empty).
  rates$.hammel_group <- index
  missing <- unique(
    unlist(
      lapply(
        calendar_rate_slices(rates, start, end, call = call),
        function(x) {
          if (nrow(x)) {
            setdiff(
              seq_len(count) - 1L, x$.hammel_group[!is.na(x$.hammel_group)]
            )
          } else {
            integer()
          }
        }
      )
    )
  )
  rates <- rates[keep, , drop = FALSE]
  if (coverage && length(missing)) {
    subject <- upper_first(process)
    absent <- vapply(sort(missing), label, "")
    shown <- cli::cli_vec(absent, list("vec-trunc" = 5))
    abort_hammel(
      c(
        "{subject} rates are missing schedules for some groups.",
        "x" = "No rates for {shown}.",
        "i" = paste(
          "Rates stratified by a group need a schedule for every category",
          "combination of the group columns they use."
        )
      ),
      "groups",
      call = call
    )
  }
  check_strata_overlap(rates, process, call = call)
  spec <- list(
    dims = vapply(columns, `[[`, integer(1), "dim") - 1L,
    roles = vapply(columns, `[[`, integer(1), "role"),
    strides = strides,
    count = count
  )
  list(rates = rates, spec = spec)
}

# The engine's view of every dimension: each person's category codes, the
# reporting strata over all dimensions (unrealized births), and per-process
# strata.
# `match_dims` are the matching dimensions: the engine's matching category is
# their combination (first fastest), as for rate strata.
group_engine_input <- function(persons, registries, match_dims, specs) {
  n <- length(registries)
  categories <- matrix(0L, nrow(persons), n)
  for (d in seq_len(n)) {
    categories[, d] <- match(
      as.character(persons[[names(registries)[d]]]),
      as.character(registries[[d]])
    ) - 1L
  }
  counts <- vapply(registries, length, integer(1))
  report <- list(
    dims = seq_len(n) - 1L,
    roles = integer(n),
    strides = as.integer(cumprod(c(1, utils::head(counts, -1)))),
    count = as.integer(prod(counts))
  )
  c(
    list(
      dimensions = n,
      categories = categories,
      counts = unname(counts),
      labels = unname(lapply(registries, as.character)),
      names = names(registries),
      matching = strata_spec(
        match(match_dims, names(registries)) - 1L,
        counts[match_dims]
      ),
      report = report
    ), specs
  )
}

# A person-level strata specification over dimensions (0-based) of given
# category counts.
strata_spec <- function(dims, counts) {
  counts <- as.integer(counts)
  list(
    dims = as.integer(dims),
    roles = integer(length(dims)),
    strides = as.integer(cumprod(c(1, utils::head(counts, -1))))[seq_along(
      dims
    )],
    count = as.integer(prod(counts))
  )
}

# Decode reporting strata (all dimensions) into one column per dimension.
decode_report_strata <- function(index, registries) {
  counts <- vapply(registries, length, integer(1))
  strides <- cumprod(c(1, utils::head(counts, -1)))
  stats::setNames(
    lapply(
      seq_along(registries), function(d) {
        registries[[d]][(index %/% strides[d]) %% counts[d] + 1L]
      }
    ), names(registries)
  )
}

# Within a stratum, bands must not overlap. Columns that are not declared group
# dimensions do not stratify rates, so a table split only by such a column
# (validated as separate strata at construction) is rejected here.
check_strata_overlap <- function(rates, process, call = rlang::caller_env()) {
  if (!nrow(rates)) {
    return(invisible(rates))
  }
  keys <- intersect(
    c(
      ".hammel_group", "sex", "union_type", "previous_union",
      "parity_start",
      "destination", "time_start"
    ), names(rates)
  )
  shown_keys <- setdiff(keys, ".hammel_group")
  for (axis in c("age", "dur")) {
    lo <- paste0(axis, "_start")
    if (!lo %in% names(rates)) {
      next
    }
    rows <- which(!is.na(rates[[lo]]))
    key <- if (length(keys)) {
      do.call(
        paste, c(
          lapply(rates[rows, keys, drop = FALSE], as.character),
          sep = "\r"
        )
      )
    } else {
      rep("", length(rows))
    }
    for (r in split(rows, key)) {
      schedule <- if (length(shown_keys)) {
        values <- vapply(
          shown_keys,
          function(column) {
            as.character(rates[[column]][r[1]])
          },
          character(1)
        )
        paste0(" for ", paste(shown_keys, "=", values, collapse = ", "))
      } else {
        ""
      }
      with_error_context(
        check_band_overlap(
          rates, r, axis,
          schedule = schedule,
          call = call
        ),
        paste(
          "In", process, "rates. Only declared group dimensions stratify",
          "rates, so rows that differ only in other columns form one",
          "schedule."
        )
      )
    }
  }
  invisible(rates)
}

# GROUP LABELS AND SCHEDULE SPLITTING --------------------------------------

# Shared group validation and execution-only schedule indexing. Public labels
# and their types are preserved; the engine uses compact integer indices.
check_groups <- function(table, name, call = rlang::caller_env()) {
  if (!"group" %in% names(table)) {
    return(invisible(table))
  }
  check_group_labels(table$group, paste0(name, "$group"), call = call)
  invisible(table)
}

# Group labels: nonmissing, nonempty character labels or finite numeric
# codes. `name` is written as `table$column`.
check_group_labels <- function(x, name, call = rlang::caller_env()) {
  subject <- describe_column(name)
  if (!(is.character(x) || is.numeric(x)) || !is.null(dim(x))) {
    abort_hammel(
      c(
        "{subject} must contain group labels.",
        "x" = "It is {describe_value(x)}.",
        "i" = "Use character labels or finite numeric codes."
      ),
      "groups",
      argument = name,
      call = call
    )
  }
  rows <- which(
    is.na(x) | !nzchar(trimws(as.character(x))) |
      (is.numeric(x) & !is.finite(x))
  )
  if (length(rows)) {
    abort_hammel(
      c(
        "{subject} must contain nonmissing, nonempty group labels.",
        "x" = "Invalid labels at {format_rows(rows)}."
      ),
      "groups",
      argument = name,
      rows = rows,
      call = call
    )
  }
  invisible(x)
}

# Row sets of each stratum of a rate table, by the given group columns.
split_rate_groups <- function(
  rates,
  columns = "group",
  call = rlang::caller_env()
) {
  columns <- intersect(columns, names(rates))
  for (column in columns) {
    check_group_labels(rates[[column]], paste0("rates$", column), call = call)
  }
  if (!length(columns)) {
    return(list(seq_len(nrow(rates))))
  }
  key <- do.call(paste, c(lapply(rates[columns], as.character), sep = "\r"))
  split(seq_len(nrow(rates)), key, drop = TRUE)
}

# BAND VALIDATION ----------------------------------------------------------

# Display names of band axes, as in "age bands".
band_axes <- c(
  age = "age",
  dur = "duration",
  parity = "parity",
  time = "calendar"
)

# Check a table's `<axis>_start`/`<axis>_end` columns row by row: finite
# starts of at least `minimum`, and numeric, nonmissing ends (possibly Inf)
# greater than their starts.
check_bands <- function(
  table,
  axis,
  name = "rates",
  minimum = 0,
  kind = "rates",
  call = rlang::caller_env()
) {
  lo <- paste0(axis, "_start")
  hi <- paste0(axis, "_end")
  check_numeric_column(table[[lo]], paste0(name, "$", lo), call = call)
  ends <- table[[hi]]
  end_column <- describe_column(paste0(name, "$", hi))
  if (!is.numeric(ends) || !is.null(dim(ends))) {
    abort_hammel(
      c(
        "{end_column} must contain numbers.",
        "x" = "It is {describe_value(ends)}."
      ),
      "table",
      argument = paste0(name, "$", hi),
      call = call
    )
  }
  rows <- which(is.na(ends))
  if (length(rows)) {
    abort_hammel(
      c(
        "{end_column} must not be missing.",
        "x" = "Missing at {format_rows(rows)}.",
        "i" = "Use {.code Inf} for a band with no upper limit."
      ),
      "table",
      argument = paste0(name, "$", hi),
      rows = rows,
      call = call
    )
  }
  rows <- which(ends <= table[[lo]])
  if (length(rows)) {
    axis_name <- band_axes[[axis]]
    abort_hammel(
      c(
        "Each {axis_name} band in {.arg {name}} must end after it starts.",
        "x" = paste(
          "{.field {hi}} is not greater than {.field {lo}} at",
          "{format_rows(rows)}."
        )
      ),
      kind,
      argument = name,
      rows = rows,
      call = call
    )
  }
  rows <- which(table[[lo]] < minimum)
  if (length(rows)) {
    start_column <- describe_column(paste0(name, "$", lo))
    abort_hammel(
      c(
        "{start_column} must be nonnegative.",
        "x" = "Negative at {format_rows(rows)}."
      ),
      kind,
      argument = name,
      rows = rows,
      call = call
    )
  }
  invisible(table)
}

# Rows `rows` of `table` form one schedule, whose `<axis>` bands must not
# overlap. `schedule` describes it for the message, as in ' for sex "f"'.
check_band_overlap <- function(
  table,
  rows,
  axis,
  name = "rates",
  schedule = "",
  kind = "rates",
  call = rlang::caller_env()
) {
  lo <- paste0(axis, "_start")
  hi <- paste0(axis, "_end")
  rows <- rows[order(table[[lo]][rows])]
  if (length(rows) < 2L) {
    return(invisible(table))
  }
  reach <- cummax(table[[hi]][rows])
  overlaps <- which(table[[lo]][rows[-1L]] < reach[-length(rows)])
  if (length(overlaps)) {
    k <- overlaps[1]
    earlier <- rows[seq_len(k)][which.max(table[[hi]][rows[seq_len(k)]])]
    later <- rows[k + 1L]
    first <- format_band(table[[lo]][earlier], table[[hi]][earlier])
    second <- format_band(table[[lo]][later], table[[hi]][later])
    axis_name <- band_axes[[axis]]
    abort_hammel(
      c(
        "{upper_first(axis_name)} bands in {.arg {name}} must not overlap.",
        "x" = "Bands {first} and {second} overlap{schedule}."
      ),
      kind,
      argument = name,
      call = call
    )
  }
  invisible(table)
}

# MORTALITY AND FERTILITY SCHEDULE VALIDATION ------------------------------

# Validate a mortality or fertility rate table in two passes: row checks on
# the whole table, so reported rows are the input's rows, then schedule
# checks within each calendar period and group stratum.
validate_rates <- function(rates, mortality, call = rlang::caller_env()) {
  check_stratified_columns(
    rates, if (mortality) "sex" else "union_type",
    if (mortality) "Mortality" else "Fertility",
    call = call
  )
  check_rate_rows(rates, mortality, call = call)
  validate_rate_schedules(rates, mortality, call = call)
  invisible(rates)
}

check_rate_rows <- function(rates, mortality, call = rlang::caller_env()) {
  columns <- c(if (mortality) "sex", "age_start", "age_end", "rate")
  check_columns(rates, columns, "rates", call = call)
  if (!nrow(rates)) {
    abort_hammel(
      c(
        "{.arg rates} must contain at least one age band.",
        "x" = "It has no rows."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  calendar_rate_segments(rates, call = call)
  for (column in candidate_strata(rates)) {
    check_group_labels(rates[[column]], paste0("rates$", column), call = call)
  }
  check_bands(rates, "age", call = call)
  check_numeric_column(rates$rate, "rates$rate", call = call)
  check_nonnegative(rates$rate, "rates$rate", call = call)
  if (mortality) {
    check_labels(rates$sex, "rates$sex", c("m", "f"), call = call)
    return(invisible(rates))
  }
  if ("sex" %in% names(rates)) {
    abort_hammel(
      c(
        "Fertility rates can't have a {.field sex} column.",
        "i" = "Fertility rates apply to women; remove the column."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  if ("union_status" %in% names(rates)) {
    abort_hammel(
      c(
        "Fertility rates can't have a {.field union_status} column.",
        "i" = "It has been renamed {.field union_type}."
      ),
      "rates",
      argument = "rates",
      call = call
    )
  }
  if ("union_type" %in% names(rates)) {
    check_labels(
      rates$union_type, "rates$union_type",
      c("unpartnered", "partnered", "cohab", "marriage"),
      call = call
    )
  }
  if (any(c("parity_start", "parity_end") %in% names(rates))) {
    check_columns(rates, c("parity_start", "parity_end"), "rates", call = call)
    check_bands(rates, "parity", call = call)
    start <- rates$parity_start
    end <- rates$parity_end
    rows <- which(start != floor(start) | (is.finite(end) & end != floor(end)))
    if (length(rows)) {
      abort_hammel(
        c(
          "Parity bands in {.arg rates} must use whole numbers.",
          "x" = "Fractional parities at {format_rows(rows)}.",
          "i" = "Use {.code Inf} for an open-ended parity band."
        ),
        "rates",
        argument = "rates",
        rows = rows,
        call = call
      )
    }
    limit <- .Machine$integer.max
    rows <- which(start > limit | (is.finite(end) & end > limit))
    if (length(rows)) {
      abort_hammel(
        c(
          "Parity bands in {.arg rates} must be at most {limit}.",
          "x" = "Larger parities at {format_rows(rows)}."
        ),
        "rates",
        argument = "rates",
        rows = rows,
        call = call
      )
    }
  }
  invisible(rates)
}

# Schedule checks within each calendar period and group stratum, after
# check_rate_rows(): both sexes, a complete set of union types, parity bands
# covering every parity, and age bands that do not overlap.
validate_rate_schedules <- function(
  rates,
  mortality,
  call = rlang::caller_env()
) {
  if (validate_calendar_rates(
        rates,
        function(x) {
          validate_rate_schedules(x, mortality, call = call)
        },
        call = call
      )) {
    return(invisible(rates))
  }
  strata <- candidate_strata(rates)
  if (length(strata)) {
    for (rows in split_rate_groups(rates, strata, call = call)) {
      schedule <- rates[rows, setdiff(names(rates), strata), drop = FALSE]
      values <- vapply(
        strata,
        function(column) {
          as.character(rates[[column]][rows[1]])
        },
        character(1)
      )
      with_error_context(
        validate_rate_schedules(schedule, mortality, call = call),
        paste0("In group stratum ", paste(strata, "=", values, collapse = ", "),
          "."
        )
      )
    }
    return(invisible(rates))
  }
  if (mortality) {
    absent <- setdiff(c("m", "f"), rates$sex)
    if (length(absent)) {
      abort_hammel(
        c(
          "Mortality rates must include schedules for both sexes.",
          "x" = "No rates for sex {.val {absent}}."
        ),
        "rates",
        argument = "rates",
        call = call
      )
    }
    schedules <- split(seq_len(nrow(rates)), rates$sex)
    labels <- paste0(' for sex "', names(schedules), '"')
  } else if ("union_type" %in% names(rates)) {
    types <- unique(rates$union_type)
    if (!(setequal(types, c("partnered", "unpartnered")) ||
            setequal(types, c("unpartnered", "cohab", "marriage")))) {
      abort_hammel(
        c(
          paste(
            "Column {.field union_type} of {.arg rates} must contain",
            "{.val unpartnered} and {.val partnered}, or {.val unpartnered},",
            "{.val cohab} and {.val marriage}."
          ),
          "x" = "It contains {.val {types}}."
        ),
        "rates",
        argument = "rates",
        call = call
      )
    }
    schedules <- split(seq_len(nrow(rates)), rates$union_type)
    labels <- paste0(' for union type "', names(schedules), '"')
  } else {
    schedules <- list(seq_len(nrow(rates)))
    labels <- ""
  }
  if (!mortality && "parity_start" %in% names(rates)) {
    parity <- split_parity_schedules(rates, schedules, labels, call)
    schedules <- parity$schedules
    labels <- parity$labels
  }
  for (i in seq_along(schedules)) {
    check_band_overlap(
      rates, schedules[[i]], "age",
      schedule = labels[i],
      call = call
    )
  }
  invisible(rates)
}

# Split each schedule by its parity bands, which must cover every parity
# from 0 without gaps or overlaps.
split_parity_schedules <- function(rates, schedules, labels, call) {
  out <- list()
  out_labels <- character()
  for (i in seq_along(schedules)) {
    rows <- schedules[[i]]
    bands <- unique(rates[rows, c("parity_start", "parity_end"), drop = FALSE])
    bands <- bands[order(bands$parity_start), , drop = FALSE]
    if (bands$parity_start[1] != 0 || tail(bands$parity_end, 1) != Inf ||
          (nrow(bands) > 1L &&
             any(head(bands$parity_end, -1) != tail(bands$parity_start, -1)))) {
      shown <- format_band(bands$parity_start, bands$parity_end)
      schedule <- labels[i]
      abort_hammel(
        c(
          "Parity bands must cover every parity from 0 to {.code Inf}.",
          "x" = "The parity bands{schedule} are {shown}.",
          "i" = "Bands must meet without gaps or overlaps."
        ),
        "rates",
        argument = "rates",
        call = call
      )
    }
    for (j in seq_len(nrow(bands))) {
      out[[length(out) + 1L]] <- rows[
        rates$parity_start[rows] == bands$parity_start[j] &
          rates$parity_end[rows] == bands$parity_end[j]
      ]
      band <- format_band(bands$parity_start[j], bands$parity_end[j])
      out_labels <- c(
        out_labels,
        if (nzchar(labels[i])) {
          paste0(labels[i], ", parity ", band)
        } else {
          paste0(" for parity ", band)
        }
      )
    }
  }
  list(schedules = out, labels = out_labels)
}
