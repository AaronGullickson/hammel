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

calendar_rate_segments <- function(rates, start = -Inf, end = Inf) {
  if (!has_calendar_rates(rates)) {
    return(list(list(start = -Inf, end = Inf, rates = rates)))
  }
  check_columns(rates, c("time_start", "time_end"), "rates")
  check_numeric_column(rates$time_start, "rates$time_start")
  ends <- rates$time_end
  if (!nrow(rates) || !is.numeric(ends) || !is.null(dim(ends)) || anyNA(ends) ||
        any(ends <= rates$time_start)) {
    stop(
      paste0(
        "Calendar intervals need finite starts and greater ",
        "ends (possibly Inf)."
      ),
      call. = FALSE
    )
  }
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

calendar_rate_slices <- function(rates, start = -Inf, end = Inf) {
  lapply(calendar_rate_segments(rates, start, end), `[[`, "rates")
}

validate_calendar_rates <- function(rates, validate, with_interval = FALSE) {
  if (!has_calendar_rates(rates)) {
    return(FALSE)
  }
  for (segment in calendar_rate_segments(rates)) {
    if (nrow(segment$rates)) {
      if (with_interval) {
        validate(segment$rates, c(segment$start, segment$end))
      } else {
        validate(segment$rates)
      }
    }
  }
  TRUE
}

check_time_coverage <- function(starts, ends, start, end, description) {
  times <- sort(unique(c(start, starts, ends)))
  times <- times[times >= start & (times < end | times == start)]
  for (time in times) {
    if (!any(starts <= time & time < ends)) {
      stop(sprintf("%s at time %g.", description, time), call. = FALSE)
    }
  }
  invisible(NULL)
}

check_calendar_coverage <- function(rates, start, end, name) {
  if (is.null(rates) || !has_calendar_rates(rates)) {
    return(invisible(NULL))
  }
  calendar_rate_slices(rates) # validate bounds, including edited inputs
  check_time_coverage(
    rates$time_start, rates$time_end, start, end,
    sprintf("%s rates have a calendar coverage gap", name)
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
  coverage = TRUE
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
    stop(
      sprintf(
        if (union_level) {
          paste0(
            "%s rates use %s; union rates are stratified by ",
            "the wife's and husband's groups, `<group>_w` and ",
            "`<group>_h`."
          )
        } else {
          paste0(
            "%s rates use %s; these rates are stratified by a ",
            "group dimension's own name."
          )
        }, process,
        paste0("`", wrong, "`", collapse = ", ")
      ),
      call. = FALSE
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
    check_strata_overlap(rates, process)
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
      stop(
        sprintf(
          paste0(
            "%s rates use `%s` values that are not categories of ",
            "group `%s`: %s."
          ),
          process,
          x$column, names(registries)[x$dim],
          paste(
            unique(as.character(rates[[x$column]])[is.na(code)]),
            collapse = ", "
          )
        ),
        call. = FALSE
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
        calendar_rate_slices(rates, start, end),
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
    stop(
      sprintf(
        "%s rates are missing schedules for group(s): %s.", process,
        paste(vapply(sort(missing), label, ""), collapse = "; ")
      ),
      call. = FALSE
    )
  }
  check_strata_overlap(rates, process)
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
check_strata_overlap <- function(rates, process) {
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
  for (axis in c("age", "dur")) {
    lo <- paste0(axis, "_start")
    hi <- paste0(axis, "_end")
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
      r <- r[order(rates[[lo]][r])]
      if (length(r) > 1L &&
            any(
              rates[[lo]][r[-1L]] < utils::head(cummax(rates[[hi]][r]), -1L)
            )) {
        stop(
          sprintf(
            paste0(
              "%s rates have overlapping bands within a stratum; ",
              "only declared group dimensions stratify rates."
            ),
            process
          ),
          call. = FALSE
        )
      }
    }
  }
  invisible(rates)
}

# GROUP LABELS AND SCHEDULE SPLITTING --------------------------------------

# Shared group validation and execution-only schedule indexing. Public labels
# and their types are preserved; the engine uses compact integer indices.
check_groups <- function(table, name) {
  if (!"group" %in% names(table)) {
    return(invisible(table))
  }
  x <- table$group
  if (!(is.character(x) || is.numeric(x)) || !is.null(dim(x)) || anyNA(x) ||
        any(!nzchar(trimws(as.character(x)))) ||
        (is.numeric(x) && any(!is.finite(x)))) {
    stop(
      sprintf(
        "`%s$group` must contain nonmissing, nonempty group labels.",
        name
      ),
      call. = FALSE
    )
  }
  invisible(table)
}

# Row sets of each stratum of a rate table, by the given group columns.
split_rate_groups <- function(rates, columns = "group") {
  columns <- intersect(columns, names(rates))
  for (column in columns) {
    check_groups(
      tibble::tibble(group = rates[[column]]),
      paste0("rates$", column)
    )
  }
  if (!length(columns)) {
    return(list(seq_len(nrow(rates))))
  }
  key <- do.call(paste, c(lapply(rates[columns], as.character), sep = "\r"))
  split(seq_len(nrow(rates)), key, drop = TRUE)
}

# MORTALITY AND FERTILITY SCHEDULE VALIDATION ------------------------------

validate_rates <- function(rates, mortality) {
  check_stratified_columns(
    rates, if (mortality) "sex" else "union_type",
    if (mortality) "Mortality" else "Fertility"
  )
  if (validate_calendar_rates(
        rates,
        function(x) validate_rates(x, mortality)
      )) {
    return(invisible(rates))
  }
  columns <- c("age_start", "age_end", "rate")
  if (mortality) {
    columns <- c("sex", columns)
  }
  check_columns(rates, columns, "rates")
  # Possible group columns: each stratum is validated as its own schedule.
  strata <- candidate_strata(rates)
  if (length(strata)) {
    for (rows in split_rate_groups(rates, strata)) {
      schedule <- rates[rows, setdiff(names(rates), strata), drop = FALSE]
      validate_rates(schedule, mortality)
    }
    if (!nrow(rates)) {
      stop("`rates` must contain at least one age interval.", call. = FALSE)
    }
    return(invisible(rates))
  }
  if (nrow(rates) == 0L) {
    stop("`rates` must contain at least one age interval.", call. = FALSE)
  }
  check_numeric_column(rates$age_start, "rates$age_start")
  check_numeric_column(rates$rate, "rates$rate")
  ends <- rates$age_end
  if (!is.numeric(ends) || !is.null(dim(ends)) || anyNA(ends) ||
        any(ends <= rates$age_start) ||
        any(rates$age_start < 0)) {
    stop(
      paste0(
        "Age intervals must have nonnegative starts and ",
        "greater ends (possibly Inf)."
      ),
      call. = FALSE
    )
  }
  if (any(rates$rate < 0)) {
    stop("`rates$rate` must be nonnegative.", call. = FALSE)
  }
  if (mortality) {
    if (!is.character(rates$sex) || !is.null(dim(rates$sex)) ||
          anyNA(rates$sex) ||
          any(!rates$sex %in% c("m", "f")) ||
          !all(c("m", "f") %in% rates$sex)) {
      stop(
        'Mortality rates must include both "m" and "f" sex schedules.',
        call. = FALSE
      )
    }
    schedules <- split(seq_len(nrow(rates)), rates$sex)
  } else {
    if ("sex" %in% names(rates)) {
      stop(
        paste0(
          "Fertility schedules do not take a sex column in ",
          "this initial interface."
        ),
        call. = FALSE
      )
    }
    if ("union_status" %in% names(rates)) {
      stop(
        "Fertility input `union_status` has been replaced by `union_type`.",
        call. = FALSE
      )
    }
    if ("union_type" %in% names(rates)) {
      status <- rates$union_type
      if (!is.character(status) || !is.null(dim(status)) || anyNA(status) ||
            !(setequal(status, c("partnered", "unpartnered")) ||
                setequal(status, c("unpartnered", "cohab", "marriage")))) {
        stop(
          paste0(
            "Fertility `union_type` must contain ",
            "unpartnered/partnered or unpartnered/cohab/marriage, ",
            "with no other or missing values."
          ),
          call. = FALSE
        )
      }
      schedules <- split(seq_len(nrow(rates)), status)
    } else {
      schedules <- list(seq_len(nrow(rates)))
    }
  }
  if (!mortality && any(c("parity_start", "parity_end") %in% names(rates))) {
    check_columns(rates, c("parity_start", "parity_end"), "rates")
    check_numeric_column(rates$parity_start, "rates$parity_start")
    upper <- rates$parity_end
    if (!is.numeric(upper) || !is.null(dim(upper)) || anyNA(upper) ||
          any(
            rates$parity_start < 0 |
              rates$parity_start != floor(rates$parity_start)
          ) ||
          any(upper <= rates$parity_start) ||
          any(is.finite(upper) & upper != floor(upper)) ||
          any(rates$parity_start > .Machine$integer.max) ||
          any(is.finite(upper) & upper > .Machine$integer.max)) {
      stop(
        paste0(
          "Parity bands need nonnegative integer starts and ",
          "greater integer ends or Inf within engine count capacity."
        ),
        call. = FALSE
      )
    }
    schedules <- unlist(
      lapply(
        schedules, function(rows) {
          bands <- unique(
            rates[rows, c("parity_start", "parity_end"), drop = FALSE]
          )
          bands <- bands[order(bands$parity_start), , drop = FALSE]
          if (bands$parity_start[1] != 0 || tail(bands$parity_end, 1) != Inf ||
                (nrow(bands) > 1L &&
                   any(
                     head(bands$parity_end, -1) !=
                       tail(bands$parity_start, -1)
                   ))) {
            stop(
              paste0(
                "Parity bands must cover [0, Inf) without gaps or ",
                "overlaps in every group/union-type schedule."
              ),
              call. = FALSE
            )
          }
          lapply(
            seq_len(nrow(bands)), function(i) {
              rows[
                rates$parity_start[rows] == bands$parity_start[i] &
                  rates$parity_end[rows] == bands$parity_end[i]
              ]
            }
          )
        }
      ),
      recursive = FALSE
    )
  }
  for (rows in schedules) {
    rows <- rows[order(rates$age_start[rows])]
    if (length(rows) > 1L &&
          any(
            rates$age_start[rows[-1L]] <
              cummax(rates$age_end[rows])[-length(rows)]
          )) {
      stop(
        "Age intervals must not overlap within a rate schedule.",
        call. = FALSE
      )
    }
  }
  invisible(rates)
}
