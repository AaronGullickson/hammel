# Prepare population histories, strata, and output tables for demographic
# analysis.

# HELPER FUNCTIONS ---------------------------------------------------------

# Shared history preparation and stratification for demographic analysis.
# Used by rate calculations, population counts, and relationship extraction.

# Population validation already checks label types and missing values.
observation_groups <- function(persons, by_group) {
  if (!is.logical(by_group) || length(by_group) != 1L ||
        !is.null(dim(by_group)) ||
        is.na(by_group)) {
    stop("`by_group` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!by_group) {
    return(NULL)
  }
  if (!"group" %in% names(persons)) {
    stop("`by_group = TRUE` requires a persons `group` column.", call. = FALSE)
  }
  groups <- attr(persons, ".hammel_group_values")
  if (is.null(groups)) unique(persons$group) else groups
}

# Private observation copies carry indexed histories to all native counters.
prepare_group_observations <- function(persons, population, by_group) {
  if (!isTRUE(by_group)) {
    return(persons)
  }
  # A combined view of several dimensions carries indexed codes directly.
  combined <- population[[".hammel_combined"]]
  if (!is.null(combined)) {
    attr(persons, ".hammel_group_values") <- combined$values
    h <- combined$history
    if (nrow(h)) {
      attr(persons, ".hammel_group_history") <- list(
        row = match(h$person_id, persons$person_id) - 1L,
        time = h$time,
        from = h$from,
        to = h$to
      )
    }
    return(persons)
  }
  tx <- if (inherits(population, "hammel_population")) {
    population$transitions
  } else {
    NULL
  }
  if (is.null(tx) || !nrow(tx)) {
    return(persons)
  }
  groups <- extend_identity_groups(
    unique(persons$group), c(
      as.character(tx$origin), as.character(tx$destination)
    )
  )
  attr(persons, ".hammel_group_values") <- groups
  attr(persons, ".hammel_group_history") <- list(
    row = match(tx$person_id, persons$person_id) - 1L,
    time = tx$time,
    from = match(as.character(tx$origin), as.character(groups)) - 1L,
    to = match(as.character(tx$destination), as.character(groups)) - 1L
  )
  persons
}

# Shared support for the `by` argument of rate calculators (D133).

# Validate `by` against the variables a calculator accepts.
# `dims` (declared group dimensions) replaces the internal `group` columns in
# the list of accepted names.
check_by <- function(by, allowed, fun, dims = NULL) {
  if (is.null(by)) {
    return(character())
  }
  if (!is.character(by) || !is.null(dim(by)) || anyNA(by) ||
        anyDuplicated(by)) {
    stop(
      "`by` must be NULL or a character vector of distinct variable names.",
      call. = FALSE
    )
  }
  unknown <- setdiff(by, allowed)
  if (length(unknown) && length(dims)) {
    shown <- unlist(
      lapply(
        allowed, function(a) {
          suffix <- sub("^group", "", a)
          if (a %in% c("group", "group_w", "group_h")) {
            paste0(dims, suffix)
          } else {
            a
          }
        }
      )
    )
    allowed <- shown
  }
  if (length(unknown)) {
    stop(
      sprintf(
        "`%s()` cannot split by %s; `by` accepts %s.", fun,
        paste0("\"", unknown, "\"", collapse = ", "),
        paste0("\"", allowed, "\"", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  by
}

# Validate the previous-union detail flag. Pooled categories (FALSE) only
# apply when `previous_union` is in `by`.
check_detailed_previous_union <- function(detailed, by) {
  if (!is.logical(detailed) || length(detailed) != 1L ||
        !is.null(dim(detailed)) ||
        is.na(detailed)) {
    stop("`use_detailed_previous_union` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!detailed && !"previous_union" %in% by) {
    stop(
      paste0(
        "`use_detailed_previous_union = FALSE` requires ",
        "\"previous_union\" in `by`."
      ),
      call. = FALSE
    )
  }
  invisible(detailed)
}

previous_labels <- function(use_detailed_previous_union) {
  if (use_detailed_previous_union) {
    c("none", "separated", "widowed")
  } else {
    c("none", "dissolved")
  }
}

# Persons, unions and group indices for the person-level engine. `population`
# may be a persons data frame when no union variable is needed.
person_observation <- function(
  population,
  by,
  end,
  requires_typed = FALSE,
  grouped = "group" %in% by,
  validate = TRUE
) {
  if (!inherits(population, "hammel_population")) {
    stop(
      paste0(
        "`population` must be a hammel_population; wrap a persons table ",
        "with hammel_population(persons, end_time = ...)."
      ),
      call. = FALSE
    )
  }
  if (validate) {
    validate_population(population)
  }
  if (is.null(population$end_time) || end > population$end_time) {
    stop(
      "Observation requires a known endpoint and cannot extend beyond it.",
      call. = FALSE
    )
  }
  persons <- population$persons
  unions <- population$unions
  if (grouped && !declared_group(population)) {
    stop(
      paste0(
        "Splitting by \"group\" requires a declared group ",
        "dimension; see set_groups() or ",
        "hammel_population(groups = )."
      ),
      call. = FALSE
    )
  }
  persons <- prepare_group_observations(persons, population, grouped)
  groups <- if (grouped) observation_groups(persons, TRUE) else NULL
  typed <- !is.null(unions) && "union_type" %in% names(unions)
  if (requires_typed && !typed) {
    stop(
      "Splitting by \"union_type\" requires typed union histories.",
      call. = FALSE
    )
  }
  if (is.null(unions)) {
    unions <- data.frame(
      union_id = numeric(),
      wife_id = numeric(),
      husband_id = numeric(),
      start_time = numeric(),
      end_time = numeric(),
      end_reason = character()
    )
  }
  list(
    persons = persons,
    unions = unions,
    typed = typed,
    groups = groups,
    group_indices = if (grouped) {
      match(persons$group, groups) - 1L
    } else {
      rep(0L, nrow(persons))
    },
    counts_cohab = isTRUE(population$simulation$cohab_as_previous_union)
  )
}

union_type_labels <- function(typed) {
  if (typed) {
    c("unpartnered", "cohab", "marriage")
  } else {
    c("unpartnered", "partnered")
  }
}

# Output rows for dimensions listed fastest first. Each dimension is a list
# with `columns` (a data frame of output columns, one row per value), `code`
# (0-based native codes) and `stride` (multiplier in the native array), and
# optionally `key` (global sort codes, for combining blocks). Returns the
# label columns plus `.index` and one `.key` column per dimension.
build_rows <- function(dims) {
  sizes <- vapply(dims, function(d) length(d$code), integer(1))
  grid <- expand.grid(lapply(sizes, seq_len), KEEP.OUT.ATTRS = FALSE)
  parts <- lapply(
    seq_along(dims), function(k) {
      d <- dims[[k]]
      x <- d$columns[grid[[k]], , drop = FALSE]
      rownames(x) <- NULL
      x
    }
  )
  out <- if (length(parts)) {
    do.call(cbind, parts)
  } else {
    data.frame(row.names = seq_len(nrow(grid)))
  }
  index <- numeric(nrow(grid))
  for (k in seq_along(dims)) {
    d <- dims[[k]]
    index <- index + d$code[grid[[k]]] * d$stride
    out[[paste0(".key", k)]] <- (
      if (is.null(d$key)) {
        seq_along(d$code)
      } else {
        d$key
      }
    )[grid[[k]]]
  }
  out$.index <- index + 1
  out
}

dim_of <- function(
  columns,
  code = seq_len(nrow(columns)) - 1L,
  stride = 0,
  key = NULL
) {
  list(
    columns = as.data.frame(columns, stringsAsFactors = FALSE),
    code = code,
    stride = stride,
    key = key
  )
}

label_dim <- function(
  name,
  labels,
  stride,
  code = seq_along(labels) - 1L,
  key = NULL
) {
  columns <- data.frame(labels, stringsAsFactors = FALSE)
  names(columns) <- name
  dim_of(columns, code, stride, key)
}

# Dimensions for the person-level `by` variables of one engine run.
person_by_dims <- function(
  by,
  obs,
  counts,
  use_detailed_previous_union,
  previous_codes = NULL
) {
  lapply(
    by, function(name) {
      switch(
        name,
        sex = label_dim("sex", c("m", "f"), counts$strides[["sex"]]),
        group = label_dim("group", obs$groups, counts$strides[["group"]]),
        union_type = label_dim(
          "union_type", union_type_labels(obs$typed),
          counts$strides[["union_type"]]
        ),
        previous_union = {
          labels <- previous_labels(use_detailed_previous_union)
          codes <- if (is.null(previous_codes)) {
            seq_along(labels) - 1L
          } else {
            previous_codes
          }
          label_dim(
            "previous_union",
            labels[codes + 1L],
            counts$strides[["previous_union"]],
            code = codes,
            key = codes + 1L
          )
        }
      )
    }
  )
}

finish_rows <- function(rows) {
  rows <- rows[setdiff(names(rows), ".index")]
  rows <- tibble::as_tibble(rows)
  rows
}

# The calculators split by one internal `group` column (`group_w`/`group_h`
# for union rates). `by` names declared dimensions (`"race"`, or
# `"race_w"`/`"race_h"`); with one dimension it maps directly to `group`. With
# several requested, `group` is a combined code over them (declaration order,
# first fastest), with a merged history of their moves, and finish_strata()
# splits it back into one column per dimension. `extra` adds dimensions to
# the combination (a transition's own dimension, first).
observation_alias <- function(
  population,
  by,
  suffixes = "",
  extra = character()
) {
  none <- list(population = population, by = by, rename = character())
  if (!inherits(population, "hammel_population") ||
        !length(population$groups)) {
    return(none)
  }
  dims <- names(population$groups)
  none$dims <- dims
  external <- as.vector(outer(dims, suffixes, paste0))
  if (length(dims) == 1L && !length(extra)) {
    dim <- dims
    if (dim != "group" && any(paste0("group", suffixes) %in% by)) {
      stop(
        sprintf(
          "Split by the declared group dimension `%s`, not `group`.",
          dim
        ),
        call. = FALSE
      )
    }
    population <- dimension_view(population, dim)
    hit <- by %in% external
    if (any(hit)) {
      by[hit] <- paste0("group", suffixes)[match(
        by[hit], paste0(dim, suffixes)
      )]
    }
    return(
      list(
        population = population,
        by = by,
        dims = dims,
        rename = stats::setNames(
          paste0(dim, suffixes), paste0("group", suffixes)
        )
      )
    )
  }
  if (!"group" %in% dims && any(paste0("group", suffixes) %in% by)) {
    stop(
      sprintf(
        "Split by declared group dimensions (%s), not `group`.",
        paste0("`", dims, "`", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  requests <- lapply(suffixes, function(x) dims[paste0(dims, x) %in% by])
  used <- unique(c(extra, dims[dims %in% unlist(requests)]))
  if (!length(used)) {
    return(none)
  }
  combined <- combine_dimensions(population, used)
  view <- population
  view$persons$group <- combined$current
  view$groups <- list(group = group())
  view$transitions <- empty_transitions()
  view["ancestry"] <- list(NULL)
  view[[".hammel_combined"]] <- list(
    values = seq_len(combined$count) - 1L,
    history = combined$history
  )
  # Internal `by`: one combined column per suffix, where its first dimension
  # was.
  # Positions, not names, so a dimension named `group` is replaced correctly.
  internal <- by
  drop <- rep(FALSE, length(by))
  for (k in seq_along(suffixes)) {
    hits <- which(by %in% paste0(requests[[k]], suffixes[k]))
    if (!length(hits)) {
      next
    }
    internal[hits[1]] <- paste0("group", suffixes[k])
    drop[hits[-1]] <- TRUE
  }
  internal <- internal[!drop]
  list(
    population = view,
    by = internal,
    rename = character(),
    external = by,
    dims = dims,
    suffixes = suffixes,
    requests = requests,
    combined = combined,
    extra = extra
  )
}

# A combined code over dimensions `dims` (the first fastest): each person's
# current code, and a merged history with one row per move in any of them,
# moves at one instant in declaration order.
combine_dimensions <- function(population, dims) {
  persons <- population$persons
  tx <- population$transitions
  tx <- if (is.null(tx)) {
    empty_transitions()
  } else {
    tx[tx$group %in% dims, , drop = FALSE]
  }
  categories <- lapply(
    stats::setNames(dims, dims), function(d) {
      x <- tx[tx$group == d, , drop = FALSE]
      extend_identity_groups(
        unique(persons[[d]]), c(
          as.character(x$origin),
          as.character(x$destination)
        )
      )
    }
  )
  sizes <- vapply(categories, length, integer(1))
  strides <- as.integer(cumprod(c(1, utils::head(sizes, -1))))
  code <- function(d, values) {
    match(as.character(values), as.character(categories[[d]])) - 1L
  }
  current <- Reduce(
    `+`, lapply(
      seq_along(dims), function(k) {
        code(dims[k], persons[[dims[k]]]) * strides[k]
      }
    )
  )
  history <- data.frame(
    person_id = numeric(),
    time = numeric(),
    from = integer(),
    to = integer()
  )
  if (nrow(tx)) {
    tx <- tx[order(tx$person_id, tx$time, match(tx$group, dims)), ,
      drop = FALSE
    ]
    n <- nrow(tx)
    first_row <- !duplicated(tx$person_id)
    # each row's person's first row
    start <- cummax(ifelse(first_row, seq_len(n), 0L))
    row_person <- match(tx$person_id, persons$person_id)
    after <- before <- 0L
    for (k in seq_along(dims)) {
      d <- dims[k]
      mine <- tx$group == d
      # The person's category in d before any move: the origin of their first
      # move in d, or their current category.
      first_d <- mine & !duplicated(ifelse(mine, tx$person_id, NA)) & mine
      initial <- code(d, persons[[d]])[row_person]
      first_origin <- tapply(
        code(d, tx$origin[first_d]), tx$person_id[first_d],
        function(x) x[1]
      )
      known <- match(as.character(tx$person_id), names(first_origin))
      initial[!is.na(known)] <- as.integer(first_origin[known[!is.na(known)]])
      # Category after each row: the latest destination in d so far.
      last <- cummax(ifelse(mine, seq_len(n), 0L))
      value_after <- ifelse(
        last >= start & last > 0L,
        code(d, tx$destination)[pmax(last, 1L)],
        initial
      )
      previous <- c(NA, value_after[-n])
      value_before <- ifelse(first_row, initial, previous)
      after <- after + value_after * strides[k]
      before <- before + value_before * strides[k]
    }
    history <- data.frame(
      person_id = tx$person_id,
      time = tx$time,
      from = as.integer(before),
      to = as.integer(after)
    )
  }
  decode <- as.data.frame(
    lapply(
      seq_along(dims), function(k) {
        categories[[k]][
          (seq_len(prod(sizes)) - 1L) %/% strides[k] %% sizes[k] + 1L
        ]
      }
    ),
    stringsAsFactors = FALSE,
    col.names = dims
  )
  list(
    current = as.integer(current),
    history = history,
    count = as.integer(prod(sizes)),
    categories = categories,
    decode = decode,
    dims = dims
  )
}

# Split combined columns back into one column per requested dimension, sum
# over dimensions combined but not requested (union rates), and order rows and
# columns by the user's `by`. `events`, `exposure` and `rate` name the count,
# exposure and rate columns to recompute after summing.
finish_strata <- function(
  out,
  alias,
  events = NULL,
  exposure = NULL,
  rate = NULL
) {
  if (is.null(alias$combined)) {
    return(rename_group_columns(out, alias$rename))
  }
  internal <- intersect(alias$by, names(out))
  block <- match(internal, names(out))
  before <- names(out)[seq_len(min(block) - 1L)]
  after <- setdiff(names(out)[-seq_len(max(block))], internal)
  decode <- alias$combined$decode
  for (k in seq_along(alias$suffixes)) {
    column <- paste0("group", alias$suffixes[k])
    if (!column %in% names(out)) {
      next
    }
    rows <- match(out[[column]], seq_len(nrow(decode)) - 1L)
    for (d in alias$requests[[k]]) {
      out[[paste0(d, alias$suffixes[k])]] <- decode[[d]][rows]
    }
  }
  keys <- c(before, alias$external)
  out <- out[c(keys, after)]
  dropped <- any(
    vapply(
      seq_along(alias$suffixes), function(k) {
        paste0("group", alias$suffixes[k]) %in% internal &&
          !setequal(alias$requests[[k]], alias$combined$dims)
      }, logical(1)
    )
  )
  if (dropped) {
    key <- do.call(paste, c(lapply(out[keys], as.character), sep = "\r"))
    first <- !duplicated(key)
    group <- match(key, key[first])
    summed <- out[first, , drop = FALSE]
    summed[[events]] <- as.vector(tapply(out[[events]], group, sum))
    summed[[exposure]] <- as.vector(tapply(out[[exposure]], group, sum))
    summed[[rate]] <- as.numeric(
      ifelse(
        summed[[exposure]] > 0, summed[[events]] / summed[[exposure]], NA_real_
      )
    )
    out <- summed
  }
  # The first `by` variable varies fastest; ties keep the calculator's order.
  order_keys <- lapply(
    rev(alias$external), function(b) {
      x <- as.character(out[[b]])
      dim <- sub("_[wh]$", "", b)
      levels <- if (dim %in% names(alias$combined$categories)) {
        as.character(alias$combined$categories[[dim]])
      } else {
        unique(x)
      }
      match(x, levels)
    }
  )
  out <- out[do.call(order, c(order_keys, list(seq_len(nrow(out))))), ,
    drop = FALSE
  ]
  tibble::as_tibble(out)
}

# Rename internal group columns to the dimension's own names.
rename_group_columns <- function(out, rename) {
  for (from in names(rename)) {
    names(out)[names(out) == from] <- rename[[from]]
  }
  out
}
