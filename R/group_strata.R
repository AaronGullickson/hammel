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
  registries <- stats::setNames(lapply(names(declaration), function(name) {
    categories <- declaration[[name]]$categories
    if (!is.null(categories)) return(categories)
    values <- unique(persons[[name]])
    # Transition histories and this dimension's transition rates.
    tx <- population$transitions
    if (!is.null(tx)) tx <- tx[as.character(tx$group) == name, , drop = FALSE]
    rates <- transition_rates[[name]]
    values <- extend_identity_groups(values, c(as.character(tx$origin),
      as.character(tx$destination), as.character(rates$origin), as.character(rates$destination)))
    shares <- population[["ancestry"]][[name]]
    if (!is.null(shares)) {
      known <- setdiff(names(shares), c("person_id", "unknown"))
      positive <- known[vapply(shares[known], function(x) any(x > 0), logical(1))]
      values <- extend_identity_groups(values, positive)
    }
    values
  }), names(declaration))
  attr(registries, "declared") <- vapply(declaration, function(g) !is.null(g$categories), logical(1))
  registries
}

# The group columns of a rate table: list of (dimension, role) pairs, with role
# 0 for a person or wife and 1 for a husband.
rate_strata_columns <- function(rates, registries, union_level = FALSE) {
  out <- list()
  for (d in seq_along(registries)) {
    name <- names(registries)[d]
    if (!union_level) {
      if (name %in% names(rates)) out[[length(out) + 1L]] <- list(column = name, dim = d, role = 0L)
    } else {
      for (role in 0:1) {
        column <- paste0(name, if (role) "_h" else "_w")
        if (column %in% names(rates)) out[[length(out) + 1L]] <- list(column = column, dim = d, role = role)
      }
    }
  }
  out
}

# Index each rate row by its stratum and check that every stratum is covered
# in each active calendar slice. Returns the rates (with `.hammel_group` when
# stratified) and the engine specification of the process's strata.
# Without `coverage`, uncovered strata have zero risk (group transitions).
prepare_strata_rates <- function(rates, process, registries, start = -Inf, end = Inf,
                                 union_level = FALSE, coverage = TRUE) {
  rates$.hammel_group <- NULL
  # Person-level tables use the dimension's name; union-level tables the
  # wife's and husband's columns.
  dims <- names(registries)
  wrong <- if (union_level) c(intersect(c(dims, "group_pair"), names(rates))) else
    intersect(c(paste0(dims, "_w"), paste0(dims, "_h")), names(rates))
  if (length(wrong))
    stop(sprintf(if (union_level)
      "%s rates use %s; union rates are stratified by the wife's and husband's groups, `<group>_w` and `<group>_h`." else
      "%s rates use %s; these rates are stratified by a group dimension's own name.", process,
      paste0("`", wrong, "`", collapse = ", ")), call. = FALSE)
  columns <- rate_strata_columns(rates, registries, union_level)
  spec <- list(dims = integer(), roles = integer(), strides = integer(), count = 1L)
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
    code <- match(as.character(rates[[x$column]]), as.character(registries[[x$dim]])) - 1L
    # Rows for categories absent from an inferred registry are unused; with
    # declared categories they are errors.
    if (anyNA(code) && isTRUE(declared[[x$dim]]))
      stop(sprintf("%s rates use `%s` values that are not categories of group `%s`: %s.", process,
                   x$column, names(registries)[x$dim],
                   paste(unique(as.character(rates[[x$column]])[is.na(code)]), collapse = ", ")),
           call. = FALSE)
    keep <- keep & !is.na(code)
    index <- index + code * strides[k]
  }
  rates <- rates[keep, , drop = FALSE]
  index <- index[keep]
  count <- as.integer(prod(counts))
  label <- function(k) paste(vapply(seq_along(columns), function(j) {
    value <- as.character(registries[[columns[[j]]$dim]])[(k %/% strides[j]) %% counts[j] + 1L]
    if (length(columns) == 1L) value else paste0(columns[[j]]$column, "=", value)
  }, ""), collapse = ", ")
  rates$.hammel_group <- index
  missing <- unique(unlist(lapply(calendar_rate_slices(rates, start, end), function(x)
    if (nrow(x)) setdiff(seq_len(count) - 1L, x$.hammel_group) else integer())))
  if (coverage && length(missing))
    stop(sprintf("%s rates are missing schedules for group(s): %s.", process,
                 paste(vapply(sort(missing), label, ""), collapse = "; ")), call. = FALSE)
  check_strata_overlap(rates, process)
  spec <- list(dims = vapply(columns, `[[`, integer(1), "dim") - 1L,
               roles = vapply(columns, `[[`, integer(1), "role"),
               strides = strides, count = count)
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
  for (d in seq_len(n))
    categories[, d] <- match(as.character(persons[[names(registries)[d]]]),
                             as.character(registries[[d]])) - 1L
  counts <- vapply(registries, length, integer(1))
  report <- list(dims = seq_len(n) - 1L, roles = integer(n),
                 strides = as.integer(cumprod(c(1, utils::head(counts, -1)))),
                 count = as.integer(prod(counts)))
  c(list(dimensions = n, categories = categories, counts = unname(counts),
         labels = unname(lapply(registries, as.character)), names = names(registries),
         matching = strata_spec(match(match_dims, names(registries)) - 1L, counts[match_dims]),
         report = report), specs)
}

# A person-level strata specification over dimensions (0-based) of given
# category counts.
strata_spec <- function(dims, counts) {
  counts <- as.integer(counts)
  list(dims = as.integer(dims), roles = integer(length(dims)),
       strides = as.integer(cumprod(c(1, utils::head(counts, -1))))[seq_along(dims)],
       count = as.integer(prod(counts)))
}

# Decode reporting strata (all dimensions) into one column per dimension.
decode_report_strata <- function(index, registries) {
  counts <- vapply(registries, length, integer(1))
  strides <- cumprod(c(1, utils::head(counts, -1)))
  stats::setNames(lapply(seq_along(registries), function(d)
    registries[[d]][(index %/% strides[d]) %% counts[d] + 1L]), names(registries))
}

# Within a stratum, bands must not overlap. Columns that are not declared group
# dimensions do not stratify rates, so a table split only by such a column
# (validated as separate strata at construction) is rejected here.
check_strata_overlap <- function(rates, process) {
  if (!nrow(rates)) return(invisible(rates))
  keys <- intersect(c(".hammel_group", "sex", "union_type", "previous_union", "parity_start",
                      "destination", "time_start"), names(rates))
  for (axis in c("age", "dur")) {
    lo <- paste0(axis, "_start"); hi <- paste0(axis, "_end")
    if (!lo %in% names(rates)) next
    rows <- which(!is.na(rates[[lo]]))
    key <- if (length(keys)) do.call(paste, c(lapply(rates[rows, keys, drop = FALSE], as.character), sep = "\r")) else
      rep("", length(rows))
    for (r in split(rows, key)) {
      r <- r[order(rates[[lo]][r])]
      if (length(r) > 1L && any(rates[[lo]][r[-1L]] < utils::head(cummax(rates[[hi]][r]), -1L)))
        stop(sprintf("%s rates have overlapping bands within a stratum; only declared group dimensions stratify rates.",
                     process), call. = FALSE)
    }
  }
  invisible(rates)
}
