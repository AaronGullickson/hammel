# Shared support for the `by` argument of rate calculators (D133).

# Validate `by` against the variables a calculator accepts.
check_by <- function(by, allowed, fun) {
  if (is.null(by)) return(character())
  if (!is.character(by) || !is.null(dim(by)) || anyNA(by) || anyDuplicated(by))
    stop("`by` must be NULL or a character vector of distinct variable names.", call. = FALSE)
  unknown <- setdiff(by, allowed)
  if (length(unknown))
    stop(sprintf("`%s()` cannot split by %s; `by` accepts %s.", fun,
                 paste0("\"", unknown, "\"", collapse = ", "),
                 paste0("\"", allowed, "\"", collapse = ", ")), call. = FALSE)
  by
}

check_pool_dissolved <- function(pool_dissolved, by) {
  if (!is.logical(pool_dissolved) || length(pool_dissolved) != 1L ||
      !is.null(dim(pool_dissolved)) || is.na(pool_dissolved))
    stop("`pool_dissolved` must be TRUE or FALSE.", call. = FALSE)
  if (pool_dissolved && !"previous_union" %in% by)
    stop("`pool_dissolved = TRUE` requires \"previous_union\" in `by`.", call. = FALSE)
  invisible(pool_dissolved)
}

previous_labels <- function(pool_dissolved) {
  if (pool_dissolved) c("none", "dissolved") else c("none", "separated", "widowed")
}

# Persons, unions and group indices for the person-level engine. `population`
# may be a persons data frame when no union variable is needed.
person_observation <- function(population, by, end, needs_unions = FALSE,
                               requires_typed = FALSE, grouped = "group" %in% by,
                               validate = TRUE) {
  union_vars <- intersect(by, c("union_type", "previous_union"))
  if (inherits(population, "hammel_population")) {
    if (validate) validate_population(population)
    if (is.null(population$end_time) || end > population$end_time)
      stop("Observation requires a known endpoint and cannot extend beyond it.", call. = FALSE)
    persons <- population$persons
    unions <- population$unions
  } else {
    if (needs_unions || length(union_vars))
      stop(sprintf("Splitting by %s requires a hammel_population with union histories.",
                   paste0("\"", union_vars, "\"", collapse = ", ")), call. = FALSE)
    if (!is.data.frame(population))
      stop("`population` must be a hammel_population or persons data frame.", call. = FALSE)
    persons <- population
    hammel_population(persons)
    unions <- NULL
  }
  if (grouped && !declared_group(population))
    stop("Splitting by \"group\" requires a declared group dimension; see set_groups() or hammel_population(groups = ).",
         call. = FALSE)
  persons <- prepare_group_observations(persons, population, grouped)
  groups <- if (grouped) observation_groups(persons, TRUE) else NULL
  typed <- !is.null(unions) && "union_type" %in% names(unions)
  if (requires_typed && !typed)
    stop("Splitting by \"union_type\" requires typed union histories.", call. = FALSE)
  if (is.null(unions)) unions <- data.frame(union_id = numeric(), wife_id = numeric(),
    husband_id = numeric(), start_time = numeric(), end_time = numeric(),
    end_reason = character())
  list(persons = persons, unions = unions, typed = typed, groups = groups,
       group_indices = if (grouped) match(persons$group, groups) - 1L else rep(0L, nrow(persons)),
       counts_cohab = inherits(population, "hammel_population") &&
         isTRUE(population$simulation$cohabitation_counts_as_previous_union))
}

union_type_labels <- function(typed) {
  if (typed) c("unpartnered", "cohab", "marriage") else c("unpartnered", "partnered")
}

# Run the person-level engine. Returns the native counts and the strides of
# each dimension in its arrays, for mapping to output rows.
run_person_lexis <- function(obs, by, time_breaks, axis_breaks, axis = 0L,
                             parity_breaks = numeric(), sex = FALSE, group = "group" %in% by,
                             utype_filter = 0L, previous_mask = c(TRUE, TRUE, TRUE),
                             pool_dissolved = FALSE, exclusive = NULL) {
  utype <- "union_type" %in% by
  previous <- "previous_union" %in% by
  if (is.null(exclusive)) exclusive <- obs$typed || previous || utype_filter > 0L
  map <- if (pool_dissolved) c(0L, 1L, 1L) else 0:2
  options <- list(typed = obs$typed, exclusive = exclusive, counts_cohab = obs$counts_cohab,
    axis = as.integer(axis), utype_filter = as.integer(utype_filter),
    previous_mask = previous_mask, sex = sex, group = group, union_type = utype,
    previous_union = previous, previous_map = map)
  ng <- if (is.null(obs$groups)) 1L else length(obs$groups)
  counts <- person_lexis_cpp(obs$persons, obs$unions, obs$group_indices, ng,
    time_breaks, axis_breaks, if (length(parity_breaks)) parity_breaks else numeric(1), options)
  nt <- length(time_breaks) - 1L
  na <- length(axis_breaks) - 1L
  sizes <- c(t = nt, a = na, sex = if (sex) 2L else 1L,
    parity = if (length(parity_breaks) > 1L) length(parity_breaks) - 1L else 1L,
    group = if (group) ng else 1L, union_type = if (utype) length(union_type_labels(obs$typed)) else 1L,
    previous_union = if (previous) length(unique(map)) else 1L)
  strides <- cumprod(c(1, head(sizes, -1)))
  names(strides) <- names(sizes)
  counts$strides <- strides
  counts$cells <- prod(sizes)
  counts
}

# Output rows for dimensions listed fastest first. Each dimension is a list
# with `columns` (a data frame of output columns, one row per value), `code`
# (0-based native codes) and `stride` (multiplier in the native array), and
# optionally `key` (global sort codes, for combining blocks). Returns the
# label columns plus `.index` and one `.key` column per dimension.
build_rows <- function(dims) {
  sizes <- vapply(dims, function(d) length(d$code), integer(1))
  grid <- expand.grid(lapply(sizes, seq_len), KEEP.OUT.ATTRS = FALSE)
  parts <- lapply(seq_along(dims), function(k) {
    d <- dims[[k]]
    x <- d$columns[grid[[k]], , drop = FALSE]
    rownames(x) <- NULL
    x
  })
  out <- if (length(parts)) do.call(cbind, parts) else data.frame(row.names = seq_len(nrow(grid)))
  index <- numeric(nrow(grid))
  for (k in seq_along(dims)) {
    d <- dims[[k]]
    index <- index + d$code[grid[[k]]] * d$stride
    out[[paste0(".key", k)]] <- (if (is.null(d$key)) seq_along(d$code) else d$key)[grid[[k]]]
  }
  out$.index <- index + 1
  out
}

# Order rows combined from several blocks by their keys, the last dimension
# varying slowest.
order_rows_keep <- function(rows) {
  keys <- rev(grep("^\\.key", names(rows), value = TRUE))
  rows[do.call(order, unname(as.list(rows[keys]))), , drop = FALSE]
}

dim_of <- function(columns, code = seq_len(nrow(columns)) - 1L, stride = 0, key = NULL) {
  list(columns = as.data.frame(columns, stringsAsFactors = FALSE), code = code,
       stride = stride, key = key)
}

band_dim <- function(breaks, prefix, stride) {
  n <- length(breaks) - 1L
  columns <- data.frame(breaks[seq_len(n)], breaks[seq_len(n) + 1L])
  names(columns) <- paste0(prefix, c("_start", "_end"))
  dim_of(columns, stride = stride)
}

label_dim <- function(name, labels, stride, code = seq_along(labels) - 1L, key = NULL) {
  columns <- data.frame(labels, stringsAsFactors = FALSE)
  names(columns) <- name
  dim_of(columns, code, stride, key)
}

# Dimensions for the person-level `by` variables of one engine run.
person_by_dims <- function(by, obs, counts, pool_dissolved, previous_codes = NULL) {
  lapply(by, function(name) switch(name,
    sex = label_dim("sex", c("m", "f"), counts$strides[["sex"]]),
    group = label_dim("group", obs$groups, counts$strides[["group"]]),
    union_type = label_dim("union_type", union_type_labels(obs$typed), counts$strides[["union_type"]]),
    previous_union = {
      labels <- previous_labels(pool_dissolved)
      codes <- if (is.null(previous_codes)) seq_along(labels) - 1L else previous_codes
      label_dim("previous_union", labels[codes + 1L], counts$strides[["previous_union"]],
                code = codes, key = codes + 1L)
    }))
}

finish_rows <- function(rows) {
  rows <- rows[setdiff(names(rows), ".index")]
  rows <- tibble::as_tibble(rows)
  rows
}

# Group strata for union-level rates (dissolution, conversion): the wife's and
# husband's groups and, optionally, union type. Returns native indices, the
# native group count and the output dimensions for the `by` variables.
union_strata <- function(population, unions, by, stride, typed_split = "union_type" %in% by) {
  wives <- "group_w" %in% by
  husbands <- "group_h" %in% by
  mode <- if (wives && husbands) 2L else if (husbands) 1L else 0L
  identities <- NULL
  w <- h <- rep(0L, nrow(unions))
  if (wives || husbands) {
    if (inherits(population, "hammel_population")) {
      if (!declared_group(population))
        stop("Splitting by \"group_w\" or \"group_h\" requires a declared group dimension.", call. = FALSE)
      persons <- prepare_group_observations(population$persons, population, TRUE)
      identities <- observation_groups(persons, TRUE)
      w <- match(persons$group[match(unions$wife_id, persons$person_id)], identities) - 1L
      h <- match(persons$group[match(unions$husband_id, persons$person_id)], identities) - 1L
    } else {
      needed <- c(if (wives) "group_w", if (husbands) "group_h")
      missing <- setdiff(needed, names(unions))
      if (length(missing))
        stop(sprintf("Splitting a union table by %s requires union column(s) %s.",
                     paste0("\"", needed, "\"", collapse = " and "),
                     paste0("`", missing, "`", collapse = ", ")), call. = FALSE)
      values <- c(if (wives) as.character(unions$group_w), if (husbands) as.character(unions$group_h))
      if (anyNA(values) || any(!nzchar(trimws(values))))
        stop("Union `group_w` and `group_h` must contain nonmissing, nonempty labels.", call. = FALSE)
      identities <- sort(unique(values), method = "radix")
      if (wives) w <- match(as.character(unions$group_w), identities) - 1L
      if (husbands) h <- match(as.character(unions$group_h), identities) - 1L
    }
  }
  ng <- if (is.null(identities)) 1L else length(identities)
  combos <- if (mode == 2L) ng * ng else if (wives || husbands) ng else 1L
  index <- if (mode == 2L) w + ng * h else if (husbands) h else w
  types <- c("cohab", "marriage")
  if (typed_split) {
    if (!"union_type" %in% names(unions))
      stop("Splitting by \"union_type\" requires a union `union_type` column.", call. = FALSE)
    if (!is.character(unions$union_type) || !is.null(dim(unions$union_type)) ||
        anyNA(unions$union_type) || any(!unions$union_type %in% types))
      stop("`union_type` must contain only cohab or marriage without missing values.", call. = FALSE)
    index <- index + combos * (match(unions$union_type, types) - 1L)
  }
  if ((wives || husbands) && inherits(population, "hammel_population"))
    unions <- prepare_union_group_observations(unions, population, mode, typed_split, ng, combos)
  dims <- lapply(by, function(name) switch(name,
    group_w = label_dim("group_w", identities, stride),
    group_h = label_dim("group_h", identities, stride * if (mode == 2L) ng else 1L),
    union_type = label_dim("union_type", types, stride * combos)))
  list(unions = unions, index = as.integer(index),
       count = combos * if (typed_split) 2L else 1L, dims = dims)
}
