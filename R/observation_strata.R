# Shared support for the `by` argument of rate calculators (D133).

# Validate `by` against the variables a calculator accepts.
# `dims` (declared group dimensions) replaces the internal `group` columns in
# the list of accepted names.
check_by <- function(by, allowed, fun, dims = NULL) {
  if (is.null(by)) return(character())
  if (!is.character(by) || !is.null(dim(by)) || anyNA(by) || anyDuplicated(by))
    stop("`by` must be NULL or a character vector of distinct variable names.", call. = FALSE)
  unknown <- setdiff(by, allowed)
  if (length(unknown) && length(dims)) {
    shown <- unlist(lapply(allowed, function(a) {
      suffix <- sub("^group", "", a)
      if (a %in% c("group", "group_w", "group_h")) paste0(dims, suffix) else a
    }))
    allowed <- shown
  }
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
                             pool_dissolved = FALSE, exclusive = NULL, group_target = NULL) {
  utype <- "union_type" %in% by
  previous <- "previous_union" %in% by
  if (is.null(exclusive)) exclusive <- obs$typed || previous || utype_filter > 0L
  map <- if (pool_dissolved) c(0L, 1L, 1L) else 0:2
  options <- list(typed = obs$typed, exclusive = exclusive, counts_cohab = obs$counts_cohab,
    axis = as.integer(axis), utype_filter = as.integer(utype_filter),
    previous_mask = previous_mask, sex = sex, group = group, union_type = utype,
    previous_union = previous, previous_map = map)
  if (!is.null(group_target)) options$group_target <- as.integer(group_target)
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

# The calculators split by one internal `group` column (`group_w`/`group_h`
# for union rates). `by` names declared dimensions (`"race"`, or
# `"race_w"`/`"race_h"`); with one dimension it maps directly to `group`. With
# several requested, `group` is a combined code over them (declaration order,
# first fastest), with a merged history of their moves, and finish_strata()
# splits it back into one column per dimension. `extra` adds dimensions to
# the combination (a transition's own dimension, first).
observation_alias <- function(population, by, suffixes = "", extra = character()) {
  none <- list(population = population, by = by, rename = character())
  if (!inherits(population, "hammel_population") || !length(population$groups)) return(none)
  dims <- names(population$groups)
  none$dims <- dims
  external <- as.vector(outer(dims, suffixes, paste0))
  if (length(dims) == 1L && !length(extra)) {
    dim <- dims
    if (dim != "group" && any(paste0("group", suffixes) %in% by))
      stop(sprintf("Split by the declared group dimension `%s`, not `group`.", dim), call. = FALSE)
    population <- dimension_view(population, dim)
    hit <- by %in% external
    if (any(hit)) by[hit] <- paste0("group", suffixes)[match(by[hit], paste0(dim, suffixes))]
    return(list(population = population, by = by, dims = dims,
                rename = stats::setNames(paste0(dim, suffixes), paste0("group", suffixes))))
  }
  if (!"group" %in% dims && any(paste0("group", suffixes) %in% by))
    stop(sprintf("Split by declared group dimensions (%s), not `group`.",
                 paste0("`", dims, "`", collapse = ", ")), call. = FALSE)
  requests <- lapply(suffixes, function(x) dims[paste0(dims, x) %in% by])
  used <- unique(c(extra, dims[dims %in% unlist(requests)]))
  if (!length(used)) return(none)
  combined <- combine_dimensions(population, used)
  view <- population
  view$persons$group <- combined$current
  view$groups <- list(group = group())
  view$transitions <- empty_transitions()
  view["ancestry"] <- list(NULL)
  view[[".hammel_combined"]] <- list(values = seq_len(combined$count) - 1L, history = combined$history)
  # Internal `by`: one combined column per suffix, where its first dimension was.
  internal <- by
  for (k in seq_along(suffixes)) {
    names_k <- paste0(requests[[k]], suffixes[k])
    if (!length(names_k)) next
    first <- min(match(names_k, internal))
    internal[first] <- paste0("group", suffixes[k])
    internal <- internal[!internal %in% names_k]
  }
  list(population = view, by = internal, rename = character(), external = by, dims = dims,
       suffixes = suffixes, requests = requests,
       combined = combined, extra = extra)
}

# A combined code over dimensions `dims` (the first fastest): each person's
# current code, and a merged history with one row per move in any of them,
# moves at one instant in declaration order.
combine_dimensions <- function(population, dims) {
  persons <- population$persons
  tx <- population$transitions
  tx <- if (is.null(tx)) empty_transitions() else tx[tx$group %in% dims, , drop = FALSE]
  categories <- lapply(stats::setNames(dims, dims), function(d) {
    x <- tx[tx$group == d, , drop = FALSE]
    extend_identity_groups(unique(persons[[d]]), c(as.character(x$origin), as.character(x$destination)))
  })
  sizes <- vapply(categories, length, integer(1))
  strides <- as.integer(cumprod(c(1, utils::head(sizes, -1))))
  code <- function(d, values) match(as.character(values), as.character(categories[[d]])) - 1L
  current <- Reduce(`+`, lapply(seq_along(dims), function(k) code(dims[k], persons[[dims[k]]]) * strides[k]))
  history <- data.frame(person_id = numeric(), time = numeric(), from = integer(), to = integer())
  if (nrow(tx)) {
    tx <- tx[order(tx$person_id, tx$time, match(tx$group, dims)), , drop = FALSE]
    n <- nrow(tx)
    first_row <- !duplicated(tx$person_id)
    start <- cummax(ifelse(first_row, seq_len(n), 0L))   # each row's person's first row
    row_person <- match(tx$person_id, persons$person_id)
    after <- before <- 0L
    for (k in seq_along(dims)) {
      d <- dims[k]
      mine <- tx$group == d
      # The person's category in d before any move: the origin of their first
      # move in d, or their current category.
      first_d <- mine & !duplicated(ifelse(mine, tx$person_id, NA)) & mine
      initial <- code(d, persons[[d]])[row_person]
      first_origin <- tapply(code(d, tx$origin[first_d]), tx$person_id[first_d], function(x) x[1])
      known <- match(as.character(tx$person_id), names(first_origin))
      initial[!is.na(known)] <- as.integer(first_origin[known[!is.na(known)]])
      # Category after each row: the latest destination in d so far.
      last <- cummax(ifelse(mine, seq_len(n), 0L))
      value_after <- ifelse(last >= start & last > 0L, code(d, tx$destination)[pmax(last, 1L)], initial)
      previous <- c(NA, value_after[-n])
      value_before <- ifelse(first_row, initial, previous)
      after <- after + value_after * strides[k]
      before <- before + value_before * strides[k]
    }
    history <- data.frame(person_id = tx$person_id, time = tx$time, from = as.integer(before),
                          to = as.integer(after))
  }
  decode <- as.data.frame(lapply(seq_along(dims), function(k)
    categories[[k]][(seq_len(prod(sizes)) - 1L) %/% strides[k] %% sizes[k] + 1L]),
    stringsAsFactors = FALSE, col.names = dims)
  list(current = as.integer(current), history = history, count = as.integer(prod(sizes)),
       categories = categories, decode = decode, dims = dims)
}

# Encode per-dimension labels (columns named `<dim><suffix>`) as combined codes.
encode_combined <- function(alias, data, suffix = "") {
  combined <- alias$combined
  sizes <- vapply(combined$categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  code <- 0
  for (k in seq_along(combined$dims)) {
    column <- paste0(combined$dims[k], suffix)
    if (!column %in% names(data)) return(rep(NA_integer_, nrow(data)))
    code <- code + (match(as.character(data[[column]]), as.character(combined$categories[[k]])) - 1L) *
      strides[k]
  }
  as.integer(code)
}

# Split combined columns back into one column per requested dimension, sum
# over dimensions combined but not requested (union rates), and order rows and
# columns by the user's `by`. `events`, `exposure` and `rate` name the count,
# exposure and rate columns to recompute after summing.
finish_strata <- function(out, alias, events = NULL, exposure = NULL, rate = NULL) {
  if (is.null(alias$combined)) return(rename_group_columns(out, alias$rename))
  internal <- intersect(alias$by, names(out))
  block <- match(internal, names(out))
  before <- names(out)[seq_len(min(block) - 1L)]
  after <- setdiff(names(out)[-seq_len(max(block))], internal)
  decode <- alias$combined$decode
  for (k in seq_along(alias$suffixes)) {
    column <- paste0("group", alias$suffixes[k])
    if (!column %in% names(out)) next
    rows <- match(out[[column]], seq_len(nrow(decode)) - 1L)
    for (d in alias$requests[[k]]) out[[paste0(d, alias$suffixes[k])]] <- decode[[d]][rows]
  }
  keys <- c(before, alias$external)
  out <- out[c(keys, after)]
  dropped <- any(vapply(seq_along(alias$suffixes), function(k)
    paste0("group", alias$suffixes[k]) %in% internal &&
      !setequal(alias$requests[[k]], alias$combined$dims), logical(1)))
  if (dropped) {
    key <- do.call(paste, c(lapply(out[keys], as.character), sep = "\r"))
    first <- !duplicated(key)
    group <- match(key, key[first])
    summed <- out[first, , drop = FALSE]
    summed[[events]] <- as.vector(tapply(out[[events]], group, sum))
    summed[[exposure]] <- as.vector(tapply(out[[exposure]], group, sum))
    summed[[rate]] <- as.numeric(ifelse(summed[[exposure]] > 0, summed[[events]] / summed[[exposure]], NA_real_))
    out <- summed
  }
  # The first `by` variable varies fastest; ties keep the calculator's order.
  order_keys <- lapply(rev(alias$external), function(b) {
    x <- as.character(out[[b]])
    dim <- sub("_[wh]$", "", b)
    levels <- if (dim %in% names(alias$combined$categories)) as.character(alias$combined$categories[[dim]]) else unique(x)
    match(x, levels)
  })
  out <- out[do.call(order, c(order_keys, list(seq_len(nrow(out))))), , drop = FALSE]
  tibble::as_tibble(out)
}

# Rename internal group columns to the dimension's own names.
rename_group_columns <- function(out, rename) {
  for (from in names(rename)) names(out)[names(out) == from] <- rename[[from]]
  out
}
