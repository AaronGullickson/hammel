matching_group_labels <- function(x, name) {
  check_groups(tibble::tibble(group = x), name)
  as.character(x)
}

expand_matching_pairs <- function(effects, symmetric = FALSE) {
  empty <- tibble::tibble(female_group = character(), male_group = character(), log_odds = numeric())
  if (is.null(effects)) return(empty)
  if (is.data.frame(effects)) {
    check_columns(effects, c("female_group", "male_group", "log_odds"), "group_exogamy")
    table <- tibble::tibble(
      female_group = matching_group_labels(effects$female_group, "group_exogamy$female_group"),
      male_group = matching_group_labels(effects$male_group, "group_exogamy$male_group"),
      log_odds = effects$log_odds)
  } else {
    labels <- names(effects)
    if (!is.numeric(effects) || is.complex(effects) || !is.null(dim(effects)) ||
        any(!is.finite(effects)) || (length(effects) &&
        (is.null(labels) || anyNA(labels) || any(!nzchar(labels)))))
      stop("`group_exogamy` must be a finite named numeric vector or a pair table.", call. = FALSE)
    tables <- lapply(seq_along(effects), function(i) {
      key <- labels[i]
      pattern <- if (symmetric) "^[^:]+:[^:]+$" else "^[^:]+:{1,2}[^:]+$"
      if (!grepl(pattern, key))
        stop(if (symmetric) "Use single ':' unordered pair names in `group_exogamy`." else
          "Use 'female:male' or symmetric 'group::group' names in `prevent_exogamy`.", call. = FALSE)
      pair <- strsplit(key, ":{1,2}")[[1]]
      if (any(!nzchar(trimws(pair)))) stop("Pair group labels must be nonempty.", call. = FALSE)
      table <- tibble::tibble(female_group = pair[1], male_group = pair[2],
                              log_odds = unname(effects[i]))
      if (!symmetric && grepl("::", key, fixed = TRUE) && pair[1] != pair[2])
        table <- rbind(table, tibble::tibble(female_group = pair[2], male_group = pair[1],
                                            log_odds = unname(effects[i])))
      table
    })
    table <- if (length(tables)) do.call(rbind, tables) else empty
  }
  if (!is.numeric(table$log_odds) || is.complex(table$log_odds) ||
      !is.null(dim(table$log_odds)) || any(!is.finite(table$log_odds)))
    stop("`group_exogamy` log odds must be finite numeric values.", call. = FALSE)
  if (symmetric && nrow(table)) {
    reversed <- table[table$female_group != table$male_group, , drop = FALSE]
    reversed[c("female_group", "male_group")] <- reversed[c("male_group", "female_group")]
    table <- rbind(table, reversed)
  }
  if (anyDuplicated(table[c("female_group", "male_group")]))
    stop(if (symmetric) "Conflicting or overlapping unordered pairs in `group_exogamy`." else
      "Conflicting or overlapping directed pairs in `prevent_exogamy`.", call. = FALSE)
  table
}

expand_group_exogamy <- function(effects) expand_matching_pairs(effects, symmetric = TRUE)

matching_has_group_formula <- function(matching) {
  any(grepl("^(female|male|same)_.+", all.vars(matching$formula)))
}

matching_predictor_data <- function(matching, gaps, female_group = NULL, male_group = NULL) {
  data <- data.frame(age_difference = gaps)
  if (!matching_has_group_formula(matching)) return(data)
  if (is.null(female_group) || is.null(male_group))
    stop("Matching group predictors require female and male group data.", call. = FALSE)
  female <- matching_group_labels(female_group, "female_group")
  male <- matching_group_labels(male_group, "male_group")
  if (!length(gaps)) return(data)
  if (!length(female) %in% c(1L, length(gaps)) || length(male) != length(gaps))
    stop("Supply one female group per search (or candidate) and one male group per candidate.", call. = FALSE)
  female <- rep(female, length.out = length(gaps))
  # A combined matching category (several dimensions): decode each
  # dimension's predictors from the partners' combinations.
  if (!is.null(matching$decode)) {
    f <- match(female, matching$group_levels)
    m <- match(male, matching$group_levels)
    if (anyNA(f) || anyNA(m)) stop("Matching group data contain unknown combinations.", call. = FALSE)
    for (d in names(matching$categories)) {
      levels <- as.character(matching$categories[[d]])
      data[[paste0("female_", d)]] <- factor(matching$decode[[d]][f], levels = levels)
      data[[paste0("male_", d)]] <- factor(matching$decode[[d]][m], levels = levels)
      data[[paste0("same_", d)]] <- as.numeric(matching$decode[[d]][f] == matching$decode[[d]][m])
    }
    return(data)
  }
  data$same_group <- as.numeric(female == male)
  if (any(all.vars(matching$formula) %in% c("female_group", "male_group"))) {
    levels <- as.character(matching$group_levels)
    if (any(!c(female, male) %in% levels)) stop("Matching group data contain undeclared `group_levels`.", call. = FALSE)
    data$female_group <- factor(female, levels = levels)
    data$male_group <- factor(male, levels = levels)
  }
  data
}

matching_pair_offsets <- function(matching, female_group, male_group) {
  pairs <- expand_group_exogamy(matching$group_exogamy)
  if (!nrow(pairs)) return(rep(0, length(male_group)))
  female <- matching_group_labels(female_group, "female_group")
  male <- matching_group_labels(male_group, "male_group")
  if (!length(female) %in% c(1L, length(male)))
    stop("Supply one female group per search (or candidate).", call. = FALSE)
  female <- rep(female, length.out = length(male))
  result <- numeric(length(male))
  for (i in seq_len(nrow(pairs))) {
    rows <- female == pairs$female_group[i] & male == pairs$male_group[i]
    result[rows] <- pairs$log_odds[i]
  }
  result
}

validate_matching_groups <- function(matching, population_groups) {
  if (is_model_schedule(matching)) {
    for (m in model_entries(matching)) validate_matching_groups(m, population_groups)
    return(invisible(matching))
  }
  if (length(matching$delta) && (is.null(population_groups) ||
      any(!names(matching$delta) %in% as.character(population_groups))))
    stop("`delta` names must belong to the population group identity registry.", call. = FALSE)
  pairs <- expand_group_exogamy(matching$group_exogamy)
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  if ((nrow(pairs) || matching_has_group_formula(matching) || restrictions$all ||
       nrow(restrictions$pairs)) && is.null(population_groups))
    stop("Group matching effects require population group identities.", call. = FALSE)
  if (any(all.vars(matching$formula) %in% c("female_group", "male_group")) &&
      any(!as.character(population_groups) %in% as.character(matching$group_levels)))
    stop("Population groups must be covered by matching `group_levels`.", call. = FALSE)
  invisible(matching)
}

# The internal model's restriction field (FALSE, TRUE or a directed pair
# table), built from NA entries of the public `group_exogamy` (D135).
expand_prevent_exogamy <- function(x) {
  empty <- expand_matching_pairs(NULL)[c("female_group", "male_group")]
  if (is.null(x)) x <- FALSE
  if (is.logical(x) && length(x) == 1L && is.null(dim(x)) && !is.na(x))
    return(list(all = x, pairs = empty))
  pairs <- tryCatch({
    if (is.character(x) && is.null(dim(x)) && !anyNA(x)) {
      expand_matching_pairs(stats::setNames(rep(0, length(x)), x))
    } else if (is.data.frame(x)) {
      check_columns(x, c("female_group", "male_group"), "prevent_exogamy")
      x$log_odds <- rep(0, nrow(x))
      expand_matching_pairs(x)
    } else stop("Use FALSE, TRUE, a character vector of pairs, or a pair table.")
  }, error = function(e) stop("Invalid `prevent_exogamy`: ", conditionMessage(e), call. = FALSE))
  if (any(pairs$female_group == pairs$male_group))
    stop("`prevent_exogamy` pairs must join different groups.", call. = FALSE)
  list(all = FALSE, pairs = pairs[c("female_group", "male_group")])
}

matching_group_allowed <- function(matching, female_group, male_group) {
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  female <- rep(as.character(female_group), length.out = length(male_group))
  male <- as.character(male_group)
  if (restrictions$all) return(female == male)
  allowed <- rep(TRUE, length(male))
  for (i in seq_len(nrow(restrictions$pairs))) {
    pair <- restrictions$pairs[i, ]
    allowed[female == pair$female_group & male == pair$male_group] <- FALSE
  }
  allowed
}

# Matching models name group predictors and settings by dimension (D134):
# `female_race`, `male_race`, `same_race`, and `group_exogamy`,
# and `delta` as lists keyed by dimension (a bare value
# means the only declared dimension). At simulation each model is translated
# to the internal one-dimension form used by the engine.
is_dimension_list <- function(x) is.list(x) && !is.data.frame(x)

matching_settings <- c("group_exogamy", "delta")

# Split one dimension's `group_exogamy` (D135) into finite pair effects and
# prohibitions: NA prohibits an unordered pair, and a bare NA every union
# between different categories. Returns the effects (NULL when none) and the
# internal restriction: FALSE, TRUE, or a directed table of prohibited pairs.
split_group_exogamy <- function(x) {
  out <- list(effects = NULL, prevent = FALSE)
  if (is.null(x)) return(out)
  if (!is.data.frame(x) && length(x) == 1L && is.null(names(x)) && is.atomic(x) && is.na(x)) {
    out$prevent <- TRUE
    return(out)
  }
  if (is.data.frame(x)) {
    check_columns(x, c("female_group", "male_group", "log_odds"), "group_exogamy")
    blocked <- is.na(x$log_odds)
    effects <- x[!blocked, , drop = FALSE]
    banned <- x[blocked, , drop = FALSE]
    banned$log_odds <- rep(0, nrow(banned))
    if (!nrow(effects)) effects <- NULL
  } else {
    if (!(is.numeric(x) || is.logical(x)) || is.complex(x) || !is.null(dim(x)) ||
        (length(x) && (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))))))
      stop(paste("`group_exogamy` must be a named vector of finite effects and NA prohibitions,",
                 "NA, or a pair table."), call. = FALSE)
    blocked <- is.na(x)
    effects <- if (any(!blocked)) stats::setNames(as.numeric(x[!blocked]), names(x)[!blocked])
    banned <- stats::setNames(rep(0, sum(blocked)), names(x)[blocked])
  }
  pairs <- expand_group_exogamy(effects)
  prohibited <- expand_group_exogamy(if (length(banned) || is.data.frame(banned)) banned)
  if (any(prohibited$female_group == prohibited$male_group))
    stop("`group_exogamy` cannot prohibit unions within a category (such as `\"A:A\" = NA`).", call. = FALSE)
  key <- function(t) paste(t$female_group, t$male_group, sep = "\r")
  if (any(key(prohibited) %in% key(pairs)))
    stop("Conflicting or overlapping unordered pairs in `group_exogamy`.", call. = FALSE)
  out$effects <- effects
  if (nrow(prohibited)) out$prevent <- as.data.frame(prohibited[c("female_group", "male_group")])
  out
}

# The dimensions a matching model's group effects use.
matching_dimensions <- function(model, dims) {
  vars <- grep("^(female|male|same)_.+", all.vars(model$formula), value = TRUE)
  used <- sub("^(female|male|same)_", "", vars)
  bare <- FALSE
  for (field in matching_settings) {
    x <- model[[field]]
    if (is_dimension_list(x)) used <- c(used, names(x))
    else if (!is.null(x) && !identical(x, FALSE) && length(x)) bare <- TRUE
  }
  if (bare) {
    if (length(dims) != 1L)
      stop(if (!length(dims)) "Group effects in matching require a declared group dimension; use set_groups()." else
        sprintf(paste("Name the dimension of matching group settings with a list, such as",
                      "`group_exogamy = list(%s = ...)`; this model declares %s."),
                dims[1], paste0("`", dims, "`", collapse = ", ")), call. = FALSE)
    used <- c(used, dims)
  }
  used <- unique(used)
  if (length(used) && !length(dims))
    stop("Group effects in matching require a declared group dimension; use set_groups().", call. = FALSE)
  unknown <- setdiff(used, dims)
  if (length(unknown))
    stop(sprintf("Matching uses %s, which is not a declared group dimension.",
                 paste0("`", unknown, "`", collapse = ", ")), call. = FALSE)
  used
}

# The internal model for dimension `dim` (NULL: no group effects), with
# categories as `group_levels`.
internal_matching_model <- function(model, dim, categories) {
  pick <- function(x) if (is_dimension_list(x)) x[[dim]] else x
  if (is.null(dim)) {
    internal <- structure(list(formula = model$formula, coefficients = model$coefficients,
      group_exogamy = NULL, group_levels = NULL, prevent_exogamy = FALSE, delta = NULL),
      class = "hammel_matching_model")
    validate_internal_matching(internal)
    return(internal)
  }
  public <- paste0(c("female_", "male_", "same_"), dim)
  private <- c("female_group", "male_group", "same_group")
  formula <- model$formula
  renamed <- stats::as.formula(call("~", do.call(substitute, list(formula[[2L]],
    stats::setNames(lapply(private, as.name), public)))), env = environment(formula))
  coefficients <- model$coefficients
  vars <- intersect(public, all.vars(formula))
  if (any(public[1:2] %in% vars) && length(categories) < 2L)
    stop("Categorical matching predictors require a group dimension with at least two categories.", call. = FALSE)
  if (length(vars)) {
    # Map coefficient names through the two model matrices' columns.
    labels <- as.character(categories)
    n <- max(2L, length(labels))
    data <- data.frame(age_difference = seq_len(n))
    data[[public[1]]] <- factor(rep(labels, length.out = n), levels = labels)
    data[[public[2]]] <- factor(rep(rev(labels), length.out = n), levels = labels)
    data[[public[3]]] <- rep(c(0, 1), length.out = n)
    original <- colnames(matching_design_matrix(formula, data))
    names(data)[match(public, names(data))] <- private
    translated <- colnames(matching_design_matrix(renamed, data))
    if (!setequal(original, names(coefficients)))
      stop(sprintf("Matching coefficient names must match the model-matrix columns: %s.",
                   paste(original, collapse = ", ")), call. = FALSE)
    names(coefficients) <- translated[match(names(coefficients), original)]
  }
  exogamy <- split_group_exogamy(pick(model$group_exogamy))
  internal <- structure(list(formula = renamed, coefficients = coefficients,
    group_exogamy = exogamy$effects,
    group_levels = if (any(c("female_group", "male_group") %in% all.vars(renamed))) categories,
    prevent_exogamy = exogamy$prevent,
    delta = pick(model$delta)), class = "hammel_matching_model")
  validate_internal_matching(internal)
  internal
}

# The internal model over combinations of several matching dimensions (step 2
# of D134): one table of summed pair effects over combinations, prohibited
# pairs (any dimension prohibiting them), `delta` summed per combination, and
# a decode table so formula predictors are evaluated from combinations.
# `categories` lists every matching dimension's categories (declaration
# order, first varying fastest); `labels` names the combinations.
combined_matching_model <- function(model, categories, labels) {
  dims <- names(categories)
  sizes <- vapply(categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  k <- length(labels)
  decode <- stats::setNames(lapply(seq_along(dims), function(j)
    as.character(categories[[j]])[(seq_len(k) - 1L) %/% strides[j] %% sizes[j] + 1L]), dims)
  own <- matching_dimensions(model, dims)
  pick <- function(x, d) if (is_dimension_list(x)) x[[d]] else NULL
  effects <- matrix(0, k, k)
  blocked <- matrix(FALSE, k, k)
  for (d in own) {
    exogamy <- split_group_exogamy(pick(model$group_exogamy, d))
    pairs <- expand_group_exogamy(exogamy$effects)
    for (r in seq_len(nrow(pairs))) {
      f <- decode[[d]] == as.character(pairs$female_group[r])
      m <- decode[[d]] == as.character(pairs$male_group[r])
      effects[f, m] <- effects[f, m] + pairs$log_odds[r]
    }
    if (isTRUE(exogamy$prevent)) blocked <- blocked | outer(decode[[d]], decode[[d]], `!=`)
    else if (is.data.frame(exogamy$prevent)) for (r in seq_len(nrow(exogamy$prevent))) {
      f <- decode[[d]] == as.character(exogamy$prevent$female_group[r])
      m <- decode[[d]] == as.character(exogamy$prevent$male_group[r])
      blocked[f, m] <- TRUE
    }
  }
  # Unordered pairs with an effect (the table is symmetric).
  upper <- which(upper.tri(effects, diag = TRUE) & effects != 0, arr.ind = TRUE)
  table <- if (nrow(upper)) data.frame(female_group = labels[upper[, 1]], male_group = labels[upper[, 2]],
                                       log_odds = effects[upper], stringsAsFactors = FALSE)
  off <- row(blocked) != col(blocked)
  prevent <- if (all(blocked[off])) TRUE else if (any(blocked)) {
    pairs <- which(blocked, arr.ind = TRUE)
    data.frame(female_group = labels[pairs[, 1]], male_group = labels[pairs[, 2]], stringsAsFactors = FALSE)
  } else FALSE
  delta <- numeric(k)
  for (d in own) {
    x <- pick(model$delta, d)
    if (is.null(x)) next
    unknown <- setdiff(names(x), as.character(categories[[d]]))
    if (length(unknown))
      stop(sprintf("`delta` names must be categories of `%s`; unknown: %s.", d,
                   paste(unknown, collapse = ", ")), call. = FALSE)
    delta <- delta + unname(x[decode[[d]]] %|na|% 0)
  }
  internal <- structure(list(formula = model$formula, coefficients = model$coefficients,
    group_exogamy = table, group_levels = labels, prevent_exogamy = prevent,
    delta = if (any(delta != 0)) stats::setNames(delta, labels),
    decode = decode, categories = categories), class = "hammel_matching_model")
  vars <- grep("^(female|male|same)_.+", all.vars(model$formula), value = TRUE)
  if (length(vars)) {
    categorical <- vars[grepl("^(female|male)_", vars)]
    small <- sub("^(female|male)_", "", categorical)
    if (any(sizes[small] < 2L))
      stop("Categorical matching predictors require a group dimension with at least two categories.", call. = FALSE)
    n <- max(2L, k)
    female <- labels[(seq_len(n) - 1L) %% k + 1L]
    male <- labels[(seq_len(n)) %% k + 1L]
    male[seq_len(n) %% 2L == 0L] <- female[seq_len(n) %% 2L == 0L]
    data <- matching_predictor_data(internal, seq_len(n), female, male)
    columns <- colnames(matching_design_matrix(model$formula, data))
    if (!setequal(columns, names(model$coefficients)))
      stop(sprintf("Matching coefficient names must match the model-matrix columns: %s.",
                   paste(columns, collapse = ", ")), call. = FALSE)
  }
  internal
}

`%|na|%` <- function(x, y) ifelse(is.na(x), y, x)

# Labels of combinations: categories joined with `|` (first dimension first).
combination_labels <- function(categories) {
  sizes <- vapply(categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  k <- prod(sizes)
  parts <- lapply(seq_along(categories), function(j)
    as.character(categories[[j]])[(seq_len(k) - 1L) %/% strides[j] %% sizes[j] + 1L])
  labels <- do.call(paste, c(parts, sep = "|"))
  if (anyDuplicated(labels)) labels <- paste0("combination_", seq_len(k))
  labels
}

# Translate a union process's matching models (and their schedules). The
# matching dimensions are every dimension any model uses (declaration order).
# With one, models use the one-dimension form; with several, the combined
# form over their combinations. Returns the process, the matching dimensions
# and the category labels the engine uses (NULL when there are none).
resolve_matching_models <- function(process, declaration, registries) {
  dims <- names(declaration)
  fields <- c("matching", "father_matching", "cohab_matching")
  used <- unique(unlist(lapply(fields, function(f)
    lapply(model_entries(process[[f]]), matching_dimensions, dims = dims))))
  used <- dims[dims %in% used]
  if (length(used) > 1L) {
    categories <- lapply(stats::setNames(used, used), function(d) registries[[d]])
    labels <- combination_labels(categories)
    convert <- function(model) combined_matching_model(model, categories, labels)
  } else {
    dim <- if (length(used)) used else NULL
    categories <- if (is.null(dim)) NULL else registries[[dim]]
    labels <- categories
    convert <- function(model) internal_matching_model(model,
      if (length(matching_dimensions(model, dims))) dim, categories)
  }
  for (f in fields) {
    x <- process[[f]]
    if (is.null(x)) next
    if (is_model_schedule(x)) x$model <- lapply(x$model, convert) else x <- convert(x)
    process[f] <- list(x)
  }
  list(process = process, dims = used, values = labels)
}
