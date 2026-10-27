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
  any(all.vars(matching$formula) %in% c("female_group", "male_group", "same_group"))
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
