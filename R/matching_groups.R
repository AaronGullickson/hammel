# Resolve group effects and category combinations for partner matching.

# INTERNAL CONSTANTS -------------------------------------------------------

matching_settings <- c("group_exogamy", "delta")

# HELPER FUNCTIONS ---------------------------------------------------------

matching_group_labels <- function(x, name, call = rlang::caller_env()) {
  check_group_labels(x, name, call = call)
  as.character(x)
}

expand_matching_pairs <- function(
  effects,
  symmetric = FALSE,
  call = rlang::caller_env()
) {
  setting <- if (symmetric) "group_exogamy" else "prevent_exogamy"
  empty <- tibble::tibble(
    female_group = character(),
    male_group = character(),
    log_odds = numeric()
  )
  if (is.null(effects)) {
    return(empty)
  }
  if (is.data.frame(effects)) {
    check_columns(
      effects, c("female_group", "male_group", "log_odds"), setting,
      call = call
    )
    table <- tibble::tibble(
      female_group = matching_group_labels(
        effects$female_group, paste0(setting, "$female_group"),
        call = call
      ),
      male_group = matching_group_labels(
        effects$male_group, paste0(setting, "$male_group"),
        call = call
      ),
      log_odds = effects$log_odds
    )
  } else {
    labels <- names(effects)
    if (!is.numeric(effects) || is.complex(effects) || !is.null(dim(effects)) ||
          any(!is.finite(effects)) ||
          (length(effects) &&
             (is.null(labels) || anyNA(labels) || any(!nzchar(labels))))) {
      abort_hammel(
        c(
          "{.arg {setting}} must be a named vector of finite effects.",
          "x" = "It is {describe_value(effects)}.",
          "i" = paste(
            "Name each effect by a pair of categories, such as",
            "{.code \"A:B\"}."
          )
        ),
        "specification",
        argument = setting,
        call = call
      )
    }
    tables <- lapply(
      seq_along(effects), function(i) {
        key <- labels[i]
        pattern <- if (symmetric) "^[^:]+:[^:]+$" else "^[^:]+:{1,2}[^:]+$"
        if (!grepl(pattern, key)) {
          abort_hammel(
            c(
              "{.arg {setting}} names must be pairs of categories.",
              "x" = "Invalid name {.val {key}}.",
              "i" = if (symmetric) {
                "Join the two categories with one colon, as in {.val A:B}."
              } else {
                "Use {.val female:male} or symmetric {.val A::B} names."
              }
            ),
            "specification",
            argument = setting,
            call = call
          )
        }
        pair <- strsplit(key, ":{1,2}")[[1]]
        if (any(!nzchar(trimws(pair)))) {
          abort_hammel(
            c(
              "{.arg {setting}} pairs must name two nonempty categories.",
              "x" = "Invalid name {.val {key}}."
            ),
            "specification",
            argument = setting,
            call = call
          )
        }
        table <- tibble::tibble(
          female_group = pair[1],
          male_group = pair[2],
          log_odds = unname(effects[i])
        )
        if (!symmetric && grepl("::", key, fixed = TRUE) &&
              pair[1] != pair[2]) {
          table <- rbind(
            table, tibble::tibble(
              female_group = pair[2],
              male_group = pair[1],
              log_odds = unname(effects[i])
            )
          )
        }
        table
      }
    )
    table <- if (length(tables)) do.call(rbind, tables) else empty
  }
  if (!is.numeric(table$log_odds) || is.complex(table$log_odds) ||
        !is.null(dim(table$log_odds)) ||
        any(!is.finite(table$log_odds))) {
    abort_hammel(
      c(
        "{.arg {setting}} log odds must be finite numbers.",
        "x" = "They are {describe_value(table$log_odds)}."
      ),
      "specification",
      argument = setting,
      call = call
    )
  }
  if (symmetric && nrow(table)) {
    reversed <- table[table$female_group != table$male_group, , drop = FALSE]
    reversed[c("female_group", "male_group")] <- reversed[c(
      "male_group",
      "female_group"
    )]
    table <- rbind(table, reversed)
  }
  duplicated_pairs <- duplicated(table[c("female_group", "male_group")])
  if (any(duplicated_pairs)) {
    repeated <- unique(
      paste0(
        table$female_group[duplicated_pairs], ":",
        table$male_group[duplicated_pairs]
      )
    )
    abort_hammel(
      c(
        "{.arg {setting}} gives some pairs more than once.",
        "x" = "Repeated: {.val {repeated}}.",
        "i" = if (symmetric) {
          "Pairs are unordered, so {.val A:B} and {.val B:A} are the same pair."
        }
      ),
      "specification",
      argument = setting,
      call = call
    )
  }
  table
}

expand_group_exogamy <- function(effects, call = rlang::caller_env()) {
  expand_matching_pairs(effects, symmetric = TRUE, call = call)
}

matching_has_group_formula <- function(matching) {
  any(grepl("^(female|male|same)_.+", all.vars(matching$formula)))
}

matching_predictor_data <- function(
  matching,
  gaps,
  female_group = NULL,
  male_group = NULL
) {
  data <- data.frame(age_difference = gaps)
  if (!matching_has_group_formula(matching)) {
    return(data)
  }
  if (is.null(female_group) || is.null(male_group)) {
    abort_internal("Matching group predictors need female and male groups.")
  }
  female <- matching_group_labels(female_group, "female_group")
  male <- matching_group_labels(male_group, "male_group")
  if (!length(gaps)) {
    return(data)
  }
  if (!length(female) %in% c(1L, length(gaps)) ||
        length(male) != length(gaps)) {
    abort_internal(
      "Matching needs one female group per search and one male per candidate."
    )
  }
  female <- rep(female, length.out = length(gaps))
  # A combined matching category (several dimensions): decode each
  # dimension's predictors from the partners' combinations.
  if (!is.null(matching$decode)) {
    f <- match(female, matching$group_levels)
    m <- match(male, matching$group_levels)
    if (anyNA(f) || anyNA(m)) {
      abort_internal("Matching group data contain unknown combinations.")
    }
    for (d in names(matching$categories)) {
      levels <- as.character(matching$categories[[d]])
      data[[paste0("female_", d)]] <- factor(
        matching$decode[[d]][f],
        levels = levels
      )
      data[[paste0("male_", d)]] <- factor(
        matching$decode[[d]][m],
        levels = levels
      )
      data[[paste0("same_", d)]] <- as.numeric(
        matching$decode[[d]][f] == matching$decode[[d]][m]
      )
    }
    return(data)
  }
  data$same_group <- as.numeric(female == male)
  if (any(all.vars(matching$formula) %in% c("female_group", "male_group"))) {
    levels <- as.character(matching$group_levels)
    if (any(!c(female, male) %in% levels)) {
      abort_internal("Matching group data contain undeclared group levels.")
    }
    data$female_group <- factor(female, levels = levels)
    data$male_group <- factor(male, levels = levels)
  }
  data
}

matching_pair_offsets <- function(matching, female_group, male_group) {
  pairs <- expand_group_exogamy(matching$group_exogamy)
  if (!nrow(pairs)) {
    return(rep(0, length(male_group)))
  }
  female <- matching_group_labels(female_group, "female_group")
  male <- matching_group_labels(male_group, "male_group")
  if (!length(female) %in% c(1L, length(male))) {
    abort_internal("Matching needs one female group per search.")
  }
  female <- rep(female, length.out = length(male))
  result <- numeric(length(male))
  for (i in seq_len(nrow(pairs))) {
    rows <- female == pairs$female_group[i] & male == pairs$male_group[i]
    result[rows] <- pairs$log_odds[i]
  }
  result
}

validate_matching_groups <- function(
  matching,
  population_groups,
  call = rlang::caller_env()
) {
  if (is_model_schedule(matching)) {
    for (m in model_entries(matching)) {
      validate_matching_groups(m, population_groups, call = call)
    }
    return(invisible(matching))
  }
  unknown <- setdiff(names(matching$delta), as.character(population_groups))
  if (length(matching$delta) && length(unknown)) {
    abort_hammel(
      c(
        "{.arg delta} names must be categories of the matching dimension.",
        "x" = "Unknown: {.val {unknown}}."
      ),
      "groups",
      categories = unknown,
      call = call
    )
  }
  pairs <- expand_group_exogamy(matching$group_exogamy, call = call)
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy, call = call)
  if ((nrow(pairs) || matching_has_group_formula(matching) ||
         restrictions$all ||
         nrow(restrictions$pairs)) &&
        is.null(population_groups)) {
    abort_hammel(
      c(
        "Group effects in matching need a declared group dimension.",
        "i" = "Declare dimensions with {.fn set_groups}."
      ),
      "groups",
      call = call
    )
  }
  if (any(all.vars(matching$formula) %in% c("female_group", "male_group")) &&
        any(!as.character(population_groups) %in%
            as.character(matching$group_levels)
        )) {
    abort_internal("Population groups are missing from matching group levels.")
  }
  invisible(matching)
}

# The internal model's restriction field (FALSE, TRUE or a directed pair
# table), built from NA entries of the public `group_exogamy` (D135).
expand_prevent_exogamy <- function(x, call = rlang::caller_env()) {
  empty <- expand_matching_pairs(NULL)[c("female_group", "male_group")]
  if (is.null(x)) {
    x <- FALSE
  }
  if (is.logical(x) && length(x) == 1L && is.null(dim(x)) && !is.na(x)) {
    return(list(all = x, pairs = empty))
  }
  pairs <- if (is.character(x) && is.null(dim(x)) && !anyNA(x)) {
    expand_matching_pairs(stats::setNames(rep(0, length(x)), x), call = call)
  } else if (is.data.frame(x)) {
    check_columns(
      x, c("female_group", "male_group"), "prevent_exogamy",
      call = call
    )
    x$log_odds <- rep(0, nrow(x))
    expand_matching_pairs(x, call = call)
  } else {
    abort_internal("Invalid internal exogamy restriction.")
  }
  if (any(pairs$female_group == pairs$male_group)) {
    abort_internal("Exogamy restrictions must join different groups.")
  }
  list(all = FALSE, pairs = pairs[c("female_group", "male_group")])
}

matching_group_allowed <- function(matching, female_group, male_group) {
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  female <- rep(as.character(female_group), length.out = length(male_group))
  male <- as.character(male_group)
  if (restrictions$all) {
    return(female == male)
  }
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

# Split one dimension's `group_exogamy` (D135) into finite pair effects and
# prohibitions: NA prohibits an unordered pair, and a bare NA every union
# between different categories. Returns the effects (NULL when none) and the
# internal restriction: FALSE, TRUE, or a directed table of prohibited pairs.
split_group_exogamy <- function(x, call = rlang::caller_env()) {
  out <- list(effects = NULL, prevent = FALSE)
  if (is.null(x)) {
    return(out)
  }
  if (!is.data.frame(x) && length(x) == 1L && is.null(names(x)) &&
        is.atomic(x) &&
        is.na(x)) {
    out$prevent <- TRUE
    return(out)
  }
  if (is.data.frame(x)) {
    if (!all(c("female_group", "male_group", "log_odds") %in% names(x))) {
      abort_internal(
        "Pair tables reach split_group_exogamy() only from tables."
      )
    }
    blocked <- is.na(x$log_odds)
    effects <- x[!blocked, , drop = FALSE]
    banned <- x[blocked, , drop = FALSE]
    banned$log_odds <- rep(0, nrow(banned))
    if (!nrow(effects)) {
      effects <- NULL
    }
  } else {
    if (!(is.numeric(x) || is.logical(x)) || is.complex(x) ||
          !is.null(dim(x)) ||
          (length(x) &&
             (is.null(names(x)) || anyNA(names(x)) ||
                any(!nzchar(names(x)))))) {
      abort_hammel(
        c(
          paste(
            "{.arg group_exogamy} must be a named vector of effects,",
            "{.code NA}, or a pair table."
          ),
          "x" = "It is {describe_value(x)}.",
          "i" = paste(
            "Name each finite effect or {.code NA} prohibition by a pair, as",
            "in {.code c(\"A:B\" = -2, \"A:C\" = NA)}."
          )
        ),
        "specification",
        argument = "group_exogamy",
        call = call
      )
    }
    blocked <- is.na(x)
    effects <- if (any(!blocked)) {
      stats::setNames(as.numeric(x[!blocked]), names(x)[!blocked])
    }
    banned <- stats::setNames(rep(0, sum(blocked)), names(x)[blocked])
  }
  pairs <- expand_group_exogamy(effects, call = call)
  prohibited <- expand_group_exogamy(
    if (length(banned) || is.data.frame(banned)) {
      banned
    },
    call = call
  )
  within <- prohibited$female_group == prohibited$male_group
  if (any(within)) {
    categories <- unique(prohibited$female_group[within])
    abort_hammel(
      c(
        "{.arg group_exogamy} can't prohibit unions within a category.",
        "x" = "Prohibited within {.val {categories}}."
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  key <- function(t) paste(t$female_group, t$male_group, sep = ":")
  both <- intersect(key(prohibited), key(pairs))
  if (length(both)) {
    abort_hammel(
      c(
        "{.arg group_exogamy} gives some pairs both an effect and {.code NA}.",
        "x" = "Pairs: {.val {both}}."
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  out$effects <- effects
  if (nrow(prohibited)) {
    out$prevent <- as.data.frame(prohibited[c("female_group", "male_group")])
  }
  out
}

# The table form of `group_exogamy` (D139): `<dim>_x`/`<dim>_y` columns and
# `log_odds`, one unordered pair of combinations per row, NA prohibiting.
# Returns the table's dimensions (in column order).
exogamy_table_dimensions <- function(x, call = rlang::caller_env()) {
  columns <- setdiff(names(x), "log_odds")
  hint <- paste(
    "A {.arg group_exogamy} table has {.field <group>_x} and",
    "{.field <group>_y} columns for each dimension, such as {.field race_x}",
    "and {.field race_y}, and a {.field log_odds} column."
  )
  old <- intersect(c("female_group", "male_group"), names(x))
  if (length(old)) {
    abort_hammel(
      c(
        "{.arg group_exogamy} tables no longer use {.field {old}} columns.",
        "i" = hint
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  bad <- columns[!grepl("_[xy]$", columns)]
  if (!"log_odds" %in% names(x) || !length(columns) || length(bad)) {
    abort_hammel(
      c(
        "{.arg group_exogamy} has the wrong columns for a pair table.",
        "x" = "Its columns are {.field {names(x)}}.",
        "i" = hint
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  dims <- unique(sub("_[xy]$", "", columns))
  expected <- c(paste0(dims, "_x"), paste0(dims, "_y"))
  if (!setequal(columns, expected)) {
    absent <- setdiff(expected, columns)
    abort_hammel(
      c(
        "Each dimension in a {.arg group_exogamy} table needs both columns.",
        "x" = "Missing: {.field {absent}}.",
        "i" = hint
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  dims
}

# Validate a table; returns it with character categories.
validate_exogamy_table <- function(x, call = rlang::caller_env()) {
  dims <- exogamy_table_dimensions(x, call = call)
  for (column in setdiff(names(x), "log_odds")) {
    check_group_labels(
      x[[column]], paste0("group_exogamy$", column),
      call = call
    )
  }
  v <- x$log_odds
  if (!(is.numeric(v) || is.logical(v)) || is.complex(v) || !is.null(dim(v)) ||
        any(is.infinite(v))) {
    abort_hammel(
      c(
        paste(
          "Column {.field log_odds} of {.arg group_exogamy} must contain",
          "finite effects or {.code NA} prohibitions."
        ),
        "x" = "It is {describe_value(v)}."
      ),
      "specification",
      argument = "group_exogamy",
      call = call
    )
  }
  side <- function(s) {
    do.call(
      paste, c(
        lapply(dims, function(d) as.character(x[[paste0(d, "_", s)]])),
        sep = "|"
      )
    )
  }
  a <- side("x")
  b <- side("y")
  key <- ifelse(a <= b, paste(a, b, sep = ":"), paste(b, a, sep = ":"))
  rows <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  if (length(rows)) {
    abort_hammel(
      c(
        "{.arg group_exogamy} lists some pairs more than once.",
        "x" = "Repeated pairs at {format_rows(rows)}.",
        "i" = paste(
          "Pairs are unordered, so each pair may appear once, in either",
          "order."
        )
      ),
      "specification",
      argument = "group_exogamy",
      rows = rows,
      call = call
    )
  }
  rows <- which(is.na(v) & a == b)
  if (length(rows)) {
    abort_hammel(
      c(
        "{.arg group_exogamy} can't prohibit unions within a combination.",
        "x" = "Prohibited within a combination at {format_rows(rows)}."
      ),
      "specification",
      argument = "group_exogamy",
      rows = rows,
      call = call
    )
  }
  invisible(dims)
}

# The dimensions a matching model's group effects use.
matching_dimensions <- function(model, dims, call = rlang::caller_env()) {
  vars <- grep("^(female|male|same)_.+", all.vars(model$formula), value = TRUE)
  used <- sub("^(female|male|same)_", "", vars)
  bare <- FALSE
  for (field in matching_settings) {
    x <- model[[field]]
    if (field == "group_exogamy" && is.data.frame(x)) {
      used <- c(used, exogamy_table_dimensions(x, call = call))
    } else if (is_dimension_list(x)) {
      used <- c(used, names(x))
    } else if (!is.null(x) && !identical(x, FALSE) && length(x)) {
      bare <- TRUE
    }
  }
  if (bare) {
    if (!length(dims)) {
      abort_hammel(
        c(
          "Group effects in matching need a declared group dimension.",
          "i" = "Declare dimensions with {.fn set_groups}."
        ),
        "groups",
        call = call
      )
    }
    if (length(dims) != 1L) {
      abort_hammel(
        c(
          "Matching group settings must name their dimension.",
          "x" = "The model declares several dimensions: {.val {dims}}.",
          "i" = paste0(
            "Give settings as a list by dimension, such as ",
            "{.code group_exogamy = list(", dims[1], " = ...)}."
          )
        ),
        "groups",
        call = call
      )
    }
    used <- c(used, dims)
  }
  used <- unique(used)
  if (length(used) && !length(dims)) {
    abort_hammel(
      c(
        "Group effects in matching need a declared group dimension.",
        "x" = "Matching uses {.field {used}}.",
        "i" = "Declare dimensions with {.fn set_groups}."
      ),
      "groups",
      call = call
    )
  }
  unknown <- setdiff(used, dims)
  if (length(unknown)) {
    abort_hammel(
      c(
        "Matching uses {.field {unknown}}, which {?isn't/aren't} declared.",
        "i" = "Declared dimensions: {.val {dims}}."
      ),
      "groups",
      group = unknown,
      call = call
    )
  }
  used
}

# Report coefficient names that do not match a matching formula's model-matrix
# columns.
abort_matching_coefficients <- function(
  coefficients,
  columns,
  call = rlang::caller_env()
) {
  missing <- setdiff(columns, names(coefficients))
  extra <- setdiff(names(coefficients), columns)
  abort_hammel(
    c(
      "Matching {.arg coefficients} must be named for each model column.",
      "x" = if (length(missing)) "Missing: {.val {missing}}.",
      "x" = if (length(extra)) "Unknown: {.val {extra}}.",
      "i" = if (length(columns)) {
        "The model columns are {.val {columns}}."
      } else {
        "The formula has no model columns, so give no coefficients."
      }
    ),
    "specification",
    argument = "coefficients",
    call = call
  )
}

# The internal model for dimension `dim` (NULL: no group effects), with
# categories as `group_levels`.
internal_matching_model <- function(
  model,
  dim,
  categories,
  call = rlang::caller_env()
) {
  pick <- function(x) if (is_dimension_list(x)) x[[dim]] else x
  if (is.null(dim)) {
    internal <- structure(
      list(
        formula = model$formula,
        coefficients = model$coefficients,
        group_exogamy = NULL,
        group_levels = NULL,
        prevent_exogamy = FALSE,
        delta = NULL
      ),
      class = "hammel_matching_model"
    )
    validate_internal_matching(internal, call = call)
    return(internal)
  }
  public <- paste0(c("female_", "male_", "same_"), dim)
  private <- c("female_group", "male_group", "same_group")
  formula <- model$formula
  renamed <- stats::as.formula(
    call(
      "~", do.call(
        substitute, list(
          formula[[2L]],
          stats::setNames(lapply(private, as.name), public)
        )
      )
    ),
    env = environment(formula)
  )
  coefficients <- model$coefficients
  vars <- intersect(public, all.vars(formula))
  if (any(public[1:2] %in% vars) && length(categories) < 2L) {
    abort_hammel(
      c(
        paste(
          "{.code female_} and {.code male_} matching predictors need a",
          "dimension with at least two categories."
        ),
        "x" = "{.field {dim}} has {length(categories)} category.",
        "i" = "Use a {.code same_} predictor instead."
      ),
      "groups",
      group = dim,
      call = call
    )
  }
  if (length(vars)) {
    # Map coefficient names through the two model matrices' columns.
    labels <- as.character(categories)
    n <- max(2L, length(labels))
    data <- data.frame(age_difference = seq_len(n))
    data[[public[1]]] <- factor(rep(labels, length.out = n), levels = labels)
    data[[public[2]]] <- factor(
      rep(rev(labels), length.out = n),
      levels = labels
    )
    data[[public[3]]] <- rep(c(0, 1), length.out = n)
    original <- colnames(matching_design_matrix(formula, data))
    names(data)[match(public, names(data))] <- private
    translated <- colnames(matching_design_matrix(renamed, data))
    if (!setequal(original, names(coefficients))) {
      abort_matching_coefficients(coefficients, original, call = call)
    }
    names(coefficients) <- translated[match(names(coefficients), original)]
  }
  exogamy <- if (is.data.frame(model$group_exogamy)) {
    table <- model$group_exogamy
    split_group_exogamy(
      data.frame(
        female_group = as.character(table[[paste0(dim, "_x")]]),
        male_group = as.character(table[[paste0(dim, "_y")]]),
        log_odds = as.numeric(table$log_odds),
        stringsAsFactors = FALSE
      ),
      call = call
    )
  } else {
    split_group_exogamy(pick(model$group_exogamy), call = call)
  }
  internal <- structure(
    list(
      formula = renamed,
      coefficients = coefficients,
      group_exogamy = exogamy$effects,
      group_levels = if (
        any(c("female_group", "male_group") %in% all.vars(renamed))
      ) {
        categories
      },
      prevent_exogamy = exogamy$prevent,
      delta = pick(model$delta)
    ),
    class = "hammel_matching_model"
  )
  validate_internal_matching(internal, call = call)
  internal
}

# The internal model over combinations of several matching dimensions (step 2
# of D134): one table of summed pair effects over combinations, prohibited
# pairs (any dimension prohibiting them), `delta` summed per combination, and
# a decode table so formula predictors are evaluated from combinations.
# `categories` lists every matching dimension's categories (declaration
# order, first varying fastest); `labels` names the combinations.
combined_matching_model <- function(
  model,
  categories,
  labels,
  call = rlang::caller_env()
) {
  dims <- names(categories)
  sizes <- vapply(categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  k <- length(labels)
  decode <- stats::setNames(
    lapply(
      seq_along(dims), function(j) {
        as.character(categories[[j]])[
          (seq_len(k) - 1L) %/% strides[j] %% sizes[j] + 1L
        ]
      }
    ), dims
  )
  own <- matching_dimensions(model, dims)
  pick <- function(x, d) if (is_dimension_list(x)) x[[d]] else NULL
  effects <- matrix(0, k, k)
  blocked <- matrix(FALSE, k, k)
  if (is.data.frame(model$group_exogamy)) {
    # A full table over combinations of the table's dimensions; dimensions it
    # does not name do not matter. Rows naming absent categories are unused.
    table <- model$group_exogamy
    table_dims <- exogamy_table_dimensions(table)
    for (r in seq_len(nrow(table))) {
      f <- m <- rep(TRUE, k)
      for (d in table_dims) {
        f <- f & decode[[d]] == as.character(table[[paste0(d, "_x")]][r])
        m <- m & decode[[d]] == as.character(table[[paste0(d, "_y")]][r])
      }
      if (!any(f) || !any(m)) {
        next
      }
      if (is.na(table$log_odds[r])) {
        blocked[f, m] <- TRUE
        blocked[m, f] <- TRUE
      } else {
        effects[f, m] <- effects[f, m] + table$log_odds[r]
        if (!identical(f, m)) {
          effects[m, f] <- effects[m, f] + table$log_odds[r]
        }
      }
    }
  }
  for (d in if (is.data.frame(model$group_exogamy)) character() else own) {
    exogamy <- split_group_exogamy(pick(model$group_exogamy, d), call = call)
    pairs <- expand_group_exogamy(exogamy$effects, call = call)
    for (r in seq_len(nrow(pairs))) {
      f <- decode[[d]] == as.character(pairs$female_group[r])
      m <- decode[[d]] == as.character(pairs$male_group[r])
      effects[f, m] <- effects[f, m] + pairs$log_odds[r]
    }
    if (isTRUE(exogamy$prevent)) {
      blocked <- blocked | outer(decode[[d]], decode[[d]], `!=`)
    } else if (is.data.frame(exogamy$prevent)) {
      for (r in seq_len(nrow(exogamy$prevent))) {
        f <- decode[[d]] == as.character(exogamy$prevent$female_group[r])
        m <- decode[[d]] == as.character(exogamy$prevent$male_group[r])
        blocked[f, m] <- TRUE
      }
    }
  }
  # Unordered pairs with an effect (the table is symmetric).
  upper <- which(upper.tri(effects, diag = TRUE) & effects != 0, arr.ind = TRUE)
  table <- if (nrow(upper)) {
    data.frame(
      female_group = labels[upper[, 1]],
      male_group = labels[upper[, 2]],
      log_odds = effects[upper],
      stringsAsFactors = FALSE
    )
  }
  off <- row(blocked) != col(blocked)
  prevent <- if (all(blocked[off])) {
    TRUE
  } else if (any(blocked)) {
    pairs <- which(blocked, arr.ind = TRUE)
    data.frame(
      female_group = labels[pairs[, 1]],
      male_group = labels[pairs[, 2]],
      stringsAsFactors = FALSE
    )
  } else {
    FALSE
  }
  delta <- numeric(k)
  for (d in own) {
    x <- pick(model$delta, d)
    if (is.null(x)) {
      next
    }
    unknown <- setdiff(names(x), as.character(categories[[d]]))
    if (length(unknown)) {
      abort_hammel(
        c(
          "{.arg delta} names must be categories of {.field {d}}.",
          "x" = "Unknown: {.val {unknown}}."
        ),
        "groups",
        group = d,
        categories = unknown,
        call = call
      )
    }
    adjustment <- unname(x[decode[[d]]])
    delta <- delta + ifelse(is.na(adjustment), 0, adjustment)
  }
  internal <- structure(
    list(
      formula = model$formula,
      coefficients = model$coefficients,
      group_exogamy = table,
      group_levels = labels,
      prevent_exogamy = prevent,
      delta = if (any(delta != 0)) stats::setNames(delta, labels),
      decode = decode,
      categories = categories
    ),
    class = "hammel_matching_model"
  )
  vars <- grep("^(female|male|same)_.+", all.vars(model$formula), value = TRUE)
  if (length(vars)) {
    categorical <- vars[grepl("^(female|male)_", vars)]
    small <- sub("^(female|male)_", "", categorical)
    if (any(sizes[small] < 2L)) {
      abort_hammel(
        c(
          paste(
            "{.code female_} and {.code male_} matching predictors need a",
            "dimension with at least two categories."
          ),
          "i" = paste(
            "Use {.code same_} predictors for a dimension with one",
            "category."
          )
        ),
        "groups",
        call = call
      )
    }
    n <- max(2L, k)
    female <- labels[(seq_len(n) - 1L) %% k + 1L]
    male <- labels[(seq_len(n)) %% k + 1L]
    male[seq_len(n) %% 2L == 0L] <- female[seq_len(n) %% 2L == 0L]
    data <- matching_predictor_data(internal, seq_len(n), female, male)
    columns <- colnames(matching_design_matrix(model$formula, data))
    if (!setequal(columns, names(model$coefficients))) {
      abort_matching_coefficients(model$coefficients, columns, call = call)
    }
  }
  internal
}

# Labels of combinations: categories joined with `|` (first dimension first).
combination_labels <- function(categories) {
  sizes <- vapply(categories, length, integer(1))
  strides <- cumprod(c(1, utils::head(sizes, -1)))
  k <- prod(sizes)
  parts <- lapply(
    seq_along(categories), function(j) {
      as.character(categories[[j]])[
        (seq_len(k) - 1L) %/% strides[j] %% sizes[j] + 1L
      ]
    }
  )
  labels <- do.call(paste, c(parts, sep = "|"))
  if (anyDuplicated(labels)) {
    labels <- paste0("combination_", seq_len(k))
  }
  labels
}

# Translate a union process's matching models (and their schedules). The
# matching dimensions are every dimension any model uses (declaration order).
# With one, models use the one-dimension form; with several, the combined
# form over their combinations. Returns the process, the matching dimensions
# and the category labels the engine uses (NULL when there are none).
resolve_matching_models <- function(
  process,
  declaration,
  registries,
  call = rlang::caller_env()
) {
  dims <- names(declaration)
  fields <- c("matching", "father_matching", "cohab_matching")
  used <- unique(
    unlist(
      lapply(
        fields, function(f) {
          lapply(
            model_entries(process[[f]]), matching_dimensions,
            dims = dims,
            call = call
          )
        }
      )
    )
  )
  used <- dims[dims %in% used]
  if (length(used) > 1L) {
    categories <- lapply(
      stats::setNames(used, used), function(d) registries[[d]]
    )
    labels <- combination_labels(categories)
    convert <- function(model) {
      combined_matching_model(model, categories, labels, call = call)
    }
  } else {
    dim <- if (length(used)) used else NULL
    categories <- if (is.null(dim)) NULL else registries[[dim]]
    labels <- categories
    convert <- function(model) {
      internal_matching_model(
        model,
        if (length(matching_dimensions(model, dims, call = call))) dim,
        categories,
        call = call
      )
    }
  }
  for (f in fields) {
    x <- process[[f]]
    if (is.null(x)) {
      next
    }
    x <- with_error_context(
      if (is_model_schedule(x)) {
        x$model <- lapply(x$model, convert)
        x
      } else {
        convert(x)
      },
      paste0("In the ", f, " model.")
    )
    process[f] <- list(x)
  }
  list(process = process, dims = used, values = labels)
}

# Check every matching model's ancestry endogamy (D199) against the dimensions
# that track ancestry. Returns whether any model uses it.
check_matching_ancestry <- function(
  process,
  declaration,
  population,
  registries,
  call = rlang::caller_env()
) {
  tracked <- tracked_ancestry(declaration)
  categories <- lapply(
    stats::setNames(tracked, tracked), function(d) {
      columns <- names(population$ancestry[[d]])
      union(
        as.character(registries[[d]]),
        setdiff(columns, c("person_id", "unknown"))
      )
    }
  )
  used <- FALSE
  for (f in c("matching", "father_matching", "cohab_matching")) {
    for (model in model_entries(process[[f]])) {
      weights <- with_error_context(
        resolve_ancestry_endogamy(model, tracked, categories, call = call),
        paste0("In the ", f, " model.")
      )
      used <- used || length(weights) > 0L
    }
  }
  used
}

# A model's ancestry endogamy weights as a list keyed by dimension: a bare
# vector names the only dimension that tracks ancestry, and weights name that
# dimension's categories or `unknown`.
resolve_ancestry_endogamy <- function(
  model,
  tracked,
  categories,
  call = rlang::caller_env()
) {
  x <- model$ancestry_endogamy
  if (is.null(x)) {
    return(list())
  }
  if (!is_dimension_list(x)) {
    if (!length(tracked)) {
      abort_hammel(
        c(
          "{.arg ancestry_endogamy} needs a dimension that tracks ancestry.",
          "i" = "Track ancestry with {.code group(ancestry = TRUE)}."
        ),
        "groups",
        call = call
      )
    }
    if (length(tracked) > 1L) {
      abort_hammel(
        c(
          "{.arg ancestry_endogamy} must name its dimension.",
          "x" = "Several dimensions track ancestry: {.val {tracked}}.",
          "i" = paste0(
            "Give it as a list by dimension, such as ",
            "{.code ancestry_endogamy = list(", tracked[1], " = ...)}."
          )
        ),
        "groups",
        call = call
      )
    }
    x <- stats::setNames(list(x), tracked)
  }
  untracked <- setdiff(names(x), tracked)
  if (length(untracked)) {
    abort_hammel(
      c(
        paste(
          "{.arg ancestry_endogamy} names {.field {untracked}}, which",
          "{?doesn't/don't} track ancestry."
        ),
        "i" = if (length(tracked)) {
          "Dimensions tracking ancestry: {.val {tracked}}."
        } else {
          "Track ancestry with {.code group(ancestry = TRUE)}."
        }
      ),
      "groups",
      group = untracked,
      call = call
    )
  }
  for (d in names(x)) {
    allowed <- c(categories[[d]], "unknown")
    unknown <- setdiff(names(x[[d]]), allowed)
    if (length(unknown)) {
      abort_hammel(
        c(
          paste(
            "{.arg ancestry_endogamy} names categories {.field {d}} doesn't",
            "have."
          ),
          "x" = "Unknown: {.val {unknown}}.",
          "i" = "Ancestry categories: {.val {allowed}}."
        ),
        "groups",
        group = d,
        call = call
      )
    }
  }
  x
}
