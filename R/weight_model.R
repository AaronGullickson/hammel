#' Create a weight model for selection within rate cells
#'
#' A weight model decides which member of a rate cell experiences an event,
#' without changing the cell's input intensity. The cell produces events at its
#' input rate times its members; each event goes to member `i` with probability
#' `w_i / sum(w)`, where `w_i = exp(eta_i)`. With a formula, `eta_i` is the
#' linear predictor `x_i %*% beta`, so `exp(beta)` is the hazard ratio between
#' two members of the same cell that differ by one unit in a term. Population
#' level effects can be smaller where cells are small or homogeneous.
#'
#' Attach a weight model with the `weights` argument of [set_mortality()],
#' [set_fertility()], [set_union_dissolution()], [set_group_transition()] and
#' [set_union_formation()] (searches), or `weights_conversion` of
#' [set_union_formation()].
#'
#' **Variables.** Mortality, fertility and group transitions are person
#' processes and may use `sex` (`"f"`, `"m"`), `group`, `union_type`
#' (`"unpartnered"`, `"cohab"`, `"marriage"`), `partnered` (0/1), `parity`
#' (births; 0 for men), `previous_union` (`"none"`, `"separated"`,
#' `"widowed"`), `partner_group` (NA when unpartnered) and ancestry shares
#' `ancestry_<group>` and `ancestry_unknown`. Group transitions may also use
#' `group_dest`, the destination of the risk. Dissolution and conversion are
#' union processes and may use `group_w`, `group_h`, `parity_w`,
#' `previous_union_w`, `previous_union_h`, `ancestry_<group>_w`,
#' `ancestry_<group>_h`, `same_group` (0/1), `age_difference` (husband's age
#' minus wife's age) and `union_type` (`"cohab"`, `"marriage"`). All variables
#' change only at events. The first listed level of each categorical variable
#' is its reference.
#'
#' **Formulas.** Terms are main effects of these variables and interactions
#' (`:` or `*`); transformations such as `I()` are not supported (use `fun`).
#' Categorical variables use treatment coding, and `coefficients` must be
#' named exactly as the corresponding model-matrix columns, for example
#' `groupB` or `union_typecohab`. One `group_levels` applies to every
#' group-valued variable (`group`, `partner_group`, `group_dest`, `group_w`,
#' `group_h`), and must include every population group. The intercept has no
#' effect. A term whose variable is missing (such as `partner_group` for an
#' unpartnered person) contributes zero, with a warning.
#'
#' **Custom functions.** `fun(members, bin)` receives a tibble with one row per
#' eligible member of the chosen cell, holding the variables above (with labels
#' for categorical variables) and identifiers (`person_id`, or `union_id`,
#' `wife_id` and `husband_id`), and a list `bin` describing the cell: its
#' `process`, strata, band (`age_start`/`age_end` or `dur_start`/`dur_end`) and
#' the active rate period `time_start`/`time_end` (NA for static rates). Exact
#' ages and times are not passed. It must return one finite log weight per row;
#' NA gives zero, with a warning. It shares the simulation's random-number
#' stream and cannot change simulation state. In fertility with weighted
#' matching calibration, the function is evaluated at every refresh of every
#' woman, which can be slow.
#'
#' @param formula A one-sided formula of weight variables.
#' @param coefficients A finite named numeric vector of log hazard ratios.
#' @param group_levels Ordered group labels, required when the formula uses a
#'   group-valued variable; the first is the reference.
#' @param fun A custom function returning log weights, instead of `formula`.
#'
#' **Calendar schedules.** A process's `weights` may instead be a tibble with
#' numeric `time_start`, `time_end` and a list column `model` of weight models,
#' as for matching and inheritance model schedules. Half-open intervals must
#' not overlap and must cover the run; all models in a schedule must be
#' formulas or all custom functions. The model active at an instant applies;
#' switching models changes only who is chosen within cells, not event clocks
#' or cell intensities. A formula with zero coefficients turns weighting off
#' for an interval.
#' @return A `hammel_weight_model`. This specifies a model; it does not estimate one.
#' @export
#' @examples
#' # Illustrative coefficients, not demographic defaults.
#' partnered_mortality <- create_weight_model(~ partnered, c(partnered = log(0.8)))
#' gap_dissolution <- create_weight_model(~ age_difference,
#'   c(age_difference = log(1.05)))
create_weight_model <- function(formula = NULL, coefficients = numeric(),
                                group_levels = NULL, fun = NULL) {
  if (is.null(formula) == is.null(fun))
    stop("Supply exactly one of `formula` or `fun`.", call. = FALSE)
  model <- structure(list(formula = formula, coefficients = coefficients,
    group_levels = group_levels, fun = fun), class = "hammel_weight_model")
  validate_weight_model(model)
  model
}

# Weight variables by process type, with engine codes and categorical levels.
weight_variable_codes <- c(sex = 0L, group = 1L, union_type = 2L, partnered = 3L,
  parity = 4L, previous_union = 5L, partner_group = 6L, ancestry = 7L, group_dest = 8L,
  group_w = 9L, group_h = 10L, parity_w = 11L, previous_union_w = 12L,
  previous_union_h = 13L, ancestry_w = 14L, ancestry_h = 15L, same_group = 16L,
  age_difference = 17L)
weight_person_variables <- c("sex", "group", "union_type", "partnered", "parity",
  "previous_union", "partner_group")
weight_union_variables <- c("group_w", "group_h", "parity_w", "previous_union_w",
  "previous_union_h", "same_group", "age_difference", "union_type")
weight_group_variables <- c("group", "partner_group", "group_dest", "group_w", "group_h")
weight_fixed_levels <- function(type) list(
  sex = c("f", "m"),
  union_type = if (type == "union") c("cohab", "marriage") else c("unpartnered", "cohab", "marriage"),
  previous_union = c("none", "separated", "widowed"),
  previous_union_w = c("none", "separated", "widowed"),
  previous_union_h = c("none", "separated", "widowed"))
# Engine codes of categorical levels: sex is 1 for female; union types 0
# unpartnered, 1 cohab, 2 marriage; previous categories 0 none, 1 separated,
# 2 widowed. Group levels map to engine group indices at simulation.
weight_level_codes <- list(sex = c(f = 1L, m = 0L),
  union_type = c(unpartnered = 0L, cohab = 1L, marriage = 2L),
  previous_union = c(none = 0L, separated = 1L, widowed = 2L))

weight_ancestry_variable <- function(name, type) {
  if (type == "union") grepl("^ancestry_.+_[wh]$", name) else grepl("^ancestry_.+", name)
}

validate_weight_model <- function(model) {
  if (!inherits(model, "hammel_weight_model"))
    stop("`weights` must be created with create_weight_model().", call. = FALSE)
  if (!is.null(model$fun)) {
    if (!is.function(model$fun)) stop("`fun` must be a function.", call. = FALSE)
    return(invisible(model))
  }
  formula <- model$formula
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("Use a one-sided weight formula.", call. = FALSE)
  tt <- stats::terms(formula)
  if (attr(tt, "intercept") != 1L)
    stop("Weight formulas cannot remove the intercept; it has no effect.", call. = FALSE)
  if (!is.null(attr(tt, "offset")))
    stop("Weight formula offsets are not supported.", call. = FALSE)
  factors <- attr(tt, "factors")
  if (length(factors) && any(!rownames(factors) %in% all.vars(formula)))
    stop(paste("Weight formulas use variables and their interactions only;",
               "use a custom function (`fun`) for transformations such as I()."), call. = FALSE)
  known <- c(weight_person_variables, weight_union_variables, "group_dest")
  bad <- all.vars(formula)[!all.vars(formula) %in% known & !grepl("^ancestry_.+", all.vars(formula))]
  if (length(bad))
    stop("Unknown weight variables: ", paste(bad, collapse = ", "), ".", call. = FALSE)
  coefficients <- model$coefficients
  if (!is.numeric(coefficients) || any(!is.finite(coefficients)) ||
      (length(coefficients) && (is.null(names(coefficients)) || anyNA(names(coefficients)) ||
       any(!nzchar(names(coefficients))) || anyDuplicated(names(coefficients)))))
    stop("`coefficients` must be a finite numeric vector with unique nonempty names.", call. = FALSE)
  if (any(all.vars(formula) %in% weight_group_variables)) {
    levels <- model$group_levels
    if (is.null(levels) || length(levels) < 2L || anyNA(levels) ||
        anyDuplicated(as.character(levels)) || any(!nzchar(as.character(levels))))
      stop("Weight formulas with group variables require at least two unique `group_levels`.",
           call. = FALSE)
  }
  invisible(model)
}

# Variables a formula may use for a process type ("person", "transition", "union").
weight_allowed <- function(name, type) {
  if (type == "union") name %in% weight_union_variables || weight_ancestry_variable(name, "union")
  else name %in% c(weight_person_variables, if (type == "transition") "group_dest") ||
    (weight_ancestry_variable(name, "person") && !grepl("_[wh]$", name))
}

# Decode a formula into model-matrix columns. Each column is a product of
# factors: a numeric variable, or an indicator of one level of a categorical
# variable. Follows R's coding (attr(terms, "factors"): 1 = contrasts,
# 2 = all levels) and checks the names against model.matrix().
weight_design <- function(model, type) {
  formula <- model$formula
  vars <- all.vars(formula)
  bad <- vars[!vapply(vars, weight_allowed, logical(1), type = type)]
  if (length(bad))
    stop("Weight variables not available for this process: ", paste(bad, collapse = ", "),
         ".", call. = FALSE)
  levels <- weight_fixed_levels(type)
  for (v in intersect(vars, weight_group_variables)) levels[[v]] <- as.character(model$group_levels)
  tt <- stats::terms(formula)
  labels <- attr(tt, "term.labels")
  factors <- attr(tt, "factors")
  columns <- list()
  for (j in seq_along(labels)) {
    in_term <- rownames(factors)[factors[, j] > 0]
    parts <- lapply(in_term, function(v) {
      if (!v %in% names(levels)) return(list(list(var = v, level = NA_character_, name = v)))
      lv <- levels[[v]]
      if (factors[v, j] == 1) lv <- lv[-1]
      lapply(lv, function(l) list(var = v, level = l, name = paste0(v, l)))
    })
    grid <- expand.grid(lapply(parts, seq_along), KEEP.OUT.ATTRS = FALSE)
    for (r in seq_len(nrow(grid))) {
      chosen <- lapply(seq_along(parts), function(k) parts[[k]][[grid[r, k]]])
      columns[[length(columns) + 1L]] <- list(
        name = paste(vapply(chosen, `[[`, "", "name"), collapse = ":"), factors = chosen)
    }
  }
  # Check against R's own column names.
  data <- lapply(vars, function(v) if (v %in% names(levels))
    factor(levels[[v]][1], levels = levels[[v]]) else 0)
  names(data) <- vars
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  categorical <- intersect(vars, names(levels))
  contrasts <- if (length(categorical))
    stats::setNames(rep(list("contr.treatment"), length(categorical)), categorical) else NULL
  expected <- colnames(stats::model.matrix(formula, data, contrasts.arg = contrasts))[-1]
  names_found <- vapply(columns, `[[`, "", "name")
  if (!identical(names_found, expected))
    stop("This weight formula structure is not supported; use a custom function (`fun`).",
         call. = FALSE)
  coefficients <- model$coefficients
  missing <- setdiff(expected, names(coefficients))
  extra <- setdiff(names(coefficients), expected)
  if (length(missing) || length(extra))
    stop("Weight coefficients must be named for each model column: ",
         paste(expected, collapse = ", "),
         if (length(missing)) paste0("; missing ", paste(missing, collapse = ", ")),
         if (length(extra)) paste0("; unknown ", paste(extra, collapse = ", ")), ".", call. = FALSE)
  for (k in seq_along(columns)) columns[[k]]$coefficient <- unname(coefficients[columns[[k]]$name])
  columns
}

# Validate a weight model, or a schedule of weight models, for a process type.
# A schedule is a tibble with time_start, time_end and a model list column (as
# for matching and inheritance models); its models must all be formulas or all
# custom functions.
validate_process_weights <- function(weights, type) {
  if (is.null(weights)) return(invisible(NULL))
  one <- function(model) {
    validate_weight_model(model)
    if (is.null(model$fun)) weight_design(model, type)
    invisible(model)
  }
  validate_model_argument(weights, one, "weights")
  if (is_model_schedule(weights)) {
    kinds <- vapply(weights$model, function(m) is.null(m$fun), logical(1))
    if (length(unique(kinds)) > 1L)
      stop("A weight schedule's models must all be formulas or all be custom functions.", call. = FALSE)
  }
  invisible(weights)
}
weight_models_of <- function(weights) {
  if (is.null(weights)) list() else if (is_model_schedule(weights)) weights$model else list(weights)
}

# Warn once for each formula term whose variables are all rate strata of the
# process: the term is constant within every cell and has no effect.
warn_cancelled_weight_terms <- function(weights, strata, process) {
  if (is_model_schedule(weights)) {
    for (model in weights$model) warn_cancelled_weight_terms(model, strata, process)
    return(invisible(NULL))
  }
  if (is.null(weights) || !is.null(weights$fun)) return(invisible(FALSE))
  tt <- stats::terms(weights$formula)
  factors <- attr(tt, "factors")
  cancelled <- character()
  for (j in seq_along(attr(tt, "term.labels"))) {
    in_term <- rownames(factors)[factors[, j] > 0]
    if (all(in_term %in% strata)) cancelled <- c(cancelled, attr(tt, "term.labels")[j])
  }
  for (term in cancelled)
    warning(sprintf(paste("The %s weight term `%s` has no effect: the input rates already",
                          "stratify by it, and weights are normalized within each cell."), process, term),
            call. = FALSE)
  invisible(length(cancelled) > 0L)
}

# Rate strata of each process, as weight variables constant within a cell.
# Cells stratify by group only when rates or group multipliers do.
weight_strata <- function(process, rates, multipliers = NULL) {
  n <- names(rates)
  grouped <- "group" %in% n || (!is.null(multipliers) && any(unlist(multipliers) != 1))
  switch(process,
    mortality = c("sex", if (grouped) "group"),
    fertility = c(if (grouped) "group",
                  if ("union_type" %in% n) c("union_type", "partnered"),
                  if ("parity_start" %in% n) "parity"),
    dissolution = c(if ("union_type" %in% n) "union_type", if ("group_pair" %in% n) "same_group"),
    conversion = c("union_type", if ("group" %in% n) "group_w"),
    search = c("sex", "partnered", "union_type", if (grouped) "group",
               if ("previous_union" %in% n) "previous_union", "partner_group"),
    transition = c("group", "group_dest", if ("sex" %in% n) "sex",
                   if ("union_type" %in% n) c("union_type", "partnered")))
}

# Engine specification of a weight model: NULL, a custom function, or columns
# of (variable code, level code, ancestry component) factors with coefficients.
prepare_weight_model <- function(weights, type, process, group_values, ancestry_labels,
                                 group_labels = group_values) {
  if (is.null(weights)) return(NULL)
  if (is_model_schedule(weights)) {
    rows <- order(weights$time_start)
    return(list(schedule = TRUE, process = process,
      time_start = weights$time_start[rows], time_end = weights$time_end[rows],
      models = lapply(weights$model[rows], prepare_weight_model, type = type, process = process,
                      group_values = group_values, ancestry_labels = ancestry_labels)))
  }
  if (!is.null(weights$fun)) return(list(custom = TRUE, fun = weights$fun, process = process))
  columns <- weight_design(weights, type)
  if (any(vapply(columns, function(x) any(vapply(x$factors, function(f)
      f$var %in% weight_group_variables, logical(1))), logical(1)))) {
    if (is.null(group_values))
      stop("Weight formulas with group variables require population groups.", call. = FALSE)
    missing <- setdiff(as.character(group_values), as.character(weights$group_levels))
    if (length(missing))
      stop("Weight `group_levels` must include every population group; missing ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  active <- Filter(function(x) x$coefficient != 0, columns)
  encode <- function(f) {
    v <- f$var
    if (grepl("^ancestry_", v)) {
      suffix <- if (type == "union") sub("^.*_([wh])$", "\\1", v) else ""
      label <- if (type == "union") sub("^ancestry_(.*)_[wh]$", "\\1", v) else sub("^ancestry_", "", v)
      column <- if (label == "unknown") "unknown" else paste0("group_", label)
      component <- match(column, ancestry_labels) - 1L
      if (is.na(component))
        stop("No ancestry share is recorded for weight variable `", v, "`.", call. = FALSE)
      code <- weight_variable_codes[[if (suffix == "w") "ancestry_w" else if (suffix == "h")
        "ancestry_h" else "ancestry"]]
      return(c(code, -1L, component))
    }
    level <- -1L
    if (!is.na(f$level)) {
      level <- if (v %in% weight_group_variables) match(f$level, as.character(group_values)) - 1L else
        weight_level_codes[[sub("_[wh]$", "", v)]][[f$level]]
      # A group level absent from the population never matches anyone.
      if (is.na(level)) level <- -2L
    }
    c(weight_variable_codes[[v]], level, 0L)
  }
  list(custom = FALSE, process = process,
       coefficients = vapply(active, `[[`, 0, "coefficient"),
       factors = lapply(active, function(x) do.call(rbind, lapply(x$factors, encode))),
       group_labels = as.character(group_labels))
}
