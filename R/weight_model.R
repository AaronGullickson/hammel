# Construct, validate, and prepare models for selection within event cells.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Create a weight model for selection within rate cells
#'
#' @description
#'
#' Create a `weight_model` that will be used to generate weights for determining
#' which individuals/units experience an event.
#'
#' @details
#'
#' A weight model decides which member of a rate cell experiences an event,
#' without changing the cell's input intensity and thus the input rates. The
#' cell produces events at its input rate times its members; each event goes to
#' member `i` with probability `w_i / sum(w)`, where `w_i = exp(eta_i)`.
#'
#' With a formula, `eta_i` is the linear predictor `x_i %*% beta`, so
#' `exp(beta)` is the hazard ratio between two members of the same cell that
#' differ by one unit in a term. Users can therefore create proportional hazards
#' models, without affecting overall input rates. Population level effects can
#' be smaller where cells are small or homogeneous.
#'
#' Attach a weight model with the `weights` argument of [set_mortality()],
#' [set_fertility()], [set_union_dissolution()], [set_group_transition()] and
#' [set_union_formation()], or `weights_conversion` of [set_union_formation()].
#'
#' @section Variables: Group variables are named by the declared group
#'   dimensions (see [set_groups()]); for a dimension `race` with categories A
#'   and B.
#'
#'   Mortality, fertility and group transitions are person processes and may use
#'   `sex` (`"f"`, `"m"`), `race`, `union_type` (`"unpartnered"`, `"cohab"`,
#'   `"marriage"`), `partnered` (0/1), `parity` (births; 0 for men),
#'   `previous_union` (`"none"`, `"separated"`, `"widowed"`), `partner_race` (NA
#'   when unpartnered) and, when the dimension tracks ancestry, shares
#'   `ancestry_race_A` and `ancestry_race_unknown`. Group transitions may also
#'   use `destination`, the destination of the risk in the transition's
#'   dimension.
#'
#'   Dissolution and conversion are union processes and may use `race_w`,
#'   `race_h`, `parity_w`, `previous_union_w`, `previous_union_h`,
#'   `ancestry_race_A_w`, `ancestry_race_A_h`, `same_race` (0/1),
#'   `age_difference` (husband's age minus wife's age) and `union_type`
#'   (`"cohab"`, `"marriage"`).
#'
#'   The first listed level of each categorical variable is its reference; for a
#'   group dimension it is its first category (declared, or first in the
#'   population).
#'
#' @section Formulas: Terms are main effects of these variables and interactions
#'   (`:` or `*`); transformations such as `I()` are not supported (use `fun`).
#'   Categorical variables use treatment coding, and `coefficients` must be
#'   named exactly as the corresponding model-matrix columns, for example
#'   `raceB` or `union_typecohab`.
#'
#' @section Custom functions: Instead of using a formula, users can create a
#'   custom function with `fun`. `fun(members, bin)` receives a tibble with one
#'   row per eligible member of the chosen cell, holding the variables above for
#'   every declared dimension (with labels for categorical variables) and
#'   identifiers (`person_id`, or `union_id`, `wife_id` and `husband_id`), and a
#'   list `bin` describing the cell: its `process`, strata (each dimension's
#'   category, or `<dim>_w`/`<dim>_h` for union processes, NA when the rates do
#'   not stratify by it; `origin` and `destination` for transitions), band
#'   (`age_start`/`age_end` or `dur_start`/`dur_end`) and the active rate period
#'   `time_start`/`time_end` (NA for static rates). Exact ages and times are not
#'   passed. In fertility with weighted matching calibration, the function is
#'   evaluated at every refresh of every woman, which can be slow.
#'
#' @section Time-varying weights: A process's `weights` may instead be a tibble
#'   with numeric `time_start`, `time_end` and a list column `model` of weight
#'   models, as for matching and inheritance model schedules.
#'
#' @param formula A one-sided formula of weight variables.
#' @param coefficients A finite named numeric vector of log hazard ratios.
#' @param fun A custom function returning log weights, instead of `formula`.
#' @return A `hammel_weight_model`. This specifies a model; it does not estimate
#'   one.
#' @export
#' @examples
#' # allow mortality to be 10% lower for individuals with partners
#' partnered_mortality <- weight_model(~partnered, c(partnered = log(0.9)))
#'
#' # allow dissolution rates to change by 5% per year increase in partner
#' # age differences - not in absolute terms which would require the custom
#' # function
#' gap_dissolution <- weight_model(
#'   ~age_difference,
#'   c(age_difference = log(1.05))
#' )
#' # like previous example but use custom function to get absolute age
#' # difference
#' gap_dissolution <- weight_model(
#'   fun = function(members, bin) {
#'     log(1.05) * abs(members$age_difference)
#'   }
#' )
weight_model <- function(formula = NULL, coefficients = numeric(), fun = NULL) {
  if (is.null(formula) == is.null(fun)) {
    stop("Supply exactly one of `formula` or `fun`.", call. = FALSE)
  }
  model <- structure(
    list(formula = formula, coefficients = coefficients, fun = fun),
    class = "hammel_weight_model"
  )
  validate_weight_model(model)
  model
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Engine codes of weight variables. Group variables (codes 1, 6, 9, 10, 16)
# carry their dimension's index as the factor component; ancestry variables
# (7, 14, 15) carry the ancestry column.
weight_variable_codes <- c(
  sex = 0L,
  group = 1L,
  union_type = 2L,
  partnered = 3L,
  parity = 4L,
  previous_union = 5L,
  partner_group = 6L,
  ancestry = 7L,
  destination = 8L,
  group_w = 9L,
  group_h = 10L,
  parity_w = 11L,
  previous_union_w = 12L,
  previous_union_h = 13L,
  ancestry_w = 14L,
  ancestry_h = 15L,
  same_group = 16L,
  age_difference = 17L
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Built-in weight variables of a process type, with engine codes and, for
# categorical variables, their levels (the first is the reference) and the
# engine code of each level.
weight_builtins <- function(type) {
  previous <- list(
    levels = c("none", "separated", "widowed"),
    level_codes = c(none = 0L, separated = 1L, widowed = 2L)
  )
  if (type == "union") {
    return(
      list(
        union_type = list(
          code = 2L,
          levels = c("cohab", "marriage"),
          level_codes = c(cohab = 1L, marriage = 2L)
        ),
        parity_w = list(code = 11L),
        previous_union_w = c(list(code = 12L), previous),
        previous_union_h = c(list(code = 13L), previous),
        age_difference = list(code = 17L)
      )
    )
  }
  list(
    sex = list(
      code = 0L,
      levels = c("f", "m"),
      level_codes = c(f = 1L, m = 0L)
    ),
    union_type = list(
      code = 2L,
      levels = c("unpartnered", "cohab", "marriage"),
      level_codes = c(unpartnered = 0L, cohab = 1L, marriage = 2L)
    ),
    partnered = list(code = 3L),
    parity = list(code = 4L),
    previous_union = c(list(code = 5L), previous)
  )
}

# Resolve a weight variable for a process type ("person", "transition",
# "union"): its engine code, component and levels, or NULL when unavailable.
# Group variables are named by declared dimension (D134), so they resolve only
# with `context` (at simulation): `categories` (each dimension's categories,
# named by dimension), `ancestry` (ancestry variable names, in engine column
# order) and `transition_dim`.
weight_variable <- function(name, type, context = NULL) {
  builtin <- weight_builtins(type)[[name]]
  if (!is.null(builtin)) {
    return(builtin)
  }
  if (is.null(context)) {
    return(NULL)
  }
  dims <- names(context$categories)
  group_info <- function(code, d) {
    list(
      code = code,
      component = match(d, dims) - 1L,
      levels = as.character(context$categories[[d]])
    )
  }
  if (type != "union") {
    if (name %in% dims) {
      return(group_info(1L, name))
    }
    partner <- sub("^partner_", "", name)
    if (startsWith(name, "partner_") && partner %in% dims) {
      return(group_info(6L, partner))
    }
    if (type == "transition" && name == "destination") {
      return(group_info(8L, context$transition_dim))
    }
    k <- match(name, context$ancestry)
    if (!is.na(k)) {
      return(list(code = 7L, component = k - 1L))
    }
    return(NULL)
  }
  same <- sub("^same_", "", name)
  if (startsWith(name, "same_") && same %in% dims) {
    return(list(code = 16L, component = match(same, dims) - 1L))
  }
  if (!grepl("_[wh]$", name)) {
    return(NULL)
  }
  base <- sub("_[wh]$", "", name)
  wife <- endsWith(name, "_w")
  if (base %in% dims) {
    return(group_info(if (wife) 9L else 10L, base))
  }
  k <- match(base, context$ancestry)
  if (!is.na(k)) {
    return(list(code = if (wife) 14L else 15L, component = k - 1L))
  }
  NULL
}

# Could `name` be a group variable for this process type, for some declared
# dimension? Dimension names follow check_group_names().
weight_group_name_possible <- function(name, type) {
  dim_ok <- function(x) {
    !x %in% reserved_group_names &&
      !grepl("^(partner|same|female|male|ancestry)_", x) &&
      !grepl("_(w|h)$", x) &&
      grepl("^[A-Za-z][A-Za-z0-9_.]*$", x)
  }
  if (type == "union") {
    base <- sub("_[wh]$", "", name)
    return(
      (grepl("_[wh]$", name) && (dim_ok(base) || grepl(
        "^ancestry_.+_.+",
        base
      ))) ||
        (startsWith(name, "same_") && dim_ok(sub("^same_", "", name)))
    )
  }
  dim_ok(name) || (startsWith(name, "partner_") && dim_ok(
    sub("^partner_", "", name)
  )) ||
    (grepl("^ancestry_.+_.+", name)) || 
      (type == "transition" && name == "destination")
}

validate_weight_model <- function(model) {
  if (!inherits(model, "hammel_weight_model")) {
    stop("`weights` must be created with weight_model().", call. = FALSE)
  }
  if (!is.null(model$fun)) {
    if (!is.function(model$fun)) {
      stop("`fun` must be a function.", call. = FALSE)
    }
    return(invisible(model))
  }
  formula <- model$formula
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("Use a one-sided weight formula.", call. = FALSE)
  }
  tt <- stats::terms(formula)
  if (attr(tt, "intercept") != 1L) {
    stop(
      "Weight formulas cannot remove the intercept; it has no effect.",
      call. = FALSE
    )
  }
  if (!is.null(attr(tt, "offset"))) {
    stop("Weight formula offsets are not supported.", call. = FALSE)
  }
  factors <- attr(tt, "factors")
  if (length(factors) && any(!rownames(factors) %in% all.vars(formula))) {
    stop(
      paste(
        "Weight formulas use variables and their interactions only;",
        paste0(
          "use a custom function (`fun`) for transformations ",
          "such as I()."
        )
      ),
      call. = FALSE
    )
  }
  coefficients <- model$coefficients
  if (!is.numeric(coefficients) || any(!is.finite(coefficients)) ||
        (length(coefficients) &&
           (is.null(names(coefficients)) || anyNA(names(coefficients)) ||
              any(!nzchar(names(coefficients))) ||
              anyDuplicated(names(coefficients))))) {
    stop(
      paste0(
        "`coefficients` must be a finite numeric vector ",
        "with unique nonempty names."
      ),
      call. = FALSE
    )
  }
  invisible(model)
}

# Decode a formula into model-matrix columns. Each column is a product of
# factors: a numeric variable, or an indicator of one level of a categorical
# variable. Follows R's coding (attr(terms, "factors"): 1 = contrasts,
# 2 = all levels) and checks the names against model.matrix().
# Without `context` (before the declaration is known), returns NULL when the
# formula uses variables that may be group variables, deferring the check.
weight_design <- function(model, type, context = NULL) {
  formula <- model$formula
  vars <- all.vars(formula)
  info <- lapply(
    stats::setNames(vars, vars), weight_variable,
    type = type,
    context = context
  )
  bad <- vars[vapply(info, is.null, logical(1))]
  # Before the declaration is known, defer names that could be group variables.
  if (length(bad) && is.null(context) &&
        all(vapply(bad, weight_group_name_possible, logical(1), type = type))) {
    return(NULL)
  }
  if (length(bad)) {
    stop(
      "Weight variables not available for this process: ", paste(
        bad,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  levels <- Filter(Negate(is.null), lapply(info, `[[`, "levels"))
  tt <- stats::terms(formula)
  labels <- attr(tt, "term.labels")
  factors <- attr(tt, "factors")
  columns <- list()
  for (j in seq_along(labels)) {
    in_term <- rownames(factors)[factors[, j] > 0]
    parts <- lapply(
      in_term, function(v) {
        if (!v %in% names(levels)) {
          return(list(list(var = v, level = NA_character_, name = v)))
        }
        lv <- levels[[v]]
        if (factors[v, j] == 1) {
          lv <- lv[-1]
        }
        lapply(lv, function(l) list(var = v, level = l, name = paste0(v, l)))
      }
    )
    grid <- expand.grid(lapply(parts, seq_along), KEEP.OUT.ATTRS = FALSE)
    for (r in seq_len(nrow(grid))) {
      chosen <- lapply(seq_along(parts), function(k) parts[[k]][[grid[r, k]]])
      columns[[length(columns) + 1L]] <- list(
        name = paste(vapply(chosen, `[[`, "", "name"), collapse = ":"),
        factors = chosen
      )
    }
  }
  # Check against R's own column names.
  data <- lapply(
    vars, function(v) {
      if (v %in% names(levels)) {
        factor(levels[[v]][1], levels = levels[[v]])
      } else {
        0
      }
    }
  )
  names(data) <- vars
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  categorical <- intersect(vars, names(levels))
  contrasts <- if (length(categorical)) {
    stats::setNames(
      rep(list("contr.treatment"), length(categorical)),
      categorical
    )
  } else {
    NULL
  }
  expected <- colnames(
    stats::model.matrix(formula, data, contrasts.arg = contrasts)
  )[-1]
  names_found <- vapply(columns, `[[`, "", "name")
  # Defensive: the decoding follows R's coding rules, so names should agree.
  if (!identical(names_found, expected)) {
    # nocov start
    stop(
      paste0(
        "This weight formula structure is not supported; ",
        "use a custom function (`fun`)."
      ),
      call. = FALSE
    )
  } # nocov end
  coefficients <- model$coefficients
  missing <- setdiff(expected, names(coefficients))
  extra <- setdiff(names(coefficients), expected)
  if (length(missing) || length(extra)) {
    stop(
      "Weight coefficients must be named for each model column: ",
      paste(expected, collapse = ", "),
      if (length(missing)) {
        paste0("; missing ", paste(missing, collapse = ", "))
      },
      if (length(extra)) paste0("; unknown ", paste(extra, collapse = ", ")),
      ".",
      call. = FALSE
    )
  }
  for (k in seq_along(columns)) {
    columns[[k]]$coefficient <- unname(coefficients[columns[[k]]$name])
  }
  columns
}

# Validate a weight model, or a schedule of weight models, for a process type.
# A schedule is a tibble with time_start, time_end and a model list column (as
# for matching and inheritance models); its models must all be formulas or all
# custom functions.
validate_process_weights <- function(weights, type) {
  if (is.null(weights)) {
    return(invisible(NULL))
  }
  one <- function(model) {
    validate_weight_model(model)
    if (is.null(model$fun)) {
      weight_design(model, type)
    }
    invisible(model)
  }
  validate_model_argument(weights, one, "weights")
  if (is_model_schedule(weights)) {
    kinds <- vapply(weights$model, function(m) is.null(m$fun), logical(1))
    if (length(unique(kinds)) > 1L) {
      stop(
        paste0(
          "A weight schedule's models must all be formulas ",
          "or all be custom functions."
        ),
        call. = FALSE
      )
    }
  }
  invisible(weights)
}

# Warn once for each formula term whose variables are all rate strata of the
# process: the term is constant within every cell and has no effect.
warn_cancelled_weight_terms <- function(weights, strata, process) {
  if (is_model_schedule(weights)) {
    for (model in weights$model) {
      warn_cancelled_weight_terms(model, strata, process)
    }
    return(invisible(NULL))
  }
  if (is.null(weights) || !is.null(weights$fun)) {
    return(invisible(FALSE))
  }
  tt <- stats::terms(weights$formula)
  factors <- attr(tt, "factors")
  cancelled <- character()
  for (j in seq_along(attr(tt, "term.labels"))) {
    in_term <- rownames(factors)[factors[, j] > 0]
    if (all(in_term %in% strata)) {
      cancelled <- c(cancelled, attr(tt, "term.labels")[j])
    }
  }
  for (term in cancelled) {
    warning(
      sprintf(
        paste(
          "The %s weight term `%s` has no effect: the input rates already",
          paste0(
            paste0("stratify by it, and weights are ", "normalized within "),
            "each cell."
          )
        ), process, term
      ),
      call. = FALSE
    )
  }
  invisible(length(cancelled) > 0L)
}

# Rate strata of each process, as weight variables constant within a cell.
# Group columns are known only from the declaration; at construction any
# possible group column is taken to be a dimension. A transition's own
# dimension (`group`) is constant within its cells.
weight_strata <- function(process, rates, group = NULL) {
  n <- names(rates)
  strata <- candidate_strata(rates)
  wives <- sub("_w$", "", grep("_w$", strata, value = TRUE))
  husbands <- sub("_h$", "", grep("_h$", strata, value = TRUE))
  union <- c(
    paste0(wives, if (length(wives)) "_w"), paste0(
      husbands, if (
        length(husbands)) {
        "_h"
      }
    ),
    paste0("same_", intersect(wives, husbands))
  )
  switch(
    process,
    mortality = c("sex", strata),
    fertility = c(
      strata, if ("union_type" %in% n) c("union_type", "partnered"),
      if ("parity_start" %in% n) "parity"
    ),
    dissolution = c(if ("union_type" %in% n) "union_type", union),
    conversion = c("union_type", union),
    search = c(
      "sex", "partnered", "union_type", strata,
      if ("previous_union" %in% n) "previous_union", paste0("partner_", strata)
    ),
    transition = c(
      group, "destination", strata, if ("sex" %in% n) "sex",
      if ("union_type" %in% n) c("union_type", "partnered")
    )
  )
}

# Engine specification of a weight model: NULL, a custom function, or columns
# of (variable code, level code, component) factors with coefficients.
# `context` is described at weight_variable().
prepare_weight_model <- function(weights, type, process, context) {
  if (is_model_schedule(weights)) {
    rows <- order(weights$time_start)
    return(
      list(
        schedule = TRUE,
        process = process,
        time_start = weights$time_start[rows],
        time_end = weights$time_end[rows],
        models = lapply(
          weights$model[rows], prepare_weight_model,
          type = type,
          process = process,
          context = context
        )
      )
    )
  }
  if (!is.null(weights$fun)) {
    return(list(custom = TRUE, fun = weights$fun, process = process))
  }
  columns <- weight_design(weights, type, context)
  active <- Filter(function(x) x$coefficient != 0, columns)
  encode <- function(f) {
    info <- weight_variable(f$var, type, context)
    level <- -1L
    if (!is.na(f$level)) {
      level <- if (!is.null(info$level_codes)) {
        info$level_codes[[f$level]]
      } else {
        match(f$level, info$levels) - 1L
      }
    }
    c(info$code, level, if (is.null(info$component)) 0L else info$component)
  }
  list(
    custom = FALSE,
    process = process,
    coefficients = vapply(active, `[[`, 0, "coefficient"),
    factors = lapply(
      active, function(x) do.call(rbind, lapply(x$factors, encode))
    )
  )
}
