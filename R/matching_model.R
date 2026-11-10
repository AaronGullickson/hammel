# Construct, validate, evaluate, and prepare partner-matching models.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#'Specify a conditional-logit matching model
#'
#'@description
#'
#'Specify conditional-logit model parameters for matching partners.
#'
#'
#'@details
#'
#'The `formula` argument can be used to define a formula for the utility
#'\eqn{U_{ij}} of the match between woman i and man j. This utility is defined
#'as the log of the unnormalized choice weight, such that the probability of
#'woman i choosing man j is given by:
#'
#' \deqn{
#'   P(j \mid i) = \frac{\exp(U_{ij})}{\sum_k \exp(U_{ik})}
#' }
#'
#'where the denominator sums across all eligible male partners. This equation is
#'equivalent to a conditional-logit model specification.
#'
#'In practice, when groups are in use, this formula needs to be modified to
#'adjust for expected male group shares, because the model uses women only to
#'make the choice:
#'
#' \deqn{
#'  P(j\mid i)= \frac{\omega_{g(j)}\exp(U_{ij})}{\sum_k\omega_{g(k)}\exp(U_{ik})}
#' }
#'
#'where \eqn{\omega_g} is a group-specific weight to adjust for male group
#'shares. More details are provided on the package website.
#'
#'When `delta = 0`, the \eqn{\omega} weights will adjust male group shares to
#'produce symmetry by gender in group exogamous unions. To allow for gender
#'asymmetry, users can adjust `delta` values. The `delta` value will typically
#'need to be tuned.
#'
#'Users can define formulas based on:
#'
#' - `age_difference`: male age minus female age
#' - `female_<group>`: group membership for `<group>` dimension of the woman.
#' - `male_<group>`: group membership for `<group>` dimension of the man.
#' - `same_<group>`: 0/1 numeric indicating same <group> membership for man and
#'woman.
#'
#'Formula transformations such as  `I(age_difference^2)` and
#'`abs(age_difference)` are supported through [stats::model.matrix()]. If users
#'limit the formula to `age_difference` and `I(age_difference^2)`, utilities
#'will be estimated in C++ leading to faster simulation runs.
#'
#'Coefficient names must match the non-intercept model-matrix columns exactly,
#'for example `female_raceB`. Names involving group predictors are checked when
#'the model meets the group declaration at simulation. Supply all coefficients
#'explicitly; no demographic coefficients are defaulted.
#'
#'Users may find it easier to define group exogamy by the use of the
#'`group_exogamy` argument than the use of a formula. This argument indicates
#'the likelihood of group exogamy across multiple group dimensions. It can take
#'one of two forms
#'
#' - **List form**: a named list keyed by dimension, such as
#'`list(race = c("A:B" = -2, "B:C" = NA), region = NA)`. Each value is a named
#'vector over unordered pairs (`"A:B"` applies to female A/male B and female
#'B/male A alike): a finite value is an additive log-weight effect and NA
#'prohibits the pair; `NA` alone prohibits every union between different
#'categories of the dimension. Effects from different dimensions add. Omitted
#'pairs are treated as implicitly zero. Cross-dimension effects are not possible
#'in this format.
#' - **Table form**, for cross-dimension effects: a data frame with `<group>_x`
#'and `<group>_y` columns for each dimension it uses and `log_odds`, such as
#'`data.frame(race_x = "W", region_x = "S", race_y = "B", region_y = "S",
#'log_odds = -5)`. Each row is one unordered pair of category combinations,
#'applied in both directions; NA prohibits it; unlisted pairs have an implicit
#'log odds of zero.
#'
#'Users may find [age_gap_matching_model()] more intuitive than specifying
#'a direct formula.
#'
#'@param formula A one-sided R formula using `age_difference` and group
#'  predictors `female_<group>`, `male_<group>` and `same_<group>`. Use `~ 1`
#'  and `numeric()` coefficients for pair-effects-only matching.
#'@param coefficients A finite named numeric vector of utility coefficients.
#'@param group_exogamy Either a list or data.frame indicating the log odds of
#'  group exogamous pairings. See details for syntax.
#'@param delta Optional male-target adjustments, as a named list keyed by
#'  dimension or a bare value for the only declared dimension: a finite numeric
#'  vector with unique category names.
#'@return An inspectable `hammel_matching_model` list containing `formula` and
#'  `coefficients`.
#'@export
#'@examples
#' # maximize age difference utility when husband is 5 years older. Double
#' # utility for same race unions
#' matching <- matching_model(
#'   ~ age_difference + I(age_difference^2) + same_race,
#'   c(
#'     age_difference = 0.1,
#'     "I(age_difference^2)" = -0.02,
#'     same_race = 0.69
#'   )
#' )
matching_model <- function(
  formula,
  coefficients,
  group_exogamy = NULL,
  delta = NULL
) {
  model <- structure(
    list(
      formula = formula,
      coefficients = coefficients,
      group_exogamy = group_exogamy,
      delta = delta
    ),
    class = "hammel_matching_model"
  )
  validate_matching_model(model)
  warn_absorbed_group_formula(model)
  model
}

#' Create a quadratic matching model from interpretable age-gap preferences
#'
#' @description
#'
#' Construct a matching model using age gap preferences rather than a formula.
#' Provides an intuitive front-end for [matching_model()].
#'
#' @details
#'
#' For a matching model that uses a quadratic formulation of age difference and
#' its square, the formula can be re-specified by the maximum value of the
#' parabola which gives the "preferred" age gap, and how quickly the utility
#' declines from that peak. If \eqn{x} is the age difference and \eqn{U} is the
#' utility defined as:
#'
#' \deqn{
#'  U = \beta_l x + \beta_q x^2
#' }
#'
#' then the value \eqn{x} where \eqn{U} is maximized is given by:
#'
#' \deqn{
#'  x_{max} = \frac{\beta_l}{-2\beta_q}
#' }
#'
#' Solving for \eqn{\beta_l} gives:
#'
#' \deqn{
#'  \beta_l = -2 * \beta_q x_{max}
#' }
#'
#' The degree of decline from the peak can be defined by how many years of age
#' difference \eqn{d} are necessary to reduce \eqn{U} by half. The formula for
#' \eqn{\beta_q} is given by:
#'
#' \deqn{
#'  \beta_q = \frac{-\log(2)}{d^2}
#' }
#'
#' Thus, the values of the linear and quadratic terms can be determined by the
#' desired values of \eqn{x_{max}} and \eqn{d}.
#'
#' @param preferred_age_gap Finite numeric scalar in years providing the
#'   preferred age gap where the age gap is defined by husband's age minus
#'   wife's age. Negative values allowed.
#' @param group_exogamy Optional group exogamy parameters passed to
#'   [matching_model()].
#' @param half_weight_distance Finite positive numeric scalar in years from the
#'   preferred age gap at which utility is halved.
#' @inheritParams matching_model
#' @return A `hammel_matching_model` with linear and quadratic age-difference
#'   terms.
#' @export
#' @examples
#' # preferred age gap of 3 years with a halving at 5 years, plus group
#' # exogamy terms. Exogamy allowed by race but not region.
#' matching <- age_gap_matching_model(
#'   preferred_age_gap = 3,
#'   half_weight_distance = 5,
#'   group_exogamy = list(
#'     race = c("W:B" = -3, "W:O" = -1, "B:O" = -2),
#'     region = NA
#'   )
#' )
age_gap_matching_model <- function(
  preferred_age_gap,
  half_weight_distance,
  group_exogamy = NULL,
  delta = NULL
) {
  for (name in c("preferred_age_gap", "half_weight_distance")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.null(dim(value)) ||
          !is.finite(value)) {
      stop("`", name, "` must be a finite numeric scalar.", call. = FALSE)
    }
  }
  if (half_weight_distance <= 0) {
    stop("`half_weight_distance` must be positive.", call. = FALSE)
  }
  quadratic <- -(sqrt(log(2)) / half_weight_distance)^2
  linear <- -2 * (preferred_age_gap * quadratic)
  if (!is.finite(quadratic) || quadratic >= 0 || !is.finite(linear)) {
    stop(
      paste0(
        "Inputs must yield finite coefficients with a negative ",
        "quadratic term."
      ),
      call. = FALSE
    )
  }
  matching_model(
    ~ age_difference + I(age_difference^2),
    c(
      age_difference = unname(linear),
      "I(age_difference^2)" = unname(quadratic)
    ),
    group_exogamy = group_exogamy,
    delta = delta
  )
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Male-group formula effects can be absorbed by calibrated target weights.
warn_absorbed_group_formula <- function(model) {
  male <- grep("^male_.+", all.vars(model$formula), value = TRUE)
  if (!length(male)) {
    return(invisible(FALSE))
  }
  warning(
    sprintf(
      paste(
        "Formula terms using %s are absorbed by calibrated group weights:",
        paste0(
          "main effects have no effect and interactions are ",
          "only partly reflected."
        ),
        "Use delta to adjust male-group targets."
      ), paste0("`", male, "`", collapse = ", ")
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

validate_delta <- function(x) {
  if (is.null(x)) {
    return(invisible(x))
  }
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
        any(!is.finite(x)) ||
        (length(x) &&
           (is.null(names(x)) || anyNA(names(x)) ||
              any(!nzchar(trimws(names(x)))) ||
              anyDuplicated(names(x))))) {
    stop(
      paste0(
        "`delta` must be a finite numeric vector with ",
        "unique nonempty group names."
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

# A public matching model: predictors and settings named by dimension. Checks
# that need the declaration (coefficient names of group predictors, setting
# categories) happen when it is translated at simulation.
validate_matching_model <- function(model) {
  if (!inherits(model, "hammel_matching_model") || !is.list(model)) {
    stop("`matching` must be created with matching_model().", call. = FALSE)
  }
  formula <- model$formula
  vars <- if (inherits(formula, "formula")) all.vars(formula) else character()
  if (!inherits(formula, "formula") || length(formula) != 2L ||
        any(
          vars != "age_difference" & !grepl("^(female|male|same)_.+", vars)
        )) {
    stop(
      paste(
        paste0(
          "Use a one-sided matching formula with ",
          "`age_difference` and group predictors"
        ),
        "`female_<group>`, `male_<group>` or `same_<group>`."
      ),
      call. = FALSE
    )
  }
  if (length(attr(stats::terms(formula), "offset"))) {
    stop("Matching formula offsets are not supported.", call. = FALSE)
  }
  beta <- model$coefficients
  if (!is.numeric(beta) || is.complex(beta) || !is.null(dim(beta)) ||
        any(!is.finite(beta)) ||
        (length(beta) && is.null(names(beta))) ||
        anyNA(names(beta)) ||
        any(!nzchar(names(beta))) ||
        anyDuplicated(names(beta))) {
    stop(
      paste0(
        "`coefficients` must be a finite numeric vector ",
        "with unique nonempty names."
      ),
      call. = FALSE
    )
  }
  if (is.data.frame(model$group_exogamy)) {
    validate_exogamy_table(model$group_exogamy)
  }
  for (field in c("group_exogamy", "delta")) {
    x <- model[[field]]
    if (field == "group_exogamy" && is.data.frame(x)) {
      next
    }
    values <- if (is_dimension_list(x)) {
      if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
            anyDuplicated(names(x))) {
        stop(
          sprintf(
            "A list `%s` must be named by group dimension, once each.",
            field
          ),
          call. = FALSE
        )
      }
      x
    } else {
      list(x)
    }
    for (v in values) {
      if (is.data.frame(v)) {
        stop(
          paste(
            paste0(
              "A dimension's `group_exogamy` entry is a named ",
              "vector or NA; give pair tables"
            ),
            paste0(
              "as the whole `group_exogamy` (with ",
              "`<group>_x`/`<group>_y` columns)."
            )
          ),
          call. = FALSE
        )
      }
      switch(
        field,
        group_exogamy = split_group_exogamy(v),
        delta = validate_delta(v)
      )
    }
  }
  # Without group predictors, coefficient names can be checked now.
  if (!any(grepl("^(female|male|same)_.+", vars))) {
    design <- matching_design_matrix(
      formula, data.frame(age_difference = seq(- 50, 50, length.out = 101))
    )
    if (!setequal(colnames(design), names(beta))) {
      stop(
        paste0(
          "Coefficient names must match all non-intercept ",
          "model-matrix columns."
        ),
        call. = FALSE
      )
    }
  }
  invisible(model)
}

# The internal one-dimension form (predictors `female_group`, `male_group`,
# `same_group` and explicit `group_levels`), built at simulation.
validate_internal_matching <- function(model) {
  validate_delta(model$delta)
  formula <- model$formula
  if (!inherits(formula, "formula") || length(formula) != 2L ||
        any(
          !all.vars(formula) %in%
             c("age_difference", "female_group", "male_group", "same_group")
        )) {
    stop(
      paste0(
        "Use a one-sided matching formula with ",
        "`age_difference`, `female_group`, `male_group`, ",
        "or `same_group`."
      ),
      call. = FALSE
    )
  }
  if (length(attr(stats::terms(formula), "offset"))) {
    stop("Matching formula offsets are not supported.", call. = FALSE)
  }
  beta <- model$coefficients
  if (!is.numeric(beta) || is.complex(beta) || !is.null(dim(beta)) ||
        any(!is.finite(beta)) ||
        (length(beta) && is.null(names(beta))) ||
        anyNA(names(beta)) ||
        any(!nzchar(names(beta))) ||
        anyDuplicated(names(beta))) {
    stop(
      paste0(
        "`coefficients` must be a finite numeric vector ",
        "with unique nonempty names."
      ),
      call. = FALSE
    )
  }
  expand_group_exogamy(model$group_exogamy)
  expand_prevent_exogamy(model$prevent_exogamy)
  levels <- model$group_levels
  if (!is.null(levels)) {
    labels <- matching_group_labels(levels, "group_levels")
    if (!length(labels) || anyDuplicated(labels)) {
      stop(
        "`group_levels` must be unique nonempty group labels.",
        call. = FALSE
      )
    }
  }
  categorical <- any(all.vars(formula) %in% c("female_group", "male_group"))
  if (categorical && length(levels) < 2L) {
    stop(
      paste0(
        "Categorical matching predictors require a group ",
        "dimension with at least two categories."
      ),
      call. = FALSE
    )
  }
  # A fixed validation grid covers both sexes' group levels, equality, and gaps.
  data <- data.frame(age_difference = seq(-50, 50, length.out = 101))
  if (matching_has_group_formula(model)) {
    labels <- if (categorical) as.character(levels) else c("A", "B")
    size <- max(101L, 2L * length(labels))
    female <- rep(labels, length.out = size)
    male <- rep(c(labels[-1L], labels[1L]), length.out = size)
    male[seq_len(size) %% 2L == 0L] <- female[seq_len(size) %% 2L == 0L]
    data <- matching_predictor_data(
      model, rep(data$age_difference, length.out = size), female, male
    )
  }
  design <- matching_design_matrix(formula, data)
  if (!setequal(colnames(design), names(beta))) {
    stop(
      "Coefficient names must match all non-intercept model-matrix columns.",
      call. = FALSE
    )
  }
  invisible(model)
}

matching_design_matrix <- function(formula, data) {
  design <- withCallingHandlers(
    {
      frame <- stats::model.frame(
        formula,
        data = data,
        na.action = stats::na.pass,
        drop.unused.levels = FALSE
      )
      factors <- frame[vapply(frame, is.factor, logical(1))]
      contrasts <- if (length(factors)) {
        lapply(
          factors, function(x)
            stats::contr.treatment(levels(x))
        )
      } else {
        NULL
      }
      stats::model.matrix(formula, data = frame, contrasts.arg = contrasts)
    },
    warning = function(w) {
      stop("Matching formula evaluation: ", conditionMessage(w), call. = FALSE)
    }
  )
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (nrow(design) != nrow(data) || any(!is.finite(design))) {
    stop(
      "Matching formula must produce finite values for every candidate.",
      call. = FALSE
    )
  }
  design
}

# Evaluate the same formula for R diagnostics and native fallback scoring.
calculate_matching_utilities <- function(
  matching,
  age_difference,
  female_group = NULL,
  male_group = NULL,
  add_group_exogamy = TRUE
) {
  check_numeric_column(age_difference, "age_difference")
  if (!length(age_difference)) {
    return(numeric())
  }
  data <- matching_predictor_data(
    matching, age_difference, female_group,
    male_group
  )
  design <- matching_design_matrix(matching$formula, data)
  if (!setequal(colnames(design), names(matching$coefficients))) {
    stop(
      "Matching formula columns changed for the supplied candidates.",
      call. = FALSE
    )
  }
  utility <- drop(design %*% matching$coefficients[colnames(design)])
  if (add_group_exogamy && nrow(expand_group_exogamy(matching$group_exogamy))) {
    if (is.null(female_group) || is.null(male_group)) {
      stop("Pair effects require female and male group data.", call. = FALSE)
    }
    utility <- utility + matching_pair_offsets(
      matching, female_group,
      male_group
    )
  }
  if (any(!is.finite(utility))) {
    stop("Matching utilities must be finite.", call. = FALSE)
  }
  unname(utility)
}

calculate_matching_probabilities <- function(
  matching,
  age_difference,
  female_group = NULL,
  male_group = NULL
) {
  utility <- calculate_matching_utilities(
    matching, age_difference,
    female_group, male_group
  )
  if (!length(utility)) {
    return(numeric())
  }
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  allowed <- rep(TRUE, length(utility))
  if (restrictions$all || nrow(restrictions$pairs)) {
    if (is.null(female_group) || is.null(male_group)) {
      stop(
        "Group restrictions require female and male group data.",
        call. = FALSE
      )
    }
    allowed <- matching_group_allowed(matching, female_group, male_group)
  }
  weights <- numeric(length(utility))
  if (any(allowed)) {
    weights[allowed] <- exp(utility[allowed] - max(utility[allowed]))
    weights <- weights / sum(weights)
  }
  weights
}

# Use an exact native path for the common linear/raw quadratic specification.
# All other formulas retain model.matrix semantics through an R evaluator.
prepare_matching_model <- function(matching, population_groups = NULL) {
  validate_matching_groups(matching, population_groups)
  beta <- matching$coefficients
  columns <- gsub("[[:space:]]", "", names(beta))
  native <- all(columns %in% c("age_difference", "I(age_difference^2)"))
  if ("I(age_difference^2)" %in% columns) {
    native <- native && identical(
      get("I", envir = environment(matching$formula), inherits = TRUE), base::I
    )
  }
  quadratic <- numeric()
  if (native) {
    quadratic <- c(
      sum(beta[columns == "age_difference"]),
      sum(beta[columns == "I(age_difference^2)"])
    )
  }
  pairs <- expand_group_exogamy(matching$group_exogamy)
  offsets <- matrix(numeric(), 0, 0)
  labels <- as.character(population_groups)
  if (nrow(pairs)) {
    offsets <- matrix(0, length(labels), length(labels))
    for (i in seq_len(nrow(pairs))) {
      female <- match(pairs$female_group[i], labels)
      male <- match(pairs$male_group[i], labels)
      if (!is.na(female) && !is.na(male)) {
        offsets[female, male] <- pairs$log_odds[i]
      }
    }
  }
  restrictions <- expand_prevent_exogamy(matching$prevent_exogamy)
  blocked <- matrix(FALSE, 0, 0)
  if (nrow(restrictions$pairs)) {
    blocked <- matrix(FALSE, length(labels), length(labels))
    for (i in seq_len(nrow(restrictions$pairs))) {
      female <- match(restrictions$pairs$female_group[i], labels)
      male <- match(restrictions$pairs$male_group[i], labels)
      if (!is.na(female) && !is.na(male)) {
        blocked[female, male] <- TRUE
      }
    }
  }
  group_formula <- matching_has_group_formula(matching)
  list(
    delta = unname(
      vapply(
        labels, function(g) {
          if (g %in% names(matching$delta)) matching$delta[[g]] else 0
        }, numeric(1)
      )
    ),
    prevent_exogamy_all = restrictions$all,
    blocked_pairs = blocked,
    quadratic = quadratic,
    has_squared = "I(age_difference^2)" %in% columns,
    group_formula = group_formula,
    needs_groups = group_formula || nrow(pairs) > 0,
    pair_offsets = offsets,
    evaluate = function(
      age_difference, female_group = NULL,
      male_group = NULL
    ) {
      calculate_matching_utilities(
        matching, age_difference,
        if (is.null(female_group)) NULL else labels[female_group + 1L],
        if (is.null(male_group)) NULL else labels[male_group + 1L],
        add_group_exogamy = FALSE
      )
    }
  )
}
