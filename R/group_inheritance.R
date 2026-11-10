# Construct, validate, and prepare newborn group inheritance models.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Specify newborn group identity inheritance
#'
#' @description
#'
#' A model specification indicating how the group category of newborns should be
#' assigned for a particular group dimension. Used in [group()].
#'
#' @details
#'
#' Whenever groups are in use, an inheritance rule is required to determine how
#' the group membership of newborns is determined. The function provides several
#' out of the box methods, as well as the ability for users to define a custom
#' method:
#'
#' - `random_parent`: The default method. The newborn is assigned the group
#' membership of a randomly selected parent with equal probability.
#' - `mother` - The newborn is assigned the group membership of the mother.
#' - `father` - The newborn is assigned the group membership of the father.
#' - `same_sex` - The newborn is assigned the group membership of their same-sex
#' parent.
#' - `different_sex` - The newborn is assigned the group membership of their
#' opposite-sex parent.
#' - `ancestry` - The newborn is assigned a group membership randomly, weighted
#' by their ancestry share from each group. `ancestry = TRUE` must be set for
#' the group to use this method.
#' - `fixed` - All newborns are assigned a fixed group membership given by the
#' `fixed_value` argument.
#' - `hierarchy` - The newborn is assigned the group membership from the parent
#' whose membership is higher on the `order` argument. If the `passing` argument
#' is provided, it defines the probability that this ordering is inverted.
#' - `custom` - The newborn is assigned a group membership based on a custom
#' R function provided by the user with the `fun` argument. Simulations with
#' custom R functions are slower to process.
#'
#' For methods that require no additional arguments (`random_parent`, `mother`,
#' `father`, `same_sex`, `different_sex`, `ancestry`), the name of the method
#' may be passed to arguments requiring an inheritance model like [group()]
#' instead of the model itself.
#'
#' For the custom function, the following variables are exposed and can be
#' called from the function declaration:
#'
#' - `mother_group`: The mother's group membership.
#' - `father_group`: The father's group membership. NA if unknown.
#' - `mother_ancestry`: Numeric vector of ancestry share from each group for
#' mother, as well as an "unknown" value for unknown ancestry.
#' - `father_ancestry`: Numeric vector of ancestry share from each group for
#' father, as well as an "unknown" value for unknown ancestry.
#' - `child_ancestry`: Numeric vector of ancestry share from each group for
#' child, as well as an "unknown" value for unknown ancestry.
#' - `child_sex`: Either "f" or "m".
#' - `parents_in_union`: Logical indicating whether the parents were in a union,
#' including on-demand union formation.
#'
#' @param method One of `"mother"`, `"father"`, `"random_parent"`, `"same_sex"`,
#'   `"opposite_sex"`, `"ancestry"`, `"fixed"`, `"hierarchy"`, or `"custom"`.
#'   See details for further information.
#' @param fun For `"custom"`, an R function defined by the user. See details for
#' variables exposed for this function.
#' @param order For `"hierarchy"`, unique group labels, highest-ranked first.
#'   Chooses between parental group membership, not ancestry.
#' @param passing For `"hierarchy"`, probability in `[0,1]` of reversing the
#'   order for a birth. Defaults to zero.
#' @param fixed_value For `"fixed"`, one category label or finite numeric code
#'   in the dimension's categories. Every newborn receives this category.
#' @return An inspectable `hammel_group_inheritance` object.
#' @export
#' @examples
#' group_inheritance("fixed", fixed_value = "none")
#' group_inheritance("mother")
#' group_inheritance("hierarchy", order = c("B", "W"), passing = 0.05)
#' group_inheritance("custom",
#'   fun = function(mother_group, father_group, ...) {
#'     if (is.na(father_group) || mother_group == father_group) {
#'       mother_group
#'     } else {
#'       "mixed"
#'     }
#'   }
#' )
group_inheritance <- function(
  method = "random_parent",
  fun = NULL,
  order = NULL,
  passing = 0,
  fixed_value = NULL
) {
  model <- structure(
    list(
      method = method,
      fun = fun,
      order = order,
      passing = passing,
      fixed_value = fixed_value
    ),
    class = "hammel_group_inheritance"
  )
  validate_group_inheritance(model)
  model
}

# HELPER FUNCTIONS ---------------------------------------------------------

inheritance_methods <- function() {
  c(
    "mother", "father", "random_parent", "same_sex", "opposite_sex", "ancestry",
    "hierarchy", "custom", "fixed"
  )
}

normalize_group_inheritance <- function(model) {
  if (is_model_schedule(model)) {
    validate_model_argument(model, validate_group_inheritance, "inheritance")
    return(normalize_model_argument(model))
  }
  if (is.null(model)) {
    model <- group_inheritance()
  }
  if (is.character(model)) {
    model <- group_inheritance(method = model)
  }
  validate_group_inheritance(model)
  model
}

validate_group_inheritance <- function(model) {
  # Missing inheritance fields use the current default.
  if (is.null(model)) {
    return(invisible(model))
  }
  if (!inherits(model, "hammel_group_inheritance") || !is.list(model)) {
    stop(
      "`inheritance` must be created with group_inheritance().",
      call. = FALSE
    )
  }
  if (!is.character(model$method) || length(model$method) != 1L ||
        is.na(model$method) ||
        !model$method %in% inheritance_methods()) {
    stop("Invalid group inheritance method.", call. = FALSE)
  }
  for (field in "order") {
    value <- model[[field]]
    if (!is.null(value)) {
      check_groups(tibble::tibble(group = value), field)
      if (anyDuplicated(as.character(value))) {
        stop(
          sprintf("Inheritance `%s` must have unique labels.", field),
          call. = FALSE
        )
      }
    }
  }
  passing <- model$passing
  if (!is.numeric(passing) || is.complex(passing) || length(passing) != 1L ||
        !is.null(dim(passing)) ||
        !is.finite(passing) ||
        passing < 0 ||
        passing > 1) {
    stop(
      "`passing` must be a finite probability between 0 and 1.",
      call. = FALSE
    )
  }
  if (model$method == "hierarchy") {
    if (!length(model$order)) {
      stop("Hierarchy inheritance requires an `order` vector.", call. = FALSE)
    }
  } else if (!is.null(model$order) || passing != 0) {
    stop(
      "`order` and nonzero `passing` apply only to hierarchy inheritance.",
      call. = FALSE
    )
  }
  if (model$method == "custom") {
    if (!is.function(model$fun)) {
      stop("Custom inheritance requires a function in `fun`.", call. = FALSE)
    }
  } else if (!is.null(model$fun)) {
    stop("`fun` applies only to custom inheritance.", call. = FALSE)
  }
  if (model$method == "fixed") {
    value <- model$fixed_value
    if (!(is.character(value) || is.numeric(value)) || is.complex(value) ||
          !is.null(dim(value)) ||
          length(value) != 1L ||
          anyNA(value) ||
          !nzchar(trimws(as.character(value))) ||
          (is.numeric(value) && !is.finite(value))) {
      stop(
        paste0(
          "Fixed inheritance requires one nonmissing category in ",
          "`fixed_value`."
        ),
        call. = FALSE
      )
    }
  } else if (!is.null(model$fixed_value)) {
    stop("`fixed_value` applies only to fixed inheritance.", call. = FALSE)
  }
  invisible(model)
}

extend_identity_groups <- function(values, additions) {
  extra <- setdiff(as.character(additions), as.character(values))
  if (!length(extra)) {
    return(values)
  }
  if (is.numeric(values)) {
    numbers <- suppressWarnings(as.numeric(extra))
    if (all(is.finite(numbers)) && identical(as.character(numbers), extra)) {
      if (is.integer(values) &&
            all(
              numbers == trunc(numbers) & abs(numbers) <= .Machine$integer.max
            )) {
        return(c(values, as.integer(numbers)))
      }
      return(c(values, numbers))
    }
  }
  c(as.character(values), extra)
}

# Engine form of one dimension's inheritance rule. `values` are the
# dimension's categories and `components` its ancestry columns (categories
# with shares, then unknown), empty when ancestry is not tracked.
prepare_group_inheritance <- function(model, values, components = character()) {
  model <- normalize_group_inheritance(model)
  labels <- as.character(values)
  fixed <- if (model$method == "fixed") {
    match(as.character(model$fixed_value), labels) - 1L
  } else {
    -1L
  }
  if (is.na(fixed)) {
    stop(
      paste0(
        "Fixed inheritance `fixed_value` is an undeclared ",
        "group identity; declare it in group(categories = ...)."
      ),
      call. = FALSE
    )
  }
  ranks <- match(labels, as.character(model$order)) - 1L
  ranks[is.na(ranks)] <- -1L
  ancestry_groups <- match(setdiff(components, "unknown"), labels) - 1L
  ancestry_groups[is.na(ancestry_groups)] <- -1L
  custom <- function(
    mother, father, sex, maternal, paternal, child,
    parents_in_union
  ) {
    inputs <- list(
      mother_group = values[mother + 1L],
      father_group = if (father < 0L) {
        values[NA_integer_]
      } else {
        values[father + 1L]
      },
      mother_ancestry = stats::setNames(maternal, components),
      father_ancestry = if (father < 0L) {
        NULL
      } else {
        stats::setNames(paternal, components)
      },
      child_ancestry = stats::setNames(child, components),
      child_sex = sex
    )
    if (any(c("parents_in_union", "...") %in% names(formals(model$fun)))) {
      inputs$parents_in_union <- parents_in_union
    }
    result <- do.call(model$fun, inputs)
    if (!(is.character(result) || is.numeric(result)) || is.complex(result) ||
          !is.null(dim(result)) ||
          length(result) != 1L ||
          anyNA(result) ||
          (is.numeric(result) && !is.finite(result))) {
      stop(
        "Custom inheritance must return one nonmissing group identity.",
        call. = FALSE
      )
    }
    code <- match(as.character(result), labels) - 1L
    if (is.na(code)) {
      stop(
        "Custom inheritance returned an undeclared group identity.",
        call. = FALSE
      )
    }
    as.integer(code)
  }
  list(
    values = values,
    group_count = length(values),
    method = match(model$method, inheritance_methods()) - 1L,
    fixed = as.integer(fixed),
    passing = model$passing,
    ranks = ranks,
    ancestry_groups = ancestry_groups,
    custom = custom
  )
}
