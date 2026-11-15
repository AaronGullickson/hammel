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

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a group inheritance model
#'
#' @description
#'
#' Summary method for class `"hammel_group_inheritance"`
#'
#' @param object A `hammel_group_inheritance` object returned by
#'   [group_inheritance()].
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_group_inheritance`
#' containing:
#' \describe{
#'   \item{method}{The inheritance method.}
#'   \item{order}{For `"hierarchy"`, the group labels, highest-ranked first;
#'     otherwise `NULL`.}
#'   \item{passing}{For `"hierarchy"`, the probability of reversing the
#'     order; otherwise `NULL`.}
#'   \item{fixed_value}{For `"fixed"`, the category every newborn receives;
#'     otherwise `NULL`.}
#' }
#'
#' @seealso [group_inheritance()], [group()]
#' @examples
#' summary(group_inheritance("hierarchy", order = c("B", "W"), passing = 0.05))
#' @export
summary.hammel_group_inheritance <- function(object, ...) {
  hierarchy <- identical(object$method, "hierarchy")
  structure(
    list(
      method = object$method,
      order = if (hierarchy) object$order else NULL,
      passing = if (hierarchy) object$passing else NULL,
      fixed_value = object$fixed_value
    ),
    class = "summary.hammel_group_inheritance"
  )
}

#' @rdname summary.hammel_group_inheritance
#' @param x A `summary.hammel_group_inheritance` object returned by
#'   `summary.hammel_group_inheritance()`.
#' @export
print.summary.hammel_group_inheritance <- function(x, ...) {
  print_model_summary("Group inheritance model", x)
}

# HELPER FUNCTIONS ---------------------------------------------------------

# Lines describing an inheritance model summary; see model_summary_lines().
inheritance_summary_lines <- function(x) {
  lines <- c(Method = x$method)
  if (!is.null(x$order)) {
    lines[["Order"]] <- paste(x$order, collapse = " > ")
    lines[["Passing probability"]] <- format(x$passing)
  }
  if (!is.null(x$fixed_value)) {
    lines[["Fixed value"]] <- format(x$fixed_value)
  }
  lines
}

inheritance_methods <- function() {
  c(
    "mother", "father", "random_parent", "same_sex", "opposite_sex", "ancestry",
    "hierarchy", "custom", "fixed"
  )
}

normalize_group_inheritance <- function(model, call = rlang::caller_env()) {
  if (is_model_schedule(model)) {
    validate_model_argument(
      model, validate_group_inheritance, "inheritance",
      call = call
    )
    return(normalize_model_argument(model))
  }
  if (is.null(model)) {
    model <- group_inheritance()
  }
  if (is.character(model)) {
    if (length(model) != 1L || is.na(model) ||
          !model %in% inheritance_methods()) {
      abort_inheritance_method(model, call = call)
    }
    model <- group_inheritance(method = model)
  }
  validate_group_inheritance(model, call = call)
  model
}

# Report an unknown inheritance method.
abort_inheritance_method <- function(method, call = rlang::caller_env()) {
  methods <- inheritance_methods()
  abort_hammel(
    c(
      "The inheritance method must be one of {.or {.val {methods}}}.",
      "x" = "It is {describe_value(method)}."
    ),
    "specification",
    argument = "method",
    call = call
  )
}

validate_group_inheritance <- function(model, call = rlang::caller_env()) {
  # Missing inheritance fields use the current default.
  if (is.null(model)) {
    return(invisible(model))
  }
  if (!inherits(model, "hammel_group_inheritance") || !is.list(model)) {
    abort_hammel(
      c(
        "{.arg inheritance} must be an inheritance rule.",
        "x" = "It is {describe_value(model)}.",
        "i" = paste(
          "Create one with {.fn group_inheritance}, or give a method name",
          "such as {.val mother}."
        )
      ),
      "specification",
      argument = "inheritance",
      call = call
    )
  }
  method <- model$method
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% inheritance_methods()) {
    abort_inheritance_method(method, call = call)
  }
  value <- model$order
  if (!is.null(value)) {
    check_group_labels(value, "order", call = call)
    if (anyDuplicated(as.character(value))) {
      repeated <- unique(as.character(value)[duplicated(as.character(value))])
      abort_hammel(
        c(
          "{.arg order} must list each category once.",
          "x" = "Repeated: {.val {repeated}}."
        ),
        "specification",
        argument = "order",
        call = call
      )
    }
  }
  passing <- model$passing
  if (!is.numeric(passing) || is.complex(passing) || length(passing) != 1L ||
        !is.null(dim(passing)) ||
        !is.finite(passing) ||
        passing < 0 ||
        passing > 1) {
    abort_hammel(
      c(
        "{.arg passing} must be a single probability between 0 and 1.",
        "x" = "It is {describe_value(passing)}."
      ),
      "argument",
      argument = "passing",
      call = call
    )
  }
  unused <- function(argument, applies) {
    abort_hammel(
      c(
        "{.arg {argument}} applies only to {applies} inheritance.",
        "x" = "The method is {.val {method}}."
      ),
      "specification",
      argument = argument,
      call = call
    )
  }
  if (method == "hierarchy") {
    if (!length(model$order)) {
      abort_hammel(
        c(
          "Hierarchy inheritance needs an {.arg order}.",
          "i" = "List the categories from highest to lowest rank."
        ),
        "specification",
        argument = "order",
        call = call
      )
    }
  } else if (!is.null(model$order)) {
    unused("order", "hierarchy")
  } else if (passing != 0) {
    unused("passing", "hierarchy")
  }
  if (method == "custom") {
    if (!is.function(model$fun)) {
      abort_hammel(
        c(
          "Custom inheritance needs a function in {.arg fun}.",
          "x" = "It is {describe_value(model$fun)}."
        ),
        "specification",
        argument = "fun",
        call = call
      )
    }
  } else if (!is.null(model$fun)) {
    unused("fun", "custom")
  }
  if (method == "fixed") {
    value <- model$fixed_value
    if (!(is.character(value) || is.numeric(value)) || is.complex(value) ||
          !is.null(dim(value)) ||
          length(value) != 1L ||
          anyNA(value) ||
          !nzchar(trimws(as.character(value))) ||
          (is.numeric(value) && !is.finite(value))) {
      abort_hammel(
        c(
          "Fixed inheritance needs one category in {.arg fixed_value}.",
          "x" = "It is {describe_value(value)}."
        ),
        "specification",
        argument = "fixed_value",
        call = call
      )
    }
  } else if (!is.null(model$fixed_value)) {
    unused("fixed_value", "fixed")
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
prepare_group_inheritance <- function(
  model,
  values,
  components = character(),
  call = rlang::caller_env()
) {
  model <- normalize_group_inheritance(model, call = call)
  labels <- as.character(values)
  fixed <- if (model$method == "fixed") {
    match(as.character(model$fixed_value), labels) - 1L
  } else {
    -1L
  }
  if (is.na(fixed)) {
    fixed_value <- model$fixed_value
    abort_hammel(
      c(
        "Fixed inheritance gives newborns an undeclared category.",
        "x" = "{.arg fixed_value} is {.val {fixed_value}}.",
        "i" = paste(
          "Declare it with {.code group(categories = ...)}. Categories are",
          "{.val {labels}}."
        )
      ),
      "groups",
      categories = as.character(fixed_value),
      call = call
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
      abort_hammel(
        c(
          "A custom inheritance function must return one category.",
          "x" = "It returned {describe_value(result)}."
        ),
        "runtime",
        process = "inheritance",
        call = call
      )
    }
    code <- match(as.character(result), labels) - 1L
    if (is.na(code)) {
      abort_hammel(
        c(
          "A custom inheritance function returned an undeclared category.",
          "x" = "It returned {.val {result}}.",
          "i" = "Declared categories are {.val {labels}}."
        ),
        "runtime",
        process = "inheritance",
        call = call,
        categories = as.character(result)
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
