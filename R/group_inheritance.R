#' Specify newborn group identity inheritance
#'
#' Assign one group independently of genealogical ancestry shares. Attach this
#' object to a dimension with [group()]. The default selects either
#' parent with equal probability, with maternal fallback.
#' Existing people are not reassigned, and ancestry is unchanged by assignment.
#' [group()] also accepts a schedule tibble with `time_start`, `time_end`, and a
#' `model` list column containing complete inheritance objects. Half-open
#' absolute-year intervals select the model at birth. A custom rule must return
#' one of the dimension's categories; declare categories that may only appear
#' later in [group()].
#'
#' @param method One of `"mother"`, `"father"`, `"random_parent"`, `"same_sex"`,
#'   `"opposite_sex"`, `"ancestry"`, `"hierarchy"`, `"custom"`, or `"fixed"`.
#'   Random parent chooses each with probability 0.5. Same/opposite sex chooses
#'   the parent whose sex matches/differs from the child's. Ancestry samples
#'   known child shares, conditioning out unknown ancestry. Missing information
#'   required by a built-in rule falls back to the mother, including missing
#'   fathers, entirely unknown ancestry, and missing hierarchy ranks.
#' @param fun For `"custom"`, an R function receiving named arguments
#'   `mother_group`, `father_group`, `mother_ancestry`, `father_ancestry`,
#'   `child_ancestry`, `child_sex` (`"f"`/`"m"`), and `parents_in_union`
#'   (logical union status at birth, including on-demand formation). Nonmarital
#'   births and unknown fathers give FALSE. Older callbacks without this new
#'   argument or `...` remain supported. Use `...` for unused inputs.
#'   Ancestry vectors (the dimension's own shares, when it tracks ancestry)
#'   are named by category and `unknown`. A missing father has
#'   NA group and NULL ancestry. Return exactly one nonmissing identity in the
#'   possible group set. Ordinary R random draws are supported; functions must
#'   not mutate simulation/global state. Invalid returns and errors fail.
#' @param order For `"hierarchy"`, unique group labels, highest-ranked first.
#'   Choose between parental identities, not ancestry. Complex orders use `fun`.
#' @param passing For `"hierarchy"`, probability in `[0,1]` of reversing the
#'   order for a birth. Defaults to zero. Equal parental identities and missing
#'   ranks/parents retain maternal assignment without a passing draw.
#' @param fixed_value For `"fixed"`, one category label or finite numeric code in
#'   the dimension's categories. Every newborn receives this category, regardless
#'   of parents or sex, without an identity random draw or R callback.
#' @return An inspectable `hammel_group_inheritance` object.
#' @export
#' @examples
#' create_group_inheritance("fixed", fixed_value = "none")
#' create_group_inheritance("mother")
#' create_group_inheritance("hierarchy", order = c("A", "B"), passing = 0.05)
#' create_group_inheritance("custom",
#'   fun = function(mother_group, father_group, ...) {
#'     if (is.na(father_group) || mother_group == father_group) mother_group else "mixed"
#'   })
create_group_inheritance <- function(method = "random_parent", fun = NULL,
                                     order = NULL, passing = 0, fixed_value = NULL) {
  model <- structure(list(method = method, fun = fun,
                          order = order, passing = passing, fixed_value = fixed_value), class = "hammel_group_inheritance")
  validate_group_inheritance(model)
  model
}

inheritance_methods <- function() {
  c("mother", "father", "random_parent", "same_sex", "opposite_sex", "ancestry", "hierarchy", "custom", "fixed")
}

normalize_group_inheritance <- function(model) {
  if (is_model_schedule(model)) {
    validate_model_argument(model, validate_group_inheritance, "inheritance")
    return(normalize_model_argument(model))
  }
  if (is.null(model)) model <- create_group_inheritance()
  if (is.character(model)) model <- create_group_inheritance(method = model)
  validate_group_inheritance(model)
  model
}

validate_group_inheritance <- function(model) {
  if (is.null(model)) return(invisible(model)) # Missing inheritance fields use the current default.
  if (!inherits(model, "hammel_group_inheritance") || !is.list(model))
    stop("`inheritance` must be created with create_group_inheritance().", call. = FALSE)
  if (!is.character(model$method) || length(model$method) != 1L || is.na(model$method) ||
      !model$method %in% inheritance_methods()) stop("Invalid group inheritance method.", call. = FALSE)
  for (field in "order") {
    value <- model[[field]]
    if (!is.null(value)) {
      check_groups(tibble::tibble(group = value), field)
      if (anyDuplicated(as.character(value))) stop(sprintf("Inheritance `%s` must have unique labels.", field), call. = FALSE)
    }
  }
  passing <- model$passing
  if (!is.numeric(passing) || is.complex(passing) || length(passing) != 1L ||
      !is.null(dim(passing)) || !is.finite(passing) || passing < 0 || passing > 1)
    stop("`passing` must be a finite probability between 0 and 1.", call. = FALSE)
  if (model$method == "hierarchy") {
    if (!length(model$order)) stop("Hierarchy inheritance requires an `order` vector.", call. = FALSE)
  } else if (!is.null(model$order) || passing != 0)
    stop("`order` and nonzero `passing` apply only to hierarchy inheritance.", call. = FALSE)
  if (model$method == "custom") {
    if (!is.function(model$fun)) stop("Custom inheritance requires a function in `fun`.", call. = FALSE)
  } else if (!is.null(model$fun)) stop("`fun` applies only to custom inheritance.", call. = FALSE)
  if (model$method == "fixed") {
    value <- model$fixed_value
    if (!(is.character(value) || is.numeric(value)) || is.complex(value) ||
        !is.null(dim(value)) || length(value) != 1L || anyNA(value) ||
        !nzchar(trimws(as.character(value))) || (is.numeric(value) && !is.finite(value)))
      stop("Fixed inheritance requires one nonmissing category in `fixed_value`.", call. = FALSE)
  } else if (!is.null(model$fixed_value)) stop("`fixed_value` applies only to fixed inheritance.", call. = FALSE)
  invisible(model)
}

extend_identity_groups <- function(values, additions) {
  extra <- setdiff(as.character(additions), as.character(values))
  if (!length(extra)) return(values)
  if (is.numeric(values)) {
    numbers <- suppressWarnings(as.numeric(extra))
    if (all(is.finite(numbers)) && identical(as.character(numbers), extra)) {
      if (is.integer(values) && all(numbers == trunc(numbers) & abs(numbers) <= .Machine$integer.max))
        return(c(values, as.integer(numbers)))
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
  fixed <- if (model$method == "fixed") match(as.character(model$fixed_value), labels) - 1L else -1L
  if (is.na(fixed)) stop("Fixed inheritance `fixed_value` is an undeclared group identity; declare it in group(categories = ...).", call. = FALSE)
  ranks <- match(labels, as.character(model$order)) - 1L
  ranks[is.na(ranks)] <- -1L
  ancestry_groups <- match(setdiff(components, "unknown"), labels) - 1L
  ancestry_groups[is.na(ancestry_groups)] <- -1L
  custom <- function(mother, father, sex, maternal, paternal, child, parents_in_union) {
    inputs <- list(
      mother_group = values[mother + 1L],
      father_group = if (father < 0L) values[NA_integer_] else values[father + 1L],
      mother_ancestry = stats::setNames(maternal, components),
      father_ancestry = if (father < 0L) NULL else stats::setNames(paternal, components),
      child_ancestry = stats::setNames(child, components), child_sex = sex)
    if (any(c("parents_in_union", "...") %in% names(formals(model$fun))))
      inputs$parents_in_union <- parents_in_union
    result <- do.call(model$fun, inputs)
    if (!(is.character(result) || is.numeric(result)) || is.complex(result) ||
        !is.null(dim(result)) || length(result) != 1L || anyNA(result) ||
        (is.numeric(result) && !is.finite(result)))
      stop("Custom inheritance must return one nonmissing group identity.", call. = FALSE)
    code <- match(as.character(result), labels) - 1L
    if (is.na(code)) stop("Custom inheritance returned an undeclared group identity.", call. = FALSE)
    as.integer(code)
  }
  list(values = values, group_count = length(values), method = match(model$method, inheritance_methods()) - 1L,
       fixed = as.integer(fixed), passing = model$passing, ranks = ranks, ancestry_groups = ancestry_groups, custom = custom)
}
