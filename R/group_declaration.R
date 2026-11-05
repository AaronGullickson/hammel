# Group dimensions are always declared (D134). A declaration is a named list
# of `hammel_group` objects, one per dimension, keyed by the persons column.

#' Describe one group dimension
#'
#' Settings for one declared group dimension, used with [set_groups()],
#' [hammel_population()] or [create_baseline_population()].
#'
#' @param categories Optional complete, ordered vector of the dimension's
#'   categories (character labels or finite numeric codes). `NULL` (the
#'   default) infers them from the population's values in first-occurrence
#'   order, plus transition destinations and categories with ancestry shares.
#'   When given, every value in the population and every transition destination
#'   must be listed; extra categories may start empty, for example the output of
#'   a custom inheritance rule. Group-specific rate tables must cover every
#'   category.
#' @param ancestry Logical, default `FALSE`. Track ancestry shares for this
#'   dimension: the population then has an `ancestry` table, newborns average
#'   their parents' shares, and the `"ancestry"` inheritance method becomes
#'   available.
#' @param inheritance How newborns are assigned to a category: a
#'   [create_group_inheritance()] object, a method name such as `"mother"`, or
#'   a calendar schedule tibble with `time_start`, `time_end` and list-column
#'   `model`. The default selects either parent with equal probability, with
#'   maternal fallback.
#' @return A `hammel_group` object.
#' @seealso [set_groups()]
#' @export
#' @examples
#' group(categories = c("A", "B", "AB"), ancestry = TRUE,
#'       inheritance = create_group_inheritance("mother"))
group <- function(categories = NULL, ancestry = FALSE,
                  inheritance = create_group_inheritance()) {
  x <- structure(list(categories = categories, ancestry = ancestry,
                      inheritance = normalize_group_inheritance(inheritance)),
                 class = "hammel_group")
  validate_group(x)
  x
}

validate_group <- function(x, name = "group") {
  if (!inherits(x, "hammel_group") || !is.list(x))
    stop(sprintf("Group dimension `%s` must be created with group().", name), call. = FALSE)
  if (!is.null(x$categories)) {
    check_groups(tibble::tibble(group = x$categories), "categories")
    if (!length(x$categories) || anyDuplicated(as.character(x$categories)))
      stop(sprintf("`categories` of group `%s` must be a nonempty vector of unique labels.", name),
           call. = FALSE)
  }
  if (!is.logical(x$ancestry) || length(x$ancestry) != 1L || is.na(x$ancestry))
    stop(sprintf("`ancestry` of group `%s` must be TRUE or FALSE.", name), call. = FALSE)
  if (is_model_schedule(x$inheritance))
    validate_model_argument(x$inheritance, validate_group_inheritance, "inheritance")
  else validate_group_inheritance(x$inheritance)
  invisible(x)
}

# Names Hammel uses for its own variables and columns, and patterns that would
# make derived names ambiguous.
reserved_group_names <- c("person_id", "sex", "birth_time", "death_time", "mother_id",
  "father_id", "union_type", "previous_union", "parity", "partnered", "age_difference",
  "origin", "destination", "hetfert", "unknown", "time", "rate")

check_group_names <- function(names) {
  if (is.null(names) || anyNA(names) || any(!nzchar(names)) || anyDuplicated(names))
    stop("Group dimensions must have unique, nonempty names.", call. = FALSE)
  bad <- names[names %in% reserved_group_names |
    grepl("^(partner|same|female|male|ancestry)_", names) | grepl("_(w|h)$", names) |
    !grepl("^[A-Za-z][A-Za-z0-9_.]*$", names)]
  if (length(bad))
    stop(sprintf("Group dimension name(s) %s are reserved or not valid names.",
                 paste0("`", bad, "`", collapse = ", ")), call. = FALSE)
  # Stage 1 of milestone 9 supports one dimension, named `group`; other names
  # and several dimensions arrive with the rate and engine changes.
  if (length(names) > 1L || any(names != "group"))
    stop("Only one group dimension, named `group`, is supported so far.", call. = FALSE)
  invisible(names)
}

# Accept NULL, a character vector of dimension names (default settings) or a
# named list of group() objects.
normalize_group_declaration <- function(groups) {
  if (is.null(groups) || (is.list(groups) && !length(groups))) return(NULL)
  if (is.character(groups)) {
    check_group_names(groups)
    return(stats::setNames(lapply(groups, function(name) group()), groups))
  }
  if (!is.list(groups) || inherits(groups, "hammel_group"))
    stop("`groups` must be dimension names or a named list of group() objects.", call. = FALSE)
  check_group_names(names(groups))
  for (name in names(groups)) {
    validate_group(groups[[name]], name)
    groups[[name]]["inheritance"] <- list(normalize_group_inheritance(groups[[name]]$inheritance))
  }
  groups
}

#' Declare a model's group dimensions
#'
#' Name the persons columns that are group dimensions and give each its
#' settings with [group()]. Groups are always declared: a column that is not
#' declared is an ordinary custom column, kept but never used by any process
#' or analysis. The declaration replaces any declaration carried by the
#' population (for example from [create_baseline_population()]) and may be
#' given before or after [set_population()]; it is checked against the
#' population when both are present and at simulation. Simulated populations
#' carry it, so analysis and restarts need nothing further.
#'
#' Currently one dimension, named `group`, is supported.
#'
#' @param model A `hammel_model`.
#' @param ... Named [group()] objects, one per dimension, named by the persons
#'   column. No arguments removes the declaration.
#' @return The updated `hammel_model`.
#' @export
#' @examples
#' model <- hammel_model() |>
#'   set_population(create_baseline_population(100, groups = list(group = c(A = .5, B = .5)))) |>
#'   set_groups(group = group(inheritance = "mother"))
set_groups <- function(model, ...) {
  check_model(model)
  groups <- list(...)
  declaration <- if (length(groups)) normalize_group_declaration(groups) else NULL
  model["groups"] <- list(declaration)
  # Ancestry shares and empty transition histories are built at simulation.
  if (!is.null(model$population)) model$population <- apply_group_declaration(model$population, declaration)
  model
}

# Replace a population's declaration. Ancestry shares are dropped when the new
# declaration does not track them (and initialized later when it does).
apply_group_declaration <- function(population, declaration) {
  declaration <- normalize_group_declaration(declaration)
  population$groups <- declaration
  if (!isTRUE(declaration$group$ancestry)) population$ancestry <- NULL
  validate_group_declaration(declaration, population)
  population
}

# The declaration in force for a model: set_groups() wins over the population's.
model_groups <- function(model) {
  if (!is.null(model$groups)) model$groups else model$population$groups
}

# Is a dimension declared for this population (or bare persons table)?
declared_group <- function(population, name = "group") {
  inherits(population, "hammel_population") && name %in% names(population$groups)
}

# Check a declaration against a population's persons and transition history.
validate_group_declaration <- function(declaration, population) {
  persons <- population$persons
  for (name in names(declaration)) {
    if (!name %in% names(persons))
      stop(sprintf("Declared group `%s` is not a persons column.", name), call. = FALSE)
    check_groups(persons, "persons")
    categories <- declaration[[name]]$categories
    if (!is.null(categories)) {
      values <- unique(as.character(c(persons[[name]],
        population$transitions$group_origin, population$transitions$group_dest)))
      missing <- setdiff(values, as.character(categories))
      if (length(missing))
        stop(sprintf("Group `%s` has values missing from its categories: %s.", name,
                     paste(missing, collapse = ", ")), call. = FALSE)
    }
  }
  if (!"group" %in% names(declaration)) {
    if (!is.null(population$transitions) && nrow(population$transitions))
      stop("Transition histories require a declared group dimension.", call. = FALSE)
    if (!is.null(population[["ancestry"]]))
      stop("Ancestry requires a declared group dimension with `ancestry = TRUE`.", call. = FALSE)
  } else if (!is.null(population[["ancestry"]]) && !isTRUE(declaration$group$ancestry)) {
    stop("An ancestry table requires `ancestry = TRUE` for its group dimension.", call. = FALSE)
  }
  invisible(declaration)
}
