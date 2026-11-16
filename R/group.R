# Declare group dimensions and apply their definitions to models and
# populations.

# PUBLIC FUNCTIONS ---------------------------------------------------------

# Group dimensions are always declared (D134). A declaration is a named list
# of `hammel_group` objects, one per dimension, keyed by the persons column.

#' Describe and declare a group dimension
#'
#' @description
#'
#' The `group()` function describes a single group dimension, including valid
#' categories, whether to track ancestry, and inheritance rules. The
#' `set_groups()` function declares all group dimensions for a `hammel_model`.
#'
#' For the `set_groups` function, each named argument should identify a group
#' name and provide the returned value from a `group()` function. See examples
#' for syntax.
#'
#' @param categories Optional complete, ordered vector of the dimension's
#'   categories (character labels or finite numeric codes). `NULL` (the default)
#'   either keeps the categories the population already declares, if any, or
#'   they are inferred from the population's values in first-occurrence order,
#'   plus transition destinations and categories with ancestry shares. When
#'   given, every value in the population and every transition destination must
#'   be listed; categories can be declared even if starting population contains
#'   no members of that category.
#' @param ancestry Logical, default `FALSE`. Track ancestry shares for this
#'   dimension in the simulation output. The returned population then has an
#'   `ancestry` table indicating the share of ancestry from each category of the
#'   group. Required for use of the `ancestry` inheritance method.
#' @param inheritance Either a single character string, a [group_inheritance()]
#'   object, or a calendar schedule tibble with `time_start`, `time_end` and a
#'   list-column `model` providing the group inheritance model for a time
#'   period. Determines how group membership is assigned to newborns. Defaults
#'   to `random_parent`. Character sting entries will use built-in methods
#'   provided gy the `group_inheritance` model.
#' @param model A `hammel_model`, used by `set_groups()`.
#' @param ... Named `group()` objects, one per dimension, named by the persons
#'   column. No arguments removes the declaration.
#' @return `group()` returns a `hammel_group` object describing one dimension.
#'   `set_groups()` returns an updated `hammel_model`.
#' @seealso [hammel_population()], [generate_uniform_population()],
#'   [group_inheritance()], [group_transition()]
#' @rdname group
#' @order 1
#' @export
#' @examples
#' # create two groups - ethnicity and education
#' model <- hammel_model() |>
#'   set_groups(
#'     ethnicity = group(
#'       categories = c("W", "B", "O"),
#'       ancestry = TRUE,
#'       inheritance = group_inheritance("ancestry")
#'     ),
#'     education = group(
#'       categories = c("LHS", "HS", "SC", "C"),
#'       ancestry = FALSE,
#'       inheritance = group_inheritance("fixed", fixed_value = "LHS")
#'     )
#'   )
group <- function(
  categories = NULL,
  ancestry = FALSE,
  inheritance = group_inheritance()
) {
  x <- structure(
    list(
      categories = categories,
      ancestry = ancestry,
      inheritance = normalize_group_inheritance(
        inheritance,
        call = rlang::current_env()
      )
    ),
    class = "hammel_group"
  )
  validate_group(x)
  x
}

#' @rdname group
#' @order 2
#' @export
set_groups <- function(model, ...) {
  check_model(model)
  groups <- list(...)
  declaration <- if (length(groups)) {
    normalize_group_declaration(groups, name = "...")
  } else {
    NULL
  }
  model["groups"] <- list(declaration)
  # Ancestry shares and empty transition histories are built at simulation.
  if (!is.null(model$population)) {
    model$population <- apply_group_declaration(model$population, declaration)
  }
  model
}

# SUMMARY AND PRINT FUNCTIONS ----------------------------------------------

#' Summarize a group dimension
#'
#' @description
#'
#' Summary method for class `"hammel_group"`
#'
#' @param object A `hammel_group` object returned by [group()].
#' @param categories Optional categories to report instead of the declared
#'   ones, such as the categories a population resolves for the dimension.
#'   `NULL` (default) uses the declared categories.
#' @param ... Additional arguments, currently unused.
#'
#' @return
#' `summary()` returns an object of class `summary.hammel_group` containing:
#' \describe{
#'   \item{categories}{The dimension's categories, or `NULL` when they are
#'     inferred from the population.}
#'   \item{ancestry}{Whether ancestry shares are tracked.}
#'   \item{inheritance}{A summary of the newborn inheritance model: a
#'     `summary.hammel_group_inheritance`, or a list of `periods` and
#'     `models` for a calendar schedule.}
#' }
#'
#' @seealso [group()], [group_inheritance()]
#' @examples
#' ethnicity <- group(
#'   categories = c("W", "B", "O"),
#'   ancestry = TRUE,
#'   inheritance = group_inheritance("ancestry")
#' )
#' summary(ethnicity)
#' @export
summary.hammel_group <- function(object, categories = NULL, ...) {
  structure(
    list(
      categories = if (is.null(categories)) object$categories else categories,
      ancestry = isTRUE(object$ancestry),
      inheritance = summarize_model_argument(object$inheritance)
    ),
    class = "summary.hammel_group"
  )
}

#' @rdname summary.hammel_group
#' @param x A `summary.hammel_group` object returned by
#'   `summary.hammel_group()`.
#' @export
print.summary.hammel_group <- function(x, ...) {

  cli::cli_h2("Group dimension")

  list <- open_summary_list()
  print_group_details(x)
  cli::cli_end(list)
  invisible(x)
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Names Hammel uses for its own variables and columns, and patterns that would
# make derived names ambiguous.
reserved_group_names <- c(
  "person_id", "sex", "birth_time", "death_time",
  "mother_id", "father_id",
  "union_type", "previous_union", "parity", "partnered",
  "age_difference",
  "origin", "destination", "hetfert", "unknown", "time", "rate",
  "arrival_time", "stream"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Print the bullets of a group summary. Shared by the group and model
# summaries.
print_group_details <- function(x) {
  categories <- if (is.null(x$categories)) {
    "inferred from the population"
  } else {
    paste(x$categories, collapse = ", ")
  }
  print_stat("Categories", categories)
  print_stat("Ancestry tracked", if (x$ancestry) "yes" else "no")
  print_model_argument("Inheritance", x$inheritance)
  invisible(x)
}

validate_group <- function(x, name = NULL, call = rlang::caller_env()) {
  if (!inherits(x, "hammel_group") || !is.list(x)) {
    abort_hammel(
      c(
        if (is.null(name)) {
          "A group dimension must be created with {.fn group}."
        } else {
          "Group dimension {.field {name}} must be created with {.fn group}."
        },
        "x" = "It is {describe_value(x)}."
      ),
      "specification",
      group = name,
      call = call
    )
  }
  if (!is.null(x$categories)) {
    check_group_labels(x$categories, "categories", call = call)
    if (!length(x$categories) || anyDuplicated(as.character(x$categories))) {
      repeated <- unique(
        as.character(x$categories)[duplicated(as.character(x$categories))]
      )
      abort_hammel(
        c(
          "{.arg categories} must be a nonempty vector of unique labels.",
          "x" = if (length(repeated)) {
            "Repeated: {.val {repeated}}."
          } else {
            "It is empty."
          }
        ),
        "groups",
        argument = "categories",
        group = name,
        call = call
      )
    }
  }
  if (!is.logical(x$ancestry) || length(x$ancestry) != 1L ||
        is.na(x$ancestry)) {
    abort_hammel(
      c(
        "{.arg ancestry} must be {.code TRUE} or {.code FALSE}.",
        "x" = "It is {describe_value(x$ancestry)}."
      ),
      "argument",
      argument = "ancestry",
      group = name,
      call = call
    )
  }
  if (is_model_schedule(x$inheritance)) {
    validate_model_argument(
      x$inheritance, validate_group_inheritance, "inheritance",
      call = call
    )
  } else {
    validate_group_inheritance(x$inheritance, call = call)
  }
  invisible(x)
}

check_group_names <- function(names, call = rlang::caller_env()) {
  if (is.null(names) || anyNA(names) || any(!nzchar(names)) ||
        anyDuplicated(names)) {
    abort_hammel(
      c(
        "Group dimensions must have unique, nonempty names.",
        "x" = if (is.null(names)) {
          "They have no names."
        } else {
          "Their names are {.val {names}}."
        }
      ),
      "groups",
      call = call
    )
  }
  bad <- names[
    names %in% reserved_group_names |
      grepl("^(partner|same|female|male|ancestry)_", names) |
      grepl("_(w|h)$", names) |
      !grepl("^[A-Za-z][A-Za-z0-9_.]*$", names)
  ]
  if (length(bad)) {
    abort_hammel(
      c(
        "Some group dimension names are reserved or invalid.",
        "x" = "Invalid: {.val {bad}}.",
        "i" = paste(
          "Names must start with a letter and use letters, digits, {.code _}",
          "or {.code .}. They can't start with {.code partner_},",
          "{.code same_}, {.code female_}, {.code male_} or",
          "{.code ancestry_}, end in {.code _w} or {.code _h}, or be a",
          "Hammel column name such as {.val sex} or {.val parity}."
        )
      ),
      "groups",
      group = bad,
      call = call
    )
  }
  invisible(names)
}

# Accept NULL, a character vector of dimension names (default settings) or a
# named list of group() objects.
normalize_group_declaration <- function(
  groups,
  name = "groups",
  call = rlang::caller_env()
) {
  if (is.null(groups) || (is.list(groups) && !length(groups))) {
    return(NULL)
  }
  if (is.character(groups)) {
    check_group_names(groups, call = call)
    return(stats::setNames(lapply(groups, function(name) group()), groups))
  }
  if (!is.list(groups) || inherits(groups, "hammel_group")) {
    abort_hammel(
      c(
        paste(
          "{.arg {name}} must be dimension names or {.fn group} objects",
          "named by dimension."
        ),
        "x" = "It is {describe_value(groups)}.",
        "i" = "For example, {.code list(race = group(ancestry = TRUE))}."
      ),
      "groups",
      argument = name,
      call = call
    )
  }
  check_group_names(names(groups), call = call)
  for (dimension in names(groups)) {
    with_error_context(
      validate_group(groups[[dimension]], dimension, call = call),
      paste0("In group dimension ", dimension, ".")
    )
    groups[[dimension]]["inheritance"] <- list(
      normalize_group_inheritance(groups[[dimension]]$inheritance, call = call)
    )
  }
  groups
}

# Replace a population's declaration. Ancestry shares are dropped when the new
# declaration does not track them (and initialized later when it does).
apply_group_declaration <- function(
  population,
  declaration,
  call = rlang::caller_env()
) {
  declaration <- normalize_group_declaration(declaration, call = call)
  # A dimension declared without categories keeps the categories the
  # population already declares (so turning on ancestry, say, keeps their
  # order).
  for (d in names(declaration)) {
    if (is.null(declaration[[d]]$categories) &&
          !is.null(population$groups[[d]]$categories)) {
      declaration[[d]]["categories"] <- list(population$groups[[d]]$categories)
    }
  }
  population$groups <- declaration
  # Keep ancestry only for dimensions that track it.
  if (!is.null(population$ancestry)) {
    ancestry <- as_ancestry_list(population$ancestry, declaration, call = call)
    ancestry <- ancestry[intersect(
      names(ancestry), tracked_ancestry(declaration)
    )]
    population$ancestry <- if (length(ancestry)) ancestry else NULL
  }
  validate_group_declaration(declaration, population, call = call)
  population
}

# The declaration in force for a model: set_groups() wins over the population's.
model_groups <- function(model) {
  if (!is.null(model$groups)) model$groups else model$population$groups
}

# Is a group dimension (any, or the named one) declared for this population?
declared_group <- function(population, name = NULL) {
  inherits(population, "hammel_population") &&
    (if (is.null(name)) {
      length(population$groups) > 0L
    } else {
      name %in% names(population$groups)
    })
}

# The name of a population's single declared dimension, for features that use
# one dimension (transitions, matching, observation by group); NULL otherwise.
single_dimension <- function(population) {
  if (inherits(population, "hammel_population") &&
        length(population$groups) == 1L) {
    names(population$groups)
  } else {
    NULL
  }
}

# Check a declaration against a population's persons and transition history.
validate_group_declaration <- function(
  declaration,
  population,
  call = rlang::caller_env()
) {
  persons <- population$persons
  for (name in names(declaration)) {
    if (!name %in% names(persons)) {
      abort_hammel(
        c(
          "Declared group {.field {name}} isn't a column of {.field persons}.",
          "i" = paste(
            "Add the column to the persons table, or remove the",
            "declaration."
          )
        ),
        "groups",
        group = name,
        call = call
      )
    }
    check_group_labels(persons[[name]], paste0("persons$", name), call = call)
    categories <- declaration[[name]]$categories
    # Ancestry tables name their columns by category.
    reserved <- intersect(
      c("unknown", "person_id"),
      as.character(c(persons[[name]], categories))
    )
    if (isTRUE(declaration[[name]]$ancestry) && length(reserved)) {
      abort_hammel(
        c(
          paste(
            "Group {.field {name}} tracks ancestry, so it can't use some",
            "category names."
          ),
          "x" = "Reserved: {.val {reserved}}.",
          "i" = "Ancestry tables use these names for their own columns."
        ),
        "groups",
        group = name,
        categories = reserved,
        call = call
      )
    }
    if (!is.null(categories)) {
      tx <- population$transitions
      if (!is.null(tx) && "group" %in% names(tx)) {
        tx <- tx[tx$group == name, , drop = FALSE]
      }
      values <- unique(
        as.character(c(persons[[name]], tx$origin, tx$destination))
      )
      missing <- setdiff(values, as.character(categories))
      if (length(missing)) {
        abort_hammel(
          c(
            "Group {.field {name}} has values missing from its categories.",
            "x" = "Undeclared: {.val {missing}}.",
            "i" = "Declared categories: {.val {as.character(categories)}}."
          ),
          "groups",
          group = name,
          categories = missing,
          call = call
        )
      }
    }
  }
  if (!length(declaration) && !is.null(population$transitions) &&
        nrow(population$transitions)) {
    abort_hammel(
      c(
        "Transition histories need a declared group dimension.",
        "i" = "Declare dimensions with {.arg groups} or {.fn set_groups}."
      ),
      "groups",
      argument = "transitions",
      call = call
    )
  }
  if (!is.null(population[["ancestry"]])) {
    tracked <- tracked_ancestry(declaration)
    if (!length(declaration) || !length(tracked)) {
      abort_hammel(
        c(
          "Ancestry tables need a group dimension that tracks ancestry.",
          "i" = "Declare it with {.code group(ancestry = TRUE)}."
        ),
        "groups",
        argument = "ancestry",
        call = call
      )
    }
    ancestry <- as_ancestry_list(
      population[["ancestry"]], declaration,
      call = call
    )
    extra <- setdiff(names(ancestry), tracked)
    if (length(extra)) {
      abort_hammel(
        c(
          "Ancestry is supplied for dimensions that don't track ancestry.",
          "x" = "Dimensions: {.field {extra}}.",
          "i" = paste(
            "Declare them with {.code group(ancestry = TRUE)}, or remove",
            "their tables."
          )
        ),
        "groups",
        group = extra,
        call = call
      )
    }
  }
  invisible(declaration)
}
