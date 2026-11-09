# Declare group dimensions and apply their definitions to models and
# populations.

# PUBLIC FUNCTIONS ---------------------------------------------------------

# Group dimensions are always declared (D134). A declaration is a named list
# of `hammel_group` objects, one per dimension, keyed by the persons column.

#' Describe and declare group dimensions
#'
#' `group()` describes one group dimension, for use with `set_groups()`,
#' [hammel_population()] or [create_baseline_population()]. `set_groups()`
#' declares named dimensions on a model, using one `group()` object per
#' dimension. Neither function runs a simulation.
#'
#' @section Declaring model dimensions:
#' Name the persons columns that are group dimensions and give each its settings
#' with [group()]. Groups are always declared: a column that is not declared is
#' an ordinary custom column, kept but never used by any process or analysis.
#' The declaration replaces any declaration carried by the population (for
#' example from [create_baseline_population()]), except that a dimension
#' declared without `categories` keeps the population's declared categories for
#' it; it may be given before or after [set_population()]; it is checked against
#' the population when both are present and at simulation. Simulated populations
#' carry it, so analysis and restarts need nothing further.
#'
#' @param categories Optional complete, ordered vector of the dimension's
#'   categories (character labels or finite numeric codes). `NULL` (the default)
#'   keeps the categories the population already declares (such as those of
#'   [create_baseline_population()]); without any, they are inferred from the
#'   population's values in first-occurrence order, plus transition destinations
#'   and categories with ancestry shares. When given, every value in the
#'   population and every transition destination must be listed; extra
#'   categories may start empty, for example the output of a custom inheritance
#'   rule. Group-specific rate tables must cover every category.
#' @param ancestry Logical, default `FALSE`. Track ancestry shares for this
#'   dimension: the population then has an `ancestry` table, newborns average
#'   their parents' shares, and the `"ancestry"` inheritance method becomes
#'   available.
#' @param inheritance How newborns are assigned to a category: a
#'   [group_inheritance()] object, a method name such as `"mother"`, or a
#'   calendar schedule tibble with `time_start`, `time_end` and list-column
#'   `model`. The default selects either parent with equal probability, with
#'   maternal fallback.
#' @param model A `hammel_model`, used by `set_groups()`.
#' @param ... Named `group()` objects, one per dimension, named by the persons
#'   column. No arguments removes the declaration.
#' @return `group()` returns a `hammel_group` object describing one dimension.
#'   `set_groups()` returns an updated `hammel_model`.
#' @seealso [hammel_population()], [create_baseline_population()],
#'   [group_inheritance()], [group_transition()]
#' @rdname group
#' @order 1
#' @export
#' @examples
#' group(
#'   categories = c("A", "B", "AB"), ancestry = TRUE,
#'   inheritance = group_inheritance("mother")
#' )
group <- function(
  categories = NULL,
  ancestry = FALSE,
  inheritance = group_inheritance()
) {
  x <- structure(
    list(
      categories = categories, ancestry = ancestry,
      inheritance = normalize_group_inheritance(inheritance)
    ),
    class = "hammel_group"
  )
  validate_group(x)
  x
}

#' @rdname group
#' @order 2
#' @export
#' @examples
#' model <- hammel_model() |>
#'   set_population(create_baseline_population(
#'     100,
#'     groups = list(group = c(A = .5, B = .5))
#'   )) |>
#'   set_groups(group = group(inheritance = "mother"))
set_groups <- function(model, ...) {
  check_model(model)
  groups <- list(...)
  declaration <- if (length(groups)) {
    normalize_group_declaration(groups)
  } else {
    NULL
  }
  model["groups"] <- list(declaration)
  # Ancestry shares and empty transition histories are built at simulation.
  if (!is.null(model$population)) {
    model$population <- apply_group_declaration(
      model$population, declaration
    )
  }
  model
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Names Hammel uses for its own variables and columns, and patterns that would
# make derived names ambiguous.
reserved_group_names <- c(
  "person_id", "sex", "birth_time", "death_time",
  "mother_id",
  "father_id", "union_type", "previous_union", "parity", "partnered",
  "age_difference",
  "origin", "destination", "hetfert", "unknown", "time", "rate"
)

# HELPER FUNCTIONS ---------------------------------------------------------

validate_group <- function(x, name = "group") {
  if (!inherits(x, "hammel_group") || !is.list(x)) {
    stop(sprintf("Group dimension `%s` must be created with group().", name),
      call. = FALSE
    )
  }
  if (!is.null(x$categories)) {
    check_groups(tibble::tibble(group = x$categories), "categories")
    if (!length(x$categories) || anyDuplicated(as.character(x$categories))) {
      stop(
        sprintf(paste0(
          "`categories` of group `%s` must be a nonempty ",
          "vector of unique labels."
        ), name),
        call. = FALSE
      )
    }
  }
  if (!is.logical(x$ancestry) ||
    length(x$ancestry) != 1L || is.na(x$ancestry)) {
    stop(sprintf("`ancestry` of group `%s` must be TRUE or FALSE.", name),
      call. = FALSE
    )
  }
  if (is_model_schedule(x$inheritance)) {
    validate_model_argument(
      x$inheritance, validate_group_inheritance,
      "inheritance"
    )
  } else {
    validate_group_inheritance(x$inheritance)
  }
  invisible(x)
}

check_group_names <- function(names) {
  if (is.null(names) || anyNA(names) || any(!nzchar(names)) || anyDuplicated(
    names
  )) {
    stop("Group dimensions must have unique, nonempty names.", call. = FALSE)
  }
  bad <- names[names %in% reserved_group_names |
    grepl("^(partner|same|female|male|ancestry)_", names) | grepl(
    "_(w|h)$",
    names
  ) |
    !grepl("^[A-Za-z][A-Za-z0-9_.]*$", names)]
  if (length(bad)) {
    stop(sprintf(
      "Group dimension name(s) %s are reserved or not valid names.",
      paste0("`", bad, "`", collapse = ", ")
    ), call. = FALSE)
  }
  invisible(names)
}

# Accept NULL, a character vector of dimension names (default settings) or a
# named list of group() objects.
normalize_group_declaration <- function(groups) {
  if (is.null(groups) || (is.list(groups) && !length(groups))) {
    return(NULL)
  }
  if (is.character(groups)) {
    check_group_names(groups)
    return(stats::setNames(lapply(groups, function(name) group()), groups))
  }
  if (!is.list(groups) || inherits(groups, "hammel_group")) {
    stop("`groups` must be dimension names or a named list of group() objects.",
      call. = FALSE
    )
  }
  check_group_names(names(groups))
  for (name in names(groups)) {
    validate_group(groups[[name]], name)
    groups[[name]]["inheritance"] <- list(normalize_group_inheritance(groups[[
      name
    ]]$inheritance))
  }
  groups
}

# Replace a population's declaration. Ancestry shares are dropped when the new
# declaration does not track them (and initialized later when it does).
apply_group_declaration <- function(population, declaration) {
  declaration <- normalize_group_declaration(declaration)
  # A dimension declared without categories keeps the categories the
  # population already declares (so turning on ancestry, say, keeps their
  # order).
  for (d in names(declaration)) {
    if (is.null(declaration[[d]]$categories) && !is.null(population$groups[[
      d
    ]]$categories)) {
      declaration[[d]]["categories"] <- list(population$groups[[d]]$categories)
    }
  }
  population$groups <- declaration
  # Keep ancestry only for dimensions that track it.
  if (!is.null(population$ancestry)) {
    ancestry <- as_ancestry_list(population$ancestry, declaration)
    ancestry <- ancestry[intersect(names(ancestry), tracked_ancestry(
      declaration
    ))]
    population$ancestry <- if (length(ancestry)) ancestry else NULL
  }
  validate_group_declaration(declaration, population)
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
      name %in% names(
        population$groups
      )
    })
}

# The name of a population's single declared dimension, for features that use
# one dimension (transitions, matching, observation by group); NULL otherwise.
single_dimension <- function(population) {
  if (inherits(population, "hammel_population") && length(population$groups) ==
    1L) {
    names(population$groups)
  } else {
    NULL
  }
}

# Check a declaration against a population's persons and transition history.
validate_group_declaration <- function(declaration, population) {
  persons <- population$persons
  for (name in names(declaration)) {
    if (!name %in% names(persons)) {
      stop(sprintf("Declared group `%s` is not a persons column.", name),
        call. = FALSE
      )
    }
    check_groups(tibble::tibble(group = persons[[name]]), paste0(
      "persons$",
      name
    ))
    categories <- declaration[[name]]$categories
    # Ancestry tables name their columns by category.
    if (isTRUE(declaration[[name]]$ancestry) &&
      any(c("unknown", "person_id") %in% as.character(c(
        persons[[name]],
        categories
      )))) {
      stop(sprintf(
        paste0(
          "Group `%s` tracks ancestry, so its categories ",
          "cannot be named `unknown` and `person_id`."
        ),
        name
      ), call. = FALSE)
    }
    if (!is.null(categories)) {
      tx <- population$transitions
      if (!is.null(tx) && "group" %in% names(tx)) {
        tx <- tx[tx$group == name, ,
          drop = FALSE
        ]
      }
      values <- unique(as.character(c(
        persons[[name]], tx$origin,
        tx$destination
      )))
      missing <- setdiff(values, as.character(categories))
      if (length(missing)) {
        stop(sprintf(
          "Group `%s` has values missing from its categories: %s.",
          name,
          paste(missing, collapse = ", ")
        ), call. = FALSE)
      }
    }
  }
  if (!length(declaration) && !is.null(population$transitions) && nrow(
    population$transitions
  )) {
    stop("Transition histories require a declared group dimension.",
      call. =
        FALSE
    )
  }
  if (!is.null(population[["ancestry"]])) {
    tracked <- tracked_ancestry(declaration)
    if (!length(declaration)) {
      stop(
        paste0(
          "Ancestry requires a declared group dimension with ",
          "`ancestry = TRUE`."
        ),
        call. = FALSE
      )
    }
    if (!length(tracked)) {
      stop(paste0(
        "An ancestry table requires `ancestry = TRUE` for ",
        "its group dimension."
      ), call. = FALSE)
    }
    ancestry <- as_ancestry_list(population[["ancestry"]], declaration)
    extra <- setdiff(names(ancestry), tracked)
    if (length(extra)) {
      stop(sprintf(
        paste0(
          "Ancestry is supplied for %s, which does not track ",
          "ancestry (`ancestry = TRUE`)."
        ),
        paste0("`", extra, "`", collapse = ", ")
      ), call. = FALSE)
    }
  }
  invisible(declaration)
}
