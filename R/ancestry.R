# Validate, normalize, and prepare ancestry shares for populations and
# simulations.

# HELPER FUNCTIONS ---------------------------------------------------------

# Ancestry shares are independent of current group identity. Each dimension
# declared with `ancestry = TRUE` has its own table (D134): `person_id`, one
# column per category and `unknown`. Missing ancestry belongs in the explicit
# unknown component, never in NA values or missing rows.

# Dimensions whose declaration tracks ancestry, in declaration order.
tracked_ancestry <- function(declaration) {
  names(declaration)[vapply(
    declaration, function(g) isTRUE(g$ancestry),
    logical(1)
  )]
}

# Accept one bare table when exactly one dimension tracks ancestry.
as_ancestry_list <- function(ancestry, declaration) {
  if (is.null(ancestry) || !is.data.frame(ancestry)) return(ancestry)
  tracked <- tracked_ancestry(declaration)
  if (length(tracked) != 1L)
    stop(
      paste0(
        "`ancestry` must be a named list of tables, one ",
        "per group dimension tracking ancestry."
      ),
      call. = FALSE
    )
  stats::setNames(list(ancestry), tracked)
}

validate_ancestry <- function(ancestry, persons, declaration = NULL) {
  if (is.null(ancestry)) return(invisible(ancestry))
  ancestry <- as_ancestry_list(ancestry, declaration)
  if (!is.list(ancestry) || is.null(names(ancestry)) || anyNA(names(
    ancestry
  )) ||
    any(!nzchar(names(ancestry))) || anyDuplicated(names(ancestry)) ||
    !all(vapply(ancestry, is.data.frame, logical(1))))
    stop(
      paste0(
        "`ancestry` must be a named list of tables, one ",
        "per group dimension tracking ancestry."
      ),
      call. = FALSE
    )
  for (name in names(ancestry)) validate_ancestry_table(
    ancestry[[name]],
    persons
  )
  invisible(ancestry)
}

validate_ancestry_table <- function(ancestry, persons) {
  check_columns(ancestry, c("person_id", "unknown"), "ancestry")
  check_numeric_column(ancestry$person_id, "ancestry$person_id", id = TRUE)
  if (anyDuplicated(ancestry$person_id) || nrow(ancestry) != nrow(persons) ||
    !setequal(ancestry$person_id, persons$person_id)) {
    stop(paste0(
      "`ancestry` must have exactly one row per person, ",
      "with matching person IDs."
    ), call. = FALSE)
  }
  columns <- setdiff(names(ancestry), c("person_id", "unknown"))
  if (any(!nzchar(trimws(columns)))) {
    stop("Known ancestry columns must be named by nonempty category labels.",
      call. = FALSE
    )
  }
  for (column in c(columns, "unknown")) {
    values <- ancestry[[column]]
    if (!is.numeric(values) || is.complex(values) || !is.null(dim(values)) ||
      anyNA(values) || any(!is.finite(values)) || any(values < 0 | values >
      1)) {
      stop(
        paste0(
          "Ancestry proportions must be finite numeric ",
          "values between 0 and 1, without missing values."
        ),
        call. = FALSE
      )
    }
  }
  sums <- rowSums(as.matrix(ancestry[c(columns, "unknown")]))
  if (any(abs(sums - 1) > 1e-8)) {
    stop(paste0(
      "Ancestry proportions, including `unknown`, must ",
      "sum to 1 for each person."
    ), call. = FALSE)
  }
  invisible(ancestry)
}

normalize_population_ancestry <- function(population) {
  persons <- population$persons
  tracked <- tracked_ancestry(population$groups)
  if (!length(tracked)) {
    population$ancestry <- NULL
    return(population)
  }
  ancestry <- as_ancestry_list(population[["ancestry"]], population$groups)
  out <- list()
  for (name in tracked) {
    x <- ancestry[[name]]
    if (is.null(x)) {
      # Without shares, people start with all their ancestry in their category.
      labels <- as.character(persons[[name]])
      x <- tibble::tibble(person_id = persons$person_id)
      for (category in unique(labels)) x[[category]] <- as.numeric(labels ==
        category)
      x$unknown <- rep(0, nrow(persons))
    } else {
      columns <- setdiff(names(x), c("person_id", "unknown"))
      x <- x[match(persons$person_id, x$person_id), c(
        "person_id", columns,
        "unknown"
      ), drop = FALSE]
    }
    out[[name]] <- tibble::as_tibble(x)
  }
  population$ancestry <- out
  population
}

# Initial ancestry shares for the engine: each tracked dimension's block of
# category shares and unknown, in declaration order.
ancestry_matrix <- function(population, tracked) {
  n <- nrow(population$persons)
  if (!length(tracked)) return(matrix(numeric(), nrow = n, ncol = 0L))
  do.call(cbind, lapply(tracked, function(d) {
    x <- population$ancestry[[d]]
    as.matrix(x[setdiff(names(x), "person_id")])
  }))
}
