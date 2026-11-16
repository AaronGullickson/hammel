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
as_ancestry_list <- function(
  ancestry,
  declaration,
  call = rlang::caller_env()
) {
  if (is.null(ancestry) || !is.data.frame(ancestry)) {
    return(ancestry)
  }
  tracked <- tracked_ancestry(declaration)
  if (length(tracked) != 1L) {
    abort_hammel(
      c(
        "{.arg ancestry} must be a named list of tables, one per dimension.",
        "x" = paste(
          "It is a single table, but {length(tracked)} dimensions track",
          "ancestry."
        ),
        "i" = paste(
          "A single table is allowed only when exactly one dimension tracks",
          "ancestry."
        )
      ),
      "population",
      argument = "ancestry",
      call = call
    )
  }
  stats::setNames(list(ancestry), tracked)
}

validate_ancestry <- function(
  ancestry,
  persons,
  declaration = NULL,
  call = rlang::caller_env()
) {
  if (is.null(ancestry)) {
    return(invisible(ancestry))
  }
  ancestry <- as_ancestry_list(ancestry, declaration, call = call)
  if (!is.list(ancestry) || is.null(names(ancestry)) ||
        anyNA(names(ancestry)) ||
        any(!nzchar(names(ancestry))) ||
        anyDuplicated(names(ancestry)) ||
        !all(vapply(ancestry, is.data.frame, logical(1)))) {
    abort_hammel(
      c(
        paste(
          "{.arg ancestry} must be a list of tables named by the dimensions",
          "that track ancestry."
        ),
        "x" = "It is {describe_value(ancestry)}."
      ),
      "population",
      argument = "ancestry",
      call = call
    )
  }
  for (name in names(ancestry)) {
    with_error_context(
      validate_ancestry_table(ancestry[[name]], persons, call = call),
      paste0("In the ancestry table for ", name, ".")
    )
  }
  invisible(ancestry)
}

validate_ancestry_table <- function(
  ancestry,
  persons,
  call = rlang::caller_env()
) {
  check_columns(ancestry, c("person_id", "unknown"), "ancestry", call = call)
  check_numeric_column(
    ancestry$person_id, "ancestry$person_id",
    id = TRUE,
    call = call
  )
  check_unique_ids(ancestry$person_id, "ancestry$person_id", call = call)
  absent <- setdiff(persons$person_id, ancestry$person_id)
  extra <- setdiff(ancestry$person_id, persons$person_id)
  if (length(absent) || length(extra)) {
    absent <- cli::cli_vec(absent, list("vec-trunc" = 5))
    extra <- cli::cli_vec(extra, list("vec-trunc" = 5))
    abort_hammel(
      c(
        "{.arg ancestry} must have exactly one row per person.",
        "x" = if (length(absent)) "No row for persons {.val {absent}}.",
        "x" = if (length(extra)) "Rows for unknown persons {.val {extra}}."
      ),
      "population",
      argument = "ancestry",
      call = call
    )
  }
  columns <- setdiff(names(ancestry), c("person_id", "unknown"))
  if (any(!nzchar(trimws(columns)))) {
    abort_hammel(
      "Ancestry columns must be named by nonempty category labels.",
      "population",
      argument = "ancestry",
      call = call
    )
  }
  for (column in c(columns, "unknown")) {
    values <- ancestry[[column]]
    check_numeric_column(values, paste0("ancestry$", column), call = call)
    rows <- which(values < 0 | values > 1)
    if (length(rows)) {
      abort_hammel(
        c(
          "Ancestry shares must be between 0 and 1.",
          "x" = paste(
            "Column {.field {column}} is out of range at",
            "{format_rows(rows)}."
          )
        ),
        "population",
        argument = "ancestry",
        rows = rows,
        call = call
      )
    }
  }
  sums <- rowSums(as.matrix(ancestry[c(columns, "unknown")]))
  rows <- which(abs(sums - 1) > 1e-8)
  if (length(rows)) {
    abort_hammel(
      c(
        paste(
          "Each person's ancestry shares, including {.field unknown}, must",
          "sum to 1."
        ),
        "x" = "Other sums at {format_rows(rows)}."
      ),
      "population",
      argument = "ancestry",
      rows = rows,
      call = call
    )
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
      for (category in unique(labels)) {
        x[[category]] <- as.numeric(labels == category)
      }
      x$unknown <- rep(0, nrow(persons))
    } else {
      columns <- setdiff(names(x), c("person_id", "unknown"))
      x <- x[
        match(persons$person_id, x$person_id),
        c("person_id", columns, "unknown"),
        drop = FALSE
      ]
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
  if (!length(tracked)) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  do.call(
    cbind, lapply(
      tracked, function(d) {
        x <- population$ancestry[[d]]
        as.matrix(x[setdiff(names(x), "person_id")])
      }
    )
  )
}
