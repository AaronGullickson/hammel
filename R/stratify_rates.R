#' Add proportional strata to a rate table
#'
#' Copy each row of a rate table once per category of one or more new columns,
#' multiplying its `rate` by that category's factor. The result is an ordinary
#' stratified rate table for any rate process, such as group-specific mortality
#' that is 20% higher in one group than in a reference group.
#'
#' Several columns in `by` are crossed: every combination of categories gets
#' the product of its factors. All other columns, including calendar bands
#' (`time_start`, `time_end`), are copied unchanged. Every category a process
#' needs must be listed; a population group omitted from `by` is a coverage
#' error at simulation, as for any group-specific table. Factors of one still
#' create separate strata.
#'
#' A column that already exists in `rates` is an error. To give different
#' effects within existing strata, for example different group effects for
#' partnered and unpartnered fertility, stratify each subset and bind the
#' results with [rbind()].
#'
#' Columns in `by` are either built-in strata (`sex`, `union_type`,
#' `previous_union`, `origin`, `destination`) or group dimensions
#' declared with [set_groups()]: a dimension's name (such as `race`) for
#' person-level rates, or `<name>_w` and `<name>_h` (the wife's and husband's
#' categories) for dissolution and conversion rates. Band, rate and calendar
#' columns cannot be added. Process setters reject built-in strata their
#' process does not support: mortality accepts `sex`; fertility `union_type`;
#' union searches `previous_union`; dissolution `union_type`; group
#' transitions `origin`, `destination`, `sex` and `union_type`. Group
#' columns are checked against the declaration at simulation.
#'
#' @param rates A data frame with a numeric `rate` column, such as a mortality,
#'   fertility, search, conversion, dissolution or transition rate table.
#' @param by A named list. Each name is a column to add; each element is a
#'   named numeric vector of finite, nonnegative factors, one per category.
#'   Factors are not normalized, and zero suppresses that stratum's risk.
#' @return A tibble with the new columns (character labels) and multiplied
#'   rates. Other attributes of `rates`, such as provenance, are kept.
#' @export
#' @examples
#' mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0,
#'   age_end = Inf, rate = c(0.01, 0.008))
#' # Race A has 20% higher mortality than reference race B.
#' stratify_rates(mortality_rates, list(race = c(A = 1.2, B = 1)))
#'
#' # Unions of an A wife and a B husband dissolve twice as often.
#' dissolution <- data.frame(dur_start = 0, dur_end = Inf, rate = 0.02)
#' stratify_rates(dissolution, list(race_w = c(A = 1, B = 1), race_h = c(A = 2, B = 1)))
stratify_rates <- function(rates, by) {
  if (!is.data.frame(rates) || !"rate" %in% names(rates) || !is.numeric(rates$rate))
    stop("`rates` must be a data frame with a numeric `rate` column.", call. = FALSE)
  validate_strata_factors(by)
  existing <- intersect(names(by), names(rates))
  if (length(existing))
    stop(sprintf("`rates` already has column(s) %s; stratify a table without them, or stratify subsets and bind them.",
                 paste0("`", existing, "`", collapse = ", ")), call. = FALSE)
  # One row per combination of categories, the first column varying slowest.
  grid <- expand.grid(lapply(rev(by), names), KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)[rev(names(by))]
  factor <- Reduce(`*`, lapply(names(by), function(col) unname(by[[col]][grid[[col]]])))
  n <- nrow(rates)
  rows <- rep(seq_len(n), times = nrow(grid))
  out <- as.data.frame(rates)[rows, , drop = FALSE]
  combination <- rep(seq_len(nrow(grid)), each = n)
  for (col in names(by)) out[[col]] <- grid[[col]][combination]
  out$rate <- out$rate * factor[combination]
  if (any(is.infinite(out$rate)))
    stop("Stratified rates must remain finite.", call. = FALSE)
  rownames(out) <- NULL
  out <- tibble::as_tibble(out)
  # Keep provenance and other attributes of the input table.
  extra <- setdiff(names(attributes(rates)), c("names", "row.names", "class"))
  for (a in extra) attr(out, a) <- attr(rates, a)
  out
}

# Built-in categorical strata of rate tables; each process accepts some.
builtin_strata <- c("sex", "union_type", "previous_union", "origin", "destination")

# Rate-table columns that are not strata: bands, rates and calendar periods.
rate_schema_columns <- c("age_start", "age_end", "dur_start", "dur_end", "rate", "prob",
  "time_start", "time_end", "parity_start", "parity_end")

# Columns of a rate table that may be group strata: everything else. Group
# dimensions are known only once the table meets the population's declaration.
candidate_strata <- function(rates) {
  setdiff(names(rates), c(rate_schema_columns, builtin_strata, "union_status",
                          "marital_status", ".hammel_group"))
}

validate_strata_factors <- function(by) {
  columns <- names(by)
  if (!is.list(by) || is.data.frame(by) || !length(by) || is.null(columns) ||
      anyNA(columns) || any(!nzchar(trimws(columns))) || anyDuplicated(columns))
    stop("`by` must be a nonempty named list with unique column names.", call. = FALSE)
  schema <- intersect(columns, rate_schema_columns)
  if (length(schema))
    stop(sprintf("`by` cannot stratify by %s, which are band, rate or calendar columns.",
                 paste0("`", schema, "`", collapse = ", ")), call. = FALSE)
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", columns)))
    stop("`by` names must be valid column names.", call. = FALSE)
  for (col in columns) {
    x <- by[[col]]
    labels <- names(x)
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || anyNA(x) ||
        any(!is.finite(x)) || any(x < 0) || is.null(labels) || anyNA(labels) ||
        any(!nzchar(trimws(labels))) || anyDuplicated(labels))
      stop(sprintf("`by$%s` must be a named numeric vector of finite nonnegative factors with unique nonempty category names.",
                   col), call. = FALSE)
  }
  invisible(by)
}

# Setters call this with the stratum columns their process supports. It checks
# the table's actual columns, whether stratified by hand or by stratify_rates().
check_stratified_columns <- function(rates, allowed, process) {
  bad <- intersect(setdiff(builtin_strata, allowed), names(rates))
  if (length(bad))
    stop(sprintf("%s rates cannot be stratified by %s; supported columns: %s.", process,
                 paste0("`", bad, "`", collapse = ", "),
                 paste0("`", allowed, "`", collapse = ", ")), call. = FALSE)
  invisible(rates)
}
