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
#' Columns in `by` must be stratum columns used by Hammel rate tables:
#' `sex`, `group`, `union_type`, `previous_union`, `group_pair`,
#' `group_origin` or `group_dest`. Each process setter checks the table's
#' columns, however it was built, and rejects stratum columns its process does
#' not support: mortality accepts `sex` and `group`; fertility `group` and
#' `union_type`; union searches `group` and `previous_union`; conversion
#' `group`; dissolution `group_pair` and `union_type`; group transitions
#' `group_origin`, `group_dest`, `sex` and `union_type`.
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
#' # Group A has 20% higher mortality than reference group B.
#' stratify_rates(mortality_rates, list(group = c(A = 1.2, B = 1)))
#'
#' # Exogamous unions dissolve twice as often.
#' dissolution <- data.frame(dur_start = 0, dur_end = Inf, rate = 0.02)
#' stratify_rates(dissolution, list(group_pair = c("A:A" = 1, "A:B" = 2, "B:B" = 1)))
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

# Categorical stratum columns of rate tables across all processes. Other
# columns (custom notes, provenance) are carried along and ignored.
stratum_columns <- c("sex", "group", "union_type", "previous_union", "group_pair",
                     "group_origin", "group_dest")

validate_strata_factors <- function(by) {
  columns <- names(by)
  if (!is.list(by) || is.data.frame(by) || !length(by) || is.null(columns) ||
      anyNA(columns) || any(!nzchar(trimws(columns))) || anyDuplicated(columns))
    stop("`by` must be a nonempty named list with unique column names.", call. = FALSE)
  unknown <- setdiff(columns, stratum_columns)
  if (length(unknown))
    stop(sprintf("`by` can only add stratum columns (%s), not %s.",
                 paste0("`", stratum_columns, "`", collapse = ", "),
                 paste0("`", unknown, "`", collapse = ", ")), call. = FALSE)
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
  bad <- intersect(setdiff(stratum_columns, allowed), names(rates))
  if (length(bad))
    stop(sprintf("%s rates cannot be stratified by %s; supported columns: %s.", process,
                 paste0("`", bad, "`", collapse = ", "),
                 paste0("`", allowed, "`", collapse = ", ")), call. = FALSE)
  invisible(rates)
}
