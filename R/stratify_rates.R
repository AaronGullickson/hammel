# Create proportional strata in input rate tables and check supported stratum
# columns.

# PUBLIC FUNCTIONS ---------------------------------------------------------

#' Add proportional strata to a rate table
#'
#' @description
#'
#' Stratify a rate schedule by additional columns, applying proportional changes
#' to rates by stratum.
#'
#' @details
#'
#' Any of the input rate processes such as fertility, mortality, marriage,
#' transition, dissolution, and conversion can be stratified by this function.
#' Stratified rates can also be further stratified by this function.
#'
#' Columns in `by` are either built-in strata (`sex`, `union_type`,
#' `previous_union`, `origin`, `destination`) or group dimensions declared with
#' [set_groups()]: a dimension's name (such as `race`) for person-level rates,
#' or `<name>_w` and `<name>_h` (the wife's and husband's categories) for
#' dissolution and conversion rates. Band, rate and calendar columns cannot be
#' added. Process setters reject built-in strata their process does not support.
#'
#' If `by` includes multiple columns, every combination of categories gets the
#' product of its factors.
#'
#' @param rates A data frame with a numeric `rate` column, such as a mortality,
#'   fertility, search, conversion, dissolution or transition rate table.
#' @param by A named list. Each name is a column to add; each element is a named
#'   numeric vector of finite, nonnegative factors, one per category. The final
#'   rate will be determined by taking the original rate times the factor
#'   provided.
#' @return A tibble with the new columns (character labels) and modified
#'   rates.
#' @export
#' @examples
#' # stratify fertility rates by group and partnership status
#' usa_2000_fertility |> calculate_tfr()
#' fertility_rates <- stratify_rates(
#'   usa_2000_fertility,
#'   by = list(
#'     group = c(A = 1.1, B = 1),
#'     # this is not really a great idea for union type as within union
#'     # type fertility rates have a totally different age profile, but
#'     # ok for testing purposes
#'     union_type = c("unpartnered" = 0.3, "partnered" = 1.3)
#'   )
#' )
#' fertility_rates |> calculate_tfr()
#'
#' # although this works, a proportional hazards assumption by union type
#' # is unrealistic because the within union type age profile of fertility is
#' # very different from the overall age profile. So for this case, its better
#' # to use already stratified rates by union type and then further stratify
#' # by group
#' fertility_rates <- stratify_rates(
#'   usa_2000_fertility_union,
#'   by = list(group = c(A = 1.1, B = 1))
#' )
#' fertility_rates |> calculate_tfr()
stratify_rates <- function(rates, by) {
  if (!is.data.frame(rates) || !"rate" %in% names(rates) ||
        !is.numeric(rates$rate)) {
    stop(
      "`rates` must be a data frame with a numeric `rate` column.",
      call. = FALSE
    )
  }
  validate_strata_factors(by)
  existing <- intersect(names(by), names(rates))
  if (length(existing)) {
    stop(
      sprintf(
        paste0(
          "`rates` already has column(s) %s; stratify a ",
          "table without them, or stratify subsets and bind them."
        ),
        paste0("`", existing, "`", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # One row per combination of categories, the first column varying slowest.
  grid <- expand.grid(
    lapply(rev(by), names),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )[rev(names(by))]
  factor <- Reduce(
    `*`,
    lapply(names(by), function(col) unname(by[[col]][grid[[col]]]))
  )
  n <- nrow(rates)
  rows <- rep(seq_len(n), times = nrow(grid))
  out <- as.data.frame(rates)[rows, , drop = FALSE]
  combination <- rep(seq_len(nrow(grid)), each = n)
  for (col in names(by)) {
    out[[col]] <- grid[[col]][combination]
  }
  out$rate <- out$rate * factor[combination]
  if (any(is.infinite(out$rate))) {
    stop("Stratified rates must remain finite.", call. = FALSE)
  }
  rownames(out) <- NULL
  out <- tibble::as_tibble(out)
  # Keep provenance and other attributes of the input table.
  extra <- setdiff(names(attributes(rates)), c("names", "row.names", "class"))
  for (a in extra) {
    attr(out, a) <- attr(rates, a)
  }
  out
}

# INTERNAL CONSTANTS -------------------------------------------------------

# Built-in categorical strata of rate tables; each process accepts some.
builtin_strata <- c(
  "sex", "union_type", "previous_union", "origin", "destination"
)

# Rate-table columns that are not strata: bands, rates and calendar periods.
rate_schema_columns <- c(
  "age_start", "age_end", "dur_start", "dur_end", "rate",
  "prob",
  "time_start", "time_end", "parity_start", "parity_end"
)

# HELPER FUNCTIONS ---------------------------------------------------------

# Columns of a rate table that may be group strata: everything else. Group
# dimensions are known only once the table meets the population's declaration.
candidate_strata <- function(rates) {
  setdiff(
    names(rates), c(
      rate_schema_columns, builtin_strata, "union_status",
      "marital_status", ".hammel_group"
    )
  )
}

validate_strata_factors <- function(by) {
  columns <- names(by)
  if (!is.list(by) || is.data.frame(by) || !length(by) || is.null(columns) ||
        anyNA(columns) ||
        any(!nzchar(trimws(columns))) ||
        anyDuplicated(columns)) {
    stop(
      "`by` must be a nonempty named list with unique column names.",
      call. = FALSE
    )
  }
  schema <- intersect(columns, rate_schema_columns)
  if (length(schema)) {
    stop(
      sprintf(
        paste0(
          "`by` cannot stratify by %s, which are band, rate or ",
          "calendar columns."
        ),
        paste0("`", schema, "`", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", columns))) {
    stop("`by` names must be valid column names.", call. = FALSE)
  }
  for (col in columns) {
    x <- by[[col]]
    labels <- names(x)
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || anyNA(x) ||
          any(!is.finite(x)) ||
          any(x < 0) ||
          is.null(labels) ||
          anyNA(labels) ||
          any(!nzchar(trimws(labels))) ||
          anyDuplicated(labels)) {
      stop(
        sprintf(
          paste0(
            "`by$%s` must be a named numeric vector of finite ",
            "nonnegative factors with unique nonempty category names."
          ),
          col
        ),
        call. = FALSE
      )
    }
  }
  invisible(by)
}

# Setters call this with the stratum columns their process supports. It checks
# the table's actual columns, whether stratified by hand or by stratify_rates().
check_stratified_columns <- function(rates, allowed, process) {
  bad <- intersect(setdiff(builtin_strata, allowed), names(rates))
  if (length(bad)) {
    stop(
      sprintf(
        "%s rates cannot be stratified by %s; supported columns: %s.",
        process,
        paste0("`", bad, "`", collapse = ", "),
        paste0("`", allowed, "`", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(rates)
}
