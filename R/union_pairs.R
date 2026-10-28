# Symmetric labels shared by dissolution schedules and union observations.
# Pair strings follow the single-colon convention used for matching preferences.
canonical_union_pairs <- function(first, second = NULL, name = "group_pair") {
  if (is.null(second)) {
    if (!is.character(first) || !is.null(dim(first)) || anyNA(first) ||
        any(!grepl("^[^:]+:[^:]+$", first)))
      stop(sprintf("`%s` must use single-colon pair labels such as 'A:A' or 'A:B'.", name), call. = FALSE)
    parts <- strsplit(first, ":", fixed = TRUE)
    first <- vapply(parts, `[`, character(1), 1L)
    second <- vapply(parts, `[`, character(1), 2L)
  }
  first <- as.character(first)
  second <- as.character(second)
  if (length(first) != length(second) || anyNA(first) || anyNA(second) ||
      any(!nzchar(trimws(c(first, second)))) || any(grepl(":", c(first, second), fixed = TRUE)))
    stop("Union pair labels require nonempty, colon-free group labels.", call. = FALSE)
  if (!length(first)) return(character())
  labels <- sort(unique(c(first, second)), method = "radix")
  ordered <- match(first, labels) <= match(second, labels)
  ifelse(ordered, paste0(first, ":", second), paste0(second, ":", first))
}

possible_union_pairs <- function(groups) {
  labels <- as.character(groups)
  pairs <- expand.grid(first = labels, second = labels,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sort(unique(canonical_union_pairs(pairs$first, pairs$second)), method = "radix")
}

validate_dissolution_multipliers <- function(factors) {
  validate_group_prop_hazard(factors)
  if (is.null(factors)) return(invisible(factors))
  pairs <- canonical_union_pairs(names(factors), name = "group_prop_hazard")
  if (anyDuplicated(pairs))
    stop("`group_prop_hazard` must specify each unordered pair once, not both A:B and B:A.", call. = FALSE)
  invisible(factors)
}
