# Functions useful for printing methods

# HELPER FUNCTIONS ------------------------------------------------------------

# Open the bullet list used by summary print methods. Wrapped lines are
# indented under their bullet text. Returns the id of the enclosing div;
# cli::cli_end(id) closes the list and its theme.
open_summary_list <- function() {
  div <- cli::cli_div(
    theme = list(li = list("text-exdent" = 2)),
    .auto_close = FALSE
  )
  cli::cli_ul(.auto_close = FALSE)
  div
}

# Print one statistic as a cli bullet in the current list. Values are
# interpolated, not parsed, so braces in category names print literally.
print_stat <- function(label, value) {
  cli::cli_li("{label}: {value}", .auto_close = FALSE)
}

# Open a bullet whose following bullets are nested beneath it. Returns the
# item id; close it, and everything inside it, with cli::cli_end(id).
open_nested_stat <- function(label, value = NULL) {
  item <- if (is.null(value)) {
    cli::cli_li("{label}", .auto_close = FALSE)
  } else {
    cli::cli_li("{label}: {value}", .auto_close = FALSE)
  }
  cli::cli_ul(.auto_close = FALSE)
  item
}

# Print metadata returned by summarize_rate_strata() as bullets, invisibly
# returning it. Band axes are printed separately; exact combinations remain
# in the metadata. NULL represents an unconfigured schedule and prints
# nothing.
print_rate_strata <- function(x) {
  if (is.null(x)) {
    return(invisible(x))
  }

  print_stat("Stratified by", format_rate_strata(x))

  if (!nrow(x$bands)) {
    print_stat("Bands", "none supplied")
    return(invisible(x))
  }
  if (!x$bands_by_stratum) {
    print_rate_bands(x$bands)
    return(invisible(x))
  }

  by_stratum <- open_nested_stat("Bands by stratum")
  for (i in seq_len(nrow(x$strata))) {
    rows <- rep(TRUE, nrow(x$bands))
    labels <- character(length(x$stratified_by))
    for (j in seq_along(x$stratified_by)) {
      column <- x$stratified_by[j]
      value <- x$strata[[column]][i]
      same <- if (is.na(value)) {
        is.na(x$bands[[column]])
      } else {
        !is.na(x$bands[[column]]) & x$bands[[column]] == value
      }
      rows <- rows & same
      labels[j] <- paste0(column, " = ", value)
    }
    stratum <- open_nested_stat(paste(labels, collapse = "; "))
    print_rate_bands(x$bands[rows, ])
    cli::cli_end(stratum)
  }
  cli::cli_end(by_stratum)
  invisible(x)
}

# Describe stratum columns and their observed values on one line.
format_rate_strata <- function(x) {
  if (!length(x$stratified_by)) {
    return("none")
  }
  columns <- vapply(
    x$stratified_by,
    function(column) {
      values <- paste(unique(x$strata[[column]]), collapse = ", ")
      paste0(column, " (", values, ")")
    },
    character(1)
  )
  paste(columns, collapse = "; ")
}

# Group calendar intervals by their exact non-calendar band combinations.
# Nonadjacent periods can share a structure; their intervals remain separate.
print_rate_bands <- function(bands) {
  time_columns <- c("time_start", "time_end")
  if (all(time_columns %in% names(bands)) && nrow(bands)) {
    columns <- intersect(
      c(
        "age_start", "age_end", "dur_start", "dur_end",
        "parity_start", "parity_end"
      ),
      names(bands)
    )
    periods <- unique(bands[time_columns])
    periods <- periods[order(periods$time_start, periods$time_end), ]
    structures <- list()
    period_groups <- list()
    for (i in seq_len(nrow(periods))) {
      rows <- bands$time_start == periods$time_start[i] &
        bands$time_end == periods$time_end[i]
      structure <- unique(as.data.frame(bands[rows, columns]))
      if (length(columns)) {
        structure <- structure[do.call(order, structure), , drop = FALSE]
      }
      rownames(structure) <- NULL
      matches <- vapply(structures, identical, logical(1), structure)
      if (any(matches)) {
        group <- which(matches)[1]
        period_groups[[group]] <- c(period_groups[[group]], i)
      } else {
        structures[[length(structures) + 1L]] <- structure
        period_groups[[length(period_groups) + 1L]] <- i
      }
    }
    if (length(structures) > 1L) {
      by_period <- open_nested_stat("Bands by calendar period")
      for (i in seq_along(structures)) {
        calendar <- format_rate_band_axes(periods[period_groups[[i]], ])
        period <- open_nested_stat(names(calendar), calendar[[1]])
        print_rate_band_axes(structures[[i]])
        cli::cli_end(period)
      }
      cli::cli_end(by_period)
      return(invisible(bands))
    }
  }
  print_rate_band_axes(bands)
  invisible(bands)
}

# Print marginal band sets, rather than expanding age by calendar time.
# These lists describe bands present, not every possible joint combination.
print_rate_band_axes <- function(bands) {
  axes <- format_rate_band_axes(bands)
  if (!length(axes)) {
    print_stat("Bands", "none supplied")
  }
  for (axis in names(axes)) {
    print_stat(axis, axes[[axis]])
  }
  invisible(bands)
}

# Named character vector of formatted marginal band sets, one per axis
# present, named by its display label.
format_rate_band_axes <- function(bands) {
  axes <- c(
    age = "Age bands (years)",
    dur = "Duration bands (years)",
    parity = "Parity bands",
    time = "Calendar bands (years)"
  )
  out <- character()
  for (axis in names(axes)) {
    columns <- paste0(axis, c("_start", "_end"))
    if (!all(columns %in% names(bands))) {
      next
    }
    pairs <- unique(bands[columns])
    pairs <- pairs[stats::complete.cases(pairs), ]
    if (!nrow(pairs)) {
      next
    }
    pairs <- pairs[order(pairs[[1]], pairs[[2]]), ]
    out[[axes[[axis]]]] <- format_rate_bands(pairs[[1]], pairs[[2]])
  }
  out
}

# Abbreviate only contiguous, equal-width finite bands as regular sequences.
# Other long sets show the first three and last two bands, with a total.
# A non-breaking space inside each band keeps cli from wrapping within it.
format_rate_bands <- function(start, end) {
  labels <- paste0("[", start, ",\u00a0", end, ")")
  n <- length(labels)
  if (n <= 8L) {
    return(paste(labels, collapse = ", "))
  }
  widths <- end - start
  regular <- all(is.finite(c(start, end))) &&
    isTRUE(all.equal(widths, rep(widths[1], n))) &&
    isTRUE(all.equal(head(end, -1), tail(start, -1)))
  if (regular) {
    shown <- c(labels[1:2], "...", labels[n])
    detail <- paste0(n, " bands; width ", widths[1])
  } else {
    shown <- c(head(labels, 3), "...", tail(labels, 2))
    detail <- paste0(n, " bands; showing 5")
  }
  paste0(paste(shown, collapse = ", "), " (", detail, ")")
}
