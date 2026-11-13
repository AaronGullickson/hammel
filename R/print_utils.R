# Functions useful for printing methods

# HELPER FUNCTIONS ------------------------------------------------------------

# format the printing of a stat
print_stat <- function(label, value) {
  cat(sprintf("%-28s %12s\n", label, value))
}

# Print metadata returned by summarize_rate_strata(), invisibly returning it.
# Band axes are printed separately; exact combinations remain in the metadata.
# NULL represents an unconfigured schedule and produces no output.
print_rate_strata <- function(x) {
  if (is.null(x)) {
    return(invisible(x))
  }

  cat("\n", cli::style_bold("Stratified by"), "\n", sep = "")
  if (length(x$stratified_by)) {
    width <- max(nchar(x$stratified_by, type = "width"))
    for (column in x$stratified_by) {
      values <- paste(unique(x$strata[[column]]), collapse = ", ")
      label <- format(column, width = width, justify = "left")
      cat("  ", label, "  ", values, "\n", sep = "")
    }
  } else {
    cat("  None\n")
  }

  if (!nrow(x$bands)) {
    cat("\nNo bands supplied.\n")
    return(invisible(x))
  }
  if (!x$bands_by_stratum) {
    print_rate_bands(x$bands)
    return(invisible(x))
  }

  cat("\n", cli::style_bold("Bands by stratum"), "\n", sep = "")
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
    cat("\n  ", cli::style_bold(paste(labels, collapse = "; ")),
      "\n", sep = ""
    )
    print_rate_bands(x$bands[rows, ], indent = "    ")
  }
  invisible(x)
}

# Group calendar intervals by their exact non-calendar band combinations.
# Nonadjacent periods can share a structure; their intervals remain separate.
print_rate_bands <- function(bands, indent = "  ") {
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
      cat("\n", indent, cli::style_bold("Bands by calendar period"),
        "\n", sep = ""
      )
      for (i in seq_along(structures)) {
        print_rate_band_axes(periods[period_groups[[i]], ], indent)
        print_rate_band_axes(structures[[i]], paste0(indent, "  "))
      }
      return(invisible(bands))
    }
  }
  print_rate_band_axes(bands, indent)
  invisible(bands)
}

# Print marginal band sets, rather than expanding age by calendar time.
# These lists describe bands present, not every possible joint combination.
print_rate_band_axes <- function(bands, indent = "  ") {
  axes <- c(
    age = "Age bands (years)",
    dur = "Duration bands (years)",
    parity = "Parity bands",
    time = "Calendar bands (years)"
  )
  printed <- FALSE
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
    text <- format_rate_bands(pairs[[1]], pairs[[2]])
    cat("\n", indent, cli::style_bold(axes[[axis]]), "\n", sep = "")
    cat(paste0(indent, "  ", text), "\n", sep = "")
    printed <- TRUE
  }
  if (!printed) {
    cat(indent, "No bands supplied.\n", sep = "")
  }
  invisible(bands)
}

# Abbreviate only contiguous, equal-width finite bands as regular sequences.
# Other long sets show the first three and last two bands, with a total.
format_rate_bands <- function(start, end) {
  labels <- paste0("[", start, ", ", end, ")")
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
