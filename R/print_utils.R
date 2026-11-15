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

# Print a component model summary on its own: a level-2 heading and one
# bullet per line from model_summary_lines().
print_model_summary <- function(title, x) {
  cli::cli_h2(title)
  list <- open_summary_list()
  lines <- model_summary_lines(x)
  for (i in seq_along(lines)) {
    print_stat(names(lines)[i], lines[[i]])
  }
  cli::cli_end(list)
  invisible(x)
}

# Summarize a process's model argument: NULL, one model, or a calendar
# schedule (a tibble with time_start, time_end and a list column of models).
summarize_model_argument <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is_model_schedule(x)) {
    return(summary(x))
  }
  structure(
    list(
      periods = tibble::as_tibble(x[c("time_start", "time_end")]),
      models = lapply(x$model, summary)
    ),
    class = "summary.hammel_model_schedule"
  )
}

# Print a summarized model argument as a bullet in a process summary. A model
# with one line prints on one bullet; otherwise its first line is the
# bullet's value and the others are nested beneath it. Each period of a
# schedule is nested in turn. `none` describes a NULL argument.
print_model_argument <- function(label, x, none = "none") {
  if (is.null(x)) {
    print_stat(label, none)
    return(invisible(x))
  }
  if (!inherits(x, "summary.hammel_model_schedule")) {
    print_model_lines(label, model_summary_lines(x))
    return(invisible(x))
  }
  schedule <- open_nested_stat(label, "calendar schedule")
  for (i in seq_along(x$models)) {
    period <- format_rate_bands(
      x$periods$time_start[i],
      x$periods$time_end[i]
    )
    print_model_lines(
      paste("Years", period),
      model_summary_lines(x$models[[i]])
    )
  }
  cli::cli_end(schedule)
  invisible(x)
}

# Print one model's lines under `label`; see print_model_argument().
print_model_lines <- function(label, lines) {
  if (length(lines) == 1L) {
    print_stat(label, lines[[1]])
    return(invisible(lines))
  }
  item <- open_nested_stat(label, lines[[1]])
  for (i in seq_along(lines)[-1]) {
    print_stat(names(lines)[i], lines[[i]])
  }
  cli::cli_end(item)
  invisible(lines)
}

# Named character vector describing a component model summary, one element
# per line, labels as names; the first element is its headline.
model_summary_lines <- function(x) {
  if (inherits(x, "summary.hammel_matching_model")) {
    return(matching_summary_lines(x))
  }
  if (inherits(x, "summary.hammel_weight_model")) {
    return(weight_summary_lines(x))
  }
  if (inherits(x, "summary.hammel_group_inheritance")) {
    return(inheritance_summary_lines(x))
  }
  abort_internal("Unsupported model summary.")
}

# Format named numbers as "name = value", with optional parenthetical notes.
format_named_values <- function(values, notes = NULL) {
  if (!length(values)) {
    return("none")
  }
  text <- paste0(names(values), " = ", signif(values, 4))
  if (!is.null(notes)) {
    text <- paste0(text, " (", notes, ")")
  }
  paste(text, collapse = ", ")
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
