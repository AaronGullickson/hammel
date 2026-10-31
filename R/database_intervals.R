# Source intervals for native HMD products; no interpolation.
check_database_years <- function(year) {
  if (!is.numeric(year) || !is.null(dim(year)) || !length(year) ||
      anyNA(year) || any(!is.finite(year)) || any(year != trunc(year)) ||
      anyDuplicated(year))
    stop("`year` must contain distinct whole-number calendar years.", call. = FALSE)
  invisible(year)
}

check_hmd_shape <- function(shape) {
  shapes <- c("1x1", "1x5", "1x10", "5x1", "5x5", "5x10")
  if (!is.character(shape) || length(shape) != 1L || is.na(shape) ||
      !is.null(dim(shape)) || !shape %in% shapes)
    stop(sprintf("`shape` must be one of: %s.", paste(shapes, collapse = ", ")),
         call. = FALSE)
  invisible(shape)
}

hmd_period_intervals <- function(data, shape = "1x1", coverage_years = NULL) {
  check_hmd_shape(shape)
  check_columns(data, "Year", "HMD data")
  labels <- trimws(as.character(data$Year))
  if (!length(labels) || anyNA(labels) ||
      any(!grepl("^[0-9]+(-[0-9]+)?$", labels)))
    stop("Expected HMD calendar years or inclusive source year ranges.", call. = FALSE)
  width <- as.numeric(strsplit(shape, "x", fixed = TRUE)[[1L]][2L])
  starts <- as.numeric(sub("-.*$", "", labels))
  ranged <- grepl("-", labels, fixed = TRUE)
  ends <- starts + width
  ends[ranged] <- as.numeric(sub("^.*-", "", labels[ranged])) + 1
  if (width > 1 && any(!ranged)) {
    if (is.null(coverage_years))
      stop("Parsed grouped HMD periods require annual coverage to recover source bounds.", call. = FALSE)
    check_database_years(coverage_years)
    periods <- sort(unique(starts))
    following <- c(tail(periods, -1L), Inf)[match(starts, periods)]
    ends[!ranged] <- pmin(ends[!ranged], following[!ranged], max(coverage_years) + 1)
  }
  if (any(!is.finite(starts)) || any(!is.finite(ends)) ||
      any(ends <= starts | ends - starts > width))
    stop("HMD source periods are inconsistent with the requested shape.", call. = FALSE)
  periods <- unique(tibble::tibble(time_start = starts, time_end = ends))
  periods <- periods[order(periods$time_start), ]
  if (anyDuplicated(periods$time_start) || (nrow(periods) > 1L &&
      any(tail(periods$time_start, -1L) < head(periods$time_end, -1L))))
    stop("HMD source periods must not overlap or have conflicting bounds.", call. = FALSE)
  if (!is.null(coverage_years)) {
    for (i in seq_len(nrow(periods))) {
      years <- seq(periods$time_start[i], periods$time_end[i] - 1)
      if (any(!years %in% coverage_years))
        stop("HMD source periods are inconsistent with annual coverage.", call. = FALSE)
    }
  }
  tibble::tibble(time_start = starts, time_end = ends)
}

hmd_age_intervals <- function(data, shape = "1x1") {
  check_hmd_shape(shape)
  ages <- trimws(as.character(data$Age))
  width <- as.numeric(strsplit(shape, "x", fixed = TRUE)[[1L]][1L])
  pattern <- if (width == 1) "^[0-9]+[+]?$" else "^[0-9]+(-[0-9]+|[+])?$"
  if (anyNA(ages) || any(!grepl(pattern, ages)))
    stop("Expected HMD age labels consistent with the requested shape.", call. = FALSE)
  starts <- as.numeric(sub("[-+].*$", "", ages))
  open <- grepl("[+]$", ages)
  # readHMDweb() 2.0.8 always parses ages, even with fixup = FALSE.
  if ("OpenInterval" %in% names(data)) {
    if (!is.logical(data$OpenInterval) || anyNA(data$OpenInterval))
      stop("HMD OpenInterval must be a nonmissing logical indicator.", call. = FALSE)
    open <- open | data$OpenInterval
  }
  ends <- if (width == 1) starts + 1 else
    ifelse(starts == 0, 1, ifelse(starts == 1, 5, starts + 5))
  ranged <- grepl("-", ages, fixed = TRUE)
  if (any(ranged)) {
    upper <- as.numeric(sub("^.*-", "", ages[ranged])) + 1
    if (any(upper != ends[ranged]))
      stop("HMD source age ranges are inconsistent with the requested shape.", call. = FALSE)
  }
  if (width == 5 && any(!open & starts > 1 & starts %% 5 != 0))
    stop("HMD five-year ages must use the native infant/child/five-year layout.", call. = FALSE)
  ends[open] <- Inf
  tibble::tibble(age_start = starts, age_end = ends)
}

# Retain source row ordering within each year, return years chronologically.
combine_database_years <- function(year, convert) {
  check_database_years(year)
  years <- sort(year)
  rates <- lapply(years, function(y) {
    rows <- convert(y)
    rows$time_start <- y
    rows$time_end <- y + 1
    rows
  })
  tibble::as_tibble(do.call(rbind, rates))
}
