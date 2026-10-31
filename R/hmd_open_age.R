# Explicit source-count pooling for a user-selected open final age category.
check_hmd_open_age <- function(open_age, shape) {
  if (is.null(open_age)) return(invisible(open_age))
  if (!is.numeric(open_age) || length(open_age) != 1L || !is.null(dim(open_age)) ||
      !is.finite(open_age) || open_age < 0 || open_age != trunc(open_age))
    stop("`open_age` must be NULL or a finite nonnegative whole age.", call. = FALSE)
  if (startsWith(shape, "5x") && !open_age %in% c(0, 1) && open_age %% 5 != 0)
    stop("`open_age` must align with a native five-year age-band start: 0, 1, 5, 10, etc.",
         call. = FALSE)
  invisible(open_age)
}

# Diagnose invalid retained Mx with source context, rather than hiding its age
# behind the generic rate validator. Missing replaced-tail Mx never enter here.
check_hmd_retained_rates <- function(values, ages, year, sex) {
  if (!(is.numeric(values) || (is.logical(values) && all(is.na(values)))))
    stop(sprintf("HMD mortality rates must be numeric (year %g, sex %s).", year, sex),
         call. = FALSE)
  bad <- which(!is.finite(values))
  if (length(bad)) {
    i <- bad[1L]
    stop(sprintf(paste0("HMD mortality rates contain missing or non-finite values: ",
      "year %g, sex %s, age [%g, %g)."), year, sex, ages$age_start[i], ages$age_end[i]),
      call. = FALSE)
  }
  invisible(values)
}

index_hmd_pool_source <- function(data, year, shape, coverage_years, name) {
  check_columns(data, c("Year", "Age", "Male", "Female"), name)
  periods <- hmd_period_intervals(data, shape, coverage_years)
  if (any(!year %in% periods$time_start))
    stop(sprintf("%s is missing requested source periods.", name), call. = FALSE)
  rows <- periods$time_start %in% year
  ages <- hmd_age_intervals(data[rows, , drop = FALSE], shape)
  tibble::tibble(time_start = periods$time_start[rows], time_end = periods$time_end[rows],
    age_start = ages$age_start, age_end = ages$age_end,
    Male = data$Male[rows], Female = data$Female[rows])
}

pool_hmd_open_age <- function(ages, year, period_end, open_age, deaths, exposures) {
  bands <- ages[ages$age_start >= open_age, , drop = FALSE]
  bands <- bands[order(bands$age_start), ]
  if (!nrow(bands) || bands$age_start[1L] != open_age ||
      tail(bands$age_end, 1L) != Inf || anyDuplicated(bands$age_start) ||
      (nrow(bands) > 1L && any(head(bands$age_end, -1L) != tail(bands$age_start, -1L))))
    stop(sprintf("HMD rates need a complete source tail from open_age %g through Inf in year %g.",
                 open_age, year), call. = FALSE)
  expected <- paste(bands$age_start, bands$age_end, sep = ":")
  align <- function(data, name) {
    rows <- data$time_start == year & data$age_start >= open_age
    out <- data[rows, ]
    keys <- paste(out$age_start, out$age_end, sep = ":")
    if (anyDuplicated(keys) || !setequal(keys, expected) || any(out$time_end != period_end))
      stop(sprintf("HMD %s tail age/period bands do not match mortality in year %g.", name, year),
           call. = FALSE)
    out[match(expected, keys), ]
  }
  d <- align(deaths, "deaths")
  e <- align(exposures, "exposure")
  totals <- lapply(c(m = "Male", f = "Female"), function(column) {
    sex <- if (column == "Male") "m" else "f"
    for (name in c("deaths", "exposure")) {
      values <- if (name == "deaths") d[[column]] else e[[column]]
      bad <- if (!is.numeric(values)) seq_along(values) else
        which(!is.finite(values) | values < 0)
      if (length(bad)) {
        i <- bad[1L]
        stop(sprintf(paste0("HMD pooled-tail %s must be finite nonnegative numeric values: ",
          "year %g, sex %s, age [%g, %g)."),
          name, year, sex, bands$age_start[i], bands$age_end[i]), call. = FALSE)
      }
    }
    count <- sum(d[[column]])
    exposure <- sum(e[[column]])
    if (!is.finite(count) || !is.finite(exposure) || exposure <= 0)
      stop(sprintf("HMD pooled tail needs finite deaths and positive finite exposure: year %g, sex %s.",
                   year, sex), call. = FALSE)
    rate <- count / exposure
    if (!is.finite(rate))
      stop(sprintf("HMD pooled-tail rate is not finite: year %g, sex %s.", year, sex), call. = FALSE)
    tibble::tibble(time_start = year, time_end = period_end, sex = sex,
      age_start = open_age, age_end = Inf, deaths = count, person_years = exposure, rate = rate)
  })
  tibble::as_tibble(do.call(rbind, totals))
}
