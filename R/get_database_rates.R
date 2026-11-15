# Retrieve HMD and HFD rates, convert their tables, and preserve source
# provenance.

# PUBLIC FUNCTIONS ---------------------------------------------------------

# HMD MORTALITY IMPORT -----------------------------------

#' Download prepared mortality rates from the HMD
#'
#' @description
#'
#' Download mortality rate schedules from the Human Mortality Database in the
#' format that [mortality()] expects. Requires user name and password and the
#' HMDHFDplus package.
#'
#' @details
#'
#' To use this function, users will need to install the HMDHFDplus package with:
#'
#' ```r
#' install.packages("HMDHFDplus")
#' ```
#'
#' Users will also need to either pass a user name and password to the function
#' arguments or set the these values in the `.Renviron` as follows:
#'
#' ```bash
#' HMD_USERNAME="username here"
#' HMD_PASSWORD="password here"
#' ```
#'
#' A single requested year returns a static schedule without calendar columns.
#' Multiple years return `time_start`/`time_end` for the selected source
#' periods, in chronological order.
#'
#' HMD rate schedules include an open-ended category at age 110+, but with many
#' historical sources this will lead to missing values at some earlier ages.
#' Users can avoid this by using the optional `open_age` argument to specify an
#' open-ended age. The function will pool source death and exposure at this open
#' age band from HMD to construct the open-age death rate.
#'
#' @param country A single HMD population code, such as `"USA"`.
#' @param year One whole-number calendar year or a nonempty vector of distinct
#'   years. For grouped HMD periods, supply their starting years. All requested
#'   years must exist in the product.
#' @param username Database username. Defaults to the `HMD_USERNAME` environment
#'   variable.
#' @param password Database password. Defaults to the `HMD_PASSWORD` environment
#'   variable.
#' @param shape Native HMD age-by-period shape: `"1x1"` (default), `"1x5"`,
#'   `"1x10"`, `"5x1"`, `"5x5"`, or `"5x10"`.
#' @param open_age Optional nonnegative whole age at a native age-band start.
#'   Default NULL keeps the source age layout. Otherwise retain rates below this
#'   age and pool `Deaths_<shape>` and `Exposures_<shape>` into one `[open_age,
#'   Inf)` rate for each sex and source period.
#' @return A validated mortality rate tibble with a `source` attribute recording
#'   years, shape and represented source periods. Multiple-year requests include
#'   `time_start`/`time_end`. With `open_age`, source metadata also records the
#'   pooling method, products and tail totals.
#' @export
#' @examples
#' \dontrun{
#' mort_rates_usa <- get_mortality_hmd("USA", year = 2000)
#' mort_rates_swe <- get_mortality_hmd("SWE", year = 1891:2000, open_age = 95)
#' }
get_mortality_hmd <- function(
  country,
  year,
  username = Sys.getenv("HMD_USERNAME"),
  password = Sys.getenv("HMD_PASSWORD"),
  shape = "1x1",
  open_age = NULL
) {
  check_hmd_shape(shape)
  check_hmd_open_age(open_age, shape)
  check_database_request(country, year, username, password)
  item <- paste0("Mx_", shape)
  data <- read_database(
    "HMD", country, item,
    username = username,
    password = password,
    fixup = FALSE
  )
  coverage_years <- NULL
  # The installed reader discards raw upper period bounds. Recover them using
  # native starts and annual coverage instead of inventing years at the tail.
  grouped_time <- !endsWith(shape, "x1")
  raw_ranges <- grepl("^[0-9]+-[0-9]+$", trimws(as.character(data$Year)))
  if (grouped_time && any(!raw_ranges)) {
    annual <- read_database(
      "HMD", country, "Mx_1x1",
      username = username,
      password = password,
      fixup = FALSE
    )
    check_database_columns(annual, "Year", "HMD annual coverage")
    coverage_years <- sort(unique(annual$Year))
    check_database_years(coverage_years)
  }
  deaths <- exposures <- NULL
  if (!is.null(open_age)) {
    deaths <- read_database(
      "HMD", country, paste0("Deaths_", shape),
      username = username,
      password = password,
      fixup = FALSE
    )
    exposures <- read_database(
      "HMD", country, paste0("Exposures_", shape),
      username = username,
      password = password,
      fixup = FALSE
    )
  }
  process <- mortality_from_hmd(
    data, year, shape, coverage_years, open_age, deaths, exposures
  )
  rates <- process$rates
  periods <- unique(hmd_period_intervals(data, shape, coverage_years))
  periods <- periods[periods$time_start %in% year, ]
  periods <- periods[order(periods$time_start), ]
  attr(rates, "source") <- database_source("HMD", country, sort(year), item)
  attr(rates, "source")$shape <- shape
  attr(rates, "source")$periods <- periods
  if (!is.null(coverage_years)) {
    attr(rates, "source")$coverage_item <- "Mx_1x1"
  }
  if (!is.null(open_age)) {
    attr(rates, "source")$open_age <- open_age
    attr(rates, "source")$open_age_method <- "pooled_deaths_exposure"
    attr(rates, "source")$pooling_items <- c(
      deaths = paste0("Deaths_", shape),
      exposure = paste0("Exposures_", shape)
    )
    attr(rates, "source")$pooled_tail <- process$pooled_tail
  }
  rates
}

# HFD FERTILITY IMPORT -----------------------------------

#' Download prepared fertility rates from the HFD
#'
#' @description
#'
#' Download fertility rate schedules from the Human Fertility Database in the
#' format that [fertility()] expects. Requires user name and password and the
#' HMDHFDplus package.
#'
#' @details
#'
#' To use this function, users will need to install the HMDHFDplus package with:
#'
#' ```r
#' install.packages("HMDHFDplus")
#' ```
#'
#' Users will also need to either pass a user name and password to the function
#' arguments or set the these values in the `.Renviron` as follows:
#'
#' ```bash
#' HFD_USERNAME="username here"
#' HFD_PASSWORD="password here"
#' ```
#'
#' A single requested year returns a static schedule without calendar columns.
#' Multiple years return `time_start`/`time_end` for the selected source
#' periods, in chronological order.
#'
#' HFD endpoint rows retained: `12-` is represented as `[12,13)` using its
#' supplied rate, and `55+` as `[55,Inf)`.
#'
#' @param country A single HFD population code, such as `"USA"`.
#' @param year One whole-number calendar year or a nonempty vector of distinct
#'   years. All requested years must exist.
#' @param username HFD username, defaulting to `HFD_USERNAME` environmental
#'   variable.
#' @param password HFD password, defaulting to `HFD_PASSWORD` environmental
#'   variable.
#' @param update Optional HFD revision date as `"YYYYMMDD"`, passed as `Update`
#'   to HMDHFDplus. `NULL` requests its current revision.
#' @param by_parity Logical, default FALSE. TRUE returns conditional rates with
#'   `parity_start`/`parity_end` in addition to the age bounds and `rate`.
#'   Requires the HFD conditional-rate product.
#' @return A fertility rate tibble with a `source` attribute recording selected
#'   years and source periods. Multiple years include `time_start`/`time_end`;
#'   length-one requests omit them.
#' @export
#' @examples
#' \dontrun{
#' fert_rates_usa <- get_fertility_hfd("USA", 2000)
#' fert_rates_swe <- get_fertility_hfd("SWE", 1891:2000)
#' # Conditional age/parity rates
#' fert_rates_usa_parity <- get_fertility_hfd("USA", 2000, by_parity = TRUE)
#' }
get_fertility_hfd <- function(
  country,
  year,
  username = Sys.getenv("HFD_USERNAME"),
  password = Sys.getenv("HFD_PASSWORD"),
  update = NULL,
  by_parity = FALSE
) {
  check_flag(by_parity, "by_parity")
  check_database_request(country, year, username, password)
  if (!is.null(update) &&
        (!is.character(update) || length(update) != 1L || is.na(update) ||
           !grepl("^[0-9]{8}$", update))) {
    abort_hammel(
      c(
        paste(
          "{.arg update} must be {.code NULL} or a revision date as",
          "{.val YYYYMMDD}."
        ),
        "x" = "It is {describe_value(update)}."
      ),
      "argument",
      argument = "update"
    )
  }
  item <- if (by_parity) "mi" else "asfrRR"
  data <- read_database(
    "HFD", country, item,
    username = username,
    password = password,
    fixup = FALSE,
    Update = update
  )
  rates <- if (by_parity) {
    fertility_by_parity_from_hfd(data, year)
  } else {
    fertility_from_hfd(data, year)$rates
  }
  attr(rates, "source") <- database_source(
    "HFD", country, sort(year), item, update
  )
  attr(rates, "source")$periods <- tibble::tibble(
    time_start = sort(year),
    time_end = sort(year) + 1
  )
  if (by_parity) {
    attr(rates, "source")$rate_type <- "conditional_age_parity"
    attr(rates, "source")$parity_mapping <-
      "birth orders 1,2,3,4,5+ -> current parity 0,1,2,3,4+"
    attr(rates, "source")$missing_rates <- "zero"
    columns <- c(paste0("m", 1:4, "x"), "m5px")
    attr(rates, "source")$missing_rates_replaced <- sum(
      is.na(data[data$Year %in% year, columns, drop = FALSE])
    )
  }
  attr(rates, "source")$open_ages <- "lower_as_single_year_upper_open"
  rates
}

# HELPER FUNCTIONS ---------------------------------------------------------

# HMD/HFD retrieval, conversion, and source provenance.

# SHARED REQUEST VALIDATION AND PROVENANCE ---------------------------------

select_database_year <- function(
  data,
  year,
  columns,
  call = rlang::caller_env()
) {
  check_database_years(year, call = call)
  check_database_columns(
    data, c("Year", "Age", columns), "database",
    call = call
  )
  check_database_numeric(data$Year, "Year", call = call)
  missing <- setdiff(year, data$Year)
  if (length(missing)) {
    abort_database(
      sprintf(
        "The requested year(s) are not present in the downloaded data: %s.",
        paste(missing, collapse = ", ")
      ),
      call = call
    )
  }
  selected <- data[data$Year %in% year, , drop = FALSE]
  selected
}

database_source <- function(database, country, year, item, update = NULL) {
  list(
    database = database,
    country = country,
    year = year,
    item = item,
    update = update,
    retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    reader_version = as.character(utils::packageVersion("HMDHFDplus"))
  )
}

check_database_request <- function(
  country,
  year,
  username,
  password,
  call = rlang::caller_env()
) {
  if (!is.character(country) || length(country) != 1L || is.na(country) ||
        !nzchar(trimws(country))) {
    abort_hammel(
      c(
        "{.arg country} must be a single population code, such as {.val USA}.",
        "x" = "It is {describe_value(country)}."
      ),
      "argument",
      argument = "country",
      call = call
    )
  }
  check_database_years(year, call = call)
  credentials <- list(username = username, password = password)
  for (name in names(credentials)) {
    credential <- credentials[[name]]
    if (!is.character(credential) || length(credential) != 1L ||
          is.na(credential) ||
          !nzchar(credential)) {
      abort_hammel(
        c(
          "A database {.arg {name}} is required.",
          "i" = paste(
            "Pass {.arg username} and {.arg password}, or set the",
            "{.envvar HMD_USERNAME}/{.envvar HMD_PASSWORD} or",
            "{.envvar HFD_USERNAME}/{.envvar HFD_PASSWORD} environment",
            "variables."
          )
        ),
        "argument",
        argument = name,
        call = call
      )
    }
  }
  if (!requireNamespace("HMDHFDplus", quietly = TRUE)) {
    abort_hammel(
      c(
        "Downloading database rates needs the {.pkg HMDHFDplus} package.",
        "i" = "Install it with {.code install.packages(\"HMDHFDplus\")}."
      ),
      "database",
      call = call
    )
  }
}

# Download one database item, reporting a failure (credentials, network or an
# unavailable item) as a database error with the reader's error attached.
read_database <- function(
  database,
  country,
  item,
  ...,
  call = rlang::caller_env()
) {
  tryCatch(
    if (database == "HMD") {
      HMDHFDplus::readHMDweb(CNTRY = country, item = item, ...)
    } else {
      HMDHFDplus::readHFDweb(CNTRY = country, item = item, ...)
    },
    error = function(e) {
      abort_hammel(
        c(
          paste(
            "Couldn't download {.val {item}} for {.val {country}} from",
            "{database}."
          ),
          "i" = "Check the country code, item, credentials and connection."
        ),
        "database",
        database = database,
        parent = e,
        call = call
      )
    }
  )
}

# Downloaded tables must have the expected columns.
check_database_columns <- function(
  data,
  columns,
  source,
  call = rlang::caller_env()
) {
  absent <- if (is.data.frame(data)) setdiff(columns, names(data)) else columns
  if (length(absent)) {
    abort_hammel(
      c(
        "The downloaded {source} data doesn't have the expected columns.",
        "x" = "Missing: {.field {absent}}."
      ),
      "database",
      columns = absent,
      call = call
    )
  }
  invisible(data)
}

# Downloaded values must be finite numbers (or NA when `missing`).
check_database_numeric <- function(
  x,
  label,
  missing = FALSE,
  call = rlang::caller_env()
) {
  observed <- x[!is.na(x)]
  if (!is.numeric(x) || !is.null(dim(x)) || any(is.nan(x)) ||
        (!missing && anyNA(x)) ||
        any(!is.finite(observed))) {
    abort_hammel(
      c(
        "The downloaded {label} values must be finite numbers.",
        "x" = "They are {describe_value(x)}."
      ),
      "database",
      call = call
    )
  }
  invisible(x)
}

# Raise a database error whose message is preformatted text; braces in it
# are printed literally.
abort_database <- function(text, ..., call = rlang::caller_env()) {
  text <- gsub("}", "}}", gsub("{", "{{", text, fixed = TRUE), fixed = TRUE)
  abort_hammel(text, "database", ..., call = call)
}

# Source intervals for native HMD products; no interpolation.
check_database_years <- function(year, call = rlang::caller_env()) {
  if (!is.numeric(year) || !is.null(dim(year)) || !length(year) ||
        anyNA(year) ||
        any(!is.finite(year)) ||
        any(year != trunc(year)) ||
        anyDuplicated(year)) {
    abort_hammel(
      c(
        "{.arg year} must contain distinct whole-number calendar years.",
        "x" = "It is {describe_value(year)}."
      ),
      "argument",
      argument = "year",
      call = call
    )
  }
  invisible(year)
}

# Retain source row ordering within each year, return years chronologically.
combine_database_years <- function(year, convert, call = rlang::caller_env()) {
  check_database_years(year, call = call)
  years <- sort(year)
  rates <- lapply(
    years, function(y) {
      rows <- convert(y)
      rows$time_start <- y
      rows$time_end <- y + 1
      rows
    }
  )
  tibble::as_tibble(do.call(rbind, rates))
}

# HMD CONVERSION HELPERS ---------------------------------------------------

mortality_from_hmd <- function(
  data,
  year,
  shape = "1x1",
  coverage_years = NULL,
  open_age = NULL,
  deaths = NULL,
  exposures = NULL,
  call = rlang::caller_env()
) {
  check_hmd_shape(shape, call = call)
  check_hmd_open_age(open_age, shape, call = call)
  check_database_years(year, call = call)
  check_database_columns(
    data, c("Year", "Age", "Male", "Female"), "HMD",
    call = call
  )
  periods <- hmd_period_intervals(data, shape, coverage_years, call = call)
  missing <- setdiff(year, periods$time_start)
  if (length(missing)) {
    abort_database(
      sprintf(
        "The requested year(s) are not present as source period starts: %s.",
        paste(missing, collapse = ", ")
      ),
      call = call
    )
  }
  if (!is.null(open_age)) {
    deaths <- index_hmd_pool_source(
      deaths, year, shape, coverage_years, "HMD deaths",
      call = call
    )
    exposures <- index_hmd_pool_source(
      exposures, year, shape, coverage_years, "HMD exposure",
      call = call
    )
  }
  converted <- lapply(
    sort(year), function(y) {
      rows <- periods$time_start == y
      selected <- data[rows, , drop = FALSE]
      ages <- hmd_age_intervals(selected, shape, call = call)
      period_end <- unique(periods$time_end[rows])
      keep <- if (is.null(open_age)) {
        rep(TRUE, nrow(ages))
      } else {
        ages$age_start < open_age
      }
      totals <- if (is.null(open_age)) {
        NULL
      } else {
        pool_hmd_open_age(
          ages, y, period_end, open_age, deaths, exposures,
          call = call
        )
      }
      out <- lapply(
        c(m = "Male", f = "Female"), function(column) {
          sex <- if (column == "Male") "m" else "f"
          check_hmd_retained_rates(
            selected[[column]][keep], ages[keep, ], y, sex,
      call = call
    )
          below <- data.frame(
            sex = rep(sex, sum(keep)),
            age_start = ages$age_start[keep],
            age_end = ages$age_end[keep],
            rate = selected[[column]][keep]
          )
          if (!is.null(totals)) {
            final <- totals[totals$sex == sex, c(
              "sex", "age_start", "age_end", "rate"
            )]
            below <- rbind(below, as.data.frame(final))
          }
          below
        }
      )
      out <- do.call(rbind, out)
      if (length(year) > 1L) {
        out$time_start <- y
        out$time_end <- unique(periods$time_end[rows])
      }
      list(rates = out, pooled_tail = totals)
    }
  )
  process <- mortality(do.call(rbind, lapply(converted, `[[`, "rates")))
  if (!is.null(open_age)) {
    process$pooled_tail <- tibble::as_tibble(
      do.call(rbind, lapply(converted, `[[`, "pooled_tail"))
    )
  }
  process
}

check_hmd_shape <- function(shape, call = rlang::caller_env()) {
  shapes <- c("1x1", "1x5", "1x10", "5x1", "5x5", "5x10")
  if (!is.character(shape) || length(shape) != 1L || is.na(shape) ||
        !is.null(dim(shape)) ||
        !shape %in% shapes) {
    abort_hammel(
      c(
        "{.arg shape} must be one of {.or {.val {shapes}}}.",
        "x" = "It is {describe_value(shape)}.",
        "i" = "The first number is the age-band width, the second the period."
      ),
      "argument",
      argument = "shape",
      call = call
    )
  }
  invisible(shape)
}

hmd_period_intervals <- function(
  data,
  shape = "1x1",
  coverage_years = NULL,
  call = rlang::caller_env()
) {
  check_hmd_shape(shape, call = call)
  check_database_columns(data, "Year", "HMD", call = call)
  labels <- trimws(as.character(data$Year))
  if (!length(labels) || anyNA(labels) ||
        any(!grepl("^[0-9]+(-[0-9]+)?$", labels))) {
    abort_database(
      "Expected HMD calendar years or inclusive source year ranges.",
      call = call
    )
  }
  width <- as.numeric(strsplit(shape, "x", fixed = TRUE)[[1L]][2L])
  starts <- as.numeric(sub("-.*$", "", labels))
  ranged <- grepl("-", labels, fixed = TRUE)
  ends <- starts + width
  ends[ranged] <- as.numeric(sub("^.*-", "", labels[ranged])) + 1
  if (width > 1 && any(!ranged)) {
    if (is.null(coverage_years)) {
      abort_database(
        paste0(
          "Parsed grouped HMD periods require annual ",
          "coverage to recover source bounds."
        ),
        call = call
      )
    }
    check_database_years(coverage_years, call = call)
    periods <- sort(unique(starts))
    following <- c(tail(periods, -1L), Inf)[match(starts, periods)]
    ends[!ranged] <- pmin(
      ends[!ranged], following[!ranged], max(coverage_years) + 1
    )
  }
  if (any(!is.finite(starts)) || any(!is.finite(ends)) ||
        any(ends <= starts | ends - starts > width)) {
    abort_database(
      "HMD source periods are inconsistent with the requested shape.",
      call = call
    )
  }
  periods <- unique(tibble::tibble(time_start = starts, time_end = ends))
  periods <- periods[order(periods$time_start), ]
  if (anyDuplicated(periods$time_start) ||
        (nrow(periods) > 1L &&
           any(tail(periods$time_start, -1L) < head(periods$time_end, -1L)))) {
    abort_database(
      "HMD source periods must not overlap or have conflicting bounds.",
      call = call
    )
  }
  if (!is.null(coverage_years)) {
    for (i in seq_len(nrow(periods))) {
      years <- seq(periods$time_start[i], periods$time_end[i] - 1)
      if (any(!years %in% coverage_years)) {
        abort_database(
          "HMD source periods are inconsistent with annual coverage.",
          call = call
        )
      }
    }
  }
  tibble::tibble(time_start = starts, time_end = ends)
}

hmd_age_intervals <- function(data, shape = "1x1", call = rlang::caller_env()) {
  check_hmd_shape(shape, call = call)
  ages <- trimws(as.character(data$Age))
  width <- as.numeric(strsplit(shape, "x", fixed = TRUE)[[1L]][1L])
  pattern <- if (width == 1) "^[0-9]+[+]?$" else "^[0-9]+(-[0-9]+|[+])?$"
  if (anyNA(ages) || any(!grepl(pattern, ages))) {
    abort_database(
      "Expected HMD age labels consistent with the requested shape.",
      call = call
    )
  }
  starts <- as.numeric(sub("[-+].*$", "", ages))
  open <- grepl("[+]$", ages)
  # readHMDweb() 2.0.8 always parses ages, even with fixup = FALSE.
  if ("OpenInterval" %in% names(data)) {
    if (!is.logical(data$OpenInterval) || anyNA(data$OpenInterval)) {
      abort_database(
        "HMD OpenInterval must be a nonmissing logical indicator.",
        call = call
      )
    }
    open <- open | data$OpenInterval
  }
  ends <- if (width == 1) {
    starts + 1
  } else {
    ifelse(starts == 0, 1, ifelse(starts == 1, 5, starts + 5))
  }
  ranged <- grepl("-", ages, fixed = TRUE)
  if (any(ranged)) {
    upper <- as.numeric(sub("^.*-", "", ages[ranged])) + 1
    if (any(upper != ends[ranged])) {
      abort_database(
        "HMD source age ranges are inconsistent with the requested shape.",
        call = call
      )
    }
  }
  if (width == 5 && any(!open & starts > 1 & starts %% 5 != 0)) {
    abort_database(
      paste0(
        "HMD five-year ages must use the native ",
        "infant/child/five-year layout."
      ),
      call = call
    )
  }
  ends[open] <- Inf
  tibble::tibble(age_start = starts, age_end = ends)
}

# Explicit source-count pooling for a user-selected open final age category.
check_hmd_open_age <- function(open_age, shape, call = rlang::caller_env()) {
  if (is.null(open_age)) {
    return(invisible(open_age))
  }
  if (!is.numeric(open_age) || length(open_age) != 1L ||
        !is.null(dim(open_age)) ||
        !is.finite(open_age) ||
        open_age < 0 ||
        open_age != trunc(open_age)) {
    abort_hammel(
      c(
        paste(
          "{.arg open_age} must be {.code NULL} or a single nonnegative",
          "whole age."
        ),
        "x" = "It is {describe_value(open_age)}."
      ),
      "argument",
      argument = "open_age",
      call = call
    )
  }
  if (startsWith(shape, "5x") && !open_age %in% c(0, 1) && open_age %% 5 != 0) {
    abort_hammel(
      c(
        "{.arg open_age} must start one of the five-year age bands.",
        "x" = "It is {open_age}.",
        "i" = "With {.code shape = \"{shape}\"}, use 0, 1, 5, 10 and so on."
      ),
      "argument",
      argument = "open_age",
      call = call
    )
  }
  invisible(open_age)
}

# Diagnose invalid retained Mx with source context, rather than hiding its age
# behind the generic rate validator. Missing replaced-tail Mx never enter here.
check_hmd_retained_rates <- function(
  values,
  ages,
  year,
  sex,
  call = rlang::caller_env()
) {
  if (!(is.numeric(values) || (is.logical(values) && all(is.na(values))))) {
    abort_database(
      sprintf(
        "HMD mortality rates must be numeric (year %g, sex %s).", year, sex
      ),
      call = call
    )
  }
  bad <- which(!is.finite(values))
  if (length(bad)) {
    i <- bad[1L]
    abort_database(
      sprintf(
        paste0(
          "HMD mortality rates contain missing or non-finite values: ",
          "year %g, sex %s, age [%g, %g)."
        ), year, sex, ages$age_start[i],
        ages$age_end[i]
      ),
      call = call
    )
  }
  invisible(values)
}

index_hmd_pool_source <- function(
  data,
  year,
  shape,
  coverage_years,
  name,
  call = rlang::caller_env()
) {
  check_database_columns(
    data, c("Year", "Age", "Male", "Female"), name,
    call = call
  )
  periods <- hmd_period_intervals(data, shape, coverage_years, call = call)
  if (any(!year %in% periods$time_start)) {
    abort_database(
      sprintf("%s is missing requested source periods.", name),
      call = call
    )
  }
  rows <- periods$time_start %in% year
  ages <- hmd_age_intervals(data[rows, , drop = FALSE], shape, call = call)
  tibble::tibble(
    time_start = periods$time_start[rows],
    time_end = periods$time_end[rows],
    age_start = ages$age_start,
    age_end = ages$age_end,
    Male = data$Male[rows],
    Female = data$Female[rows]
  )
}

pool_hmd_open_age <- function(
  ages,
  year,
  period_end,
  open_age,
  deaths,
  exposures,
  call = rlang::caller_env()
) {
  bands <- ages[ages$age_start >= open_age, , drop = FALSE]
  bands <- bands[order(bands$age_start), ]
  if (!nrow(bands) || bands$age_start[1L] != open_age ||
        tail(bands$age_end, 1L) != Inf ||
        anyDuplicated(bands$age_start) ||
        (nrow(bands) > 1L &&
           any(head(bands$age_end, -1L) != tail(bands$age_start, -1L)))) {
    abort_database(
      sprintf(
        paste0(
          "HMD rates need a complete source tail from ",
          "open_age %g through Inf in year %g."
        ),
        open_age, year
      ),
      call = call
    )
  }
  expected <- paste(bands$age_start, bands$age_end, sep = ":")
  align <- function(data, name) {
    rows <- data$time_start == year & data$age_start >= open_age
    out <- data[rows, ]
    keys <- paste(out$age_start, out$age_end, sep = ":")
    if (anyDuplicated(keys) || !setequal(keys, expected) ||
          any(out$time_end != period_end)) {
      abort_database(
        sprintf(
          "HMD %s tail age/period bands do not match mortality in year %g.",
          name, year
        ),
        call = call
      )
    }
    out[match(expected, keys), ]
  }
  d <- align(deaths, "deaths")
  e <- align(exposures, "exposure")
  totals <- lapply(
    c(m = "Male", f = "Female"), function(column) {
      sex <- if (column == "Male") "m" else "f"
      for (name in c("deaths", "exposure")) {
        values <- if (name == "deaths") d[[column]] else e[[column]]
        bad <- if (!is.numeric(values)) {
          seq_along(values)
        } else {
          which(!is.finite(values) | values < 0)
        }
        if (length(bad)) {
          i <- bad[1L]
          abort_database(
            sprintf(
              paste0(
                "HMD pooled-tail %s must be finite ",
                "nonnegative numeric values: ",
                "year %g, sex %s, age [%g, %g)."
              ),
              name, year, sex, bands$age_start[i], bands$age_end[i]
            ),
            call = call
          )
        }
      }
      count <- sum(d[[column]])
      exposure <- sum(e[[column]])
      if (!is.finite(count) || !is.finite(exposure) || exposure <= 0) {
        abort_database(
          sprintf(
            paste0(
              "HMD pooled tail needs finite deaths and positive ",
              "finite exposure: year %g, sex %s."
            ),
            year, sex
          ),
          call = call
        )
      }
      rate <- count / exposure
      if (!is.finite(rate)) {
        abort_database(
          sprintf(
            "HMD pooled-tail rate is not finite: year %g, sex %s.", year, sex
          ),
          call = call
        )
      }
      tibble::tibble(
        time_start = year,
        time_end = period_end,
        sex = sex,
        age_start = open_age,
        age_end = Inf,
        deaths = count,
        person_years = exposure,
        rate = rate
      )
    }
  )
  tibble::as_tibble(do.call(rbind, totals))
}

# HFD CONVERSION HELPERS ---------------------------------------------------

hfd_age_intervals <- function(data, call = rlang::caller_env()) {
  ages <- trimws(as.character(data$Age))
  if (anyNA(ages) || any(!grepl("^[0-9]+[-+]?$", ages))) {
    abort_database(
      "Expected single-year HFD ages with optional endpoint signs.",
      call = call
    )
  }
  starts <- as.numeric(sub("[-+]$", "", ages))
  if (anyDuplicated(starts)) {
    abort_database("HFD rates have duplicate age bands.", call = call)
  }
  tibble::tibble(
    age_start = starts,
    age_end = ifelse(grepl("[+]$", ages), Inf, starts + 1)
  )
}

fertility_from_hfd <- function(data, year, call = rlang::caller_env()) {
  check_database_years(year, call = call)
  if (length(year) > 1L) {
    return(
      fertility(
        combine_database_years(
          year, function(y) {
            fertility_from_hfd(data, y, call = call)$rates
          },
      call = call
    )
      )
    )
  }
  data <- select_database_year(data, year, "ASFR", call = call)
  ages <- hfd_age_intervals(data, call = call)
  fertility(
    data.frame(
      age_start = ages$age_start,
      age_end = ages$age_end,
      rate = data$ASFR
    )
  )
}

fertility_by_parity_from_hfd <- function(
  data,
  year,
  call = rlang::caller_env()
) {
  check_database_years(year, call = call)
  if (length(year) > 1L) {
    rates <- combine_database_years(
      year, function(y) fertility_by_parity_from_hfd(data, y, call = call),
      call = call
    )
    validate_rates(rates, mortality = FALSE)
    return(rates)
  }
  columns <- c(paste0("m", 1:4, "x"), "m5px")
  data <- select_database_year(data, year, columns, call = call)
  ages <- hfd_age_intervals(data, call = call)
  starts <- ages$age_start
  for (column in columns) {
    values <- data[[column]]
    if (is.logical(values) && all(is.na(values))) {
      values <- as.numeric(values)
    }
    check_database_numeric(
      values, paste("HFD", column),
      missing = TRUE,
      call = call
    )
    if (any(data[[column]] < 0, na.rm = TRUE)) {
      abort_database(
        "HFD conditional fertility rates must be nonnegative.",
        call = call
      )
    }
  }
  rates <- tibble::tibble(
    age_start = rep(starts, times = 5L),
    age_end = rep(ages$age_end, times = 5L),
    parity_start = rep(0:4, each = length(starts)),
    parity_end = rep(c(1:4, Inf), each = length(starts)),
    rate = as.numeric(unlist(data[columns], use.names = FALSE))
  )
  rates$rate[is.na(rates$rate)] <- 0
  rates
}
