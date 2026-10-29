#' Download prepared mortality rates from the Human Mortality Database
#'
#' Download the HMD `Mx_1x1` product through HMDHFDplus and select one calendar
#' year. Male and female central death rates are retained in person-year units
#' and used under Hammel's piecewise constant hazard interpretation. The selected
#' year's schedule is constant over simulation time, regardless of the model's
#' calendar origin. No probabilities or monthly rates are calculated.
#'
#' Single-year ages become `[age, age + 1)`; an upper open age such as `110+`
#' becomes `[110, Inf)`. Missing rates are rejected, not interpolated or dropped.
#' Credentials are passed to the reader and are not stored in the table.
#'
#' @param country A single HMD population code, such as `"USA"`.
#' @param year A single whole-number calendar year to select.
#' @param username Database username. Defaults to the `HMD_USERNAME`
#'   environment variable.
#' @param password Database password. Defaults to the `HMD_PASSWORD`
#'   environment variable.
#' @return A validated mortality rate tibble with a `source` attribute.
#' @export
#' @examples
#' \dontrun{
#' rates <- get_mortality_hmd("USA", year = 2000)
#' saveRDS(rates, "rates.rds")
#' model <- hammel_model() |> set_mortality(readRDS("rates.rds"))
#' }
get_mortality_hmd <- function(country, year, username = Sys.getenv("HMD_USERNAME"),
                          password = Sys.getenv("HMD_PASSWORD")) {
  check_database_request(country, year, username, password)
  data <- HMDHFDplus::readHMDweb(CNTRY = country, item = "Mx_1x1",
                               username = username, password = password,
                               fixup = FALSE)
  rates <- mortality_from_hmd(data, year)$rates
  attr(rates, "source") <- database_source("HMD", country, year, "Mx_1x1")
  rates
}

select_database_year <- function(data, year, columns) {
  check_columns(data, c("Year", "Age", columns), "database data")
  check_numeric_column(data$Year, "database data$Year")
  selected <- data[data$Year == year, , drop = FALSE]
  if (nrow(selected) == 0L) {
    stop("The requested year is not present in the downloaded data.", call. = FALSE)
  }
  selected
}

mortality_from_hmd <- function(data, year) {
  data <- select_database_year(data, year, c("Male", "Female"))
  ages <- trimws(as.character(data$Age))
  if (anyNA(ages) || any(!grepl("^[0-9]+[+]?$", ages))) {
    stop("Expected single-year HMD ages with an optional upper '+' interval.",
         call. = FALSE)
  }
  starts <- as.numeric(sub("[+]$", "", ages))
  open <- grepl("[+]$", ages)
  # HMDHFDplus 2.0.8 readHMDweb always parses ages, even with fixup = FALSE.
  # Preserve its explicit upper-open indicator as well as raw '+' labels.
  if ("OpenInterval" %in% names(data)) {
    if (!is.logical(data$OpenInterval) || anyNA(data$OpenInterval)) {
      stop("HMD OpenInterval must be a nonmissing logical indicator.", call. = FALSE)
    }
    open <- open | data$OpenInterval
  }
  ends <- ifelse(open, Inf, starts + 1)
  mortality(rbind(
    data.frame(sex = "m", age_start = starts, age_end = ends, rate = data$Male),
    data.frame(sex = "f", age_start = starts, age_end = ends, rate = data$Female)
  ))
}

database_source <- function(database, country, year, item, update = NULL) {
  list(database = database, country = country, year = year, item = item,
       update = update, retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
       reader_version = as.character(utils::packageVersion("HMDHFDplus")))
}

#' Download prepared fertility rates from the Human Fertility Database
#'
#' Download total period age-specific fertility rates (`asfrRR`) through
#' HMDHFDplus and select one explicit year. Rates remain births per woman-year
#' and the schedule is constant over simulation time. Birth-order-specific
#' columns are not used by default. Missing total rates are rejected.
#' With `by_parity = TRUE`, request conditional age/parity rates (`mi`) instead:
#' first births map to parity [0,1), through fifth-or-higher births at [4,Inf).
#' Missing conditional cells are replaced with zero; source metadata records
#' that policy and the number of replacements in the selected year.
#' Parity tables prepare data for milestone 5; parity-dependent execution is
#' not yet implemented and fertility constructors reject parity columns.
#'
#' HFD endpoint rows are retained: `12-` is represented as `[12,13)` using its
#' supplied rate, and `55+` as `[55,Inf)`. TFR omits the upper-open band and warns
#' if its rate is nonzero. No rates are supplied outside the represented schedule.
#' Credentials are never retained in the table.
#'
#' @param country A single HFD population code.
#' @param year A single whole-number calendar year.
#' @param username HFD username, defaulting to `HFD_USERNAME`.
#' @param password HFD password, defaulting to `HFD_PASSWORD`.
#' @param update Optional HFD revision date as `"YYYYMMDD"`, passed as `Update`
#'   to HMDHFDplus. `NULL` requests its current revision.
#' @param by_parity Logical, default FALSE. TRUE returns conditional rates with
#'   `parity_start`/`parity_end` in addition to the age bounds and `rate`.
#'   Requires the HFD conditional-rate product; unavailable products/years error
#'   rather than falling back to unconditional rates.
#' @return A fertility rate tibble with a `source` attribute. Total-rate tables
#'   are validated for execution. Parity tables currently support preparation/
#'   inspection only; unavailable conditional rates are zero-filled.
#' @export
#' @examples
#' \dontrun{
#' rates <- get_fertility_hfd("USA", 2000)
#' saveRDS(rates, "rates.rds")
#' model <- hammel_model() |> set_fertility(readRDS("rates.rds"))
#' # Data preparation for future parity-dependent execution:
#' parity_rates <- get_fertility_hfd("USA", 2000, by_parity = TRUE)
#' }
get_fertility_hfd <- function(country, year, username = Sys.getenv("HFD_USERNAME"),
                          password = Sys.getenv("HFD_PASSWORD"), update = NULL,
                          by_parity = FALSE) {
  if (!is.logical(by_parity) || length(by_parity) != 1L ||
      !is.null(dim(by_parity)) || is.na(by_parity))
    stop("`by_parity` must be TRUE or FALSE.", call. = FALSE)
  check_database_request(country, year, username, password)
  if (!is.null(update) &&
      (!is.character(update) || length(update) != 1L || is.na(update) ||
       !grepl("^[0-9]{8}$", update))) {
    stop("`update` must be NULL or a YYYYMMDD revision date.", call. = FALSE)
  }
  item <- if (by_parity) "mi" else "asfrRR"
  data <- HMDHFDplus::readHFDweb(CNTRY = country, item = item,
                               username = username, password = password,
                               fixup = FALSE, Update = update)
  rates <- if (by_parity) fertility_by_parity_from_hfd(data, year) else
    fertility_from_hfd(data, year)$rates
  attr(rates, "source") <- database_source("HFD", country, year, item, update)
  if (by_parity) {
    attr(rates, "source")$rate_type <- "conditional_age_parity"
    attr(rates, "source")$parity_mapping <- "birth orders 1,2,3,4,5+ -> current parity 0,1,2,3,4+"
    attr(rates, "source")$missing_rates <- "zero"
    columns <- c(paste0("m", 1:4, "x"), "m5px")
    attr(rates, "source")$missing_rates_replaced <- sum(is.na(data[data$Year == year, columns, drop = FALSE]))
  }
  attr(rates, "source")$open_ages <- "lower_as_single_year_upper_open"
  rates
}

hfd_age_intervals <- function(data) {
  ages <- trimws(as.character(data$Age))
  if (anyNA(ages) || any(!grepl("^[0-9]+[-+]?$", ages)))
    stop("Expected single-year HFD ages with optional endpoint signs.", call. = FALSE)
  starts <- as.numeric(sub("[-+]$", "", ages))
  if (anyDuplicated(starts)) stop("HFD rates have duplicate age bands.", call. = FALSE)
  tibble::tibble(age_start = starts,
    age_end = ifelse(grepl("[+]$", ages), Inf, starts + 1))
}

fertility_from_hfd <- function(data, year) {
  data <- select_database_year(data, year, "ASFR")
  ages <- hfd_age_intervals(data)
  fertility(data.frame(age_start = ages$age_start, age_end = ages$age_end,
                       rate = data$ASFR))
}

fertility_by_parity_from_hfd <- function(data, year) {
  columns <- c(paste0("m", 1:4, "x"), "m5px")
  data <- select_database_year(data, year, columns)
  ages <- hfd_age_intervals(data)
  starts <- ages$age_start
  for (column in columns) {
    values <- data[[column]]
    if (is.logical(values) && all(is.na(values))) values <- as.numeric(values)
    check_numeric_column(values, paste0("HFD ", column), missing = TRUE)
    if (any(data[[column]] < 0, na.rm = TRUE))
      stop("HFD conditional fertility rates must be nonnegative.", call. = FALSE)
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

check_database_request <- function(country, year, username, password) {
  if (!is.character(country) || length(country) != 1L || is.na(country) ||
      !nzchar(trimws(country))) {
    stop("`country` must be a single nonempty population code.", call. = FALSE)
  }
  if (!is.numeric(year) || length(year) != 1L || !is.finite(year) ||
      year != trunc(year)) {
    stop("`year` must be a single whole-number calendar year.", call. = FALSE)
  }
  for (credential in list(username, password)) {
    if (!is.character(credential) || length(credential) != 1L ||
        is.na(credential) || !nzchar(credential)) {
      stop("Supply database credentials as arguments or environment variables.",
           call. = FALSE)
    }
  }
  if (!requireNamespace("HMDHFDplus", quietly = TRUE)) {
    stop("Install HMDHFDplus to download database rates.", call. = FALSE)
  }
}
