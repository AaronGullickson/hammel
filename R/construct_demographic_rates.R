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
#' columns are not used. Missing rates are rejected.
#'
#' HFD open-age rows such as `12-` and `55+` are omitted. The resulting
#' schedule covers only closed single-year age bands, so its TFR excludes those
#' endpoint contributions. No rates are supplied outside the retained schedule.
#' Credentials are never retained in the table.
#'
#' @param country A single HFD population code.
#' @param year A single whole-number calendar year.
#' @param username HFD username, defaulting to `HFD_USERNAME`.
#' @param password HFD password, defaulting to `HFD_PASSWORD`.
#' @param update Optional HFD revision date as `"YYYYMMDD"`, passed as `Update`
#'   to HMDHFDplus. `NULL` requests its current revision.
#' @return A validated fertility rate tibble with a `source` attribute.
#' @export
#' @examples
#' \dontrun{
#' rates <- get_fertility_hfd("USA", 2000)
#' saveRDS(rates, "rates.rds")
#' model <- hammel_model() |> set_fertility(readRDS("rates.rds"))
#' }
get_fertility_hfd <- function(country, year, username = Sys.getenv("HFD_USERNAME"),
                          password = Sys.getenv("HFD_PASSWORD"), update = NULL) {
  check_database_request(country, year, username, password)
  if (!is.null(update) &&
      (!is.character(update) || length(update) != 1L || is.na(update) ||
       !grepl("^[0-9]{8}$", update))) {
    stop("`update` must be NULL or a YYYYMMDD revision date.", call. = FALSE)
  }
  data <- HMDHFDplus::readHFDweb(CNTRY = country, item = "asfrRR",
                               username = username, password = password,
                               fixup = FALSE, Update = update)
  rates <- fertility_from_hfd(data, year)$rates
  attr(rates, "source") <- database_source("HFD", country, year, "asfrRR", update)
  attr(rates, "source")$open_ages <- "omit"
  rates
}

fertility_from_hfd <- function(data, year) {
  data <- select_database_year(data, year, "ASFR")
  ages <- trimws(as.character(data$Age))
  if (anyNA(ages) || any(!grepl("^[0-9]+[-+]?$", ages))) {
    stop("Expected single-year HFD ages with optional endpoint signs.",
         call. = FALSE)
  }
  open <- grepl("[-+]$", ages)
  data <- data[!open, , drop = FALSE]
  ages <- ages[!open]
  starts <- as.numeric(sub("[-+]$", "", ages))
  fertility(data.frame(age_start = starts, age_end = starts + 1,
                       rate = data$ASFR))
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
