# Download prepared mortality rates from the HMD

Download mortality rate schedules from the Human Mortality Database in
the format that
[`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
expects. Requires user name and password and the HMDHFDplus package.

## Usage

``` r
get_mortality_hmd(
  country,
  year,
  username = Sys.getenv("HMD_USERNAME"),
  password = Sys.getenv("HMD_PASSWORD"),
  shape = "1x1",
  open_age = NULL
)
```

## Arguments

- country:

  A single HMD population code, such as `"USA"`.

- year:

  One whole-number calendar year or a nonempty vector of distinct years.
  For grouped HMD periods, supply their starting years. All requested
  years must exist in the product.

- username:

  Database username. Defaults to the `HMD_USERNAME` environment
  variable.

- password:

  Database password. Defaults to the `HMD_PASSWORD` environment
  variable.

- shape:

  Native HMD age-by-period shape: `"1x1"` (default), `"1x5"`, `"1x10"`,
  `"5x1"`, `"5x5"`, or `"5x10"`.

- open_age:

  Optional nonnegative whole age at a native age-band start. Default
  NULL keeps the source age layout. Otherwise retain rates below this
  age and pool `Deaths_<shape>` and `Exposures_<shape>` into one
  `[open_age, Inf)` rate for each sex and source period.

## Value

A validated mortality rate tibble with a `source` attribute recording
years, shape and represented source periods. Multiple-year requests
include `time_start`/`time_end`. With `open_age`, source metadata also
records the pooling method, products and tail totals.

## Details

To use this function, users will need to install the HMDHFDplus package
with:

    install.packages("HMDHFDplus")

Users will also need to either pass a user name and password to the
function arguments or set the these values in the `.Renviron` as
follows:

    HMD_USERNAME="username here"
    HMD_PASSWORD="password here"

A single requested year returns a static schedule without calendar
columns. Multiple years return `time_start`/`time_end` for the selected
source periods, in chronological order.

HMD rate schedules include an open-ended category at age 110+, but with
many historical sources this will lead to missing values at some earlier
ages. Users can avoid this by using the optional `open_age` argument to
specify an open-ended age. The function will pool source death and
exposure at this open age band from HMD to construct the open-age death
rate.

## Examples

``` r
if (FALSE) { # \dontrun{
mort_rates_usa <- get_mortality_hmd("USA", year = 2000)
mort_rates_swe <- get_mortality_hmd("SWE", year = 1891:2000, open_age = 95)
} # }
```
