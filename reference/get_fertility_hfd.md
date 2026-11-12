# Download prepared fertility rates from the HFD

Download fertility rate schedules from the Human Fertility Database in
the format that
[`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md)
expects. Requires user name and password and the HMDHFDplus package.

## Usage

``` r
get_fertility_hfd(
  country,
  year,
  username = Sys.getenv("HFD_USERNAME"),
  password = Sys.getenv("HFD_PASSWORD"),
  update = NULL,
  by_parity = FALSE
)
```

## Arguments

- country:

  A single HFD population code, such as `"USA"`.

- year:

  One whole-number calendar year or a nonempty vector of distinct years.
  All requested years must exist.

- username:

  HFD username, defaulting to `HFD_USERNAME` environmental variable.

- password:

  HFD password, defaulting to `HFD_PASSWORD` environmental variable.

- update:

  Optional HFD revision date as `"YYYYMMDD"`, passed as `Update` to
  HMDHFDplus. `NULL` requests its current revision.

- by_parity:

  Logical, default FALSE. TRUE returns conditional rates with
  `parity_start`/`parity_end` in addition to the age bounds and `rate`.
  Requires the HFD conditional-rate product.

## Value

A fertility rate tibble with a `source` attribute recording selected
years and source periods. Multiple years include
`time_start`/`time_end`; length-one requests omit them.

## Details

To use this function, users will need to install the HMDHFDplus package
with:

    install.packages("HMDHFDplus")

Users will also need to either pass a user name and password to the
function arguments or set the these values in the `.Renviron` as
follows:

    HFD_USERNAME="username here"
    HFD_PASSWORD="password here"

A single requested year returns a static schedule without calendar
columns. Multiple years return `time_start`/`time_end` for the selected
source periods, in chronological order.

HFD endpoint rows retained: `12-` is represented as `[12,13)` using its
supplied rate, and `55+` as `[55,Inf)`.

## Examples

``` r
if (FALSE) { # \dontrun{
fert_rates_usa <- get_fertility_hfd("USA", 2000)
fert_rates_swe <- get_fertility_hfd("SWE", 1891:2000)
# Conditional age/parity rates
fert_rates_usa_parity <- get_fertility_hfd("USA", 2000, by_parity = TRUE)
} # }
```
