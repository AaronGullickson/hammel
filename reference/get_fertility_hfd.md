# Download prepared fertility rates from the Human Fertility Database

Download total period age-specific fertility rates (`asfrRR`) through
HMDHFDplus and select one year or a vector of years from one download.
Rates remain births per woman-year. One year gives a static schedule
without calendar columns; multiple years give chronological
`time_start`/`time_end` bands \[year, year + 1). Nonconsecutive years
retain gaps. Birth-order-specific columns are not used by default.
Missing total rates are rejected. These HFD products use their native
annual, single-age layout; HMD-style `shape` options are not available.
With `by_parity = TRUE`, request conditional age/parity rates (`mi`)
instead: first births map to parity \[0,1), through fifth-or-higher
births at \[4,Inf). Missing conditional cells are replaced with zero;
source metadata records that policy and the number of replacements
across the selected years. Parity tables can be supplied to fertility
constructors for parity-dependent execution. Parity initializes from
recorded maternal birth histories.

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

  A single HFD population code.

- year:

  One whole-number calendar year or a nonempty vector of distinct years.
  All requested years must exist; duplicates are rejected.

- username:

  HFD username, defaulting to `HFD_USERNAME`.

- password:

  HFD password, defaulting to `HFD_PASSWORD`.

- update:

  Optional HFD revision date as `"YYYYMMDD"`, passed as `Update` to
  HMDHFDplus. `NULL` requests its current revision.

- by_parity:

  Logical, default FALSE. TRUE returns conditional rates with
  `parity_start`/`parity_end` in addition to the age bounds and `rate`.
  Requires the HFD conditional-rate product; unavailable products/years
  error rather than falling back to unconditional rates.

## Value

A fertility rate tibble with a `source` attribute recording selected
years and source periods. Multiple years include
`time_start`/`time_end`; length-one requests omit them. Total-rate
tables are validated for execution. Parity tables use half-open
age/parity bounds; unavailable conditional rates are zero-filled.

## Details

HFD endpoint rows are retained: `12-` is represented as `[12,13)` using
its supplied rate, and `55+` as `[55,Inf)`. TFR omits the upper-open
band and warns if its rate is nonzero. No rates are supplied outside the
represented schedule. Credentials are never retained in the table.

## Examples

``` r
if (FALSE) { # \dontrun{
rates <- get_fertility_hfd("USA", 2000)
saveRDS(rates, "rates.rds")
model <- hammel_model() |> set_fertility(readRDS("rates.rds"))
history <- get_fertility_hfd("SWE", 1891:2000)
# Conditional age/parity rates for simulation:
parity_rates <- get_fertility_hfd("USA", 2000, by_parity = TRUE)
} # }
```
