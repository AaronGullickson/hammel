# Download prepared mortality rates from the Human Mortality Database

Download a native HMD mortality product through HMDHFDplus. Male and
female central death rates retain their person-year units and Hammel's
piecewise constant hazard interpretation. Rates are not rescaled by
age/period width. Missing retained rates are rejected. Optional
`open_age` pools source deaths and exposure into a final open age band,
permitting missing Mx in that tail.

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
  years must exist in the product; duplicates are rejected.

- username:

  Database username. Defaults to the `HMD_USERNAME` environment
  variable.

- password:

  Database password. Defaults to the `HMD_PASSWORD` environment
  variable.

- shape:

  Native HMD age-by-period shape: `"1x1"` (default), `"1x5"`, `"1x10"`,
  `"5x1"`, `"5x5"`, or `"5x10"`. This selects `Mx_<shape>`; `open_age`
  can separately pool the final age band.

- open_age:

  Optional nonnegative whole age at a native age-band start. Default
  NULL keeps the source age layout. Otherwise retain rates below this
  age and pool `Deaths_<shape>` and `Exposures_<shape>` into one
  `[open_age, Inf)` rate for each sex and source period: total deaths
  divided by total person-years. This requires two additional product
  downloads. Counts/exposures must be finite and nonnegative with
  positive pooled exposure; individual zero-exposure cells may
  contribute zero. No missing counts/exposures are imputed. Cutoffs must
  align with source bands (0, 1, 5, 10, etc. for five-year ages). The
  cutoff is not a maximum lifespan.

## Value

A validated mortality rate tibble with a `source` attribute recording
years, shape and represented source periods. Multiple-year requests
include `time_start`/`time_end`; length-one requests omit them. With
`open_age`, source metadata also records the pooling method, products
and tail totals.

## Details

A single requested year returns a static schedule without calendar
columns. Multiple years return `time_start`/`time_end` for the selected
source periods, in chronological order. For grouped-time shapes, `year`
selects period starts: for example, 2000 selects the source interval
2000-2004 in a five-year product. Nonconsecutive requests retain
calendar gaps; no rates are carried forward.

Native five-year ages are 0, 1-4, 5-9, etc.; returned bands are \[0,1),
\[1,5), \[5,10), etc. Upper-open ages end at Inf. Source period ends are
inclusive and map to exclusive next-year endpoints. Shortened source
periods are retained. Each required product is downloaded once per call.
For grouped-time products, readers that discard period end labels also
require one `Mx_1x1` request to establish annual coverage; this is
recorded in source metadata. Credentials are passed to the reader and
never stored in the returned table.

## Examples

``` r
if (FALSE) { # \dontrun{
rates <- get_mortality_hmd("USA", year = 2000)
saveRDS(rates, "rates.rds")
model <- hammel_model() |> set_mortality(readRDS("rates.rds"))
history <- get_mortality_hmd("SWE", year = 1891:2000, open_age = 95)
grouped <- get_mortality_hmd("SWE", year = seq(1900, 1995, 5),
  shape = "5x5", open_age = 95)
} # }
```
