# Calculate marriage entry rates including cohabitation conversions

Count direct marriages and conversions to marriage over nonmarried
person-years: unpartnered and cohabiting time both contribute exposure.
Married time contributes none. Each marriage episode counts once per
spouse, including a marriage beginning at conversion. Untyped unions are
interpreted as marriages. Men describe realized outcomes, not male
search-rate inputs.

## Usage

``` r
calculate_marriage_rates(
  population,
  time_breaks,
  age_breaks,
  by_group = FALSE,
  by_previous_union = FALSE,
  conf.level = NULL,
  pool_dissolved = FALSE
)
```

## Arguments

- population:

  A `hammel_population` with person and union histories and a known
  endpoint. Observation cannot extend beyond that endpoint.

- time_breaks:

  Strictly increasing finite calendar boundaries in years.

- age_breaks:

  Strictly increasing nonnegative age boundaries in years. A finite
  final boundary gains an open-ended band through Inf.

- by_group:

  Logical, default FALSE. With TRUE, require person groups and return
  separate rows by each individual's group. No pooled row is added.

- by_previous_union:

  Logical, default FALSE. With TRUE, separate exposure and starts by
  previous union category and include `previous_union`. Otherwise pool
  counts and exposure before calculating rates, and omit that column.
  Supplying `dur_breaks` automatically enables this option.

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

- pool_dissolved:

  Logical, default FALSE. With TRUE, return `"none"` and `"dissolved"`,
  pooling separated and widowed events and exposure before calculating
  rates and confidence intervals. Requires `by_previous_union` (also
  enabled by `dur_breaks`). With FALSE, stratified output uses the
  detailed `"none"`, `"separated"`, and `"widowed"` categories.

## Value

A tibble with `time_start`, `time_end`, `sex`, `age_start`, `age_end`,
`person_years`, `marriages`, and `rate_marriage`, plus `group` and
`previous_union` when requested. Empty cells are retained. Person-years
include cohabiting time. With non-NULL `conf.level`, also includes
`conf_low` and `conf_high`.

## Details

Complete, nonoverlapping histories are required. Cells are half-open;
starts on boundaries enter the cell beginning there. Exposure is
censored at death and the observation endpoint. Rates are NA when
exposure is zero. This age-specific measure does not provide duration
strata. With `by_previous_union = TRUE`, cohabiting exposure and
conversion retain the category held immediately before that cohabitation
began. Ending categories follow the recorded cohabitation policy, as in
[`calculate_union_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_union_rates.md).
For entry into any partnership over unpartnered exposure, use
[`calculate_union_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_union_rates.md).
