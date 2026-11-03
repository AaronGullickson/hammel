# Calculate cohabitation-to-marriage conversion rates

Count conversions over cohabitation union-years, once per episode.
Marriage episodes contribute neither exposure nor events. Conversion,
dissolution and partner death end cohabitation exposure; only
`end_reason = "conversion"` counts as an event. Ongoing episodes are
censored at the observation endpoint. Cells are half-open; boundary
events enter the cell beginning there, and events at the final time
boundary are excluded. Zero-duration episodes can contribute events
without exposure. Rates are NA when exposure is zero.

## Usage

``` r
calculate_cohab_convert_rate(
  population,
  time_breaks,
  dur_breaks = NULL,
  age_breaks = NULL,
  by_group = FALSE,
  conf.level = NULL
)
```

## Arguments

- population:

  A `hammel_population` with typed union histories and a known endpoint.
  Observation cannot extend beyond it.

- time_breaks:

  Strictly increasing finite calendar boundaries in years.

- dur_breaks:

  Optional strictly increasing nonnegative duration boundaries, in years
  since cohabitation began. Supply exactly one of `dur_breaks` and
  `age_breaks`. A finite final boundary gains an open-ended band through
  Inf.

- age_breaks:

  Optional strictly increasing nonnegative boundaries for the woman's
  attained age, matching the conversion schedule's age convention.

- by_group:

  Logical, default FALSE. Separate by the woman's group, matching
  conversion input schedules. Require person groups; retain their input
  type.

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

## Value

A tibble with `time_start`, `time_end`, `dur_start`, `dur_end` (or
`age_start`, `age_end`), `union_years`, `conversions`,
`rate_cohab_convert`, and optional `group`. Empty cells are retained.
Rates are conversions per cohabitation union-year, not conversion
probabilities with competing risks. With non-NULL `conf.level`, also
includes `conf_low` and `conf_high`.
