# Calculate union exposure and formation rates

Count union starts and unpartnered person-years by age, period, sex, and
optionally `previous_union`: `"none"`, `"separated"`, or `"widowed"`.
The most recent episode's ending determines the category; each start
counts in the category immediately before formation. Partnered time
contributes no exposure. Men describe realized matching outcomes, not a
male search-rate input. Typed entries pool cohabitation and direct
marriage; conversion is excluded because partners remain partnered.
Previous-category classification follows the population's recorded
cohabitation policy, default FALSE without metadata. With `dur_breaks`,
only none rows use age; separated/widowed rows use duration since the
most recent ending; `pool_dissolved` can combine them.

## Usage

``` r
calculate_union_rates(
  population,
  time_breaks,
  age_breaks,
  by_group = FALSE,
  by_previous_union = FALSE,
  dur_breaks = NULL,
  by_union_type = FALSE,
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

- dur_breaks:

  Optional nonnegative duration boundaries in years. When supplied,
  separated/widowed rows use years since the latest union ending; none
  rows still use `age_breaks`. Automatically enables `by_previous_union`
  regardless of its supplied value. A finite final boundary gains Inf.
  Adds `dur_start`/`dur_end`, with NA in the unused age/duration pair.
  Duration rows include all ages; this is not an age-by-duration
  surface.

- by_union_type:

  Logical, default FALSE. Separate cohabitation and direct marriage
  entries. Requires typed histories. Both types use the same full
  unpartnered exposure; exposure must not be summed across type rows.

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
`person_years`, `unions`, and `rate_union`, plus `group` and
`previous_union` when requested, and duration columns when `dur_breaks`
is supplied, and `union_type` when requested. Empty cells are retained.
With non-NULL `conf.level`, also includes `conf_low` and `conf_high`.

## Details

Histories must be complete for the intended interpretation. `"none"`
means no recorded earlier union; omitted history cannot establish
never-partnered status. Overlapping episodes are rejected. Touching
episodes are allowed; tied starts follow union ID order. Time and age
cells are half-open, with events on internal boundaries assigned to the
cell beginning there. Every episode counts, including zero-duration
episodes. Rates are formations per unpartnered person-year, NA at zero
exposure. Search attempts that fail to form unions are not counted, so
these rates can differ from search hazards.
