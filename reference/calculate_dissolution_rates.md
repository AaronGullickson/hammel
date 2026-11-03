# Calculate union-years, dissolutions and duration-specific rates

Split each union episode across calendar periods and elapsed-duration
bands, counting union-years once per episode, not once per spouse. Death
and dissolution both stop exposure, but only
`end_reason = "dissolution"` counts toward dissolutions. Missing ends
are censored at the observation endpoint. Repeated episodes are separate
observations; distinct overlapping records are not merged. Analyze one
independent history per call, then combine counts.

## Usage

``` r
calculate_dissolution_rates(
  population,
  time_breaks,
  dur_breaks,
  by_group = FALSE,
  by_union_type = FALSE,
  conf.level = NULL
)
```

## Arguments

- population:

  A `hammel_population`, or a union data frame containing unique
  positive whole-number `union_id`, `start_time`, `end_time`, and
  `end_reason`. Ends and reasons are NA for ongoing episodes; ended
  reasons are `"dissolution"`, `"death_husband"`, `"death_wife"`, or
  `"conversion"`. Population objects need a known endpoint; observation
  cannot extend beyond it. For bare tables, ensure complete follow-up
  through the requested endpoint.

- time_breaks:

  Strictly increasing finite calendar boundaries in years.

- dur_breaks:

  Strictly increasing nonnegative duration boundaries in years. If the
  last boundary is finite, append Inf for an open final band.

- by_group:

  Logical, default FALSE. With TRUE, use unordered group pairs and
  return `group_pair` (e.g. `"A:A"`, `"A:B"`, `"B:B"`). Both directions
  of a cross-group union share one stratum. Population inputs require
  person groups and retain every unordered pair in the full
  period/duration grid. Bare union tables require `group_pair`; reversed
  labels are normalized. Pair labels use nonempty colon-free character
  representations of groups, ordered alphabetically (radix). Empty
  grouped inputs return a typed zero-row tibble; no pooled row is added.

- by_union_type:

  Logical, default FALSE. Separate cohabitation and marriage episodes
  using each type's own exposure. Requires a valid `union_type` column.
  Conversion stops cohabitation exposure without counting as
  dissolution; duration resets at the start of the linked marriage
  episode.

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

## Value

A tibble with `time_start`, `time_end`, `dur_start`, `dur_end`,
`union_years`, `dissolutions`, and `rate_dissolution`, plus `group_pair`
and `union_type` when requested. Period varies fastest within duration
bands and group blocks. With non-NULL `conf.level`, also includes
`conf_low` and `conf_high`.

## Details

Exposure and event cells are half-open. Dissolutions at internal
boundaries belong to the cell beginning there, and events at the last
time boundary are excluded. Zero-duration episodes may contribute events
with no exposure. Rates are dissolutions per union-year: zero with
positive exposure and no events, and NA when exposure is zero. This is a
cause-specific rate, not the probability of dissolution in the presence
of competing partner deaths.

## Examples

``` r
unions <- data.frame(union_id = 1:3, start_time = c(0, 1, 2),
  end_time = c(3, 4, NA),
  end_reason = c("dissolution", "death_husband", NA))
calculate_dissolution_rates(unions, c(0, 2, 5), c(0, 2, Inf))
#> # A tibble: 4 × 7
#>   time_start time_end dur_start dur_end union_years dissolutions
#>        <dbl>    <dbl>     <dbl>   <dbl>       <dbl>        <dbl>
#> 1          0        2         0       2           3            0
#> 2          2        5         0       2           3            0
#> 3          0        2         2     Inf           0            0
#> 4          2        5         2     Inf           3            1
#> # ℹ 1 more variable: rate_dissolution <dbl>
```
