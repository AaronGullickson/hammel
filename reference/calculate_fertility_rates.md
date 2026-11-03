# Calculate age-specific fertility rates, optionally by current parity

Report maternal births per woman-year, including time excluded by
minimum birth spacing. Only female rows are returned. Current parity is
reconstructed from all recorded maternal births, including pre-window
births and deceased children; custom person parity columns are ignored.
Each birth counts at the parity immediately before it and advances
subsequent exposure. Simultaneous births progress successively in
recorded person-row order.

## Usage

``` r
calculate_fertility_rates(
  population,
  time_breaks,
  age_breaks,
  by_union_type = FALSE,
  by_group = FALSE,
  parity_breaks = NULL,
  unrealized = FALSE,
  conf.level = NULL
)
```

## Arguments

- population:

  A hammel_population or persons data frame.

- time_breaks:

  Strictly increasing finite period boundaries in years.

- age_breaks:

  Strictly increasing nonnegative age boundaries in years. If the final
  boundary is finite, Inf is appended to create an open-ended final age
  band starting at that boundary.

- by_union_type:

  Logical, default FALSE. Require population histories and return
  `union_type` rows for unpartnered and each discovered union type,
  without overall rows. Untyped histories use partnered/unpartnered.
  Births tied to dissolution or conversion retain the ending type;
  deaths use the type immediately before death. Typed overlapping
  histories are rejected; untyped overlapping coverage is merged. For
  formation rates use
  [`calculate_union_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_union_rates.md)
  or
  [`calculate_marriage_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_marriage_rates.md).

- by_group:

  Logical, default `FALSE`. With `TRUE`, require a persons `group`
  column and return the complete period/age/sex grid for each group
  appearing in the person history, in first-occurrence order. Labels
  retain their type. Exposure and deaths use the person's group; births
  use the mother's group, even when the child has a different group.
  Combines with `by_union_type`. No pooled group row is added.
  Unknown-mother births remain in the single period-level
  `unassigned_births` attribute, without a group assignment.

- parity_breaks:

  NULL to pool parities, or strictly increasing nonnegative integer
  boundaries starting at zero. Bands are half-open. A finite final
  boundary receives Inf, as for age breaks.

- unrealized:

  Logical, default `FALSE`. With `TRUE`, the population must come from
  [`simulate()`](https://aarongullickson.github.io/hammel/reference/simulate.md);
  add `unrealized_births` (expected births dropped from fertility cells
  that had members but no eligible woman, from
  `simulation$unrealized_births`) and
  `rate_fertility_adjusted = (births + unrealized_births) / person_years`.
  Values are `NA` where they cannot be attributed exactly: periods
  before the most recent simulation's start (restarts keep only the
  latest run's record), breaks that cut through a recorded one-year
  period or fertility cell, and `union_type` or parity strata the
  fertility schedule did not use.

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

## Value

A tibble with period, age, sex (female), person_years, births and
rate_fertility; optional parity_start/parity_end, group and union_type,
and with `unrealized = TRUE`, unrealized_births and
rate_fertility_adjusted. Unknown-mother births are in the period-level
unassigned_births attribute.
