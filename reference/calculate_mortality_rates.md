# Calculate age-specific mortality rates

Report deaths per person-year, with exact half-open period and age
splitting.

## Usage

``` r
calculate_mortality_rates(
  population,
  time_breaks,
  age_breaks,
  by_union_type = FALSE,
  by_group = FALSE,
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

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

## Value

A tibble containing period, age, sex, person_years, deaths and
rate_mortality, with optional group and union_type strata.
