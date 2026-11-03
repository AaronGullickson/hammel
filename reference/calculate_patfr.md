# Calculate parity- and age-adjusted total fertility rates

Progress a synthetic cohort from parity zero through conditional
age/parity fertility schedules using the HFD period fertility-table
midpoint convention. Women in the final open parity band remain there
and can have repeated births. This is a deterministic fertility table
without mortality or explicit birth spacing, not a microsimulation or an
exact continuous-time hazard solution.

## Usage

``` r
calculate_patfr(rates)
```

## Arguments

- rates:

  A fertility process or a table with `age_start`, `age_end`,
  `parity_start`, `parity_end` and `rate` or `rate_fertility` in births
  per woman-year at that current parity. Optional `time_start`,
  `time_end`, `group`, `union_type` and legacy `union_status` define
  separate strata.

## Value

A tibble with the supplied stratification columns, finite `age_start`,
`age_end`, and `patfr`. Empty results retain their column types. A
supplied `unassigned_births` attribute is preserved.

## Details

Separate the same period, group and union-type strata as
[`calculate_tfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md);
select female rows when sex is present. Parity itself is integrated
over. Each stratum has its own cohort: a union-type-specific index
assumes that type throughout, without partnership or group transitions.
Every age must have complete, nonoverlapping parity coverage from zero
to one open final category. Finite parity bounds must be nonnegative
integers. Each parity must cover the same contiguous finite age range.
Missing included rates yield NA for that stratum. No missing-rate zero
fill is performed here.

Upper-open age rows are omitted, warning on known nonzero omitted female
rates. The returned finite age range defines the calculation: all women
start at parity zero at its lower bound. Split at integer ages and input
boundaries into steps of at most one year. On each step of length d, use
progression probability `d*m / (1 + d*m/2)` and midpoint entry/exit
exposure, following HFD Methods Protocol (2015), sections 8.3-8.4.
Probabilities above one are rejected. Fractional steps generalize the
annual convention. HFD endpoint mapping and missing-cell choices can
change agreement with its published indices. See
<https://www.humanfertility.org/File/GetDocumentFree/Docs/methods.pdf>.
