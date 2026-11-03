# Calculate total fertility rates

Sum age-specific fertility rates times their actual interval widths. If
parity-bound columns are present, delegate to
[`calculate_patfr()`](https://aarongullickson.github.io/hammel/reference/calculate_patfr.md)
instead and return its `patfr` column rather than `tfr`, with the same
other strata. Open-ended age bands are omitted, warning if an omitted
female rate is known and nonzero that TFR underestimates fertility over
the full age range. Zero or missing omitted rates do not trigger that
warning. The remaining finite, contiguous age range defines the TFR
reported. Any missing ASFR in those bands makes the TFR missing. When
`union_type` (or a legacy `union_status` column) is present, calculate a
separate sum for each stratum within each period. These conditional sums
describe fertility under continuous exposure to that status; they are
not components to add into an overall TFR. A `group` column
automatically adds separate group strata, with its type and labels
preserved; groups are not pooled. Age bands are checked independently
within each period/group/status combination.

## Usage

``` r
calculate_tfr(rates)
```

## Arguments

- rates:

  A table returned by
  [`calculate_fertility_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_fertility_rates.md),
  a fertility process object, or a data frame with `age_start`,
  `age_end`, and `rate_fertility` or `rate` in births per woman-year.
  Bare `rate` tables are interpreted as fertility rates; mortality
  processes are not valid inputs.

## Value

With parity columns, the result of
[`calculate_patfr()`](https://aarongullickson.github.io/hammel/reference/calculate_patfr.md)
with `patfr`. Otherwise, a tibble with period columns, `union_type` or
legacy `union_status` when supplied, `group` when supplied, `age_start`,
`age_end`, and `tfr`. If `sex` is present, only female rows are used.
