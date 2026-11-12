# Calculate total fertility and parity-adjusted total fertility rates

Constructs total fertility rates or parity-adjusted total fertility rats
from input tables that provide age-specific fertility rates. Input data
can come from user-provided input rate schedules, fertility process
objects, or estimated rates from a simulation via
[calculate_fertility_rates](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
and may be stratified by additional variables.

## Usage

``` r
calculate_tfr(rates)

calculate_patfr(rates)
```

## Arguments

- rates:

  A table returned by
  [`calculate_fertility_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
  a fertility process object from
  [`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
  or a data frame with at least `age_start`, `age_end`, and
  `rate_fertility` or `rate` columns. Input may be further stratified by
  other variables.

## Value

A tibble with the supplied period and other stratification columns,
finite `age_start` and `age_end`, and either `tfr` or `patfr`.

## Details

When the supplied `rates` table contains columns other than `age_start`,
`age_end`, and either `rate_fertility` or `rate`, the function will
treat these additional columns as stratifying dimensions and will
stratify the output appropriately. It will also stratify output by time
when `time_start` and `time_end` are included. The function also ignores
additional calculated columns from the output of
[calculate_fertility_rates](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
so the results of that function can be directly piped into this
function.

When the input `rates` table is stratified by parity with `parity_start`
and `parity_end` columns, the function `calculate_patfr` will be
automatically used in place of `calculate_tfr` to calculate
parity-adjusted total fertility rates. See See
<https://www.humanfertility.org/File/GetDocumentFree/Docs/methods.pdf>
for details on the estimation of parity-adjusted total fertility rates.

The function will automatically remove any open-ended age bands because
TFR cannot be estimated across open-ended age bands. If this band
contains a non-zero rate, the function will produce a warning that TFR
may be underestimated.

Missing values on rates will result in missing values for the total
fertility rate.

## See also

[`calculate_fertility_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
[`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md)

## Examples

``` r
# TFR from age-specific rates
calculate_tfr(usa_2000_fertility)
#> # A tibble: 1 × 3
#>   age_start age_end   tfr
#>       <int>   <int> <dbl>
#> 1        15      45  2.05

# TFR stratified by union type
calculate_tfr(usa_2000_fertility_union)
#> # A tibble: 2 × 4
#>   union_type  age_start age_end   tfr
#>   <chr>           <int>   <int> <dbl>
#> 1 partnered          15      45  4.20
#> 2 unpartnered        15      45  1.17
```
