# Summarize a fertility process

Summary method for class `"hammel_fertility"`

## Usage

``` r
# S3 method for class 'hammel_fertility'
summary(object, ...)

# S3 method for class 'summary.hammel_fertility'
print(x, ...)
```

## Arguments

- object:

  A `hammel_fertility` object returned by
  [`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_fertility` object returned by
  `summary.hammel_fertility()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_fertility` containing:

- birth_interval:

  The minimum maternal birth-to-birth interval, in years.

- ffab:

  The probability that a newborn is female.

- rate_strata:

  A list describing the input rate schedule, containing `stratified_by`
  (stratum column names), `strata` (a tibble of observed stratum
  combinations), `bands` (a tibble of exact band combinations), and
  `bands_by_stratum` (whether band combinations differ by stratum).
  Shared bands are stored once; otherwise `bands` includes the stratum
  columns. Band endpoints define half-open intervals `[start, end)`.

- weights:

  A summary of the weight model used for selection within rate cells: a
  `summary.hammel_weight_model`, a list of `periods` and `models` for a
  calendar schedule, or `NULL` when selection is unweighted.

- open_band_non_zero:

  Whether any age band ends at `Inf` and has a positive fertility rate.

## See also

[`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
[`calculate_tfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md)

## Examples

``` r
process <- fertility(usa_2000_fertility)
summary(process)
#> 
#> ── Fertility process ──
#> 
#> • Minimum birth interval: 0.75 years
#> • Fraction female at birth: 0.4886
#> • Stratified by: none
#> • Age bands (years): [15, 20), [20, 25), [25, 30), [30, 35), [35, 40), [40, 45)
#> • Weight model: none
```
