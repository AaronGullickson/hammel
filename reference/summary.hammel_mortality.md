# Summarize a mortality process

Summary method for class `"hammel_mortality"`

## Usage

``` r
# S3 method for class 'hammel_mortality'
summary(object, ...)

# S3 method for class 'summary.hammel_mortality'
print(x, ...)
```

## Arguments

- object:

  A `hammel_mortality` object returned by
  [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_mortality` object returned by
  `summary.hammel_mortality()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_morality` containing:

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

## See also

[`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
[`calculate_life_table()`](https://aarongullickson.github.io/hammel/reference/calculate_life_table.md)

## Examples

``` r
process <- mortality(usa_2000_mortality)
summary(process)
#> 
#> ── Mortality process ──
#> 
#> • Stratified by: sex (f, m)
#> • Age bands (years): [0, 1), [1, 5), [5, 10), ..., [95, 100), [100, Inf) (22
#>   bands; showing 5)
#> • Weight model: none
```
