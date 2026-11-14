# Summarize a union dissolution process

Summary method for class `"hammel_union_dissolution"`

## Usage

``` r
# S3 method for class 'hammel_union_dissolution'
summary(object, ...)

# S3 method for class 'summary.hammel_union_dissolution'
print(x, ...)
```

## Arguments

- object:

  A `hammel_union_dissolution` object returned by
  [`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_union_dissolution` object returned by
  `[summary.hammel_union_dissolution()]`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_union_dissolution` containing:

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

[`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md)

## Examples

``` r
diss_rates <- tibble::tribble(
 ~union_type,  ~dur_start, ~dur_end, ~rate,
  "cohab",      0,          5,        0.08,
  "cohab",      5,          10,       0.04,
  "cohab",      10,         Inf,      0.02,
  "marriage",   0,          5,        0.06,
  "marriage",   5,          10,       0.02,
  "marriage",   10,         15,       0.01,
  "marriage",   15,         Inf,      0.015
)
process <- union_dissolution(diss_rates)
summary(process)
#> 
#> ── Union dissolution process ──
#> 
#> • Stratified by: union_type (cohab, marriage)
#> • Bands by stratum
#>   • union_type = cohab
#>     • Duration bands (years): [0, 5), [5, 10), [10, Inf)
#>   • union_type = marriage
#>     • Duration bands (years): [0, 5), [5, 10), [10, 15), [15, Inf)
#> • Weight model: none
```
