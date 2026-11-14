# Summarize a union formation process

Summary method for class `"hammel_union_formation"`

## Usage

``` r
# S3 method for class 'hammel_union_formation'
summary(object, ...)

# S3 method for class 'summary.hammel_union_formation'
print(x, ...)
```

## Arguments

- object:

  A `hammel_union_formation` object returned by
  [`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_union_formation` object returned by
  `summary.hammel_union_formation()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_union_formation` containing:

- method:

  The union initiation method: `"rate"` or `"on_demand"`.

- min_male_age:

  The minimum eligible male age, in years.

- kinship_generations:

  The configured number of generations used for kinship exclusion when
  selecting partners.

- cancel_birth_if_no_match:

  Whether a birth is cancelled when its on-demand partner search or
  rate-mode nonmarital father selection fails, rather than occurring
  with an unknown father.

- cohab_as_previous_union:

  Whether cohabitations count as previous unions when determining
  previous-union status.

- rate_strata:

  A list describing the union-search rate schedule, or `NULL` when no
  rate schedule is configured. The list contains `stratified_by`
  (stratum column names), `strata` (a tibble of observed stratum
  combinations), `bands` (a tibble of exact band combinations), and
  `bands_by_stratum` (whether band combinations differ by stratum).
  Shared bands are stored once; otherwise `bands` includes the stratum
  columns. Band endpoints define half-open intervals `[start, end)`.

- cohabitation:

  `NULL` when cohabitation is disabled; otherwise a list containing
  `cohab_prob_range` (the minimum and maximum supplied cohabitation
  probabilities), `cohab_rate_strata` (strata and bands of the
  cohabitation probability schedule), `convert_rate_strata` (strata and
  bands of the conversion rate schedule, or `NULL` when no conversion
  schedule is configured), `matching` (a summary of the cohabitation
  matching model, or `NULL` when the general matching model is used) and
  `conversion_weights` (a summary of the conversion weight model, or
  `NULL` when unweighted). Both schedule summaries use the same
  structure as `rate_strata`.

- matching:

  A summary of the matching model: a `summary.hammel_matching_model`, or
  a list of `periods` and `models` for a calendar schedule.

- father_matching:

  A summary of the matching model for fathers of unpartnered births, as
  for `matching`, or `NULL` when the general matching model is used.

- search_weights:

  A summary of the weight model choosing which woman searches: a
  `summary.hammel_weight_model`, a list of `periods` and `models` for a
  calendar schedule, or `NULL` when searching women are selected without
  weights.

## See also

[`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md)

## Examples

``` r
process <- union_formation("rate", usa_2000_marriage)
summary(process)
#> 
#> ── Union formation process ──
#> 
#> • Search method: rate
#> • Stratified by: none
#> • Age bands (years): [15, 20), [20, 25), [25, 30), ..., [65, 70), [70, Inf) (12
#>   bands; showing 5)
#> • Search weight model: none
#> • Matching model: ~age_difference + I(age_difference^2)
#>   • Coefficients: age_difference = 0, I(age_difference^2) = -0.02773
#>   • Preferred age gap: 0 years
#>   • Half-weight distance: 5 years
#> • Father matching model: same as matching model
#> • Births without a matched father: cancelled
#> • Male eligibility restrictions
#>   • Minimum male age: 15
#>   • Generational proximity: 2
#> • Cohabitation: disabled
```
