# Summarize a matching model

Summary method for class `"hammel_matching_model"`

## Usage

``` r
# S3 method for class 'hammel_matching_model'
summary(object, ...)

# S3 method for class 'summary.hammel_matching_model'
print(x, ...)
```

## Arguments

- object:

  A `hammel_matching_model` object returned by
  [`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md)
  or
  [`age_gap_matching_model()`](https://aarongullickson.github.io/hammel/reference/age_gap_matching_model.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_matching_model` object returned by
  `summary.hammel_matching_model()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_matching_model` containing:

- formula:

  The utility formula.

- coefficients:

  The named utility coefficients.

- age_gap:

  `NULL`, or a list with the implied `preferred_age_gap` and
  `half_weight_distance`, in years, of a quadratic age-gap model.

- group_exogamy:

  The `group_exogamy` specification, or `NULL`.

- delta:

  The `delta` specification, or `NULL`.

- ancestry_endogamy:

  The `ancestry_endogamy` specification, or `NULL`.

- ancestry_resolution:

  The ancestry rounding resolution.

## Details

When the formula is exactly `~ age_difference + I(age_difference^2)`
with a negative quadratic coefficient, as built by
[`age_gap_matching_model()`](https://aarongullickson.github.io/hammel/reference/age_gap_matching_model.md),
the summary also reports the implied preferred age gap and half-weight
distance.

## See also

[`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md),
[`age_gap_matching_model()`](https://aarongullickson.github.io/hammel/reference/age_gap_matching_model.md)

## Examples

``` r
matching <- age_gap_matching_model(
  preferred_age_gap = 3,
  half_weight_distance = 5,
  group_exogamy = list(race = c("W:B" = -3, "W:O" = NA), region = NA)
)
summary(matching)
#> 
#> ── Matching model ──
#> 
#> • Formula: ~age_difference + I(age_difference^2)
#> • Coefficients: age_difference = 0.1664, I(age_difference^2) = -0.02773
#> • Preferred age gap: 3 years
#> • Half-weight distance: 5 years
#> • Group exogamy (race): W:B = -3, W:O = prohibited
#> • Group exogamy (region): all pairs prohibited
```
