# Create a quadratic matching model from interpretable age-gap preferences

Convert a preferred age gap and a symmetric half-weight distance into a
formula and coefficients for
[`create_matching_model()`](https://aarongullickson.github.io/hammel/reference/create_matching_model.md).
Age difference is male age minus female age. Candidates at the preferred
gap plus or minus `half_weight_distance` receive half the selection
weight of a candidate at the preferred gap. Availability and eligibility
also affect realized gaps.

## Usage

``` r
create_age_gap_matching_model(
  preferred_age_gap,
  half_weight_distance,
  group_exogamy = NULL,
  prevent_exogamy = FALSE,
  delta = NULL
)
```

## Arguments

- preferred_age_gap:

  Finite numeric scalar in years; negative values allowed.

- half_weight_distance:

  Finite positive numeric scalar in years from the preferred gap at
  which selection weight halves.

- group_exogamy:

  Optional additive group-pair effects passed to
  [`create_matching_model()`](https://aarongullickson.github.io/hammel/reference/create_matching_model.md).

- prevent_exogamy:

  Optional hard group restrictions passed to
  [`create_matching_model()`](https://aarongullickson.github.io/hammel/reference/create_matching_model.md).
  Defaults to `FALSE`.

- delta:

  Optional finite numeric vector with unique group names. Multiply
  mirrored male targets by `exp(delta)` and normalize. Omitted groups
  have zero adjustment. Positive values raise the male target relative
  to baseline and other groups; negative values lower it. Only
  differences matter. Does not set pair-specific ratios or change female
  search rates. Names must belong to the population's identity registry
  (including declared future groups).

## Value

A `hammel_matching_model` with linear and quadratic age-difference
terms.

## Details

The quadratic coefficient is `-log(2) / half_weight_distance^2` and the
linear coefficient is `-2 * preferred_age_gap * quadratic`. This
specifies relative selection weights, not a halving of utility or a
realized standard deviation.

## Examples

``` r
# Illustrative preferences, not demographic defaults.
matching <- create_age_gap_matching_model(3, 5)
```
