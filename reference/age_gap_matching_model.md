# Create a quadratic matching model from interpretable age-gap preferences

Construct a matching model using age gap preferences rather than a
formula. Provides an intuitive front-end for
[`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md).

## Usage

``` r
age_gap_matching_model(
  preferred_age_gap,
  half_weight_distance,
  group_exogamy = NULL,
  delta = NULL
)
```

## Arguments

- preferred_age_gap:

  Finite numeric scalar in years providing the preferred age gap where
  the age gap is defined by husband's age minus wife's age. Negative
  values allowed.

- half_weight_distance:

  Finite positive numeric scalar in years from the preferred age gap at
  which utility is halved.

- group_exogamy:

  Optional group exogamy parameters passed to
  [`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md).

- delta:

  Optional male-target adjustments, as a named list keyed by dimension
  or a bare value for the only declared dimension: a finite numeric
  vector with unique category names.

## Value

A `hammel_matching_model` with linear and quadratic age-difference
terms.

## Details

For a matching model that uses a quadratic formulation of age difference
and its square, the formula can be re-specified by the maximum value of
the parabola which gives the "preferred" age gap, and how quickly the
utility declines from that peak. If \\x\\ is the age difference and
\\U\\ is the utility defined as:

\$\$ U = \beta_l x + \beta_q x^2 \$\$

then the value \\x\\ where \\U\\ is maximized is given by:

\$\$ x\_{max} = \frac{\beta_l}{-2\beta_q} \$\$

Solving for \\\beta_l\\ gives:

\$\$ \beta_l = -2 \* \beta_q x\_{max} \$\$

The degree of decline from the peak can be defined by how many years of
age difference \\d\\ are necessary to reduce \\U\\ by half. The formula
for \\\beta_q\\ is given by:

\$\$ \beta_q = \frac{-\log(2)}{d^2} \$\$

Thus, the values of the linear and quadratic terms can be determined by
the desired values of \\x\_{max}\\ and \\d\\.

## Examples

``` r
# preferred age gap of 3 years with a halving at 5 years, plus group
# exogamy terms. Exogamy allowed by race but not region.
matching <- age_gap_matching_model(
  preferred_age_gap = 3,
  half_weight_distance = 5,
  group_exogamy = list(
    race = c("W:B" = -3, "W:O" = -1, "B:O" = -2),
    region = NA
  )
)
```
