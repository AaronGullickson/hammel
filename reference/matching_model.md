# Specify a conditional-logit matching model

Specify conditional-logit model parameters for matching partners.

## Usage

``` r
matching_model(formula, coefficients, group_exogamy = NULL, delta = NULL)
```

## Arguments

- formula:

  A one-sided R formula using `age_difference` and group predictors
  `female_<group>`, `male_<group>` and `same_<group>`. Use `~ 1` and
  [`numeric()`](https://rdrr.io/r/base/numeric.html) coefficients for
  pair-effects-only matching.

- coefficients:

  A finite named numeric vector of utility coefficients.

- group_exogamy:

  Either a list or data.frame indicating the log odds of group exogamous
  pairings. See details for syntax.

- delta:

  Optional male-target adjustments, as a named list keyed by dimension
  or a bare value for the only declared dimension: a finite numeric
  vector with unique category names.

## Value

An inspectable `hammel_matching_model` list containing `formula` and
`coefficients`.

## Details

The `formula` argument can be used to define a formula for the utility
\\U\_{ij}\\ of the match between woman i and man j. This utility is
defined as the log of the unnormalized choice weight, such that the
probability of woman i choosing man j is given by:

\$\$ P(j \mid i) = \frac{\exp(U\_{ij})}{\sum_k \exp(U\_{ik})} \$\$

where the denominator sums across all eligible male partners. This
equation is equivalent to a conditional-logit model specification.

In practice, when groups are in use, this formula needs to be modified
to adjust for expected male group shares, because the model uses women
only to make the choice:

\$\$ P(j\mid i)=
\frac{\omega\_{g(j)}\exp(U\_{ij})}{\sum_k\omega\_{g(k)}\exp(U\_{ik})}
\$\$

where \\\omega_g\\ is a group-specific weight to adjust for male group
shares. More details are provided on the package website.

When `delta = 0`, the \\\omega\\ weights will adjust male group shares
to produce symmetry by gender in group exogamous unions. To allow for
gender asymmetry, users can adjust `delta` values. The `delta` value
will typically need to be tuned.

Users can define formulas based on:

- `age_difference`: male age minus female age

- `female_<group>`: group membership for `<group>` dimension of the
  woman.

- `male_<group>`: group membership for `<group>` dimension of the man.

- `same_<group>`: 0/1 numeric indicating same membership for man and
  woman.

Formula transformations such as `I(age_difference^2)` and
`abs(age_difference)` are supported through
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). If
users limit the formula to `age_difference` and `I(age_difference^2)`,
utilities will be estimated in C++ leading to faster simulation runs.

Coefficient names must match the non-intercept model-matrix columns
exactly, for example `female_raceB`. Names involving group predictors
are checked when the model meets the group declaration at simulation.
Supply all coefficients explicitly; no demographic coefficients are
defaulted.

Users may find it easier to define group exogamy by the use of the
`group_exogamy` argument than the use of a formula. This argument
indicates the likelihood of group exogamy across multiple group
dimensions. It can take one of two forms

- **List form**: a named list keyed by dimension, such as
  `list(race = c("A:B" = -2, "B:C" = NA), region = NA)`. Each value is a
  named vector over unordered pairs (`"A:B"` applies to female A/male B
  and female B/male A alike): a finite value is an additive log-weight
  effect and NA prohibits the pair; `NA` alone prohibits every union
  between different categories of the dimension. Effects from different
  dimensions add. Omitted pairs are treated as implicitly zero.
  Cross-dimension effects are not possible in this format.

- **Table form**, for cross-dimension effects: a data frame with
  `<group>_x` and `<group>_y` columns for each dimension it uses and
  `log_odds`, such as
  `data.frame(race_x = "W", region_x = "S", race_y = "B", region_y = "S", log_odds = -5)`.
  Each row is one unordered pair of category combinations, applied in
  both directions; NA prohibits it; unlisted pairs have an implicit log
  odds of zero.

Users may find
[`age_gap_matching_model()`](https://aarongullickson.github.io/hammel/reference/age_gap_matching_model.md)
more intuitive than specifying a direct formula.

## Examples

``` r
# maximize age difference utility when husband is 5 years older. Double
# utility for same race unions
matching <- matching_model(
  ~ age_difference + I(age_difference^2) + same_race,
  c(
    age_difference = 0.1,
    "I(age_difference^2)" = -0.02,
    same_race = 0.69
  )
)
```
