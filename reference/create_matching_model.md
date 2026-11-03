# Specify a conditional-logit matching model

Define relative partner utilities with an R formula and named
coefficients. Predictors include `age_difference`: male age minus female
age, in continuous years at the search time. Formula transformations
such as `I(age_difference^2)` and `abs(age_difference)` are supported
through
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). An
intercept is ignored because it cancels from partner-only
conditional-logit probabilities. No wait alternative is included.
`female_group` and `male_group` are factors with explicit `group_levels`
and fixed treatment contrasts (first level is reference). `same_group`
is numeric 1 for equal identities, 0 otherwise. Pair effects are added
after the formula, including custom formulas. Group effects require
population identities. These modify relative candidate weights, not
search hazards or eligibility.

## Usage

``` r
create_matching_model(
  formula,
  coefficients,
  group_exogamy = NULL,
  group_levels = NULL,
  prevent_exogamy = FALSE,
  delta = NULL
)
```

## Arguments

- formula:

  A one-sided R formula using `age_difference`, `female_group`,
  `male_group`, and/or `same_group`. Use `~ 1` and
  [`numeric()`](https://rdrr.io/r/base/numeric.html) coefficients for
  pair-effects-only matching.

- coefficients:

  A finite named numeric vector of utility coefficients.

- group_exogamy:

  Optional finite named numeric vector of symmetric group-pair effects:
  `"A:B"` applies the same log-weight effect to female A/male B and
  female B/male A. Specify each unordered pair once; reverse duplicates
  and `::` syntax are errors. Omitted pairs have zero effect;
  within-group entries are allowed. Alternatively provide a data frame
  with `female_group`, `male_group`, and `log_odds`; each row specifies
  an unordered pair and expands symmetrically, including for
  colon-containing labels. Effects add to formula scores and do not
  prohibit matches.

- group_levels:

  Optional ordered vector of group labels (character or numeric codes).
  Required with at least two levels when the formula uses `female_group`
  or `male_group`; must cover population groups. Sets fixed treatment
  reference/columns, independent of the available pool. Not required for
  `same_group` or pair effects alone.

- prevent_exogamy:

  `FALSE` (default) for no group restrictions; `TRUE` excludes all
  different-group candidates. A character vector such as
  `c("A::B", "A:C")` excludes selected pairs with female-first `:` and
  symmetric `::` notation. Alternatively supply a directional data frame
  with `female_group` and `male_group`. Pairs must join different
  groups. Restricted men are ineligible, alongside age and kinship
  exclusions. A search fails only when no eligible man exists in the
  whole pool, following existing failed-search and birth cancellation
  rules. Restrictions override finite preferences, require population
  identities, and do not invalidate existing partners.

- delta:

  Optional finite numeric vector with unique group names. Multiply
  mirrored male targets by `exp(delta)` and normalize. Omitted groups
  have zero adjustment. Positive values raise the male target relative
  to baseline and other groups; negative values lower it. Only
  differences matter. Does not set pair-specific ratios or change female
  search rates. Names must belong to the population's identity registry
  (including declared future groups).

## Value

An inspectable `hammel_matching_model` list containing `formula` and
`coefficients`. This specifies a model; it does not estimate one.

## Details

Coefficient names must match the non-intercept model-matrix columns
exactly. Supply all coefficients explicitly; no demographic coefficients
are defaulted. Formula functions must be deterministic and must not
mutate simulation state. Data-dependent transformations can change with
the candidate pool; use fixed transformations when coefficients should
have a fixed interpretation.

## Examples

``` r
# Illustrative coefficients, not demographic defaults.
matching <- create_matching_model(
  ~ age_difference + I(age_difference^2),
  c(age_difference = 0.1, "I(age_difference^2)" = -0.02)
)
```
