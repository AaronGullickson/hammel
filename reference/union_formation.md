# Specify a union-formation process

Constructs a union formation process specification for a simulation
model. The convenience function `set_union_formation` is equivalent to
`set_process(model, union_formation(...))` and allows for direct piping
from a model.

## Usage

``` r
union_formation(
  initiation,
  rates = NULL,
  matching = age_gap_matching_model(0, 5),
  father_matching = NULL,
  cohab_matching = NULL,
  cohab_prob = NULL,
  union_convert_rate = NULL,
  weights = NULL,
  weights_conversion = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  cancel_birth_if_no_match = TRUE,
  cohab_as_previous_union = FALSE,
  diagnostics = FALSE
)

set_union_formation(
  model,
  initiation,
  rates = NULL,
  matching = age_gap_matching_model(0, 5),
  father_matching = NULL,
  cohab_matching = NULL,
  cohab_prob = NULL,
  union_convert_rate = NULL,
  weights = NULL,
  weights_conversion = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  cancel_birth_if_no_match = TRUE,
  cohab_as_previous_union = FALSE,
  diagnostics = FALSE
)
```

## Arguments

- initiation:

  Either `"on_demand"` for birth-triggered searches or `"rate"` for
  scheduled searches. No default is selected.

- rates:

  For rate mode, a nonempty data frame defining the marriage rate
  schedules. Must contain at least numeric `age_start`, `age_end`, and
  `rate`. See details for description of stratification variables that
  can be applied.

- matching:

  A specification created by
  [`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md),
  or a tibble with numeric `time_start`, `time_end`, and list-column
  `model` holding complete matching objects. Defaults to
  `age_gap_matching_model(0, 5)`: preferred age gap zero years, with
  selection weight halving five years on either side.

- father_matching:

  Optional matching specification for fathers of unpartnered births in
  rate mode, created by
  [`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md),
  or a `time_start`/`time_end`/`model` schedule as for `matching`.
  Defaults to `NULL`, using the general `matching` model active at that
  birth.

- cohab_matching:

  Optional matching model or whole-model schedule for cohabitation; NULL
  shares the general `matching` model.

- cohab_prob:

  Optional age-specific probability table defining the probability that
  a union is a cohabitation. See details for specification.

- union_convert_rate:

  Optional age- or duration-specific conversion table, defining
  conversion rates of cohabitations into marriages. See details for
  specification.

- weights:

  NULL (default) for unweighted selection of searching women, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  using person variables, which chooses which unpartnered woman in a
  cell searches without changing the cell's input search rate. Rate mode
  only.

- weights_conversion:

  NULL (default) for unweighted selection, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  using union variables, which chooses which cohabitation in a cell
  converts to marriage without changing the cell's input rate.

- min_male_age:

  Minimum male age in continuous years for partner eligibility, default
  `15`. Finite nonnegative scalar; zero disables this age exclusion.

- kinship_generations:

  Generation bound on both sides of a known common ancestor defining
  partner ineligiblity, default `2`. Excludes parents, grandparents,
  siblings, aunts/uncles, nieces/nephews and first cousins, including
  half relationships, when known. Zero imposes no cross-person kin
  exclusion. See
  [`extract_kin()`](https://aarongullickson.github.io/hammel/reference/extract_kin.md)
  for the shared distance definitions.

- cancel_birth_if_no_match:

  Logical, default `TRUE`. Cancel a birth when its on-demand partner
  search or rate-mode nonmarital father selection fails. Set `FALSE` to
  allow births with unknown fathers. Historical missing father IDs are
  retained. Partner scarcity can reduce realized fertility below inputs.

- cohab_as_previous_union:

  Logical, default FALSE. FALSE only marriage defines the
  `previous_union` value; TRUE cohabitation or marriage define the
  `previous_union` value.

- diagnostics:

  Logical, default `FALSE`. When `TRUE`, adds additional technical
  variables from matching model to each union record.

- model:

  A `hammel_model` object.

## Value

`union_formation()` returns a `hammel_union_formation` process;
`set_union_formation()` returns an updated model.

## Details

The union formation process is female-centric. Women initiate searches
and select male partners. Two different approaches to union search
initiation are provided for the required `initiation` argument:

- The "on_demand" approach will initiate a union search when an
  unpartnered woman gives birth. This approach should not be used with
  fertility rates that are stratified by `union_type` because
  unpartnered fertility is what creates the search. This approach is
  useful if the user has no particular interest in union formation
  processes but wants to assign parentage correctly for kinship
  determination.

- The "rate" approach uses female marriage rates defined by the `rates`
  argument. The input should be a marriage rate schedule that may be
  stratified by group dimensions and/or `previous_union`. When
  stratified by `previous_union`, remarriages may be duration-specific
  rather than age-specific. This approach should generally be used with
  fertility rates that are stratified by `union_type`, although it is
  not required. Failure to stratify by `union_type` will lead to
  independence between fertility and union formation.

Once a search is initiated, a woman will select a male partner from
among all eligible male partners. Men related to the women are
considered ineligible out to a kinship depth defined by
`kinship_generations` and men below `min_male_age` are also ineligible.
the decision itself is determined by the provided `matching_model` which
uses a conditional logit formulation to score potential suitors by age
differences and group exogamy. See
[`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md)
for more details.

It is possible for this search to fail if no eligible men are available.
When using the "on_demand" method, the `cancel_birth_if_no_match`
argument controls whether the birth continues if no match is found. If
set to FALSE, the birth continues but with `NA` values for the father.
If set to TRUE, observed fertility rates may deviate from input rates.

When unpartnered fertility is allowed, the `matching_model` will also be
used to identify the father of a birth to an unpartnered woman. Users
can optionally specify a separate `father_matching` model for this case.

Cohabitation can be added to the union formation process by specifying
`cohab_prob` and `union_convert_rate`. The `cohab_prob` argument expects
a data frame that is similar to a rate process object. It contains
age-specific bands defined by `age_start` and `age_end`, but also
contains a `prob` column rather than a `rate` column. This `prob` column
defines the probability that a union formed by a woman of the given age
is a cohabitation instead of a marriage. Like rate processes,
`cohab_prob` can be stratified by time period using `time_start` and
`time_end`, group dimensions, or by `previous_union`.

The `union_convert_rate` is not required, but is necessary for
cohabitations to be able to convert into marriages. It follows the
standard rate process format and can be either duration or age
dependent. It can also be stratified by time period and group dimension.

By default, partner selection for cohabitation will use the
`matching_model` argument. Users can specify a matching model specific
to cohabitations with the `cohab_matching` argument.

## Examples

``` r
union_process <- union_formation(
  initiation = "rate",
  rates = usa_2000_marriage,
  matching = age_gap_matching_model(3, 5)
)

model <- hammel_model() |>
  set_fertility(usa_2000_fertility_union) |>
  set_mortality(usa_2000_mortality) |>
  set_process(union_process)
```
