# Group effects in partner matching

Requested by the user October 20, 2026. This extends matching utility, not
eligibility, search hazards, group identity, or ancestry inheritance.

Milestone 3 was accepted by the user October 27, 2026 (D077). The
[validation report](../validation/milestone3.qmd) records the assessment;
acceptance retains the documented limitations.

## Pair effects

Since D134 (stage 5), group predictors and settings name a declared group
dimension: predictors `female_<dim>`, `male_<dim>` and `same_<dim>`, and
`group_exogamy` and `delta` as lists keyed by dimension,
such as `group_exogamy = list(race = c("A:B" = -1))`; a bare value applies to
the only declared dimension. The matching dimensions are every dimension any
matching model of the run names (general, father, cohabitation, and schedule
entries); the engine's matching category is their combination, and effects
from different dimensions add (step 2 of D134, specified in
[multiple group dimensions](multiple-groups.md#step-2-specification-matching-across-several-dimensions)).
With one matching dimension, each model is translated to the internal
one-dimension form described here, with `group` standing for that dimension.
The examples below use a dimension named `group`.

`create_matching_model(..., group_exogamy = NULL)` accepts
finite additive log-selection-weight effects as a named numeric vector.
`"A:B"` supplies the same value for female A/male B and female B/male A.
For example, `c("A:B" = -1, "A:C" = -0.8)`. Directional pair effects and
`::` syntax are not supported (D075).
Omitted pairs have odds ratio 1.0 and contribute zero on the log scale;
within-group entries such as `"A:A"` are allowed.
Effects are added after the formula is evaluated. An effect of -1 multiplies
relative candidate weight by `exp(-1)`, holding other utility terms fixed.
It is not a target intergroup union proportion or a marriage hazard multiplier.

Reject repeated unordered pairs, including `A:B` with `B:A`, even with equal
values. `A:A` describes one within-group pair. Reject malformed separators,
missing labels, infinite values, and unnamed vectors. Finite values do not
prohibit matches; NA does (below).

## Hard group restrictions

The user selected hard restrictions on October 20, 2026; since D072 they apply
to the whole eligible pool rather than to a candidate sample. Since D135
(November 6, 2026) they are written in `group_exogamy` itself; the separate
`prevent_exogamy` argument and its directional syntax are removed, so every
pair, effect or prohibition, is symmetric:

```r
group_exogamy = list(race = c("A:B" = -2, "B:C" = NA), region = NA)
```

- `"A:B" = -2`: an effect of -2 for female A/male B and female B/male A.
- `"B:C" = NA`: B–C unions are prohibited in both directions.
- An unlisted pair (here A:C) has no effect.
- A whole dimension set to NA (`region = NA`) prohibits every union between
  different categories of it; within-category unions are allowed.

With one declared dimension a bare value works the same way:
`group_exogamy = c("A:B" = -2, "B:C" = NA)` or `group_exogamy = NA`. A
within-category prohibition (`"A:A" = NA`) is an error, since it is not
exogamy; a finite `"A:A"` effect remains allowed. In the pair-table form an NA
`log_odds` prohibits that unordered pair. With several dimensions, a pair of
combinations is prohibited when any dimension prohibits its categories, so in
the example above unions across regions are prohibited while race pairs carry
their effects.

Restricted men are ineligible, alongside age and kinship exclusions; the exact
draw skips restricted group cells. Group restrictions win over any finite
utility preference. Do not rescore or replace existing partners. A search fails
only when no eligible man exists anywhere in the pool, following the existing
search-retry/birth-cancellation rules, including births with unknown fathers
only when explicitly allowed. Restrictions belong to the selected matching
model: NULL `father_matching` inherits them, while an explicit father model has
its own restrictions.

With equally aged compatible A and incompatible B men, an A woman under a
prohibition always finds the A man. A finite penalty of -10 still selects B
when B men are the only ones available, whereas prohibiting A:B (`"A:B" = NA`)
makes the search fail. Without NA entries there is no restriction.

For labels containing colons, accept a data frame with `female_group`,
`male_group`, and `log_odds`. These column names are retained for compatibility,
but each row denotes an unordered pair and supplies the effect in both
directions. Duplicate unordered pairs are errors. Numeric group codes
are compared through their character representation, as in group-rate lookup.
`create_age_gap_matching_model(..., group_exogamy = NULL)` also accepts effects.
Pair effects require population group identities; missing group data errors
at execution, including zero-duration calls. Extra pair labels do not create
population groups and are harmless when absent from the population.

## Formula predictors

Matching formulas may use `age_difference`, `female_group`, `male_group`, and
`same_group`. Ages remain male minus female in continuous years. `same_group`
is numeric 1 for equal recorded identities and 0 otherwise. The two group
variables are unordered factors whose levels are the dimension's categories
(declared, or in order of first appearance); the first is the treatment
reference. Categorical group formulas require at least two categories, and
their coefficient names are checked against the categories at simulation. Fixed treatment
contrasts and retained levels prevent a candidate sample from changing model
columns or reference categories. Unknown levels error, not imputation.
Formulas using only `same_group` do not require declared categorical levels.

Coefficient names must match all non-intercept design columns exactly. An
intercept remains ignored, since it cancels in conditional selection. A formula
`~ 1` with `numeric()` coefficients permits pair-effects-only matching. Formula
functions remain deterministic; fixed group identity information is supplied
as candidate data, never mutable state. Female-only terms are constant within a search and cancel; calibrated
weights depend only on the man's group. Group
effects in both the formula and pair table intentionally add together.

## Execution and validation

Keep the native linear/quadratic age scorer; add compiled pair offsets afterward
without an R callback. Custom formulas receive the group predictors when used,
and receive the same additive pair offsets afterward. Both on-demand and rate
initiation use these scores. `father_matching = NULL` reuses the full matching
model for nonmarital father choice; an explicit father model replaces it,
including its pair effects. Current partners are retained without rematching.
Neither screening nor event timing changes.

Worked score case: female A faces equally aged candidates in A, B, and C.
Effects `c("A:B" = log(2), "A:C" = log(3))` and a zero formula give uncalibrated
relative weights `c(1,2,3)`. Female B sees A/B weights 2/1; female C sees A/B
weights 3/1. Actual selection also includes calibrated male-group weights.
Tests cover these utilities, stable factor designs and interactions, native/R
agreement, malformed/repeated unordered pairs, father overrides, reproducibility,
empty populations, and unchanged behavior without group preferences.

## Directional pairing

All matching uses calibrated male-group weights (D075), with target adjustments
from `delta` (D074); see [weighted-choice matching](weighted-choice-matching.md).
Symmetric preferences with zero delta and balanced group-specific sexes give
balanced two-group directional flows in the static illustration. Delta changes
male-group targets without changing static group-only cross-products. No
separate unweighted female-choice execution mode remains.
