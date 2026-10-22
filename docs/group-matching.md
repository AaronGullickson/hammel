# Group effects in partner matching

Requested by the user October 20, 2026. This extends matching utility, not
eligibility, search hazards, group identity, or ancestry inheritance.

## Pair effects

`create_matching_model(..., group_exogamy = NULL, group_levels = NULL)` accepts
finite additive log-selection-weight effects as a named numeric vector.
`"A:B"` means female A and male B; `"A::B"` supplies the same value in both
directions. For example, `c("A:B" = -1, "B:A" = -0.5, "A::C" = -0.8)`.
Omitted pairs have odds ratio 1.0 and contribute zero on the log scale;
within-group entries such as `"A:A"` are allowed.
Effects are added after the formula is evaluated. An effect of -1 multiplies
relative candidate weight by `exp(-1)`, holding other utility terms fixed.
It is not a target intergroup union proportion or a marriage hazard multiplier.

Reject overlaps after expanding symmetric entries, rather than applying
precedence. `A:B` and `A::B`, or `A::B` and `B::A`, conflict. `A::A` describes
only one within-group pair. Reject malformed separators, missing labels,
nonfinite values, and unnamed vectors. Finite values do not prohibit matches;
multiple identities remain outside this increment. Hard group restrictions are
specified separately below.

## Hard group restrictions

The user selected post-sampling restrictions on October 20, 2026.
`prevent_exogamy = FALSE` on a matching model imposes none. TRUE excludes
every different-group candidate. A character vector such as
`c("A::B", "A:C")` excludes selected cross-group pairs, using the same
female-first and symmetric syntax as `group_exogamy`. A directional table with
`female_group` and `male_group` is also accepted for colon-containing labels.
Listed pairs must join different groups; overlapping entries are errors.

Draw from the original living, unpartnered male pool, then apply age, kinship,
and group screening to that sample without refilling. Group restrictions win
over any finite utility preference. Do not rescore or replace existing partners.
An empty eligible sample follows the existing search-retry/birth-cancellation
rules, including births with unknown fathers only when explicitly allowed.
Restrictions belong to the selected matching model: NULL `father_matching`
inherits them, while an explicit father model has its own restrictions.
Diagnostic eligible counts and utility ranks refer to the final eligible set.

This retains the possibility that a minority woman samples no compatible men
even when some exist in the full pool. With equally aged compatible A and
incompatible B men, a one-man sample may fail; the full pool will find A.
A finite penalty of -10 still selects B from a sample containing B alone,
whereas prohibiting A:B makes that sample empty. Increasing `sample_size`
changes encounter coverage, not the restriction. No extra random draws occur
in group screening; FALSE must retain the existing histories and RNG state.

For labels containing colons, accept an alternative directional data frame
with `female_group`, `male_group`, and `log_odds`. Duplicate directed pairs
are errors; use two rows for a symmetric table effect. Numeric group codes
are compared through their character representation, as in group-rate lookup.
`create_age_gap_matching_model(..., group_exogamy = NULL)` also accepts effects.
Pair effects require population group identities; missing group data errors
at execution, including zero-duration calls. Extra pair labels do not create
population groups and are harmless when absent from the population.

## Formula predictors

Matching formulas may use `age_difference`, `female_group`, `male_group`, and
`same_group`. Ages remain male minus female in continuous years. `same_group`
is numeric 1 for equal recorded identities and 0 otherwise. The two group
variables are unordered factors with fixed levels from explicit `group_levels`;
the first level is the treatment reference. Categorical group formulas require
at least two declared levels, covering every population group. Fixed treatment
contrasts and retained levels prevent a candidate sample from changing model
columns or reference categories. Unknown levels error, not imputation.
Formulas using only `same_group` do not require declared categorical levels.

Coefficient names must match all non-intercept design columns exactly. An
intercept remains ignored, since it cancels in conditional selection. A formula
`~ 1` with `numeric()` coefficients permits pair-effects-only matching. Formula
functions remain deterministic; fixed group identity information is supplied
as candidate data, never mutable state. Under default female choice,
female-only terms are constant within a search and cancel. In the experimental
two-sided method, they can vary across a man's female alternatives and affect
his choice. Group effects in both the formula and pair table intentionally add
together.

## Execution and validation

Keep the native linear/quadratic age scorer; add compiled pair offsets afterward
without an R callback. Custom formulas receive the group predictors when used,
and receive the same additive pair offsets afterward. Both on-demand and rate
initiation use these scores. `father_matching = NULL` reuses the full matching
model for nonmarital father choice; an explicit father model replaces it,
including its pair effects. Current partners are retained without rematching.
Diagnostic utility ranks use the final combined utility. Neither sampling,
screening, event timing, nor the number of RNG draws changes.

Worked case: female A faces equally aged candidates in A, B, and C. Effects
`c("A:B" = log(2), "A::C" = log(3))` and a zero formula give probabilities
`c(1,2,3)/6`. A female B with A/B candidates has equal weights unless `B:A`
is separately specified. A female C gives male A weight 3 from the symmetric
entry. Tests cover these analytical scores/probabilities, stable factor designs
and interactions, native/R agreement, malformed or conflicting specifications,
matching modes and father overrides, diagnostic ranks, reproducibility,
empty populations, and unchanged behavior without group preferences.

## Experimental two-sided selection (October 22, 2026)

Current default selection normalizes scores within each initiating woman's
candidate set. Symmetric group coefficients do not guarantee equal directional
union or birth counts with unequal group sizes.

The user authorized an optional `matching_method = "two_sided"` experiment on
`experiment/two-sided-sequential-matching` (D061). A temporary market contains
M ego-eligible men, M - 1 sampled alternate women, and ego. D062 replaces
sequential depletion with fixed utilities. Normalize ego's row and each man's
column, add the probabilities of either choosing the other, normalize, and
draw once. Each woman's allowed male row is evaluated once, including for R
formulas; column normalization reuses those same scores. No hypothetical pairs
are removed. This integrates repeated fair-coin uniform initiation with replacement.
See [two-sided matching](two-sided-matching.md) for the full contract,
shortage and exclusion rules, worked examples, validation, and timing plan.

This differs from the earlier proposed matrix rescaling or sampling assignments
in proportion to `exp(sum of utilities)`. Neither rescaling by female demand
nor fixed-margin IPF was selected. The latter's constraints and the behavior
of directed and age effects motivated the two-sided experiments.
The optional prototype does not establish aggregate odds-ratio recovery,
preferred age-gap recovery, or balanced inheritance trajectories. Scientific
assessment remains with the user and the milestone 3 analyses.
