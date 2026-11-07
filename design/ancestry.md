# Group ancestry proportions

Requested by the user October 20, 2026. Ancestry is tracked independently of
assigned group identity, without changing rates, matching, or configured identity
assignment. It represents expected genealogical contributions, not stochastic
realized genetic inheritance.

Milestone 3 was accepted by the user October 27, 2026 (D077). The
[validation report](../validation/milestone3.qmd) records the assessment;
acceptance retains the documented limitations.

## Initialization and public table

Ancestry is tracked only for dimensions declared with `group(ancestry = TRUE)`
(D134). `hammel_population(persons, ..., ancestry = NULL)` accepts a named list
of ancestry data frames, one per tracked dimension (a bare data frame when only
one dimension tracks ancestry). The public `ancestry` component is a named list
of tibbles, each with one row per person, including deceased people, keyed by
`person_id`. Known ancestry columns are named by category; `unknown` is a
separate share. For a dimension `race` with categories A and B,
`ancestry$race` has columns `person_id`, `A`, `B`, and `unknown`. Categories of
a tracked dimension therefore cannot be named `unknown` or `person_id`.
Labels may be numeric codes.

Without supplied proportions, every initial person receives 1.0 in their
recorded category of each tracked dimension and zero elsewhere, including zero
unknown share. This applies to all starting people, regardless of known parent
histories: do not trace existing ancestors to replace the requested starting
reference. Reference group columns follow first occurrence in the input history.
No ancestry component is created automatically without population groups.

If supplied, use the ancestry proportions as given. Require exactly one row for
each person, no extra IDs, finite numeric shares in `[0,1]`, and row sums within
`1e-8` of 1.0, including `unknown`. Do not renormalize or infer residual shares.
Explicit unknown ancestry is represented by the `unknown` component, not NA or
missing rows. Reference ancestry groups need not coincide with current identity
groups; supplied columns are retained even if no person currently has that
identity. A supplied ancestry table can also be retained without identity groups.
Rows are normalized to person-table order, known columns retain input order,
and `unknown` is last. Base data frames are accepted and returned as tibbles.
An empty grouped population has an empty `person_id`/`unknown` ancestry tibble.

## Births and restarts

For each ancestry component, a newborn receives half the mother's share plus
half the father's share. A missing parent contributes its entire half to
`unknown`. Known parents with partially or entirely unknown ancestry transmit
their unknown share by the same averaging rule. Existing engine births always
have a known mother; an allowed unmatched birth may have an unknown father.
Ancestry is appended at birth before the newborn enters risk, and is retained
after death. Group identity selects either parent with equal probability by default (maternal
fallback for unknown fathers), or follows the configured
fertility inheritance model. Assignment does not change ancestry proportions;
see [group inheritance](group-inheritance.md).

Keep supplied/returned proportions across attachment, simulation, and restarts;
do not reinitialize them from identity. Individual group transitions leave
all ancestry shares unchanged; see [group transitions](group-transitions.md). If a grouped population lacks an ancestry
component when attached or simulated, initialize it by the same starting-person
rule. No ancestry draws are made and demographic histories/RNG consumption must
remain unchanged. Each simulation owns its ancestry state, with work per birth
and storage per person proportional to the number of ancestry components.

## Worked examples and validation

- Mother A=1 and father B=1: child A=0.5, B=0.5, unknown=0, identity A.
- That child and a pure B partner: next child A=0.25, B=0.75.
- Mother A=1 and missing father: child A=0.5, unknown=0.5.
- Mother A=0.5, unknown=0.5 and father B=1: child A=0.25, B=0.5,
  unknown=0.25.
- An initial person in group A with supplied A=0.25, B=0.75 keeps those
  proportions, including after a restart. Without a supplied table, the same
  initial person starts A=1 even if the input identifies a mother in B.

Validation covers these cases, row sums/bounds, multi-generation inheritance,
missing fathers, supplied reference groups, shuffled IDs and ancestry rows,
deaths, empty populations, restarts, numeric labels, seed preservation and
independence between runs. Current ancestry is inspectable state; this increment
does not add ancestry predictors to demographic processes or analysis helpers.
