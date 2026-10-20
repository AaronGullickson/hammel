# Group-specific demographic schedules

Milestone 3 initial scope, requested October 20, 2026. Group data belongs in
the population and its effects belong in the relevant processes (D051).
This specification does not select the deferred
inheritance-model API or proportional rate multipliers.

## Lookup and timing

An optional `group` column identifies each person's fixed group. Values must be
nonmissing character labels or finite numeric codes; character labels must be
nonempty. Types are preserved in population results. Comparisons use character
representations, including when numeric codes are matched to labels. Alphabetical
selection uses radix order on those representations: codes `10` and `2` select
`10` first when the population has no groups. Mortality,
fertility, and female union-search rate tables may independently include
`group`. A table without that column applies to everyone, including a grouped
population. A grouped table selects the person's own group; union-search rates
use the woman's group. Every population group must have an explicit schedule
in each grouped table, or simulation fails. Extra rate groups are permitted.

Without a population `group` column, each grouped rate table independently
uses its alphabetically first group and simulation emits a warning identifying
the process and selected group. Selection is independent of input row order.
This does not assign group identities to the population.

Age-band overlap is checked within group and sex for mortality, within group
and union status for conditional fertility, and within group for overall
fertility and searches. Each mortality group needs both sex schedules, with
complete coverage of `[0, Inf)` at execution. Each conditional fertility group
needs both statuses. Fertility and search gaps continue to mean zero risk.
Units remain events per person-year and all times remain continuous years.
Fertility midpoint compensation is computed independently for each group/status
schedule using the original finite bands and existing birth-spacing rules.

## Births and unions

For this increment, every newborn inherits the mother's group automatically.
The child's own schedule applies from birth, including neonatal mortality.
The father need not share that group. Group membership neither restricts the
candidate pool nor changes matching utilities. Other eligibility rules remain.
Group identity is fixed over life and retained across restarts. No ancestry
proportions, configurable inheritance model, or new group summaries are added.

## Baseline initialization

`create_baseline_population(group_weights = NULL)` optionally takes a named
numeric vector. Names are unique, nonmissing, nonempty group labels; weights
are finite, nonnegative, and include at least one positive value. Normalize
weights to probabilities and independently sample a group for every person,
after generating sex and age. For example, `c(A = 2, B = 1)` gives each person
probability 2/3 of A and 1/3 of B; realized counts are random, not quotas.
Zero weights are permitted. Omission adds no group column or random draws.
An empty population with weights has a zero-length character `group` column.
Sampling uses R's RNG and is reproducible with `set.seed()`. Initialization
does not associate group with age or sex.

## Worked cases and validation

- Population groups A and B with an ungrouped mortality schedule both use it.
- Group A mortality is zero and B mortality is positive: only B can die.
- Rate rows list B before A but the population has no group column: A is used,
  with a warning, and no population group column is created.
- A mother in A and father in B may form a union; their child belongs to A.
- Fertility is zero for A and positive for B: only B mothers can give birth.
- A grouped table containing A alone cannot serve a population containing B;
  validation fails rather than assuming a rate or treating B as zero risk.
- Identical group schedules must reproduce the ungrouped demographic histories
  with the same seed; adding group labels consumes no random draws.

Validation covers these cases, group/status compensation, search selection,
newborn schedule selection, empty populations, restart persistence, edited
objects, and unchanged behavior without group columns. Existing summaries
continue to pool groups.
