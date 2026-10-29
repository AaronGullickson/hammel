# Group-specific demographic schedules

Milestone 3 was accepted October 27, 2026 (D077), following the
[user-authored validation report](../validation/milestone3.qmd). Group data belongs in
the population and its effects belong in the relevant processes (D051).
Newborn identity rules are specified in [group inheritance](group-inheritance.md).

## Proportional rate multipliers

The user requested optional `group_prop_hazard` named numeric vectors for
mortality, fertility, and rate-based union formation on October 20, 2026.
Each vector name identifies a population group, using the same character
representation as group-rate lookup. Values are dimensionless, finite,
nonnegative factors; zero suppresses the affected risk and factors are not
normalized. Names must be unique, nonmissing, and nonempty. NULL means no
multipliers; a group omitted from the vector receives factor 1.0. Extra names
do not create population groups.

Each process owns its multipliers and preserves its supplied rate table and
source metadata. Construction, attachment, and simulation validate edited
vectors. Mortality factors apply to both sex schedules of a group, fertility
factors apply to both union statuses where present, and search factors use
the woman's group, across all `previous_union` categories when those are
supplied (D080). Factors may modify an existing group-specific schedule
or a shared schedule; they do not supply missing age bands or replace the
requirement for an explicit matching schedule in grouped rate tables.

The user selected multiplication of input rates before fertility birth-spacing
compensation. For a group with input ASFR 0.1 and factor 1.1, the effective ASFR
is 0.11. With spacing 0.75 in a constant-rate band, its adjusted eligible hazard
is `0.11 / (1 - 0.75 * 0.11)`, exactly as for an explicit group schedule of 0.11.
Mortality 0.02 with factor 0.9 becomes 0.018 per person-year. Search rate 0.2
with factor zero suppresses searches for that group without altering fertility
or partner matching. Unspecified groups retain their original rates.

Without a population `group` column, warn once per process and ignore its
multipliers, even for zero-duration simulation. Existing alphabetical fallback
for grouped rate tables remains independent and emits its own warning. A group
column on an empty population remains valid and requires no multiplier fallback.
On-demand formation has no search hazard to multiply; non-NULL multipliers are
rejected for that mode. The option does not change matching utility.

At execution, scale a private rate-table copy before constructing native
schedules. Reject nonfinite products rather than clipping them; effective
fertility must pass the existing midpoint-eligibility check. No new random
draws are introduced. Multipliers of one or equivalent explicit schedules must
reproduce the same histories under the same seed, including newborns and
restarts.

For conditional fertility, `group_prop_hazard` may instead be a named list of
named vectors, for example `list(partnered = c(A = 1.1, B = 0.9),
unpartnered = c(A = 1.05))`. List names must occur in the supplied `union_type` rates: either
`partnered`/`unpartnered` or `unpartnered`/`cohab`/`marriage`, without duplicates. The user selected factor 1.0 for an omitted
status, as for an omitted group within a status. Each supplied list element
must be a valid named numeric multiplier vector. A status list requires input
fertility rates with `union_type`; overall fertility uses a single vector.
Mortality and union formation accept only a single vector. A single fertility
vector still applies to every supplied type. For example, partnered A rate 0.1
becomes 0.11 and unpartnered A rate 0.02 becomes 0.021 with the list above;
unpartnered B stays unchanged. Compensation then uses each effective
group/status schedule independently. Without population groups, ignore the
entire vector/list and warn once for fertility.

Rate-table summaries (`calculate_tfr()` and `calculate_life_table()`,
including process inputs) continue to summarize the stored input rate table;
they do not expand population-dependent multiplier settings. Grouped simulation
summaries measure the resulting histories normally.

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
schedule using original bands and existing birth-spacing rules: finite
bands use midpoints, while upper-open bands use the interior limit (D095).

## Births and unions

Newborns select either parent's group with equal probability by default,
falling back to the mother when the father is unknown. Fertility may attach a
`create_group_inheritance()` object for other parental, ancestry, hierarchy,
or custom assignment rules. The child's assigned schedule applies from birth,
including neonatal mortality. Possible identities include initial identities,
positive reference ancestry groups, and explicitly declared outputs; grouped
schedules and categorical matching must cover the full set before execution.
The father need not share that group by default. Optional `prevent_exogamy`
restrictions make men of prohibited groups ineligible. Group formula terms and symmetric pair effects
can change matching utilities; see [group matching](group-matching.md).
Other eligibility rules remain.
Group identity is fixed over life and retained across restarts. Ancestry
proportions are tracked separately in the population's `ancestry` tibble;
see [ancestry](ancestry.md) for initialization, averaging, and unknown shares.
See [group inheritance](group-inheritance.md) for fallbacks, passing, and custom
assignment of mixed identities such as AB.

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
objects, and unchanged behavior without group columns.

## Group observations

`calculate_demog_rates(..., by_group = TRUE)` and
`calculate_population_size(..., by_group = TRUE)` return separate counts for
each group appearing anywhere in the supplied person history. Both default to
FALSE, preserving pooled output. TRUE requires a population/persons `group`
column; no group identity is inferred when that column is absent. Labels retain
their input type and groups appear in first-occurrence order. No extra pooled
group rows or groups absent from the history are added.

Lexis counts retain the full period/age/sex grid in every group, including
zero-exposure cells. Exposure and deaths use the person's group; births use
the mother's group, even if an input child's group differs. Existing maternal
age rules, half-open boundaries, and zero-exposure NA rates are unchanged.
Unknown-mother births stay in the single period-level `unassigned_births`
attribute and are not assigned a maternal group from the child. Combining
`by_group = TRUE` and `by_union_type = TRUE` returns unpartnered and discovered
union-type rows within each group, without an overall row. Each spouse's formation and
exposure belong to that spouse's group, including cross-group unions.

`calculate_tfr()` and `calculate_life_table()` automatically treat an input
`group` column as an additional stratum alongside periods and union status
(and sex for life tables). Age-contiguity checks, missing-rate rules, and
life-table assumptions apply independently within each stratum. The group
column is preserved, including in empty outputs; no averaging or pooling is
performed across groups.

Population sizes retain requested time order and repeated times, with groups
in first-occurrence order within each time. Each existing group gets zero
counts where no members are alive. Empty histories with group columns produce
typed empty group summaries; empty time vectors produce zero rows.

Worked observation case: woman A and man B form a union at 0.5 in a window
`[0, 1)`. If both live throughout it, each contributes one person-year in
their own group and 0.5 partnered person-years. Each contributes one formation
to their group's unpartnered row. Their birth at 0.75 counts for maternal
group A, regardless of the child's recorded group. Summing group exposure,
deaths, and maternal births recovers pooled cells; group population counts
sum to pooled counts at every requested instant. Tests use these accounting
identities, boundary examples, analytical group TFR/life expectancy, and empty
outputs, as well as unchanged default results.


Marriage observations now use `calculate_union_rates(..., by_group = TRUE)`
(D080), grouping each spouse by their own identity. Set `by_previous_union = TRUE`
to additionally split unpartnered exposure and union starts by previous category;
the default FALSE pools counts and exposure before calculating rates. Search-rate multipliers apply
across previous categories within each woman's group. See
[remarriage](remarriage.md).


Cohabitation extends fertility input strata to union_type (D086). Lists of
`group_prop_hazard` use names occurring in the rate table: unpartnered/partnered
or unpartnered/cohab/marriage. Partial lists are valid and omitted entries use
1; keys for an absent type are errors. A single group vector still applies
across every type. Type selection probabilities and conversion schedules can
stratify by the woman's group; dissolution retains symmetric group_pair strata.
