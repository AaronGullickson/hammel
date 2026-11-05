# Group-specific demographic schedules

Milestone 3 was accepted October 27, 2026 (D077), following the
[user-authored validation report](../validation/milestone3.qmd). Group data belongs in
the population and its effects belong in the relevant processes (D051).
Newborn identity rules are specified in [group inheritance](group-inheritance.md).

## Proportional stratification with `stratify_rates()` (D132)

`stratify_rates(rates, by)` adds proportional-hazard strata to a rate table
before it is attached to a process. It replaced the per-process
`group_prop_hazard` argument on November 5, 2026 (D132). `by` is a named list:
each name is a column to add, and each element a named numeric vector of
factors, one per category:

```r
mortality_rates |> stratify_rates(list(group = c(A = 1.2, B = 1)))
```

returns two copies of the table, one with `group = "A"` and rates multiplied by
1.2, and one with `group = "B"` and rates unchanged. B is then the reference.

**Rules.**

- Every row of the input is copied once per category, with the new column set
  to the category label (character) and `rate` multiplied by its factor. All
  other columns, including calendar bands, are copied unchanged.
- Several columns in `by` are crossed: each combination of categories gets the
  product of its factors. `by = list(group = c(A = 1.2, B = 1), sex = c(f = 1, m = 1.1))`
  gives four copies, with A men at 1.32 times the input.
- Factors are finite and nonnegative; zero suppresses that stratum's risk.
  They are not normalized. Category names are unique and nonempty, and each
  column has at least one category. Products that overflow are errors.
- A column that already exists in `rates` is an error, as is `rate` itself.
  To give different effects within existing strata, for example different
  group effects for partnered and unpartnered fertility, stratify each subset
  separately and bind the results; any interaction can be built this way.
- `rates` must be a data frame with a `rate` column. The result is a tibble that
  keeps the input's other attributes, such as database provenance. Empty
  tables gain the new columns with no rows.
- The result is an ordinary stratified rate table. Every category the process
  needs must appear: a group omitted from `by` is a coverage error at
  simulation, as for any group-specific table. Factors of one still create
  separate strata, so cells split by those categories (D128) and a weight-model
  term on the same variable cancels with a warning.

`by` may add any valid column name except band, rate and calendar columns.
Built-in strata are `sex`, `union_type`, `previous_union`, `group_origin` and
`group_dest`; any other name may be a declared group dimension (D134), such as
`race`, or for union-level processes its wife's and husband's columns,
`race_w` and `race_h`. `stratify_rates()` cannot know which process or
declaration a table is for, so each process setter checks the table's built-in
stratum columns, whether they came from `stratify_rates()`, an earlier
stratification or direct entry, and rejects those its process does not
support. Group columns are checked at simulation against the declaration;
there an undeclared column is an error. Other columns, such as notes, are
carried along:

| Process | Stratum columns accepted from `stratify_rates()` |
| --- | --- |
| Mortality | `sex`, `<group>` |
| Fertility | `<group>`, `union_type` |
| Union search (rate mode) | `<group>` (woman's), `previous_union` |
| Conversion (`union_convert_rate`) | `<group>_w`, `<group>_h` |
| Dissolution | `<group>_w`, `<group>_h`, `union_type` |
| Group transitions | `group_origin`, `group_dest`, `sex`, `union_type` |

Category values are then validated as for any supplied table: `union_type`
values must be valid for the process, and group values must be categories of
the dimension (rows for categories absent from an inferred registry are
unused).

**Examples.** Group A fertility 10% higher than B: an ASFR of 0.1 becomes 0.11
for A, and A's cells generate births at 0.11 times their members, exactly as an
explicit schedule would. Directional effects multiply: `stratify_rates(rates,
list(group_w = c(A = 1, B = 1.5), group_h = c(A = 1, B = 2)))`. Effects that
are not products of partner effects, such as doubling dissolution in exogamous
unions, are written as an explicit table with one row per wife-husband
combination.
Different group effects by union status:

```r
rbind(
  stratify_rates(fert[fert$union_type == "partnered", ], list(group = c(A = 1.1, B = 1))),
  stratify_rates(fert[fert$union_type == "unpartnered", ], list(group = c(A = 1.05, B = 1))))
```

**Populations without groups.** A table stratified by `group` follows the
existing fallback for group-specific rates: simulation warns and uses the
alphabetically first group's schedule, including its factor.

Rate-table summaries (`calculate_tfr()`, `calculate_life_table()`) summarize
each stratum of a stratified table like any other stratified input.

## Lookup and timing

An optional `group` column identifies each person's current group. Values must be
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
Fertility cells are built separately for each group/status schedule (D101).

## Births and unions

Newborns select either parent's group with equal probability by default,
falling back to the mother when the father is unknown. Fertility may attach a
`create_group_inheritance()` object for other parental, ancestry, hierarchy,
or custom assignment rules. The child's assigned schedule applies from birth,
including neonatal mortality. Possible identities include initial identities,
positive reference ancestry groups, transition origins/destinations, and explicitly declared outputs; grouped
schedules and categorical matching must cover the full set before execution.
The father need not share that group by default. Optional `prevent_exogamy`
restrictions make men of prohibited groups ineligible. Group formula terms and symmetric pair effects
can change matching utilities; see [group matching](group-matching.md).
Other eligibility rules remain.
Group identity changes only through the optional [group-transition process](group-transitions.md); current identity and transition history are retained across restarts. Ancestry
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

Validation covers these cases, group/status fertility cells, search selection,
newborn schedule selection, empty populations, restart persistence, edited
objects, and unchanged behavior without group columns.

## Group observations

`calculate_fertility_rates(..., by = "group")` and
`calculate_population_size(..., by = "group")` return separate counts for
each group appearing in current person identities or the retained transition history. Both default to
FALSE, preserving pooled output. TRUE requires a population/persons `group`
column; no group identity is inferred when that column is absent. Labels retain
their input type; current groups appear in first-occurrence order, followed by any historical groups absent from current identities. No extra pooled
group rows or groups absent from the history are added.

Lexis counts retain the full period/age/sex grid in every group, including
zero-exposure cells. Exposure is split at group transitions, deaths use the
person's group at death, and births use the mother's group at birth, even if
an input child's group differs. Existing maternal
age rules, half-open boundaries, and zero-exposure NA rates are unchanged.
Unknown-mother births stay in the single period-level `unassigned_births`
attribute and are not assigned a maternal group from the child. Combining
`by = c("group", "union_type")` returns unpartnered and all
union-type rows within each group, without an overall row. Each spouse's formation and
exposure belong to that spouse's group at the event or exposure time, including cross-group unions.

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


Marriage observations now use `calculate_union_rates(..., by = "group")`
(D080), grouping each spouse by their own identity. Set `by = "previous_union"`
to additionally split unpartnered exposure and union starts by previous category;
the default FALSE pools counts and exposure before calculating rates. See
[remarriage](remarriage.md).


Cohabitation extends fertility input strata to union_type (D086). Type selection probabilities and conversion schedules can
stratify by the woman's group; conversion and dissolution use directional
`<group>_w`/`<group>_h` strata (D134 stage 2).
