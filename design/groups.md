# Group-specific demographic schedules

Milestone 3 was accepted October 27, 2026 (D077), following the
[user-authored validation report](../validation/milestone3.qmd). Group data belongs in
the population and its effects belong in the relevant processes (D051).
Milestone 9 (D134) generalized groups to any number of declared dimensions,
such as `race` and `region`; this page describes the current behavior, and
[multiple group dimensions](multiple-groups.md) records the decisions and the
stages that implemented them. Newborn identity rules are specified in
[group inheritance](group-inheritance.md).

## Declaring dimensions

Group dimensions are persons columns declared with `set_groups(model, race =
group(...), region = group(...))`, or with `hammel_population(..., groups =)`
for analysis without a model. There is no implicit `group` column: a
one-dimension model also declares it (it may be named `group`). An undeclared
column is an ordinary custom column, kept unchanged and `NA` for newborns, and
any rate table, weight model, matching model or `by` that names it is an error.
`group()` sets the dimension's optional ordered `categories` (otherwise kept
from the population's existing declaration, or, without one, inferred from the
population in first-appearance order plus transition categories), whether it
tracks `ancestry`, and its `inheritance` rule. The first category is the
reference level for treatment-coded group terms.
Dimension names exclude reserved names (`sex`, `parity`, `origin`, ...) and the
patterns `partner_*`, `same_*`, `female_*`, `male_*`, `ancestry_*`, `*_w` and
`*_h`, which name derived variables.

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
Built-in strata are `sex`, `union_type`, `previous_union`, `origin` and
`destination`; any other name may be a declared group dimension (D134), such as
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
| Group transitions | `origin`, `destination`, `sex`, `union_type`, other dimensions |

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

Rate-table summaries (`calculate_tfr()`, `calculate_life_table()`) summarize
each stratum of a stratified table like any other stratified input.

## Lookup and timing

A rate table stratifies by the declared dimensions it has columns for: the
dimension's name for person-level rates (mortality, fertility, union search,
cohabitation probability, group transitions), meaning the person's (or
searching woman's) category, and `<dim>_w`/`<dim>_h` for union-level rates
(dissolution, conversion), meaning the wife's and husband's categories.
Dimensions without a column share the schedule, so a table without group
columns applies to everyone. Each process's cells are only as fine as its
rates (D128). Every combination of the categories of the dimensions a table
uses must have a schedule, including categories that are empty at the start,
or simulation fails; group transitions instead treat uncovered strata as zero
risk. With declared `categories`, rate rows for other categories are errors;
with inferred ones they are unused. Values are compared through their
character representations, and labels keep their persons-column type in
results.

Age-band overlap is checked within each stratum: sex for mortality, union
status for conditional fertility, and the group strata of every table. Each
mortality stratum needs both sex schedules, with complete coverage of
`[0, Inf)` at execution. Each conditional fertility stratum needs both
statuses. Fertility and search gaps continue to mean zero risk. Units remain
events per person-year and all times remain continuous years. Fertility
cells are built separately for each stratum and status schedule (D101).

## Births and unions

Newborns receive a category in every declared dimension, each from its own
`group(inheritance = ...)` rule (random parent with maternal fallback by
default; parental, ancestry, hierarchy, custom and calendar-scheduled rules
are available). The child's categories apply from birth, including neonatal
mortality. Possible categories include declared ones, initial ones, positive
ancestry shares and transition origins and destinations; group-specific
schedules must cover them all before execution.

Partners need not share a category. Matching group effects (formula
predictors such as `female_race`, symmetric pair effects and NA prohibitions
in `group_exogamy`, and `delta`) may use several dimensions, whose effects add
over their combinations; see
[group matching](group-matching.md). Other eligibility rules remain.
Categories change only through the optional
[group-transition processes](group-transitions.md), one per dimension; current
categories and transition histories are retained across restarts. Ancestry
shares are tracked for dimensions declared with `ancestry = TRUE`, one tibble
each in the population's `ancestry` list; see [ancestry](ancestry.md). See
[group inheritance](group-inheritance.md) for fallbacks, passing, and custom
assignment of mixed categories such as AB.

## Baseline initialization

`create_baseline_population(groups = NULL)` optionally samples group
dimensions, independently of age and sex. `groups` is either a named list of
named weight vectors, one per dimension, drawn independently, such as
`list(race = c(A = 2, B = 1), region = c(N = 1, S = 1))`, or a data frame with
one column per dimension and a `weight` column giving the joint distribution.
Weights are finite and nonnegative with a positive value, normalized to
probabilities; realized counts are random, not quotas. The population declares
each dimension with `group()` defaults, except that its categories are those
given, in order: the weights' names (zero-weight categories included, starting
empty), or each dimension's values in order of first appearance in the joint
table's rows. The reference category therefore does not depend on which person
happens to come first. `set_groups()` can replace the declaration; a dimension
it declares without `categories` keeps these (user decision, November 6, 2026).
Omission adds no groups and consumes no random draws. Sampling uses R's RNG
and is reproducible with `set.seed()`.

## Worked cases and validation

- Populations in A and B with a mortality schedule without group columns both
  use it.
- Group A mortality is zero and B mortality is positive: only B can die.
- A rate table with a column for an undeclared persons column is an error that
  asks for the dimension to be declared.
- A mother in A and father in B may form a union; their child's category comes
  from the inheritance rule (random parent by default).
- Fertility is zero for A and positive for B: only B mothers can give birth.
- A group table containing A alone cannot serve a population containing B;
  validation fails rather than assuming a rate or treating B as zero risk.
- Rates by `race` and `region` together need a schedule for every race and
  region combination; rates by `race` alone apply to every region.
- Identical group schedules must reproduce the ungrouped demographic histories
  with the same seed; adding group labels consumes no random draws.

Validation covers these cases, several dimensions in rates, inheritance,
ancestry, transitions, analysis, weights and matching, fertility cells,
search selection, newborn schedule selection, empty populations, restart
persistence, edited objects, and unchanged behavior without groups.

## Group observations

The calculators' `by` argument (D133) names declared dimensions:
`calculate_population_size(..., by = "race")`,
`calculate_fertility_rates(..., by = c("race", "region"))`, and
`<dim>_w`/`<dim>_h` for union rates. Results get one column per dimension,
covering every category in current categories or the retained transition
history (first-appearance order), with the full grid of combinations. No
category is inferred when no dimension is declared. See
[demographic analysis](demographic-analysis.md) for the shared rules.

Lexis counts retain the full period/age/sex grid in every group, including
zero-exposure cells. Exposure is split at group transitions, deaths use the
person's group at death, and births use the mother's group at birth, even if
an input child's group differs. Existing maternal
age rules, half-open boundaries, and zero-exposure NA rates are unchanged.
Unknown-mother births stay in the single period-level `unassigned_births`
attribute and are not assigned a maternal group from the child. Combining
`by = c("race", "union_type")` returns unpartnered and all
union-type rows within each group, without an overall row. Each spouse's formation and
exposure belong to that spouse's group at the event or exposure time, including cross-group unions.

`calculate_tfr()`, `calculate_patfr()` and `calculate_life_table()`
automatically treat any column that is not a band or measure, such as a group
dimension, as an additional stratum alongside periods and union status (and
sex for life tables). Age-contiguity checks, missing-rate rules, and
life-table assumptions apply independently within each stratum. Stratum
columns are preserved, including in empty outputs; no averaging or pooling is
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


Marriage observations use `calculate_union_rates(..., by = "race")`
(D080), grouping each spouse by their own category. Add `"previous_union"` to
`by` to additionally split unpartnered exposure and union starts by previous
category; without it, counts and exposure are pooled before calculating rates.
See [remarriage](remarriage.md).


Cohabitation extends fertility input strata to union_type (D086). Type selection probabilities
stratify by the woman's categories; conversion and dissolution use directional
`<group>_w`/`<group>_h` strata (D134 stage 2).
