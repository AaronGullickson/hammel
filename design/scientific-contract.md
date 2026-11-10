# Scientific contract

Updated November 5, 2026. This is the current mortality, fertility, union-formation,
dissolution, remarriage, cohabitation/conversion, group identity/ancestry, and
observation contract. Milestones 1–8 are accepted as complete
(D029, D050, D077, D093, D105, D116, D127, D130), retaining documented limitations.
Milestone 8 schedules every rate process by cells with one clock per process
(D128) and adds weight models that choose who in a cell experiences an event
(D129); see [cell scheduling](#cell-scheduling-of-rate-processes-d128) and
[weight models](#weight-models-d129).
Individual group transitions are accepted under the
[transition contract](group-transitions.md). Milestone 7, time-varying inputs,
has implemented calendar-varying rates under the user-approved
[rate contract](time-varying-rates.md) (D118), accepted complete November 2,
2026 (D127). Matching and inheritance model schedules are implemented under
[their contract](time-varying-models.md) (D123).
Earlier monthly, exact
age-varying and midpoint compensation specifications are available in Git
history; the [decision log](decisions.md) records their replacement by
cell-accounted fertility (D101). Acceptance
of a model assumption is distinct from evidence of demographic accuracy.

## Units and inputs

Public dates, ages, and durations are numeric years; rates are events per
person-year and exposure is person-years. The user supplies the time origin.
Every population requires one finite numeric `end_time`, which supplies the
next simulation's starting time (D166).
Execution is continuous: there is no month index, configurable execution step,
time truncation, or rounding to birthdays.

Mortality, fertility and union-search rate tables contain `age_start`,
`age_end`, and `rate`. Previously partnered union-search categories can instead
use `dur_start`/`dur_end` with age boundaries NA (D083), as described below. Dissolution uses `dur_start`/`dur_end` for years since
union formation, with rates per partnered woman-year (equivalently union-year).
Mortality also contains
`sex`, with separate `"m"` and `"f"` schedules. Intervals are half-open
`[age_start, age_end)`. Starts must be finite and nonnegative; ends must exceed
starts. Rates must be finite and nonnegative, and may exceed one because they
are rates rather than probabilities. Overlapping or duplicate risk bands are rejected within each active calendar
slice. Optional `time_start`/`time_end` select half-open absolute-year calendar
bands on all six annual-rate inputs (D118). Starts are finite, ends may be Inf,
and fractional years are supported. Calendar changes precede same-time events;
coverage is required during a run. See [calendar-varying rates](time-varying-rates.md). Constructors preserve
input order and extra columns; native scheduling sorts its own copy.

| Rule | Mortality | Fertility |
| --- | --- | --- |
| Input interpretation | Piecewise constant age/sex hazard | Overall or status-specific ASFR, births per woman-year in the specified population |
| Denominator for output rates | Living person-years by age/sex | Female person-years overall or in the indicated status, including spacing-excluded time |
| Open-ended input band | Allowed; required at execution | Allowed (D095) |
| Coverage at execution | `[0, Inf)` for each sex, with no gaps | Zero outside supplied finite bands, including gaps |
| Additional setting | None | `birth_interval`, finite nonnegative scalar years; default 0.75 |

Mortality constructors validate structure without requiring full coverage;
`run_simulation()` checks complete schedules. Specifying and attaching processes does
not execute a simulation. Edited built-in specifications are revalidated at
attachment and execution. Fertility schedules need no feasibility check under
cell accounting; parity schedules with unreachable cells receive a warning.

Treating empirical mortality central rates `mx` as hazards is an explicit
assumption. It is not a universal identity, particularly for infant mortality.
For a constant hazard `h` and elapsed years `d`, the event probability is
`1 - exp(-h*d)`; multiplying a rate by duration gives cumulative hazard,
not a probability. There is no requirement to convert to monthly units.

## Group-specific schedules

Populations may declare any number of group dimensions (D134), such as `race`
and `region`. Rate tables may stratify by any of them (union rates by the
wife's and husband's categories, `race_w`, `race_h`); tables without a
dimension's column apply to all its categories, and tables with it must cover
every category combination they use. A table naming an undeclared column is
an error. Group does not affect matching eligibility unless restrictions are
specified. Optional group formula predictors and symmetric pair effects in one
named dimension can change utility. Newborns receive each dimension's category
from its inheritance rule, by default either parent's with equal probability
(with maternal fallback for unknown fathers), entering the assigned
categories' mortality schedules immediately. Rate cells split by group only when a process's effective
schedules are group-specific (D128). Group values are
nonmissing character labels or finite numeric codes, with nonempty labels.
See [groups](groups.md) for full rules, ordering, examples, and validation scope.

Proportional differences between strata are written into the rate table with
`stratify_rates(rates, by)` (D132), which copies each row per category of a new
column and multiplies its rate by that category's factor, crossing several
columns. Factors are dimensionless, finite and nonnegative and are not
normalized; zero suppresses that risk. The result is an ordinary stratified
table, so every needed category must be listed, and fertility multiplication
precedes cell accounting exactly as an explicit schedule would. Setters reject
added columns that a process does not support. Fertility cells split by group
only when the schedules are group-specific (D128). See
[groups](groups.md#proportional-stratification-with-stratify_rates-d132).

## Initialization and population history

The [population schema](population-schema.md) specifies required fields and
structural validation. `end_time` is the population reference time and supplies
the next simulation start. Persons born by that time with missing `death_time`
are alive; deceased persons remain in the history and never re-enter risk.

`generate_uniform_population(size, max_age = 70, end_time = 0,
female_probability = 0.5)` generates independent uniform ages in continuous
years and independent sex assignments. Everyone is alive, parents are unknown,
and unions are empty. IDs are `1:size`; size zero is permitted. Size must be a
nonnegative integer within R's integer range; maximum age is finite and
nonnegative; the sex probability is in `[0,1]`. The generator advances R's
current RNG state; use `set.seed()` for reproducibility.

This baseline is not claimed to be historically representative, stable, or
stationary. A mother with no recorded birth starts eligible for fertility;
that does not establish that she truly had no previous children. For supplied
histories and restarts, recover the latest recorded maternal birth, including
children who have died. No initial histories are fabricated.

### Ancestry proportions

Optional ancestry, for dimensions declared with `ancestry = TRUE`, stores one
tibble per such dimension with one row per person, a share per category and an
explicit `unknown` share. Without supplied proportions, every starting person
begins at 100% of their recorded category,
irrespective of earlier parent histories. Supplied proportions are used as given
and retained across restarts. Shares must be finite, in `[0,1]`, and sum to one
within `1e-8`; missing ancestry is represented explicitly as unknown, not NA.
No table is created automatically without groups, though supplied ancestry may
be retained without current group identities.

At birth, average each parent's shares with weights 0.5 and 0.5. A missing
parent contributes its entire half to unknown; known parents' unknown shares
are averaged normally. This is deterministic expected genealogical ancestry,
independent of assigned group identity. All shares are retained after death;
the engine appends them before scheduling the newborn without extra RNG draws.
See [ancestry](ancestry.md) for the public schema, examples, and validation.

## Cell scheduling of rate processes (D128)

Every rate process (mortality, fertility, rate-mode union search, dissolution,
cohabitation-to-marriage conversion, and each group-transition slot) keeps
cells of at-risk units. A unit is a person, or for dissolution and conversion a
union. A cell is one band of a stratum's schedule on the unit's clock (age, or
duration since union start, last union end or group entry). Cells stratify
only as finely as the input schedules: by group only when a process's
effective schedules are group-specific, and by unordered group pair for
dissolution only when its rates are. Fertility splits its bands into cells of
at most one year and by actual parity, because spacing is an eligibility rule
(D101); the other processes use their input bands as cells, since without
eligibility rules every member's hazard is the cell rate. Zero-rate bands are
gaps with no cells.

A cell produces events at its input rate times its members. Each process keeps
one clock: an Exp(1) threshold from R's generator, spent at the process's total
intensity (the sum over cells of rate times members). Any change in membership
or rates first spends the threshold for the elapsed interval and then re-times
the clock, with no new random draw. When the clock fires, a uniform draw picks
a cell in proportion to its intensity, a second picks the unit within it
(uniformly, or in proportion to weights), and a new threshold is drawn. Units
move between cells at band edges through bookkeeping boundaries, which consume
no randomness and are processed before events at the same time. If nobody can
reach a positive band, no event occurs; no maximum age forces death.

Worked example: a mortality process has two cells, rate .01 with 300 members
and rate .05 with 100. Its total intensity is 3 + 5 = 8 deaths per year, so with
threshold .4 the next death is due in .05 years. If one member of the second
cell dies first by another route, the elapsed .02 years spend .16 of the
threshold and the remaining .24 is re-timed at the new total 7.95.

At a calendar rate change, a process whose new schedules keep its cells (same
strata and nonzero band edges; for fertility also the same union types, group
stratification and age splitting; for transitions the same slots) changes its
cell rates in place and its clock continues. Otherwise that part is rebuilt and
everyone re-enters with new thresholds, valid after survival to the boundary.

## Fertility and birth spacing

Without `union_type`, all living females face the supplied overall fertility
schedule. With a `union_type` column, require both `partnered` and
`unpartnered` schedules, interpreted per woman-year in that status. Alternatively supply all three labels unpartnered/cohab/marriage under D086;
partial or mixed tables are errors. Check
age overlap separately within each status; boundaries may differ and gaps are
zero risk. Use explicit zeros to suppress births for a status. On-demand
formation rejects status-specific fertility schedules, including at zero duration.
Rate formation accepts overall rates with one warning per simulation call. After a recorded birth at time `t`, a woman is
ineligible until `t + birth_interval`. Mortality remains active throughout spacing.
A birth at 1.03 with spacing 0.75 permits the next birth no earlier than 1.78.
Zero spacing makes every woman eligible at all times. Arbitrary finite
nonnegative intervals are allowed; there is no whole-month restriction.

### Newborn group identity

Each declared dimension holds a `group_inheritance()` object through
`group(inheritance = ...)`, defaulting to `random_parent`.
The `"fixed"` method (D142) assigns every newborn the explicitly specified
`fixed_value` category, regardless of parents or sex, with no identity RNG draw.
The category must belong to the registry; existing people are not reassigned.
Other built-in methods follow mother, father, same/opposite-sex parent,
known ancestry shares, or a strict hierarchy of parental identities. Required
missing inputs fall back to the mother. Ancestry sampling conditions out unknown
shares; all-unknown ancestry falls back to the mother. Hierarchy `passing` is
an independent probability of reversed order, default 0. Missing ranks or
father use the maternal fallback. Custom functions receive copied parental
identities/ancestry, child ancestry, and child sex and must return one valid
identity. Invalid returns and undeclared outputs error.

Infer possible groups from starting identities, positive ancestry reference
shares, and optional explicit declarations. Grouped rates and categorical
matching levels must cover the full set before execution. Assignment follows
sex and ancestry calculation, consumes draws only when the chosen rule requires
randomness, and finishes before scheduling the newborn. Custom R draws share
the engine RNG stream and explicit seeds restore state even on errors. Identity
does not alter ancestry: a custom AB identity may retain only A/B ancestry.
See [group inheritance](group-inheritance.md) for rules and worked examples.

### Cell-accounted fertility (D101)

Births are generated by accounting of the simulated population. Cells come
from the effective schedules: union type or overall, and group (after
multipliers) when those schedules are group-specific. Input age bands are split at whole ages and band edges into sub-cells
of at most one year; upper-open bands are split up to the oldest attainable age
and then left open, except that an upper-open band stays one cell when
`birth_interval * rate < 0.001` (a numerical tolerance). Zero-rate bands form no
cells. Parity-specific schedules have one cell per actual parity, so a pooled
band applies its rate to each contained parity; without parity columns cells
pool all parities.

A cell generates births at its input rate times all its living members,
including women inside the birth interval, and assigns each birth to a
chosen eligible member, uniformly by default or proportional to persistent
`hetfert` weights under D102. The equivalent individual hazard is
`rate * members * individual_weight / eligible_weight_sum`. Observed rates equal inputs in expectation
for any composition (initialization without histories, restarts with new
schedules or intervals, union-type transitions), at cell resolution, while a
cell has an eligible member. Within a one-year cell, selection does not vary
with exact age. A cell with members but none eligible keeps its intensity; when
the birth clock draws it, no child is born and one unrealized birth is counted
in the population's `unrealized_births` (a tibble of counts by one-year period
from the run start, group, union type, parity band and age band; D128 revises
D101's expected-value reporting). The record accumulates across restarts from
`unrealized_since` (D160). In a cell shared by groups, the draw that would pick
a mother picks the group in proportion to its members.
`calculate_fertility_rates(unrealized = TRUE)` attributes the counts to
observation rows. Cancelled attempts keep their spent intensity. Union changes move a woman between cells; restarts rebuild cells
from histories. A woman's hazard depends on recent births by others in her cell;
this is a calibration device rather than a behavioral model. Births come from
the fertility process's clock (D128). At exact ties a death precedes a birth,
which precedes the other processes. Individual weights set selection
probability within a cell, not the cell's rate. See [cell-accounted fertility](fertility-cell-accounting.md).

For example, a cell with input rate .2 has 100 members, 15 of them within the
birth interval. It produces 20 births per year, so each of the 85 eligible
women has hazard 20/85 = .235. If a restart leaves only 5 women in spacing, the
eligible hazard becomes 20/95 = .211.

For parity-specific schedules with positive spacing, `fertility()` warns about
positive rates in one-year age cells that no eligible woman can occupy. Parity
k is first eligible after k singleton births at least `birth_interval` apart,
beginning at the schedule's first positive rates. Union types are pooled;
partly reachable cells and higher parities within an open parity band are not
reported.

### Persistent fertility weights (D102)

`set_fertility(..., hetfert = het_fert(cv, alpha))` enables positive persistent
lognormal weights, with fresh mean one and requested CV. If persons lacks
`hetfert`, simulation draws values for every recorded woman under its seed;
otherwise supplied female values must be finite and strictly positive and are
retained. Male values are unused; newly generated male entries are NA. Output
always contains numeric `hetfert` when enabled, including zero-year runs.

Daughters inherit on the log scale with coefficient alpha and independent
normal innovation scaled by sqrt(1-alpha^2). Alpha zero gives independent
weights, alpha one exact maternal copies. Reproducing mothers are selected, so
marginal cohort variance is not guaranteed to remain constant. Existing weights
persist through restarts; changing CV does not rescale them. `hetfert = NULL`
retains uniform selection and removes any existing weight column from output,
including zero-year runs (D104). Re-enabling from a disabled output draws fresh
weights for all recorded women. Input populations are not mutated.
The [heterogeneity specification](fertility-heterogeneity-notes.md) defines
zero CV, missing/invalid supplied weights, RNG handling, weighted calibration
and remaining scope. Cell targets, spacing and failed-search rules are unchanged.

### Birth transition and parent assignment

Union histories use `wife_id` for the female partner and `husband_id` for the
male partner. These roles are validated for ended and ongoing unions and do
not imply legal marriage.

Each birth creates one child with a new unique ID above all supplied IDs,
records the mother, and assigns female sex with probability `ffab`, configured
on `fertility()` / `set_fertility()` (default 0.4886). `ffab` is one finite,
dimensionless numeric probability in `[0, 1]`, shared by all births in that
fertility process. Sex is drawn at each actual birth before identity inheritance
and newborn event scheduling. Each birth retains one uniform sex draw, including
at probabilities zero and one; cancelled attempts do not draw child sex.
Changing `ffab` does not change existing people's sex or adjust input fertility
rates. A restart uses the newly attached fertility setting for subsequent births.
Older fertility specifications without this field retain the default.

For example, `set_fertility(model, rates, ffab = 0)` produces only male
newborns, `ffab = 1` only female newborns, and `ffab = 0.5` gives equal
probabilities per birth, not an exactly balanced realized population. The
setting is recorded in the latest-run simulation metadata. Validation covers
invalid and edited specifications, endpoint outcomes in both initiation modes,
sex-dependent inheritance, RNG isolation, and unchanged seeded default histories.
A union-formation process is required for every simulation, including zero
duration. Missing processes are rejected rather than given an implicit default.
The former milestone 1 path without unions is no longer executable.

When matching is omitted, both initiation modes use
`age_gap_matching_model(0, 5)`: preferred male-minus-female gap zero
years, with relative weight halved five years either side. Explicit models
override that default; it does not fix the realized age-gap distribution.

Matching models optionally accept `group_exogamy`: a named vector using `A:B`
for equal effects on female A/male B and female B/male A (D075). Effects
are finite log-scale additions after the full formula score. An omitted pair
has odds ratio 1.0 (zero log effect). Within-group effects are allowed; overlaps
after symmetric expansion are errors. A pair table is also accepted
for labels containing colons; each row expands in both directions. These are selection preferences, not hazards or
prohibitions. For a declared dimension such as `race`, `female_race`,
`male_race`, and numeric 0/1 `same_race` are available in formulas;
categorical terms use the dimension's categories with treatment contrasts. Missing population identity data errors when a model
uses group effects. See [group matching](group-matching.md) for worked cases.

Matching models also prohibit pairs through NA in `group_exogamy` (D135): `NA`
for a dimension allows only same-category candidates, and `"A:B" = NA`
prohibits that unordered pair in both directions. Restricted men are ineligible, alongside the
age and kinship exclusions below. A search fails only when no eligible man
exists in the whole pool; existing retry/cancellation rules then apply.
Restrictions override finite scores, apply to the selected father or partner
model, and do not invalidate existing partners. Group restrictions require
known population identities; exact group exclusions do not use NA or infinite
scores.

### Exact full-pool partner draws (D072)

There is no candidate sample. With `union_formation(initiation = "on_demand", ...)`,
use the woman's ongoing living male partner as father. If she has no partner,
draw a man from all eligible living, unpartnered men with probability
proportional to `exp(V)`, where V is his utility: the matching formula on
male-minus-female continuous age differences, multiplied by its named
coefficients, plus symmetric pair effects, plus the calibrated log weight for his group
(specified below). Form the union at the birth time
and assign the partner as father. Nonfinite formula values/utilities are errors.

The draw is exact. Unpartnered men are kept by group and quarter-year birth
bin. A bin is chosen with probability proportional to its size times the
largest value of `exp(V)` any member could have for this woman; a member is
chosen uniformly and accepted with probability `exp(V - that maximum)`; only
that man is then screened for age and kinship, and a rejected man leads to a
fresh draw. Rejection sampling makes the accepted man's distribution the
conditional logit over the eligible pool. After 200 rejections, the eligible
pool is enumerated and drawn from directly. For the native linear/raw-quadratic
age model the within-bin maximum is exact. Custom formulas are bounded by
their values on a 1/32-year gap grid; if a drawn man's utility exceeds the
bound, the draw is counted as approximate and simulation warns once. Very
steep utilities reduce acceptance and may reach the enumeration fallback;
outcomes are unchanged but random-number use differs.

Exclude men younger than `min_male_age`
(default 15 continuous years; finite nonnegative scalar, fractional values
allowed). A man at exactly the threshold is eligible. Also exclude a man if he and the woman share a recorded
ancestor within `kinship_generations` parent steps of EACH (default 2), counting
each person themself at distance zero. This includes known parents, grandparents,
siblings, aunts/uncles and first cousins, including half kin; it does not impose
an unlimited direct-ancestor exclusion. Zero disables cross-person kin exclusion.
Only demonstrated links count: unknown parents do not imply unrelatedness and
do not hide other known relationships. Age and kin search restrictions do not invalidate existing unions,
and no extra minimum female age is imposed. Each person can have at most one ongoing union. Formation removes the man
from the pool; widowhood returns a surviving husband.
If no eligible man exists, `cancel_birth_if_no_match = TRUE` (the default)
cancels the birth. FALSE explicitly preserves it with unknown father and no
union. A cancelled attempt creates no person, parentage, ID or birth count and
does not update the last actual birth or start a spacing interval. Redraw risks
from the attempt time at attained age; retain any spacing from actual prior
births. This applies to on-demand partner selection and rate-mode nonmarital
father selection. It does not rewrite input histories with unknown fathers.
Partner scarcity can depress realized fertility below supplied rates; do not
compensate away this intentional selection effect. Rate-based initiation follows the separate rules below.

On-demand input must provide unambiguous active female/male unions with living
partners and at most one ongoing union per person; unsupported or inconsistent
active histories are rejected. Supplied ended histories remain unchanged.
Death closes all affected ongoing unions and clears both partners' current
partner immediately. Her next birth may trigger another search in on-demand mode. In that mode fertility and
mortality remain independent of union status; union changes refresh the
partners' cell memberships without random draws. Optional
`diagnostics = TRUE` records only weighted choice's applied log group weights
(`.log_weight_<group>`) on newly formed unions in grouped populations; it
defaults off and consumes no additional random draws. See [union formation](union-formation-models.md)
for worked examples and implemented milestone 2 rules.

Update the mother's latest birth and refresh the newborn's and mother's cell
memberships. Newborns enter mortality risk immediately; their fertility depends on the supplied ages.
Existing custom person/union columns are retained; new rows have missing custom
values, except that `group` is assigned by the fertility inheritance model. C++ maintains union histories and shared partner caches for all runs,
closing unions at death with `end_reason = "death_husband"` or
`"death_wife"` according to the executed death event. Initiation mode
controls search timing; it does not control whether partner state is tracked.

### Rate-based union initiation

With `initiation = "rate"`, only living unpartnered women face the supplied
age-specific search hazard. Men never initiate searches. Gaps and ages outside
search bands have zero hazard. Searches come from the search process's cells
and clock (D128). Searches compete with death; death wins exact ties. Search execution uses the same exact draw,
screening, matching and formation code as on-demand mode. A failed search
leaves the woman in her search cell.
Without `previous_union` stratification, all unpartnered women use that schedule
regardless of prior union history. Stratified schedules distinguish none,
separated and widowed, with age or permitted duration clocks (D080–D083).

Women face birth risk while unpartnered or partnered: the current `union_type` selects
fertility, or an overall schedule applies equally to both. With overall rates,
warn once per simulation call that fertility does not depend on union status.
Unpartnered births in rate mode select a father from living, currently
unpartnered men by the same exact draw, using the father matching model and its
weights, with the same minimum age and known-kin screening as union searches. This creates parentage only:
no union forms, the man remains available, and neither person's partner state
changes. No eligible man cancels the birth by default; explicitly
setting `cancel_birth_if_no_match = FALSE` retains it with unknown father. Paternal links are appended to the kinship graph before
subsequent events. Partnered births use the current husband. No special
preference for an earlier child's father is implemented; it is deferred to a
later milestone (D042).
An optional `father_matching` specification applies only to unpartnered
rate-mode births; NULL uses the current general `matching` model. Union searches
always use the general model. Both follow the same formula contract. Birth and search candidates
compete while unpartnered; death wins exact ties, followed by birth then search.

A new union or husband's death moves the woman to the cells of her new state,
preserving her recorded birth spacing. Widowhood also restores search risk.
Surviving men's death hazards do not depend on union status. A woman who changes type moves to the new
type's fertility cell with her eligibility status, so that cell's births cover
her spacing-excluded exposure (D101). The milestone 2 report's conditional
ASFR discrepancies came from the earlier midpoint compensation.
See [union fertility specification](union-formation-models.md#fertility-schedules-by-union-status-d039).

### Duration-dependent union dissolution (D078)

An optional `union_dissolution` process adds one annual hazard per union (the
union is the unit at risk, D128), indexed by elapsed years since `start_time`.
Shared rates or schedules/multipliers by unordered group pair (D079) use
half-open duration bands;
gaps and uncovered durations have zero hazard. New episodes start at duration
zero; active input unions and restarts use attained duration. Births do not
change a union's dissolution cell. Exact same-person ties
select death, birth, dissolution, conversion, search, then group transition. A tied birth occurs before
closure and retains that husband as father.

Execution ends the episode with reason `dissolution`, clears partner links,
restores male availability and female calibration counts with spacing preserved,
and moves her to her unpartnered fertility and search cells. Men's mortality
cells do not depend on union status.
Either partner's death cancels the episode's risk. Without this process,
existing death-only behavior and random-number consumption are retained.
Search schedules may distinguish `previous_union` as `none`,
`separated`, or `widowed` (D080). The most recent recorded episode determines
state; complete histories are required. Without that column, all categories
share the same rates. Under D126, search-rate and cohabitation-probability
inputs supply either the complete detailed set `none`/`separated`/`widowed`,
or pooled `none`/`dissolved`, where dissolved includes separation and widowhood.
Each group/calendar slice requires one complete alternative.
`none` must use age; `separated` and `widowed` each
select age or duration since the latest union ending (D083). Each category uses
one clock in all bands/groups. Duration columns are optional for all-age tables;
when included, each row supplies one complete boundary pair with the unused
pair NA. No simultaneous age/duration hazards are permitted. Recover the latest
ending origin on restart and reset it at each subsequent ending. Search
cells use the selected clock; matching calibration takes their intensities by
cohort bin (D128). Cohabitation is implemented under D086; richer couple covariates remain deferred. See [remarriage](remarriage.md) and the
[dissolution specification](union-dissolution.md) for worked examples and limits.

### Optional cohabitation (D085/D086)

Rate formation accepts cohab_prob age/group/previous-category probabilities,
with age gaps zero. Optional `time_start`/`time_end` select half-open calendar
bands under D122, with coverage/overlap validation as for rates. At each attempt,
the current calendar/age/group/previous-category probability selects intended
type before matching;
cohab_matching optionally overrides general preferences. Unpartnered father
matching remains separate. Calibration uses lambda*p for cohabitation,
lambda*(1-p) for direct marriage, and fertility intensity for fathers.
A probability tracker updates each woman's probability at its age boundaries. Calendar
probability changes refresh role contributions before same-time searches;
probability-only changes retain all pending event dates and existing union types.
See [the probability contract](cohabitation.md#calendar-varying-type-probabilities-d122).
NULL probability
keeps the untyped/no-extra-RNG path; on-demand rejects cohabitation inputs.

Fertility input uses union_type: either unpartnered/partnered or
unpartnered/cohab/marriage per group; old input union_status errors. Fertility
multiplier-list keys must exist in those schedules, with omitted factors 1.
Dissolution rates may supply both union types within each pair or share rates.
Conversion uses a female age/cohab-duration clock, optional group, no previous
axis; NULL warns and disables conversion. Births do not change conversion cells. Ties
select death, birth, dissolution, conversion, search. Conversion closes cohab,
opens a linked marriage with the same partners and resets dissolution duration;
neither partner enters availability. Typed histories preserve legal classification
for the default D085 policy. All genuine endings reset search-duration origin.
Conversion is excluded from unpartnered entry counts and is not dissolution;
partnered exposure remains continuous. See [cohabitation](cohabitation.md).

### Weighted choice

All partner matching uses calibrated male-group weights (D075), with no
method argument. A separate `father_matching` model retains its own preferences,
restrictions and delta; a NULL father model inherits the general model.
Ego chooses among
her screened candidates with probabilities proportional to exp(U), multiplied
by a weight for each man's group. Weights are recalculated at every search so
expected male union flows by group meet a target mirroring female per-capita
intensity, multiplied by exp(delta[g]) and normalized (D074). `delta` is a
named, dimensionless finite numeric vector; omitted groups use zero. Adjustments
can represent directional differences or offset composition imbalance, but do
not guarantee a fixed pairwise ratio over time. Event initiation, eligibility,
exclusions and timing are unchanged. Rules, examples, and consequences are in
[weighted-choice matching](weighted-choice-matching.md).

### Shared kinship observations

`extract_kin()` returns the same recorded-parent distance evidence used by
eligibility, with compact nearest common-ancestor witnesses, lineage and
within-depth ancestry completeness. Deceased ancestors remain in the graph.
Simulation and extraction reject ancestry cycles; structural population
validation alone does not perform this graph check. Newborn parent links are
appended before future events; restarts reconstruct the graph from histories.
See [kinship](kinship.md) for worked examples, extraction ordering, depth limits,
and the distinction between evidence and proof of unrelatedness.

## Execution, competition, and random numbers

`run_simulation(model, years, seed = NULL)` processes `[start, start + years)`.
Duration is finite and nonnegative; zero duration returns the population with
updated run metadata. Times that cannot be represented at the chosen numeric
scale are rejected. Candidates exactly at the run endpoint are excluded.

Each process has one clock (D128); the earliest due clock fires. Exact ties
between processes select death, birth, dissolution, conversion, search, then
group transition, so a tied dissolution can follow a birth at the same time.
Calendar changes precede, and cell boundaries are processed before, events at
the same time. Father availability changes immediately when
births/deaths execute. For example, a birth at 0.02 precedes death at 0.06 even
though both occur within the same nominal month.

All engine state belongs to the call. An explicit integer seed isolates the
call and restores the caller's RNG state afterward, including on error.
Without a seed, use and advance the current RNG. Baseline generation has its
own ordinary RNG consumption. The seed contract does not promise fixed-seed
agreement with SOCSIM or older Hammel algorithms.

The output is a `hammel_population` containing all person/union histories,
`end_time`, retained extra components, and latest-run `simulation` metadata.
Restarts recover birth spacing from those histories and rebuild every cell with new thresholds.
They are not exact checkpoints, and split runs need not reproduce uninterrupted
runs. Reattach processes explicitly in a new model, or reuse the original model
with `set_population(result)` to retain its process specifications.

## Observations and derived measures

See [demographic analysis](demographic-analysis.md) for worked workflows.
Analyze one independent history per call; duplicated IDs from stacked
replications are invalid. `by` (D133) optionally separates Lexis and
population-size counts by group, union type, previous union and (where not
always included) sex at the observation time, splitting exposure at recorded
changes. TFR and life tables automatically
preserve group strata in input tables. Exposure/deaths use the person's group,
births use the mother's group at birth, and union formations use each spouse's group at formation.
Unknown-mother births remain in the period-level attribute without assigning a
maternal group. Defaults pool groups; arbitrary custom strata are not supported.
See [groups](groups.md) for output ordering, empty cases, and worked examples.

- `calculate_transition_rates()` counts directed group changes over origin
  person-years by calendar period and either attained age or duration since
  latest group entry (age without earlier transitions). `by` can add `sex`,
  `union_type` and `previous_union`, splitting numerator and denominator. Each destination
  shares its origin's exposure; denominator copies must not be summed across
  destinations. Return `transitions`, `person_years`, and `rate_transition`,
  with NA rates at zero exposure. See [transition observations](group-transitions.md#observed-transition-rates).
- `calculate_fertility_rates()` and `calculate_mortality_rates()` split recorded lifetimes at half-open period/age
  boundaries and reports exposure, deaths, births, ASMR, and ASFR. A finite last
  age breakpoint automatically gains `Inf`; for example `0:110` ends with
  `[110, Inf)`. Zero exposure yields `NA` rates; zero events with positive
  exposure yields zero. Fertility summaries contain female rows only; mortality retains both sexes.
- With `by = "union_type"`, demographic rates partition exposure, births and
  deaths into `unpartnered` and discovered union types in `union_type`, without
  an overall row. Untyped histories use partnered/unpartnered. Components sum
  to pooled counts/exposure obtained with FALSE. Deaths use the left-limit type;
  births tied to dissolution or conversion use the ending type. Exposure is
  half-open; births at formation use the new episode. Typed overlapping
  histories are rejected; untyped coverage is merged. Formation observations
  use `calculate_union_rates()` or `calculate_marriage_rates()` separately.
  Complete histories are required; unknown history is not inferred.
- Births use maternal age and all female exposure, including time excluded by
  spacing. Unknown mothers are counted by period in the tibble attribute
  `unassigned_births`; their births cannot be assigned an age-specific rate.
- `calculate_tfr()` selects female rows, omits open age bands, and sums
  `rate_fertility * band_width` over the remaining contiguous finite bands by period.
  Omit open age bands, warning once per call if an omitted female rate is known
  and nonzero because the reported TFR underestimates full-age fertility.
  Any missing ASFR in the included finite bands makes the corresponding
  period/group/status TFR `NA`, including an all-missing set. Known zero rates
  contribute zero; missing rates in excluded open bands do not affect TFR. There is no `age_range` argument.
- `calculate_life_table()` accepts `rate_mortality` or `rate` as its input
  rate column, preferring `rate_mortality` when both are present. `mx` is
  output only and cannot substitute for or override an input rate column.
  It requires contiguous bands from zero through an
  explicit open interval and separates periods/sexes and union status when
  present. Conditional tables do not model status transitions. For finite width `n`
  and hazard `m`, `q = 1-exp(-m*n)`, survival propagates exponentially, and
  `L = l*(1-exp(-m*n))/m`, with `L=l*n` at zero hazard. In the positive-rate
  open interval, `q=1` and `L=l/m`. `T` is remaining integrated person-years
  and `ex=T/l`. Missing rates propagate uncertainty; positive survivorship
  with zero terminal hazard gives infinite life expectancy; `ex` is undefined
  at zero survivorship. This is a synthetic period table, not cohort follow-up.
- `calculate_union_rates()` returns unpartnered person-years, union starts
  and `rate_union` by period, age and sex. Optional `by = "previous_union"`
  adds previous categories; `by = "group"` independently adds individual
  groups. Both default to FALSE, pooling counts and exposure before rates. The most recent episode determines category; starts
  use the category immediately before formation. Empty categories remain, with
  NA rates at zero exposure. Complete exclusive histories are required;
  overlapping episodes are rejected. Touching and zero-duration episodes count
  separately without adding time. With common age bands, summing categories
  recovers pooled unpartnered exposure and starts. Men describe matching outcomes, not male search inputs.
  Optional `dur_breaks` switches separated/widowed rows to duration since the
  latest ending and automatically separates previous categories (D084). None
  rows retain age bands; unused boundary pairs are NA. Duration rows include
  all ages, and age/duration axes are not crossed.
  `use_detailed_previous_union = FALSE` combines separated/widowed counts and exposure before rates and confidence
  intervals when previous-union stratification is enabled (D126). NULL duration
  breaks preserve age-only output.
  See [remarriage](remarriage.md) for input, restart and boundary conventions.
- `calculate_dissolution_rates()` splits each recorded union episode at calendar
  and duration boundaries, counting union-years once per episode and dissolution
  events only for `end_reason = "dissolution"`. Death also ends exposure but is
  not a dissolution event. Missing ends censor at the window endpoint. Duration
  uses each episode's start; boundaries are half-open and finite final duration
  breaks gain `Inf`. Return `rate_dissolution = dissolutions / union_years`,
  with NA at zero exposure. `by` can split by the wife's and husband's groups
  (`group_w`, `group_h`, directional) and by union type. Input must be a
  population object with a known endpoint (D159). This cause-specific
  rate is not a cumulative probability with competing deaths.
- `calculate_population_size()` includes births and excludes deaths at each
  requested instant. It preserves requested order and repeated times.

Population objects require known endpoints for observation; do not observe
beyond them. Person-level summaries require a population object, so the
endpoint is always known (D159). Missing death means censored at the observation endpoint. Zero-duration lives
contribute zero exposure, while recorded events still count in their applicable
age/period cells.

## Database preparation and provenance

`get_mortality_hmd()` and `get_fertility_hfd()` return prepared rate tibbles,
without configuring a model or writing files. Select one population and one
explicit year or a vector of distinct years (D119). A singleton retains a static
schema; multiple years add calendar bands. Missing years error and nonconsecutive
requests retain gaps. No interpolation, averaging, latest-year selection or
cache refresh occurs. Save/load with RDS to preserve tables and their `source`
attributes. Constructors copy provenance into the process. Credentials are
never retained. HMDHFDplus is optional and needed only for downloads.

HMD defaults to `Mx_1x1`; `shape` also selects native `1x5`, `1x10`, `5x1`,
`5x5`, or `5x10` rates under the accepted hazard interpretation. Native age and
source-period widths are preserved; grouped years select period starts.
Five-year ages use [0,1), [1,5), [5,10), etc.; upper-open ages end at Inf.
Inclusive source periods become half-open intervals, retaining shortened first
or last periods. Respect both raw `+` labels and parsed `OpenInterval` indicators.
For grouped-time products whose reader discards upper bounds, an auxiliary
annual download recovers coverage (HMDHFDplus 2.0.8). Optional `open_age` pools
source deaths and exposure above an aligned cutoff into one open-band rate,
separately by sex/period (D120). Younger rates are unchanged; replaced-tail Mx
may be missing, but pooled counts/exposure must be finite and nonnegative with
positive aggregate exposure. The cutoff is not a maximum lifespan. Missing
retained rates are rejected, not filled. See [database preparation](database-rates.md).

By default HFD uses total `asfrRR` in births per woman-year without scaling.
D095 retains all endpoint rows: map 12- to [12,13) with its supplied rate,
55+ to [55,Inf), and closed ages to one-year bands. This is the selected model
representation of the source aggregates, not a re-estimation of their rates.
TFR includes the lower represented band but drops the upper-open band, warning
if its known rate is nonzero; it need not equal the full published HFD TFR.

Record database, population, requested years, product, represented source
periods, retrieval timestamp and reader version, plus HMD shape and any auxiliary
coverage product. HFD records requested revision date and endpoint mapping. A null
revision requests the current release; it does not pin a historical revision.
Synthetic/local-reader tests establish conversions; live requests are opt-in.
`get_fertility_hfd(..., by_parity = TRUE)` imports
conditional mi rates with half-open current-parity bounds (D094). Endpoint
rows use the D095 mapping above. Unavailable conditional cells are zero-filled under the user-selected D096
preparation rule, with policy/count in source metadata. Parity-dependent
simulation uses actual-parity fertility cells (D101); fertility constructors
accept complete parity partitions (D099 validation rules).
See [parity fertility](parity-fertility.md) for the birth-order/parity mapping.

Reference documentation: [HMDHFDplus source](https://github.com/timriffe/HMDHFDplus),
[HFD FAQ](https://www.humanfertility.org/Project/FAQ), and
[HFD basic indicator methods](https://www.demogr.mpg.de/papers/technicalreports/tr-2011-001.pdf).

`calculate_patfr()` constructs a synthetic cohort from conditional age/parity
rates using HFD's midpoint fertility-table convention (D097). It shares TFR's
period/group/union-type strata and female selection; calculate_tfr() delegates
when parity columns are supplied and returns patfr. Missing included rates give
NA within that stratum. Upper-open ages are omitted with nonzero-rate warning;
finite age bounds are reported. The terminal parity state admits repeated births.
This does not implement parity-dependent simulation, group/union-type transitions
or birth-spacing dynamics. See [parity fertility](parity-fertility.md).

## Known limitations and deferred issues

| Issue | Current treatment |
| --- | --- |
| Initialization transient | No correction or automatic burn-in; users choose observation windows |
| Cells without an eligible woman | Would-be births counted as unrealized (D128); concentrated in sparse young high-parity cells and initialization |
| Rate changes between runs | Histories retain spacing; cell accounting uses the current population, so no transient correction is needed |
| Sparse oldest ages | Observed missing cells remain missing; HFD conditional import alone uses the user-selected zero-fill rule (D096) |
| Union formation/dissolution | Exclusive conditional-logit on-demand or female rate-based formation, minimum male age and known-kin screening; widowhood plus optional duration-dependent dissolution; age- or duration-specific searches by previous union |
| Custom mechanisms/attributes | Read-only custom weight functions (D129) and custom inheritance; no other callbacks or newborn-attribute initializer |
| Exact continuation | No saved clock/RNG checkpoint |

Milestone 3 was accepted October 27, 2026 (D077), based on the user's
[group validation report](../validation/milestone3.qmd). Its assessment covers
group-specific rates, weighted matching, ancestry and identity inheritance,
custom mixed identities, and performance. Small or imbalanced available pools
can give unexpected directional pair ratios; delta adjusts male-group targets
without guaranteeing fixed pairwise ratios. Group-only cross-product identities
are analytical benchmarks, not exact guarantees after pooling changing age
structures and markets. Kinship is omitted from weight calibration and
quarter-year bins approximate calibration matching utilities and cutoff
eligibility, while female intensities sum actual current individual hazards
(D082). Boundary updates consume no randomness and do not alter individual
search scheduling. Actual partner draws use exact ages and screening. No marriage-squeeze feedback is implemented.

Baseline generation is another source of Monte Carlo variation. Replications
sharing one baseline condition on it; independently generated baselines include
that variation. The user's inheritance comparison changed both baseline
replication and maximum initial age, without isolating their contributions.
Completion does not make either an initialization default or guarantee.

Under the earlier midpoint compensation, recently excluded women after a rate
change or status transition reflected the previous schedule, so realized
fertility temporarily differed from the new target (milestone 1 and 2 reports;
D045). Cell accounting (D101) uses the current composition and removes that
transient in expectation; remaining losses are cells without an eligible woman.

## Validation requirements and future mechanism specifications

Current checks cover analytical survival and waiting times, cell accounting
(counts, unrealized births, recovery), cell and clock consistency against
explicit scans, weight-model ratio recovery, constant-rate fertility after burn-in, spacing through restarts,
accounting, identities, father availability, seed behavior, and empty/zero-risk
runs. [Testing](testing.md) separates these from broader demographic assessment.

Milestone 1 was accepted as complete by the user on October 15, 2026 (D029),
based on the completed [validation report](../validation/milestone1.qmd).
This acceptance retains the documented limitations and does not establish
universal accuracy or additional numerical thresholds. For further validation,
make rate/age alignment, observation windows, replications, acceptance criteria,
and sparse-cell treatment explicit. Compare ASFR/ASMR as well as aggregate
TFR/life-table measures, reporting uncertainty and systematic departures
separately. Narrow controlled SOCSIM parity/spacing and heterogeneity checks are recorded
in the audits; they do not establish broad engine equivalence. Compare with a
pinned reference when appropriate. Changed algorithms do not
warrant exact seed matching.

Every new consequential mechanism needs its own plain-language specification:
status/scope, inputs/units, eligibility, timing/competition, state transitions
and invalidation, observations, worked examples, validation, and unresolved
choices. The user decides demographic assumptions. Passing software tests is
not sufficient evidence of research readiness.


## Milestone 4 acceptance

The user accepted dissolution, age/duration remarriage, and optional rate-mode
cohabitation/conversion on October 29, 2026 (D093), based on the
[final report](../validation/milestone4.qmd) and its rendered assessment. This
retained earlier limitations, including transition-aware fertility compensation
(D045, since addressed by D101), thin-pool failures, initialization transients, and approximate calibration
age bins. Reported recovery in these scenarios is not a universal exact pooled
cross-tab constraint. Static SOCSIM audits remain distinct from executable
reference comparisons. Parity-specific execution is implemented (D099 schedules,
D101 cells); D118 implements calendar-time rate schedules; D102 implements lognormal fertility weights.


## Parity-specific fertility execution

D099 permits parity_start/parity_end partitions covering [0,Inf) within each
fertility group/union-type stratum, composed with existing age bands and
multipliers. Counts initialize from recorded maternal births under D098, not
from a custom persons parity column, and increment only on successful births.
The count is not capped at the final band's lower bound. Newborns start at zero;
restarts reconstruct retained histories. No prehistory is imputed.

Births use one fertility cell per actual parity (D101), so a pooled input band
applies its rate to each contained parity. Refresh calibration contributions at
birth. No new parity matching bins are introduced. See
[parity fertility](parity-fertility.md) and
[cell-accounted fertility](fertility-cell-accounting.md). With the user's HFD
2000 inputs, conditional rates are recovered at most ages and parities; young
high-parity shortfalls equal the reported unrealized births (about 0.1% of
births at 50,000 people) and shrink with population size.

Observed parity-specific ASFRs use `calculate_fertility_rates(parity_breaks = ...)`
(D100). Bands begin at zero and contain integer parity boundaries, with Inf
appended to a finite final boundary. All recorded maternal births determine
parity; births count before advancing parity and subsequent exposure uses the
new parity. Spacing-excluded time remains exposure. See [parity observations](parity-fertility.md#observed-fertility-by-parity).


The [spacing-calibration investigation](parity-spacing-investigation.md) records
the parity-rate failures of the earlier preceding-parity correction (D099), a
controlled executable SOCSIM comparison showing the same limitation, and the
evidence that led to cell accounting (D101).


## Milestone 5 acceptance and next scope

The user accepted parity-specific fertility, cell-accounted birth spacing and
persistent lognormal fertility weights on October 29, 2026 (D105), based on
[the milestone 5 report](../validation/milestone5.qmd). Its final render was
active at documentation closeout. Acceptance retains unrealized intensity in
empty-eligibility cells, sparse-cell and structural feasibility limits,
history-dependent parity/spacing, reproductive selection of weights and
non-checkpoint restarts. Adjusted rates that include unrealized expected births
are accounting diagnostics, not realized fertility. No universal scientific
accuracy claim or new acceptance threshold is inferred.

Individual group transitions (milestone 6) are implemented under the
[group-transition contract](group-transitions.md), following explicit user
implementation authorization on October 30. Rates are annual hazards for
directed pairs with a common age or duration axis and optional sex/union-type
strata. Group duration runs from latest entry, or birth without history;
missing pairs and gaps mean zero risk. Moves preserve ancestry and refresh
affected events, fertility cells, matching pools and calibration. Grouped
outputs retain all moves in a `transitions` tibble, even when empty. Household
moves remain deferred beyond milestone 6 (D114). Milestone 9 now
covers multiple group columns, potentially including regions/internal migration
(D121, D131); those mechanisms remain unimplemented. Calendar-time rates are
implemented under D118; matching and inheritance model schedules are implemented
under D123. The user accepted milestone 6 on
October 30, 2026 (D116) after reviewing the report and implementation. Its final
render was active at acceptance; no final-render outcome is inferred.
Acceptance retains documented limitations and does not authorize the next
mechanism or household moves.


## Weight models (D129)

Weight models choose which member of a cell experiences an event, without
changing any cell's intensity. Member i is chosen with probability
w_i / sum(w), where log w_i is a formula's linear predictor x_i * beta or the
value of a custom function; exp(beta) is the hazard ratio between members of
the same cell. A model attaches through `weights` on mortality, fertility,
rate-mode search, dissolution and group transitions, and `weights_conversion`
for conversion, either as one `weight_model()` or a calendar schedule
of them (the D123 tibble structure). Weights are normalized within cells, so
each cell keeps its input rate; population-level effects can be smaller where
cells are small, homogeneous or changing in composition, and terms that are
constant within every cell cancel with a warning. In fertility, formula weights
combine with `hetfert` on the log scale. Variables, the custom-function
contract, missing values, evaluation and matching calibration are specified in
[binned-group hazards](binned-group-hazards.md#weight-models-d129).
No eligibility rules beyond fertility spacing exist; unnormalized weights are
not offered.

## Milestone 8 acceptance

Milestone 8 was accepted complete on November 5, 2026 (D130), based on the
user's [validation report](../validation/milestone8.qmd) and its
[rendered assessment](../validation/milestone8.html). Acceptance covers the
assessed scenarios and retains the documented limitations, including
composition-dependent population-level effects and performance that depends on
model structure.

## Calendar-varying matching and inheritance (D123)

`matching`, `father_matching`, `cohab_matching` and fertility's `inheritance`
accept direct models or tibbles with `time_start`, `time_end` and list-column
`model`. Half-open absolute-year intervals select complete models at event
execution; all rows validate and execution requires coverage. Model-only
boundaries preserve pending dates and do not rewrite existing partners,
identities, ancestry or histories. NULL father/cohabitation roles dynamically
share current general matching. All possible inheritance identities across the
schedule join one registry checked against grouped rates and matching levels.
See [the model-schedule contract](time-varying-models.md) for defaults, worked
examples, calibration updates, callback constraints and validation.

## Optional uncertainty for observed rates (D124)

All event-rate calculators accept `conf_level = NULL` (unchanged default output)
or a confidence level strictly between zero and one. Non-NULL values append
`conf_low` and `conf_high`, using the row's observed count and its existing
exposure in `stats::poisson.test()`. Zero exposure yields NA bounds. Fertility
intervals describe observed births, not the rate adjusted by unrealized
births. These are marginal Poisson-model intervals, not replication-based
uncertainty or a guarantee under demographic dependence. See
[the analysis contract](demographic-analysis.md#optional-poisson-rate-intervals)
for units, examples and edge cases.
