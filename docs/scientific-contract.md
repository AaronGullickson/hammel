# Scientific contract

Updated October 29, 2026. This is the current mortality, fertility, union-formation,
dissolution, remarriage, cohabitation/conversion, group identity/ancestry, and
observation contract. Milestones 1–4 are accepted as complete
(D029, D050, D077, D093), retaining documented limitations.
Earlier monthly and exact
age-varying compensation specifications are available in Git history; the
[decision log](decisions.md) records their replacement. Acceptance
of a model assumption is distinct from evidence of demographic accuracy.

## Units and inputs

Public dates, ages, and durations are numeric years; rates are events per
person-year and exposure is person-years. The user supplies the time origin.
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
are rates rather than probabilities. Overlapping bands, duplicate bands, and
calendar columns `time_start`/`time_end` are rejected. Constructors preserve
input order and extra columns; native scheduling sorts its own copy.

| Rule | Mortality | Fertility |
| --- | --- | --- |
| Input interpretation | Piecewise constant age/sex hazard | Overall or status-specific ASFR, births per woman-year in the specified population |
| Denominator for output rates | Living person-years by age/sex | Female person-years overall or in the indicated status, including spacing-excluded time |
| Open-ended input band | Allowed; required at execution | Rejected, even with zero rates/spacing |
| Coverage at execution | `[0, Inf)` for each sex, with no gaps | Zero outside supplied finite bands, including gaps |
| Additional setting | None | `birth_interval`, finite nonnegative scalar years; default 0.75 |

Mortality constructors validate structure without requiring full coverage;
`simulate()` checks complete schedules. Specifying and attaching processes does
not execute a simulation. Edited built-in specifications are revalidated at
attachment and execution. Feasibility of fertility compensation is checked
when building the native execution schedule for a positive-duration run.

Treating empirical mortality central rates `mx` as hazards is an explicit
assumption. It is not a universal identity, particularly for infant mortality.
For a constant hazard `h` and elapsed years `d`, the event probability is
`1 - exp(-h*d)`; multiplying a rate by duration gives cumulative hazard,
not a probability. There is no requirement to convert to monthly units.

## Group-specific schedules

Population persons and mortality, fertility, and female union-search rate tables
may contain `group`. Ungrouped rate tables apply to everyone; grouped tables
must cover every population group. Without population groups, use each table's
alphabetically first group and emit a simulation warning. Group does not affect
matching eligibility. Optional group formula predictors and symmetric
pair effects can change utility. Newborns select either parent's group with equal probability by default
(with maternal fallback for unknown fathers)
or the configured inheritance rule, entering the assigned group's mortality
schedule immediately. Fertility
compensation is separate within each group/status schedule. Group values are
nonmissing character labels or finite numeric codes, with nonempty labels.
See [groups](groups.md) for full rules, ordering, examples, and validation scope.

Mortality, fertility, and rate-based union formation optionally accept
`group_prop_hazard`, a named vector of dimensionless finite nonnegative factors.
Multiply the selected input rates by the person's group factor (the woman's for
searches), using 1.0 for omitted groups. Fertility multiplication precedes
midpoint spacing compensation, making it equivalent to explicitly scaled input
ASFRs. Zero suppresses that risk; factors are not normalized. Without population
groups, warn and ignore factors. Shared and stratified schedules are supported;
factors do not supply missing grouped schedules. On-demand union formation has
no search hazard and rejects this option. Stored input tables remain unchanged;
overflow and infeasible effective fertility schedules fail explicitly.
Conditional fertility also accepts a named list of vectors for partnered and/or
unpartnered rates. Omitted statuses and omitted groups within a status use 1.0.
A single vector applies to both statuses; overall fertility requires the vector
form. Each effective group/status schedule is compensated independently.

## Initialization and population history

The [population schema](population-schema.md) specifies required fields and
structural validation. `end_time` is the population reference time and supplies
the next simulation start. Persons born by that time with missing `death_time`
are alive; deceased persons remain in the history and never re-enter risk.

`create_baseline_population(size, max_age = 70, end_time = 0,
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

An optional `ancestry` tibble stores one row per person with `group_<label>`
shares and an explicit `unknown` share. Without supplied proportions, grouped
populations initialize every starting person at 100% of their recorded group,
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

## Shared waiting-time sampling

The engine draws a positive threshold `E` from R's `Exp(1)` generator,
equivalent in distribution to `-log(U)` for a uniform `U` in `(0,1)`.
It consumes `rate * duration` across age bands. If a band reaches the remaining
threshold, the additional within-band wait is `remaining_threshold / rate`.
Zero-rate bands and uncovered fertility gaps consume time without consuming
threshold. If the schedule never reaches the threshold, no future event is
scheduled; no maximum age forces death.

Mortality and adjusted fertility both use the same piecewise constant `Band`
and `waiting_time()` implementation. Crossing a known age boundary does not
require a new random draw. See [the sampling explanation](waiting-time-sampling.md)
for the derivation.

Worked example: two months remain at hazard 0.12/year, followed by 0.24/year.
For threshold 0.05, the first band consumes 0.02; the remaining 0.03 needs
0.125 years. Total wait is 3.5 months, preserved as a fractional number of years.
Months here are explanatory durations, not an execution grid.

## Fertility and birth spacing

Without `union_type`, all living females face the supplied overall fertility
schedule. With a `union_type` column, require both `partnered` and
`unpartnered` schedules, interpreted per woman-year in that status. Alternatively supply all three labels unpartnered/cohab/marriage under D086;
partial or mixed tables are errors. Check
age overlap separately within each status; boundaries may differ and gaps are
zero risk. Use explicit zeros to suppress births for a status. On-demand
formation rejects status-specific fertility schedules, including at zero duration.
Rate formation accepts overall rates with one warning per simulation call. After a recorded birth at time `t`, fertility sampling starts
at `max(now, t + birth_interval)`. Mortality remains active throughout spacing.
A birth at 1.03 with spacing 0.75 permits the next birth no earlier than 1.78.
Zero spacing disables both the exclusion and compensation. Arbitrary finite
nonnegative intervals are allowed; there is no whole-month restriction.

### Newborn group identity

Fertility holds a `create_group_inheritance()` object, defaulting to `random_parent`.
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

### Midpoint compensation (D028)

Calculate one adjusted hazard per original input band. Finite bands retain
midpoint adjustment; an upper-open band uses the interior-limit rule (D095). For band `[L,U)`,
midpoint `m = L + (U-L)/2`, input rate `f`, and spacing `b`:

    eligibility_mid = 1 - integral from max(0, m-b) to m of input_ASFR(a) da
    adjusted_rate = f / eligibility_mid

Evaluate the integral using exact duration overlaps with the unadjusted input
bands within the same status schedule. Gaps and ages outside the schedule contribute zero. Preserve the
original boundaries and use the resulting constant rates in the shared sampler.
Compute this schedule once per simulation call, not from the actual fraction
eligible in the simulated population. Reject nonpositive midpoint eligibility
and unrepresentable adjusted rates rather than clipping or imputing them.
This checks the evaluated midpoints, not the entire age-varying eligibility curve.

At midpoint 25.5, with spacing 0.75, input 0.1 on `[24,25)` and 0.2 on
`[25,26)`, the excluded fraction is `0.25*0.1 + 0.5*0.2 = 0.125`.
Eligibility is 0.875 and the adjusted rate for `[25,26)` is `0.2/0.875`.
A longer spacing window may cross several bands and gaps; all overlaps count.

For an upper-open band `[L, Inf)`, use eligibility `1 - b*f` and adjusted
hazard `f/(1-b*f)` throughout the band. Zero rates remain zero and zero spacing
leaves rates unchanged. Reject nonpositive eligibility. This is the limiting
constant-rate adjustment once the preceding spacing window lies wholly within
the open band; the first b years may overlap earlier bands and are approximated.
There is no artificial midpoint and no finite age cap on a positive open-band rate.

When the window lies wholly in a constant-rate region, adjustment is
`h = f/(1-b*f)`. With `f=0.2` and `b=0.75`, eligibility is 0.85;
mean interval `b + 1/h` is five years, so the long-run renewal rate is 0.2
away from rate boundaries. A constant region with `b*f >= 1` at its midpoint
is infeasible and rejected. For zero spacing, the original schedule is used.

The underlying eligibility identity assumes cohorts experienced the target
schedule during the preceding window. With age-only mortality shared by
eligible and excluded women, survival cancels between recent births and living
female exposure. Fertility-dependent mortality or heterogeneity would require
another specification. Midpoint evaluation approximates the age-varying
adjustment and does not guarantee exact ASFR or TFR recovery after burn-in.

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
`create_age_gap_matching_model(0, 5)`: preferred male-minus-female gap zero
years, with relative weight halved five years either side. Explicit models
override that default; it does not fix the realized age-gap distribution.

Matching models optionally accept `group_exogamy`: a named vector using `A:B`
for equal effects on female A/male B and female B/male A (D075). Effects
are finite log-scale additions after the full formula score. An omitted pair
has odds ratio 1.0 (zero log effect). Within-group effects are allowed; overlaps
after symmetric expansion are errors. A pair table is also accepted
for labels containing colons; each row expands in both directions. These are selection preferences, not hazards or
prohibitions. `female_group`, `male_group`, and numeric 0/1 `same_group` are
available in formulas; categorical terms require fixed declared `group_levels`
with treatment contrasts. Missing population identity data errors when a model
uses group effects. See [group matching](group-matching.md) for worked cases.

Matching models also accept `prevent_exogamy = FALSE` (no restrictions), TRUE
(only same-group candidates), or selected cross-group pairs with female-first
`:` and symmetric `::` notation. Restricted men are ineligible, alongside the
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
mortality remain independent of union status, so surviving partners' pending
events need no resampling merely because a union forms or ends. Optional
`diagnostics = TRUE` records only weighted choice's applied log group weights
(`.log_weight_<group>`) on newly formed unions in grouped populations; it
defaults off and consumes no additional random draws. See [union formation](union-formation-models.md)
for worked examples and implemented milestone 2 rules.

Update the mother's latest birth and resample both newborn and mother. Newborns
enter mortality risk immediately; their fertility depends on the supplied ages.
Existing custom person/union columns are retained; new rows have missing custom
values, except that `group` is assigned by the fertility inheritance model. C++ maintains union histories and shared partner caches for all runs,
closing unions at death with `end_reason = "death_husband"` or
`"death_wife"` according to the executed death event. Initiation mode
controls search timing; it does not control whether partner state is tracked.

### Rate-based union initiation

With `initiation = "rate"`, only living unpartnered women face the supplied
age-specific search hazard. Men never initiate searches. Gaps and ages outside
search bands have zero hazard. Search waiting times use the same cumulative
hazard sampler as deaths and births. Searches compete with death; death wins
exact same-person ties. Search execution uses the same exact draw,
screening, matching and formation code as on-demand mode. A failed search
reschedules from the attempt time using the same attained-age schedule.
Without `previous_union` stratification, all unpartnered women use that schedule
regardless of prior union history. Stratified schedules distinguish none,
dissolved and widowed, with age or permitted duration clocks (D080–D083).

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

A new union or husband's death invalidates the woman's pending event and
reschedules at attained age with the new fertility schedule, preserving her
recorded birth spacing. Widowhood also restores search risk. Versioned queue
entries prevent invalidated candidates executing. Surviving men's death hazards
do not depend on union status. Midpoint compensation is computed independently
within each schedule; it is not calibrated to the woman's prior-status history.
The milestone 2 report examines conditional ASFR discrepancies; general
recovery across status transitions is not established.
See [union fertility specification](union-formation-models.md#fertility-schedules-by-union-status-d039).

### Duration-dependent union dissolution (D078)

An optional `union_dissolution` process adds one annual hazard per union,
carried by its female partner and indexed by elapsed years since `start_time`.
Shared rates or schedules/multipliers by unordered group pair (D079) use
half-open duration bands;
gaps and uncovered durations have zero hazard. New episodes start at duration
zero; active input unions and restarts use attained duration. The union's
sampled dissolution date is retained through births. Exact same-person ties
select death, birth, dissolution, then search. A tied birth occurs before
closure and retains that husband as father.

Execution ends the episode with reason `dissolution`, clears partner links,
restores male availability and female calibration counts with spacing preserved,
and reschedules her unpartnered fertility and existing searches. Men's death
hazards are status-independent, so their pending death events remain valid.
Either partner's death cancels the episode's risk. Without this process,
existing death-only behavior and random-number consumption are retained.
Search schedules may distinguish `previous_union` as `none`,
`dissolved`, or `widowed` (D080). The most recent recorded episode determines
state; complete histories are required. Without that column, all categories
share the same rates. All three categories must be supplied per scheduled
group when stratifying. `none` must use age; `dissolved` and `widowed` each
select age or duration since the latest union ending (D083). Each category uses
one clock in all bands/groups. Duration columns are optional for all-age tables;
when included, each row supplies one complete boundary pair with the unused
pair NA. No simultaneous age/duration hazards are permitted. Recover the latest
ending origin on restart and reset it at each subsequent ending. The
contribution queue handles the selected search clock independently of fertility
ages and matching cohort bins. Cohabitation is implemented under D086; richer couple covariates remain deferred. See [remarriage](remarriage.md) and the
[dissolution specification](union-dissolution.md) for worked examples and limits.

### Optional cohabitation (D085/D086)

Rate formation accepts cohab_prob age/group/previous-category probabilities,
with gaps zero. At each attempt, intended type is drawn before matching;
cohab_matching optionally overrides general preferences. Unpartnered father
matching remains separate. Calibration uses lambda*p for cohabitation,
lambda*(1-p) for direct marriage, and fertility intensity for fathers.
Probability boundaries update the common contribution queue. NULL probability
keeps the untyped/no-extra-RNG path; on-demand rejects cohabitation inputs.

Fertility input uses union_type: either unpartnered/partnered or
unpartnered/cohab/marriage per group; old input union_status errors. Fertility
multiplier-list keys must exist in those schedules, with omitted factors 1.
Dissolution rates may supply both union types within each pair or share rates.
Conversion uses a female age/cohab-duration clock, optional group, no previous
axis; NULL warns and disables conversion. Pending clocks survive births. Ties
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

`simulate(model, years, seed = NULL)` processes `[start, start + years)`.
Duration is finite and nonnegative; zero duration returns the population with
updated run metadata. Times that cannot be represented at the chosen numeric
scale are rejected. Candidates exactly at the run endpoint are excluded.

Each person retains one selected next event. Exact same-person ties select
death, birth, dissolution, then search. A union retains its sampled dissolution
date through births, so a tied dissolution can follow a birth at the same time. Across people, exact ties follow queue insertion order:
initial row order; after birth, newborn before mother. Otherwise execute the
earliest continuous candidate. Father availability changes immediately when
births/deaths execute. For example, a birth at 0.02 precedes death at 0.06 even
though both occur within the same nominal month.

All engine state belongs to the call. An explicit integer seed isolates the
call and restores the caller's RNG state afterward, including on error.
Without a seed, use and advance the current RNG. Baseline generation has its
own ordinary RNG consumption. The seed contract does not promise fixed-seed
agreement with SOCSIM or older Hammel algorithms.

The output is a `hammel_population` containing all person/union histories,
`end_time`, retained extra components, and latest-run `simulation` metadata.
Restarts recover birth spacing from those histories and resample pending events.
They are not exact checkpoints, and split runs need not reproduce uninterrupted
runs. Reattach processes explicitly in a new model, or reuse the original model
with `set_population(result)` to retain its process specifications.

## Observations and derived measures

See [demographic analysis](demographic-analysis.md) for worked workflows.
Analyze one independent history per call; duplicated IDs from stacked
replications are invalid. `by_group = TRUE` optionally separates Lexis and
population-size counts by the recorded group. TFR and life tables automatically
preserve group strata in input tables. Exposure/deaths use the person's group,
births use the mother's group, and union formations use each spouse's own group.
Unknown-mother births remain in the period-level attribute without assigning a
maternal group. Defaults pool groups; arbitrary custom strata are not supported.
See [groups](groups.md) for output ordering, empty cases, and worked examples.

- `calculate_fertility_rates()` and `calculate_mortality_rates()` split recorded lifetimes at half-open period/age
  boundaries and reports exposure, deaths, births, ASMR, and ASFR. A finite last
  age breakpoint automatically gains `Inf`; for example `0:110` ends with
  `[110, Inf)`. Zero exposure yields `NA` rates; zero events with positive
  exposure yields zero. Fertility summaries contain female rows only; mortality retains both sexes.
- With `by_union_type = TRUE`, demographic rates partition exposure, births and
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
- `calculate_life_table()` requires contiguous bands from zero through an
  explicit open interval and separates periods/sexes and union status when
  present. Conditional tables do not model status transitions. For finite width `n`
  and hazard `m`, `q = 1-exp(-m*n)`, survival propagates exponentially, and
  `L = l*(1-exp(-m*n))/m`, with `L=l*n` at zero hazard. In the positive-rate
  open interval, `q=1` and `L=l/m`. `T` is remaining integrated person-years
  and `ex=T/l`. Missing rates propagate uncertainty; positive survivorship
  with zero terminal hazard gives infinite life expectancy; `ex` is undefined
  at zero survivorship. This is a synthetic period table, not cohort follow-up.
- `calculate_union_rates()` returns unpartnered person-years, union starts
  and `rate_union` by period, age and sex. Optional `by_previous_union = TRUE`
  adds previous categories; `by_group = TRUE` independently adds individual
  groups. Both default to FALSE, pooling counts and exposure before rates. The most recent episode determines category; starts
  use the category immediately before formation. Empty categories remain, with
  NA rates at zero exposure. Complete exclusive histories are required;
  overlapping episodes are rejected. Touching and zero-duration episodes count
  separately without adding time. With common age bands, summing categories
  recovers pooled unpartnered exposure and starts. Men describe matching outcomes, not male search inputs.
  Optional `dur_breaks` switches dissolved/widowed rows to duration since the
  latest ending and automatically separates previous categories (D084). None
  rows retain age bands; unused boundary pairs are NA. Duration rows include
  all ages, and clocks are not crossed or pooled. NULL preserves age-only output.
  See [remarriage](remarriage.md) for input, restart and boundary conventions.
- `calculate_dissolution_rates()` splits each recorded union episode at calendar
  and duration boundaries, counting union-years once per episode and dissolution
  events only for `end_reason = "dissolution"`. Death also ends exposure but is
  not a dissolution event. Missing ends censor at the window endpoint. Duration
  uses each episode's start; boundaries are half-open and finite final duration
  breaks gain `Inf`. Return `rate_dissolution = dissolutions / union_years`,
  with NA at zero exposure. Optional grouping returns symmetric `group_pair`
  strata using both partners, with both cross-group directions pooled. Bare
  union tables require complete caller-verified follow-up. This cause-specific
  rate is not a cumulative probability with competing deaths.
- `calculate_population_size()` includes births and excludes deaths at each
  requested instant. It preserves requested order and repeated times.

Population objects require known endpoints for observation; do not observe
beyond them. Bare persons tables require the caller to ensure complete
follow-up. Missing death means censored at the observation endpoint. Zero-duration lives
contribute zero exposure, while recorded events still count in their applicable
age/period cells.

## Database preparation and provenance

`get_mortality_hmd()` and `get_fertility_hfd()` return prepared rate tibbles,
without configuring a model or writing files. Select one population and one
explicit year; no latest-year selection, interpolation, averaging, or automatic
cache refresh occurs. Save/load with RDS to preserve tables and their `source`
attributes. Constructors copy provenance into the process. Credentials are
never retained. HMDHFDplus is optional and needed only for downloads.

HMD uses `Mx_1x1` for both sexes and the accepted hazard interpretation. Closed
single-year ages become `[age,age+1)` and the open age becomes `[age,Inf)`.
Respect both raw `+` labels and parsed `OpenInterval` indicators, including
HMDHFDplus 2.0.8 behavior. Missing retained rates are rejected, not filled.

By default HFD uses total `asfrRR` in births per woman-year without scaling.
D095 retains all endpoint rows: map 12- to [12,13) with its supplied rate,
55+ to [55,Inf), and closed ages to one-year bands. This is the selected model
representation of the source aggregates, not a re-estimation of their rates.
TFR includes the lower represented band but drops the upper-open band, warning
if its known rate is nonzero; it need not equal the full published HFD TFR.

Record database, population, year, product, retrieval timestamp, and reader
version. HFD records requested revision date and endpoint mapping. A null
revision requests the current release; it does not pin a historical revision.
Synthetic/local-reader tests establish conversions; live requests are opt-in.
`get_fertility_hfd(..., by_parity = TRUE)` imports
conditional mi rates with half-open current-parity bounds (D094). Endpoint
rows use the D095 mapping above. Unavailable conditional cells are zero-filled under the user-selected D096
preparation rule, with policy/count in source metadata. D099 implements
parity-dependent simulation with the selected preceding-parity correction.
Fertility constructors accept complete parity partitions under D099.
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
| Midpoint approximation | Validate bias against supplied age schedules; no accuracy guarantee |
| Rate changes between runs | Histories retain spacing, but compensation uses only the new schedule |
| Sparse oldest ages | Observed missing cells remain missing; HFD conditional import alone uses the user-selected zero-fill rule (D096) |
| Union formation/dissolution | Exclusive conditional-logit on-demand or female rate-based formation, minimum male age and known-kin screening; widowhood plus optional duration-dependent dissolution; age- or duration-specific searches by previous union |
| Custom mechanisms/attributes | No arbitrary callback execution or newborn-attribute initializer |
| Exact continuation | No saved queue/RNG checkpoint |

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

For changed fertility rates, recently excluded women may reflect the previous
schedule, so realized fertility can temporarily differ from the new target.
The milestone reports examine particular rate-change and status-transition
cases, but do not establish a general magnitude or duration; there is no
guaranteed nine-month convergence deadline. A recent-schedule history and time-dependent
adjustment were discussed as possible future treatments; the user explicitly
deferred correction at that time. D045 defers transition-aware compensation again, including union-status
transitions and rate changes between runs. The existing adjustment and
actual-birth spacing are retained; see the user-authored reports for case-specific
assessment rather than a general accuracy guarantee.

## Validation requirements and future mechanism specifications

Current checks cover analytical survival and waiting times, midpoint overlap
examples, constant-rate fertility after burn-in, spacing through restarts,
accounting, identities, father availability, seed behavior, and empty/zero-risk
runs. [Testing](testing.md) separates these from broader demographic assessment.

Milestone 1 was accepted as complete by the user on October 15, 2026 (D029),
based on the completed [validation report](../validation/milestone1.qmd).
This acceptance retains the documented limitations and does not establish
universal accuracy or additional numerical thresholds. For further validation,
make rate/age alignment, observation windows, replications, acceptance criteria,
and sparse-cell treatment explicit. Compare ASFR/ASMR as well as aggregate
TFR/life-table measures, reporting uncertainty and systematic departures
separately. No controlled executable SOCSIM comparison has been performed;
compare with a pinned reference when appropriate. Changed algorithms do not
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
retains earlier limitations, including transition-aware fertility compensation
(D045), thin-pool failures, initialization transients, and approximate calibration
age bins. Reported recovery in these scenarios is not a universal exact pooled
cross-tab constraint. Static SOCSIM audits remain distinct from executable
reference comparisons. Parity-specific execution is implemented under D099; calendar-time schedules
and fertility heterogeneity remain unimplemented.


## Parity-specific fertility execution

D099 permits parity_start/parity_end partitions covering [0,Inf) within each
fertility group/union-type stratum, composed with existing age bands and
multipliers. Counts initialize from recorded maternal births under D098, not
from a custom persons parity column, and increment only on successful births.
The count is not capped at the final band's lower bound. Newborns start at zero;
restarts reconstruct retained histories. No prehistory is imputed.

Actual parity zero uses unadjusted input fertility. Higher parities use their
own rate divided by an eligibility approximation computed from the preceding
actual parity's unadjusted age schedule. Retain finite-band midpoints, upper-
open interior limits and actual-birth spacing. Current group/union type is held
for this computation; earlier transitions are not reconstructed. Refresh current
fertility/calibration contributions at birth. No new parity matching bins are
introduced. See [parity fertility](parity-fertility.md) for validation, formulae
and the highest-band distinction.

This is a SOCSIM-style approximation, not an exact exposure adjustment. In the
controlled initialized-cohort check, .75-year spacing yields mean parity-1 and
2+ rates .2665 and .3335 against inputs .3 and .4. The limitation is explicit;
software correctness does not establish conditional-rate recovery or milestone 5
acceptance. Alternative exposure-based compensation remains a future decision.

Observed parity-specific ASFRs use `calculate_fertility_rates(parity_breaks = ...)`
(D100). Bands begin at zero and contain integer parity boundaries, with Inf
appended to a finite final boundary. All recorded maternal births determine
parity; births count before advancing parity and subsequent exposure uses the
new parity. Spacing-excluded time remains exposure. See [parity observations](parity-fertility.md#observed-fertility-by-parity).
