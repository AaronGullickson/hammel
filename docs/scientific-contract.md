# Scientific contract

Updated October 15, 2026. This is the current mortality/fertility contract,
reflecting D019, D024, D025, D026, and D028. Earlier monthly and exact
age-varying compensation specifications are available in Git history; the
[decision log](decisions.md) records their replacement. Acceptance
of a model assumption is distinct from evidence of demographic accuracy.

## Units and inputs

Public dates, ages, and durations are numeric years; rates are events per
person-year and exposure is person-years. The user supplies the time origin.
Execution is continuous: there is no month index, configurable execution step,
time truncation, or rounding to birthdays.

Rate tables contain `age_start`, `age_end`, and `rate`. Mortality also contains
`sex`, with separate `"m"` and `"f"` schedules. Intervals are half-open
`[age_start, age_end)`. Starts must be finite and nonnegative; ends must exceed
starts. Rates must be finite and nonnegative, and may exceed one because they
are rates rather than probabilities. Overlapping bands, duplicate bands, and
calendar columns `time_start`/`time_end` are rejected. Constructors preserve
input order and extra columns; native scheduling sorts its own copy.

| Rule | Mortality | Fertility |
| --- | --- | --- |
| Input interpretation | Piecewise constant age/sex hazard | Population ASFR, births per woman-year |
| Denominator for output rates | Living person-years by age/sex | All female person-years, including excluded time |
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

Without `union_status`, all living females face the supplied overall fertility
schedule. With a `union_status` column, require both `partnered` and
`unpartnered` schedules, interpreted per woman-year in that status. Only those
two character labels are permitted; partial or mixed tables are errors. Check
age overlap separately within each status; boundaries may differ and gaps are
zero risk. Use explicit zeros to suppress births for a status. On-demand
formation rejects status-specific fertility schedules, including at zero duration.
Rate formation accepts overall rates with one warning per simulation call. After a recorded birth at time `t`, fertility sampling starts
at `max(now, t + birth_interval)`. Mortality remains active throughout spacing.
A birth at 1.03 with spacing 0.75 permits the next birth no earlier than 1.78.
Zero spacing disables both the exclusion and compensation. Arbitrary finite
nonnegative intervals are allowed; there is no whole-month restriction.

### Midpoint compensation (D028)

Calculate one adjusted hazard per original finite input band. For band `[L,U)`,
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
records the mother, and assigns female sex with probability 0.4886.
A union-formation process is required for every simulation, including zero
duration. Missing processes are rejected rather than given an implicit default.
The former milestone 1 path without unions is no longer executable.

With `union_formation(initiation = "on_demand", ...)`, use the woman's ongoing
living male partner as father. If she has no partner, sample up to `sample_size`
living, unpartnered men without replacement (default 1000; zero means the whole pool).
Evaluate the matching formula on male-minus-female continuous age differences
and multiply its non-intercept model-matrix columns by their named coefficients.
For candidate utilities V, sample with probabilities proportional to
`exp(V - max(V))`, form the union at the birth time, and assign the partner as
father. Nonfinite formula values/utilities are errors. Linear/raw quadratic
terms have an exact native evaluator; other formulas use the R model matrix.
After sampling and before scoring, exclude men younger than `min_male_age`
(default 15 continuous years; finite nonnegative scalar, fractional values
allowed). A man at exactly the threshold is eligible. Do not refill rejected
samples. Also exclude a sampled man if he and the woman share a recorded
ancestor within `kinship_generations` parent steps of EACH (default 2), counting
each person themself at distance zero. This includes known parents, grandparents,
siblings, aunts/uncles and first cousins, including half kin; it does not impose
an unlimited direct-ancestor exclusion. Zero disables cross-person kin exclusion.
Only demonstrated links count: unknown parents do not imply unrelatedness and
do not hide other known relationships. Age and kin search restrictions do not invalidate existing unions,
and no extra minimum female age is imposed. Each person can have at most one ongoing union. Formation removes the man
from the pool before the next sample; widowhood returns a surviving husband.
If no sampled man is eligible, `cancel_birth_if_no_match = TRUE` (the default)
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
events need no resampling merely because a union forms or ends. Candidate
sets are temporary and are not returned. See [union formation](union-formation-models.md)
for worked examples and remaining milestone 2 rules.

Update the mother's latest birth and resample both newborn and mother. Newborns
enter mortality risk immediately; their fertility depends on the supplied ages.
Existing custom person/union columns are retained; new rows have missing custom
values. C++ maintains union histories and shared partner caches for all runs,
closing unions at death with `end_reason = "death_husband"` or
`"death_wife"` according to the executed death event. Initiation mode will
control search timing; it does not control whether partner state is tracked.

### Rate-based union initiation

With `initiation = "rate"`, only living unpartnered women face the supplied
age-specific search hazard. Men never initiate searches. Gaps and ages outside
search bands have zero hazard. Search waiting times use the same cumulative
hazard sampler as deaths and births. Searches compete with death; death wins
exact same-person ties. Search execution uses the same partner sampling,
screening, matching and formation code as on-demand mode. A failed search
reschedules from the attempt time using the same attained-age schedule.
All unpartnered women use that schedule regardless of prior union history.

Women face birth risk in both statuses: the current `union_status` selects
fertility, or an overall schedule applies equally to both. With overall rates,
warn once per simulation call that fertility does not depend on union status.
Unpartnered births in rate mode select a father from living, currently
unpartnered men using the same sample size, minimum age, known-kin screening
and conditional-logit matching as union searches. This creates parentage only:
no union forms, the man remains available, and neither person's partner state
changes. No qualifying sampled candidate cancels the birth by default; explicitly
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
Recovery of conditional ASFRs around status transitions remains unverified.
See [union fertility specification](union-formation-models.md#fertility-schedules-by-union-status-d039).

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

Each person retains one selected next event. Death wins an exact same-person
death/birth tie. Across people, exact ties follow queue insertion order:
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
replications are invalid. There are no custom strata/group arguments yet.

- `calculate_lexis_counts()` splits recorded lifetimes at half-open period/age
  boundaries and reports exposure, deaths, births, ASMR, and ASFR. A finite last
  age breakpoint automatically gains `Inf`; for example `0:110` ends with
  `[110, Inf)`. Zero exposure yields `NA` rates; zero events with positive
  exposure yields zero. Male birth counts/ASFR are `NA`.
- With `include_unions = TRUE`, Lexis counts add exact partnered/unpartnered
  exposure for both sexes, union starts at each spouse's age (including repeat
  unions), and births by maternal status. Union coverage is half-open
  `[start_time, end_time)`; births triggering on-demand formation are partnered.
  Overlapping records do not duplicate exposure. Complete union histories are
  needed: absence of a supplied episode is not evidence of historical singlehood.
  See the worked examples and rate denominators in demographic analysis.
- Births use maternal age and all female exposure, including time excluded by
  spacing. Unknown mothers are counted by period in the tibble attribute
  `unassigned_births`; their births cannot be assigned an age-specific rate.
- `calculate_tfr()` selects female rows, omits open age bands, and sums
  `asfr * band_width` over the remaining contiguous finite bands by period.
  Missing ASFR contributions are omitted with `na.rm = TRUE`; an all-missing
  set sums to zero. The returned finite age range does not indicate complete
  exposure, so inspect missing cells separately. There is no `age_range` argument.
- `calculate_life_table()` requires contiguous bands from zero through an
  explicit open interval and separates periods/sexes. For finite width `n`
  and hazard `m`, `q = 1-exp(-m*n)`, survival propagates exponentially, and
  `L = l*(1-exp(-m*n))/m`, with `L=l*n` at zero hazard. In the positive-rate
  open interval, `q=1` and `L=l/m`. `T` is remaining integrated person-years
  and `ex=T/l`. Missing rates propagate uncertainty; positive survivorship
  with zero terminal hazard gives infinite life expectancy; `ex` is undefined
  at zero survivorship. This is a synthetic period table, not cohort follow-up.
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

HFD uses total `asfrRR` in births per woman-year without scaling. Omit open-age
rows such as `12-` and `55+`; map retained closed ages to one-year intervals.
No fertility is inferred outside them. Input-derived TFR excludes endpoint
contributions and need not equal HFD's full published TFR. Open endpoint rates
are aggregate categories, not evidence for an individual hazard at arbitrarily
young or old ages.

Record database, population, year, product, retrieval timestamp, and reader
version. HFD records requested revision date and endpoint omission. A null
revision requests the current release; it does not pin a historical revision.
Synthetic/local-reader tests establish conversions; live requests are opt-in.
Reference documentation: [HMDHFDplus source](https://github.com/timriffe/HMDHFDplus),
[HFD FAQ](https://www.humanfertility.org/Project/FAQ), and
[HFD basic indicator methods](https://www.demogr.mpg.de/papers/technicalreports/tr-2011-001.pdf).

## Known limitations and deferred issues

| Issue | Current treatment |
| --- | --- |
| Initialization transient | No correction or automatic burn-in; users choose observation windows |
| Midpoint approximation | Validate bias against supplied age schedules; no accuracy guarantee |
| Rate changes between runs | Histories retain spacing, but compensation uses only the new schedule |
| Sparse oldest ages | Missing cells remain missing; no automatic imputation or smoothing |
| Union formation/dissolution | Conditional-logit on-demand formation, minimum male age, known-kin distance screening, and widowhood; exclusive ongoing unions enforced; female rate initiation implemented |
| Custom mechanisms/attributes | No arbitrary callback execution or newborn-attribute initializer |
| Exact continuation | No saved queue/RNG checkpoint |

For changed fertility rates, recently excluded women may reflect the previous
schedule, so realized fertility can temporarily differ from the new target.
Magnitude and duration have not been quantified; there is no guaranteed
nine-month convergence deadline. A recent-schedule history and time-dependent
adjustment were discussed as possible future treatments; the user explicitly
deferred correction at that time. D045 defers transition-aware compensation again, including union-status
transitions and rate changes between runs. The bias remains unquantified; the
existing adjustment and actual-birth spacing are retained.

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
