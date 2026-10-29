# Design and implemented interface

Updated October 29, 2026. Milestones 1–4 are complete by user acceptance
(D029, D050, D077, D093). This document describes implemented architecture.
Future extension ideas are listed separately and do not authorize new behavior.
The [scientific contract](scientific-contract.md) defines demographic semantics.

## Responsibilities and source locations

| Responsibility | Implementation |
| --- | --- |
| Model construction and population attachment | `R/hammel_model.R` |
| Group validation and rate lookup | `R/groups.R` |
| Ancestry schema and initialization | `R/ancestry.R` |
| Population schema, validation, baseline generation | `R/hammel_population.R`, `R/checks.R` |
| Mortality/fertility specification and rate validation | `R/set_process.R` |
| Union formation/matching and dissolution specifications | `R/union_formation.R`, `R/union_dissolution.R` |
| Shared kinship traversal and extraction | `src/kinship.h`, `src/kinship.cpp`, `R/kinship.R` |
| HMD/HFD downloads and preparation | `R/construct_demographic_rates.R` |
| Simulation validation, RNG isolation, result assembly | `R/simulate.R` |
| Continuous scheduling and events | `src/simulate.cpp` |
| Demographic extraction and measures | `src/demography.cpp`, `R/calc_demographic_measures.R` |
| Software and analytical validation | `tests/testthat/` |
| Replication assessment | `validation/milestone1.qmd`, `validation/milestone2.qmd`, `validation/milestone3.qmd` |

The engine is C++17 with Rcpp. R exposes ordinary list objects and tibbles.
`Band` and `waiting_time()` are shared by mortality and adjusted fertility.
Every simulation owns its population, queue, and working state; there is no
hidden process-wide engine population.

## Construction is separate from execution

```r
set.seed(42)
population <- create_baseline_population(1000, end_time = 2000)
mortality_rates <- data.frame(
  sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01
)
fertility_rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)

model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates, birth_interval = 0.75) |>
  set_union_formation("on_demand",
    create_matching_model(~ age_difference, c(age_difference = 0)))

result <- simulate(model, years = 10, seed = 123)
```

`hammel_model()` stores `population` and a named `processes` list.
`hammel_population(persons, unions = NULL, end_time = NULL)` constructs the
population; its schema is specified in [population tables](population-schema.md).
`set_population(model, population)` accepts that class, validates it, normalizes
its core tables to tibbles, and preserves extra components. The population's
`end_time` is the next run's start; no duplicate starting time is stored in the
model. An unknown endpoint is allowed for construction, but not execution.

`mortality(rates)` and `fertility(rates, birth_interval = 0.75)` return typed
`hammel_process` objects with `name`, `rates`, and `rate_unit = "person_year"`.
Fertility also stores `birth_interval` and the female-birth probability `ffab`
(default 0.4886); source attributes are copied into the
process when present. `set_process(model, process)` stores a process under its
own name. The demographic setters construct and attach it in one operation.
Replacing a process replaces that specification and retains the others.

Mortality, fertility, and rate-based union formation also accept optional
`group_prop_hazard` vectors on constructors and setters. The process stores the
vector separately from its input rates. Execution scales a private schedule
copy, expanding shared schedules by population group where needed, before
native schedule construction and fertility compensation. Omitted group factors
are 1.0; absent population groups cause a warning and the vector is ignored.
Conditional fertility may instead store a named list of partnered/unpartnered
vectors. Omitted statuses also use 1.0; a single vector scales both statuses.
See [groups](groups.md) for effective-rate examples and validation rules.

Setters validate built-in specifications again, including edited objects.
Mortality constructors permit partial nonoverlapping schedules, but simulation
requires complete coverage. Fertility bands must be finite; gaps have zero risk
at execution. Time-varying rate columns are rejected. Arbitrary named process
specifications can be stored; execution supports mortality/fertility and the
conditional-logit on-demand and rate-based union formation plus optional
duration-dependent union dissolution. Other processes are rejected.

## Execution and output

`simulate(model, years, seed = NULL)` requires population, mortality,
fertility, and an explicit union-formation process. It uses continuous years, with no `dt`, execution step, monthly
compensation grid, or birth-spacing argument on the execution function.

C++ constructs midpoint-adjusted fertility bands, samples cumulative-hazard
thresholds, and processes the earliest event in a sparse ordered queue.
R validates inputs and assembles the returned `hammel_population`, retaining
histories and existing extra person columns. Newborn group identity follows
the fertility inheritance model (random parent by default, with maternal fallback); other custom attributes are
missing. Group-specific schedules and baseline group sampling are documented in
[groups](groups.md). The `simulation` component records `start_time`, `end_time`,
`birth_interval`, `ffab`, `births`, and `deaths` for the latest run, plus
`dissolutions` when that optional process is configured.

Grouped populations also have an `ancestry` tibble, or may supply it explicitly
to `hammel_population(ancestry = ...)`. Initialization uses recorded identities
only when proportions are absent. The engine keeps call-local row-major numeric
ancestry shares, appending parental averages at each birth; it returns the table
in person order with `group_<label>` columns and `unknown`. Missing fathers
contribute half to unknown. Death and restart retain existing shares. This
requires work per birth and storage per person proportional to the number of
components; it does not traverse the genealogy or draw random numbers. See
[ancestry](ancestry.md).

With an explicit seed, R's RNG state is restored after the call. With no seed,
the current state advances normally. See the scientific contract for ties,
run boundaries, parent assignment, and resampling.

## Union-history observations

`calculate_dissolution_rates()` returns union-years, dissolution counts and rates
by calendar period and elapsed duration, optionally by symmetric `group_pair`
strata (A:B includes both partner directions).
The native `dissolution_counts_cpp()` splits recorded episodes directly at
period/duration boundaries. This counts each episode once rather than summing
both partners' exposure. Death censors dissolution exposure without entering
its numerator. See [demographic analysis](demographic-analysis.md).

## Union dissolution

`set_union_dissolution(model, rates, group_prop_hazard = NULL)` attaches the
optional duration-hazard process (D078). Public bands use `dur_start` and
`dur_end`; preparation maps a private copy to the shared `Band` sampler. The
couple's unordered group pair selects schedules and proportional multipliers
(D079). Public `group_pair` labels are normalized; private lookup uses both
partners' identity indices.
Each active union stores one sampled dissolution date, retained through births.
Her next-event selection competes that date with death/birth; birth wins a tie,
then dissolution can execute at the same instant. Death closes the union and
invalidates that risk. Dissolution restores both partners' availability and
reschedules her status-dependent hazards; the man's death event stays valid.
The process is absent by default, consuming no additional random draws then.
See [union dissolution](union-dissolution.md) for examples and deferred scope.

## Newborn identity configuration

`create_group_inheritance()` attaches through fertility's `inheritance` argument.
R validates the mode, possible identity registry, grouped schedule coverage,
and matching levels. Native built-ins assign identity after sex and ancestry
are known, then schedule with the selected group's rates. Custom R functions
receive copied inputs, return a validated registry index, and cannot access
mutable state through the interface. RNG state is synchronized around custom
calls. See [group inheritance](group-inheritance.md).

## Reusing outputs and changing rates

```r
continued <- model |>
  set_population(result) |>
  simulate(years = 5, seed = 124)
```

Reusing the original model retains its processes. Building a new model from
`result` requires explicitly attaching processes; the output stores population
history and run metadata, not a saved model or pending queue.

Users can replace rates between runs. This supports incremental workflows,
but does not implement calendar-varying schedules inside one run. Restarting
redraws pending events and is not guaranteed to match an uninterrupted fixed-seed
run. Fertility adjustment after changing rates ignores the previous schedule;
that transition issue is [documented and deferred](scientific-contract.md#known-limitations-and-deferred-issues).

## Data and observation interfaces

`get_mortality_hmd(country, year, ...)` and `get_fertility_hfd(country, year, ...)`
return rate tibbles, without configuring a model. Users can save/load them with
RDS and reuse them without downloads. HMDHFDplus is optional for these helpers.
The HFD helper's `by_parity = TRUE` option imports conditional mi rates as
age/parity interval rows for milestone 5 preparation. Unavailable
conditional cells are zero-filled under D096, with policy/count metadata. The
engine accepts complete parity partitions under D099 and initializes counts
from recorded maternal births.
See [parity fertility](parity-fertility.md). The default remains total asfrRR. `calculate_patfr()` summarizes conditional
parity schedules with a deterministic HFD-style fertility table; calculate_tfr()
delegates when parity columns are present, returning patfr (D097). Both products preserve endpoint rows under
D095: 12- becomes [12,13), and 55+ becomes [55,Inf).

There is no automatic disk cache or rate refresh.

`calculate_fertility_rates()`, `calculate_mortality_rates()`, `calculate_tfr()`, `calculate_life_table()`, and
`calculate_population_size()` operate on inspectable outputs. See
[demographic analysis](demographic-analysis.md) for exposure, missing cells,
replications, and age alignment. `calculate_fertility_rates(by_union_type = TRUE)`
returns female unpartnered and discovered union-type rows with exposure, births
and fertility rates, without an overall row. The separate mortality helper
returns exposure, deaths and mortality rates for both sexes. Optional
`parity_breaks` on fertility splits observed current parity using birth histories. Untyped histories
use partnered/unpartnered. `calculate_union_rates()` separately returns
unpartnered exposure, union starts, and partnership-entry rates by sex,
excluding conversions.
`by_previous_union = TRUE` adds previous categories; `by_group = TRUE` independently
adds individual groups. Both default to FALSE (D080). Optional `dur_breaks`
automatically separates previous categories and measures dissolved/widowed rows
by duration, leaving none rows by age (D084). It requires exclusive union histories.
Rate-mode input schedules optionally use the same `previous_union` categories;
see [remarriage](remarriage.md) for initialization and calibration details.

`calculate_marriage_rates()` instead counts direct marriages and conversions
against unpartnered plus cohabiting years. It supports age, group and previous
categories, without duration bands. `calculate_dissolution_rates()` reports
non-death dissolution per union-year, optionally by symmetric pair and union
type. `calculate_cohab_convert_rate()` reports conversion per cohabitation
year by duration or woman's age, optionally by woman's group. These risk
sets differ and must not be treated as interchangeable search-rate estimates.

`get_unpartnered_births()` extracts known-parent births outside either union
type, adding parent ages/groups for father-matching assessment. Each birth is
retained, including repeated parental pairs. Birth-first dissolution/conversion
boundaries and half-open ordinary coverage match demographic observations.

`get_detailed_unions(population)` returns a separate tibble preserving union
records and appending each partner's birth/death time and optional group, with
`_w` and `_h` suffixes. It joins by partner IDs, preserves types and row order,
and leaves the population unchanged. See [demographic analysis](demographic-analysis.md).

`get_parent_info(population)` returns a separate tibble preserving all persons
and their columns, adding `age_mother` and `age_father` in years at each
person's birth, plus `group_mother` and `group_father` when groups are recorded.
Parents are matched by ID, including deceased parents; unknown parents yield
missing fields. Row order and group types are preserved, conflicting column
names are rejected, and no endpoint is required. For example, a person born
in 2000.5 whose mother was born in 1975.25 has `age_mother = 25.25`.
Founders with unknown parents retain their rows with missing parent fields;
an empty population retains the enriched schema. These are recorded identities,
not reconstructions of past group transitions. Validation uses shuffled IDs,
fractional ages, siblings, unknown/deceased parents, optional groups, empty
tables, and unchanged inputs.

`calculate_fertility_rates(by_group = TRUE)` and
`calculate_population_size(by_group = TRUE)` optionally split observations by
recorded group. The native count kernels retain full histories and add a group
axis to their output arrays, preserving cross-group parent and spouse links.
TFR and life-table calculations automatically retain input group strata. See
[groups](groups.md) for attribution, ordering, and empty-output rules (D053).

## Extension questions

These remain proposals or unresolved requirements:

- Custom process execution and permitted state mutations/invalidation.
- Poor-pool waiting, previous-father continuity,
  and transition-aware fertility compensation; see [union formation](union-formation-models.md)
  and the [roadmap](roadmap.md).
- Heterogeneity, parity-specific fertility, and user-defined newborn attributes.
- Observation retention controls and exact stochastic checkpoints.
- Calendar-varying rate representations and additional database adapters.

Do not infer future demographic rules from SOCSIM or implement arbitrary
callbacks that can bypass consistency checks.

## Exact full-pool partner draws

There is no candidate sample (D072). `MalePool` in `src/simulate.cpp` keeps
living, unpartnered men in cells by group and quarter-year birth bin, with
constant-time swap removal. Each matching model is prepared once as a
`DrawModel`: the native linear/raw-quadratic age model has a closed-form
maximum over any age-gap interval; custom formulas are tabulated on a
1/32-year gap grid by one R evaluation. `draw_partner()` picks a cell in
proportion to its size times its maximum weight for the woman, a member
uniformly, accepts him with probability exp(U - maximum), and screens only him
for age and kinship, enumerating the eligible pool after 200 rejections. All partner and nonmarital father searches use it. See the [scientific contract](scientific-contract.md#exact-full-pool-partner-draws-d072).

## Weighted-choice matching

Matching always uses calibrated male-group weights (D075); constructors have
no method argument. Separate father models have their own preferences,
restrictions and delta. Grouped searches build the calibrator unless strict
endogamy makes weighting irrelevant. Each man's group log weight is added to
his selection score in the exact draw.
`WeightCalibrator` in `src/simulate.cpp` keeps unpartnered counts by group and
birth-time bin and incrementally maintained pool sums for each matching model's
table of utilities over bin-centre lags, updated at events, and solves the
weights at each search from those sums. `FemaleIntensityTracker` maintains
individual current search/fertility contributions and their group/cohort totals,
with one queued next age/duration/spacing boundary per counted woman
(D082/D083). Search clock origins are separate from fertility birth times. State
changes replace contributions; searches do not scan all women or enumerate a
cross-classification of rate predictors. See
[current-hazard contributions](calibration-intensity.md). Its boundary
queue is the shared `BoundaryQueue` class, also used by `FertilityCells`
(D101, `spacing_method = "cell"`), which keeps cell counts and eligible lists,
keeps one birth clock per cell in its own indexed heap (compared with the main
event queue at each step), processes its boundaries in time order, and
supplies calibration fertility totals at search time. `BoundaryQueue` is a
lazy-deletion binary heap ordered by (time, person). See [cell-accounted fertility](fertility-cell-accounting.md). Matching models carry
`delta` (omitted groups default to zero). One group gives unit weights and female
choice exactly. See [weighted-choice matching](weighted-choice-matching.md).

## Hard group restrictions

Matching models support `prevent_exogamy`: all exogamy or selected directed/
symmetric pairs are ineligible, alongside age and kin restrictions. Compiled
specifications carry either a same-group flag or blocked-pair matrix; the
exact draw skips restricted cells. Partner/father searches use their applicable
model's restrictions. See [group matching](group-matching.md).

## Milestone 2 construction

`create_matching_model()` accepts a formula with age and group predictors,
named utility coefficients, and optional symmetric additive pair effects (`A:B` in both directions).
Categorical group formulas declare fixed `group_levels`; see
[group matching](group-matching.md). `union_formation()` and `set_union_formation()` specify initiation
and a matching specification. The C++ on-demand path creates unions at births
using conditional-logit selection over all living, unpartnered men who pass
minimum-male-age screening (default 15 continuous years) and known-kin
screening (default two generations on each side of a common ancestor). Group
preferences alter utility after screening, without imposing group barriers. See the
[union specification](union-formation-models.md) for worked examples and
deferred extensions. Rate-based initiation schedules searches for unpartnered women.

The engine always builds a call-local current-partner lookup for both sexes once
from ongoing union histories. Per-person union links allow death to close all
affected unions. Formation/death update these caches immediately; restart
reconstructs them from histories. This gives constant-time partner lookup at
birth without adding redundant public `persons` columns. Multiple ongoing
unions for either sex are rejected. Formation removes the selected man from
the pool; widowhood returns a surviving husband. Pool updates use constant-time
swap removal within each group and birth-bin cell.

There is no union-enable flag or no-union execution path. Both initiation modes share the same partner and union state and search helper.
Rate mode adds a search event type and per-person event versions to invalidate
queued births after widowhood. Rate mode selects the current status fertility schedule; an overall schedule
applies to both statuses with one warning per simulation. On-demand requires
overall fertility. Birth spacing persists across status transitions.
Failed searches reschedule at attained age; age-band gaps have zero hazard.

Matching uses the entered formula and named coefficients. R prepares an exact
native linear/raw-quadratic evaluator where applicable; other formulas retain
`model.matrix()` evaluation through an R callback receiving age differences and requested fixed group predictors.
Additive pair effects use a matrix indexed by simulation group codes and are
added after either formula evaluator. No formula callback receives mutable
population state. The native path needs no R call during simulation; custom
formulas are evaluated once on the bound grid and then once per drawn man.

`min_male_age` is stored in the union process and passed to native execution.
The exact draw checks it for each drawn man at the search time; the youngest
eligible birth bin can hold younger men, who are rejected individually.
Existing partner caches bypass new-search restrictions.

The shared `Kinship` module resolves parent IDs to positions, rejects cycles,
and lazily caches distinct ancestors to the configured generation bound.
Screening marks one woman's ancestors and checks candidate membership directly
in C++, skipping these queries at zero generations; R's `extract_kin()`
uses the same kernel and an inverted ancestor index for selected relatives.
Shortest-distance lineage and within-depth completeness accompany extraction
rows. No all-pairs kin matrix or full ancestral-path expansion is constructed.
Newborn append leaves existing ancestry caches valid; restart rebuilds them.


Nonmarital rate-mode births and union searches share a `find_candidate` helper
covering weights, the exact draw and screening. Only `form_union` updates partner
caches and removes the man from the available pool. Parentage alone leaves
both parties unpartnered; the newborn's father index is appended to kinship
independently of current-partner state. No additional persistence state is
retained. `father_matching = NULL` resolves to the general model at execution; an explicit
model supplies a separately prepared draw model and weights for
nonmarital father selection. Previous-father
preference is deferred (D042).

Omitting `matching` uses `create_age_gap_matching_model(0, 5)`. The helper
converts a preferred age gap and half-weight distance in years to the existing
quadratic utility model. `diagnostics = TRUE` adds weighted choice's applied
`.log_weight_<group>` columns to new union records in grouped populations,
without consuming random draws. Diagnostics default off; see the
[union specification](union-formation-models.md#optional-union-diagnostics).


## Typed partnerships (D086)

Optional cohab_prob on union formation enables cohab/marriage histories in rate
mode. Cohab matching and conversion are part of that process; dissolution
optionally stratifies by type. See [cohabitation](cohabitation.md) for complete
inputs, chronology, linked episodes, legal-history recovery and observations.
Fertility input strata are now union_type (the old input union_status errors).
The contribution cache has marriage-search, father and cohab-search channels;
preference tables may be shared, but role weights use each role's intensity.
