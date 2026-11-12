# Design and implemented interface

Updated November 9, 2026. Milestones 1–8 are complete by user acceptance
(D029, D050, D077, D093, D105, D116, D127, D130). Individual group transitions are accepted
under D116; milestone 7 calendar-varying rates are implemented under D118,
accepted complete November 2, 2026 (D127). D123 implements matching/inheritance model
schedules. This document describes implemented architecture.
Future extension ideas are listed separately and do not authorize new behavior.
Milestone 8 schedules every rate process by cells with one clock per process
(D128) and adds weight models (D129); see
[cell scheduling and weight models](binned-group-hazards.md).
The [scientific contract](scientific-contract.md) defines demographic semantics.

## Responsibilities and source locations

| Responsibility | Implementation |
| --- | --- |
| Model construction and validation | `R/hammel_model.R` |
| Population schema, model attachment and baseline generation | `R/hammel_population.R`, `R/generate_uniform_population.R`, `R/checks.R` |
| Group declarations and histories | `R/group.R`, `R/group_history.R` |
| Ancestry schema and initialization | `R/ancestry.R` |
| Process attachment | `R/set_process.R` |
| Mortality/fertility specifications | `R/mortality.R`, `R/fertility.R` |
| Input-rate validation, calendar coverage, group strata and engine encoding | `R/rates_util.R` |
| Proportional input-rate stratification | `R/stratify_rates.R` |
| Calendar model schedules and engine preparation | `R/model_schedules.R`; switches in `src/simulate.cpp` |
| Group transition and household specifications | `R/group_transition.R`; events in `src/simulate.cpp` |
| Persistent fertility heterogeneity (internal, D176) | `R/het_fert.R`; cell weight trees in `src/fertility_cells.h` |
| Weight models | `R/weight_model.R`; stored weights and selection in `src/simulate.cpp`, `src/rate_cells.h` and `src/fertility_cells.h` |
| Union formation, dissolution and cohabitation settings | `R/union_formation.R`, `R/union_dissolution.R`, `R/cohabitation.R` |
| Matching and inheritance models | `R/matching_model.R`, `R/matching_groups.R`, `R/group_inheritance.R` |
| Kinship traversal, extraction and selective living-kin counting | `src/kinship.h`, `src/kinship.cpp`, `R/extract_kin.R` |
| HMD/HFD downloads, conversion and provenance | `R/get_database_rates.R` |
| Simulation validation, RNG isolation and result assembly | `R/run_simulation.R` |
| Continuous scheduling and events | `src/simulate.cpp` (the `Simulation` engine); rate processes in `src/rate_cells.h` and `src/fertility_cells.h`; shared queues, sum trees and clocks in `src/event_queues.h`; rate bands in `src/rate_bands.h` |
| Partner matching, exact draws and weighted-choice calibration | `src/matching.h`, `src/matching.cpp`, `src/weight_calibration.h`, `src/weight_calibration.cpp` |
| Shared history preparation and analysis strata | `R/analysis_utils.R`; native observation in `src/observation.cpp` |
| Observed event-rate calculations | `R/calculate_rates.R`; native counts in `src/observation.cpp` and `src/demography.cpp` |
| Population counts and rate summaries | `R/calculate_population_size.R`, `R/calculate_tfr.R`, `R/calculate_life_table.R`, `R/summary_rates.R` |
| History augmentation and unpartnered-birth extraction | `R/augment_histories.R` |
| Software and analytical validation | `tests/testthat/` |
| Replication assessment | `validation/milestone1.qmd` through `validation/milestone9.qmd` |

See [R source organization](r-source-organization.md) for the complete function map and placement rules (D143). Every multiline R function definition must have a braced body, including anonymous callbacks and single-expression bodies (D171).

The engine is C++17 with Rcpp. R exposes ordinary list objects and tibbles.
Grouped dissolution and conversion observations build union-level group
histories through `union_group_history_cpp()` using the shared `GroupHistory`
index. Each union visits only its partners' transition histories; pair labels
are mapped once before native traversal.
Every rate process uses the same cell machinery (`RateCells`, or
`FertilityCells` for fertility): an `IntensityTree` over cells, one
`ProcessClock` and a `BoundaryQueue` for band-edge moves (D128). Every
simulation owns its population, cells, clocks and working state; there is no
hidden process-wide engine population.

## Bundled England and Wales example data

`england_wales_1992_onward` is a named list of tibbles for the second
demographic-transition example (D173):

| Member | Contents | Preparation |
| --- | --- | --- |
| `fertility` | 18 union-type-specific 1992 fertility rates | Standard age columns; no source-label or duplicate-rate columns |
| `marriage` | 672 first-marriage/remarriage rows, 1992–2019 | Standard calendar, age and previous-union columns; 48 missing rates retained |
| `premarital_cohabitation` | 373 female shared-address proportions at marriage | `year`, `age_groups`, `prop_cohabit`; source values retained |
| `cohabitation_outcomes` | 3 Great Britain cohabitation-cohort outcome rows | Source columns, rounding and missingness retained |
| `marital_divorce_death_2005` | 52 cumulative marital divorce/death observations | Source columns and values retained |
| `marital_divorce_death_2010` | 60 cumulative marital divorce/death observations | Source columns and values retained |

The bundle provides the local numerical inputs, not an already calibrated
simulation scenario. In particular, marriage rates have not been filled or
tuned, probabilities have not been converted to hazards, cohabitation cohorts
have not been mapped to calendar periods, and mortality is not included.
Fertility uses the vignette's denominator-age bounds of 16–45; source birth
numerators are not strictly restricted to these end-band ages. Source files
and provenance remain under `validation/data/england_wales_1992_onward/`;
`data-raw/england_wales_1992_onward.R` builds the bundle offline.

## Construction is separate from execution

```r
set.seed(42)
population <- generate_uniform_population(1000, end_time = 2000)
mortality_rates <- data.frame(
  sex = c("m", "f"), age_start = 0, age_end = Inf, rate = 0.01
)
fertility_rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)

model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates, birth_interval = 0.75) |>
  set_union_formation("on_demand",
    matching = matching_model(~ age_difference, c(age_difference = 0)))

result <- run_simulation(model, years = 10, seed = 123)
```

`hammel_model()` stores `population` and a named `processes` list.
`hammel_population(persons, unions = NULL, end_time)` constructs the
population; its schema is specified in [population tables](population-schema.md).
`set_population(model, population)` accepts that class, validates it, normalizes
its core tables to tibbles, and preserves extra components. The population's
`end_time` is the next run's start; no duplicate starting time is stored in the
model. The endpoint is required and must be a single finite numeric year at
construction, validation and attachment (D166).

`mortality(rates)` and `fertility(rates, birth_interval = 0.75)` return typed
`hammel_process` objects with `name` and `rates`; rates are always per person-year.
Fertility also stores `birth_interval` and the female-birth probability `ffab`
(default 0.4886); source attributes are copied into the
process when present. `set_process(model, process)` stores a process under its
own name. The demographic setters construct and attach it in one operation.
Replacing a process replaces that specification and retains the others.

`stratify_rates(rates, by)` in `R/stratify_rates.R` adds proportional strata to
a rate table before it is attached (D132). Each setter checks the table's
stratum columns (`stratum_columns`) against those its process supports.
Conditional fertility may instead store a named list of partnered/unpartnered
vectors. Omitted statuses also use 1.0; a single vector scales both statuses.
See [groups](groups.md) for effective-rate examples and validation rules.

Setters validate built-in specifications again, including edited objects.
Mortality constructors permit partial nonoverlapping schedules, but simulation
requires complete coverage. Fertility bands must be finite; gaps have zero risk
at execution. Optional `time_start`/`time_end` columns select half-open calendar
bands on all six annual-rate inputs; static tables remain supported (D118). Arbitrary named process
specifications can be stored; execution supports mortality/fertility and the
conditional-logit on-demand and rate-based union formation plus optional
duration-dependent union dissolution. Optional `group_transition()` is also supported. Other processes are rejected.

## Execution and output

`run_simulation(model, years, seed = NULL, progress = FALSE)` requires population, mortality,
fertility, and an explicit union-formation process. It uses continuous years, with no `dt`, execution step, monthly
compensation grid, or birth-spacing argument on the execution function.

`progress` must be a single nonmissing logical value. When TRUE, a console
progress bar reports elapsed simulation time divided by the requested duration,
not the proportion of wall-clock work completed. For example, a run starting
at 2000 for 20 years reports 25% after processing time 2005. Native progress
updates occur after chronological events or bookkeeping/calendar boundaries,
at most once per percentage point; they add no events or random draws. R uses
`utils::txtProgressBar()` and closes the bar on success, error or interruption.
The bar reaches 100% only after successful result assembly. Empty populations,
zero-risk periods and zero-year runs also finish at 100%; initialization and
result assembly have no separate percentage allocation. FALSE produces no
progress output. Progress is a display setting, excluded from result metadata.
Validation compares complete results and RNG states with progress enabled and
disabled, checks nonzero/fractional starts and durations, monotone bounded
updates, no-event runs, invalid arguments and cleanup on failure.

C++ builds cells for every rate process from the input schedules (D101,
D128) and processes the earliest of the next calendar change, the next cell
boundary and the next process clock; `refresh_person()` and `refresh_union()`
update a unit's memberships in every process after each state change.
R validates inputs and assembles the returned `hammel_population`, retaining
histories and existing extra person columns. Newborn group identity follows
the dimension's inheritance model (random parent by default, with maternal fallback); other custom attributes are
missing. Group-specific schedules and baseline group sampling are documented in
[groups](groups.md). The `simulation` component records the latest run's
`start_time`, `end_time`, `birth_interval` and `ffab`, plus `hetfert` when enabled
and `cohab_as_previous_union` when cohabitation is enabled.
Event counts are not duplicated there (D167); events remain in person, union
and transition histories. `unrealized_births` and `unrealized_since` are
population components accumulated across chained runs (D160).

Grouped populations also have an `ancestry` tibble, or may supply it explicitly
to `hammel_population(ancestry = ...)`. Initialization uses recorded identities
only when proportions are absent. The engine keeps call-local row-major numeric
ancestry shares, appending parental averages at each birth; it returns the table
in person order as one tibble per tracked dimension, with a column per category and `unknown`. Missing fathers
contribute half to unknown. Death and restart retain existing shares. This
requires work per birth and storage per person proportional to the number of
components; it does not traverse the genealogy or draw random numbers. See
[ancestry](ancestry.md).

With an explicit seed, R's RNG state is restored after the call. With no seed,
the current state advances normally. See the scientific contract for ties,
run boundaries, parent assignment, and restarts.

## Union-history observations

`calculate_dissolution_rates()` returns union-years, dissolution counts and rates
by calendar period and elapsed duration, optionally by the wife's and
husband's groups and union type (`by`, D133).
The native `dissolution_counts_cpp()` splits recorded episodes directly at
period/duration boundaries. This counts each episode once rather than summing
both partners' exposure. Death censors dissolution exposure without entering
its numerator. See [demographic analysis](demographic-analysis.md).

## Union dissolution

`set_union_dissolution(model, rates, weights = NULL)` attaches the
optional duration-hazard process (D078). Public bands use `dur_start` and
`dur_end`; preparation maps a private copy to duration cells (D128). The
wife's and husband's categories select schedules through directional
`<group>_w`/`<group>_h` columns (D134; this replaced D079's symmetric
`group_pair`).
Each active union is a unit in the dissolution cells, unaffected by births;
birth wins a tie, then dissolution can execute at the same instant. Death
closes the union and removes it from the cells. Dissolution restores both
partners' availability and moves her to her unpartnered cells.
The process is absent by default, consuming no additional random draws then.
See [union dissolution](union-dissolution.md) for examples and deferred scope.

## Newborn identity configuration

`group_inheritance()` attaches per dimension through `group(inheritance = ...)`.
The `"fixed"` method (D142) uses an explicit `fixed_value` category, validated against
the registry and assigned in C++ with no callback or identity random draw.
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
  run_simulation(years = 5, seed = 124)
```

Reusing the original model retains its processes. Building a new model from
`result` requires explicitly attaching processes; the output stores population
history and run metadata, not a saved model or pending clocks.

Users can replace rates between runs or supply calendar bands inside one run.
The engine selects active rate rows and processes calendar boundaries before
same-time events. A process whose cell layout is unchanged takes the new
rates in place and keeps its clock (`same_layout`/`set_rates`,
`FertilityCells::update_rates`); otherwise its cells are rebuilt, preserving
spacing, weights, parity and histories. Matching calibration contributions
follow the cells. Weight schedules switch the active model and recompute stored
weights without touching clocks.
This stays within one native invocation with independent call-local state;
matching and inheritance models can also vary through model-list tibbles (D123).
Model-only changes preserve pending dates and histories. Role preference/draw
structures and calibration tables are cached per model within the call, with
active calibration sums rebuilt from the current male pool when preferences
change. All inheritance models share one registry including future declared
identities. Direct static models retain the existing execution path. See
[calendar-varying rates](time-varying-rates.md) and
[calendar-varying models](time-varying-models.md).

Restarting rebuilds cells with new thresholds and is not guaranteed to match an
uninterrupted fixed-seed run. A calendar boundary that rebuilds cells consumes
random draws, so it does not promise fixed-seed equality with a static table.

## Data and observation interfaces

`get_mortality_hmd(country, year, ...)` and `get_fertility_hfd(country, year, ...)`
return rate tibbles, without configuring a model. One requested year returns
static rates; a vector returns calendar-banded rates. HMD's `shape` selects one
of six native products, defaulting to `1x1`; grouped periods select source
starting years. Each required product is downloaded once, including auxiliary
annual coverage when a parsed HMD grouped-time product needs its endpoint
recovered. Optional `open_age` retains younger rates and pools deaths/exposure
into the final open band, using two additional native products (D120).
Cutoffs align with source age bands; pooled counts and exposure are retained
in provenance. HFD retains native annual/single-age products. See
[database preparation](database-rates.md). Users can save/load tables with
RDS and reuse them without downloads. HMDHFDplus is optional for these helpers.
The HFD helper's `by_parity = TRUE` option imports conditional mi rates as
age/parity interval rows for parity-specific fertility (milestone 5, accepted under D105). Unavailable
conditional cells are zero-filled under D096, with policy/count metadata. The
engine accepts complete parity partitions under D099, initializes counts
from recorded maternal births and uses one fertility cell per actual parity (D101).
See [parity fertility](parity-fertility.md). The default remains total asfrRR. `calculate_patfr()` summarizes conditional
parity schedules with a deterministic HFD-style fertility table; calculate_tfr()
delegates when parity columns are present, returning patfr (D097). Both products preserve endpoint rows under
D095: 12- becomes [12,13), and 55+ becomes [55,Inf).

There is no automatic disk cache or rate refresh.

`calculate_fertility_rates()`, `calculate_mortality_rates()`, `calculate_tfr()`, `calculate_life_table()`, and
`calculate_population_size()` operate on inspectable outputs. See
[demographic analysis](demographic-analysis.md) for exposure, missing cells,
replications, and age alignment. `calculate_fertility_rates(by = "union_type")`
returns female unpartnered and discovered union-type rows with exposure, births
and fertility rates, without an overall row. The separate mortality helper
returns exposure, deaths and mortality rates for both sexes. Optional
`parity_breaks` on fertility splits observed current parity using birth histories. Untyped histories
use partnered/unpartnered. `calculate_union_rates()` separately returns
unpartnered exposure, union starts, and partnership-entry rates by sex,
excluding conversions.
`by = "previous_union"` adds previous categories; `by = "group"` independently
adds individual groups. Both default to FALSE (D080). Optional `dur_breaks`
automatically separates previous categories and measures separated/widowed rows
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

`extract_unpartnered_births()` extracts known-parent births outside either union
type, adding parent ages/groups for father-matching assessment. Each birth is
retained, including repeated parental pairs. Birth-first dissolution/conversion
boundaries and half-open ordinary coverage match demographic observations.

`augment_unions(population)` returns a separate tibble preserving union
records and appending each partner's birth/death time and optional group at
union formation, reconstructed from transition histories, with `_w` and `_h` suffixes. It joins by partner IDs, preserves types and row order,
and leaves the population unchanged. See [demographic analysis](demographic-analysis.md).

`augment_persons(population)` returns a separate tibble preserving all persons
and their columns, adding `age_mother` and `age_father` in years at each
person's birth, plus `<group>_mother` and `<group>_father` for each declared
group dimension (such as `race_mother`).
Parents are matched by ID, including deceased parents; unknown parents yield
missing fields. Row order and group types are preserved, conflicting column
names are rejected, and no endpoint is required. For example, a person born
in 2000.5 whose mother was born in 1975.25 has `age_mother = 25.25`.
Founders with unknown parents retain their rows with missing parent fields;
an empty population retains the enriched schema. These are recorded identities,
not reconstructions of past group transitions. Validation uses shuffled IDs,
fractional ages, siblings, unknown/deceased parents, optional groups, empty
tables, and unchanged inputs.

`calculate_fertility_rates(by = "group")` and
`calculate_population_size(by = "group")` optionally split observations by
recorded group. The native count kernels retain full histories and add a group
axis to their output arrays, preserving cross-group parent and spouse links.
TFR and life-table calculations automatically retain input group strata. See
[groups](groups.md) for attribution, ordering, and empty-output rules (D053).

## Extension questions

These remain proposals or unresolved requirements:

- Custom process execution and permitted state mutations/invalidation.
- Poor-pool waiting and previous-father continuity; see [union formation](union-formation-models.md)
  and the [roadmap](roadmap.md).
- Heterogeneity, parity-specific fertility, and user-defined newborn attributes.
- Observation retention controls and exact stochastic checkpoints.
- Calendar-varying rate representations and additional database adapters.

Do not infer future demographic rules from SOCSIM or implement arbitrary
callbacks that can bypass consistency checks.

## Exact full-pool partner draws

There is no candidate sample (D072). `MalePool` in `src/matching.h` keeps
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
`WeightCalibrator` in `src/weight_calibration.h` keeps unpartnered counts by group and
birth-time bin and incrementally maintained pool sums for each matching model's
table of utilities over bin-centre lags, updated at events, and solves the
weights at each search from those sums. Female contributions come from the
cells (D082/D083/D128): `SearchCalibration` keeps per search cell and
calibration bin the sums of w(1-p) and w p, recomputing changed cells before
each solve; `FertilityCells` keeps unpartnered eligible weight per bin; and
`ProbabilityTracker` keeps each unpartnered woman's cohabitation probability
and its age boundaries. Searches do not scan all women or enumerate a
cross-classification of rate predictors. See
[current-hazard contributions](calibration-intensity.md). `BoundaryQueue` is a
bucket queue (1/365-year buckets over a 16-year ring, plus an overflow heap)
ordered by (time, unit). See [cell-accounted fertility](fertility-cell-accounting.md). Matching models carry
`delta` (omitted groups default to zero). One group gives unit weights and female
choice exactly. See [weighted-choice matching](weighted-choice-matching.md).

## Hard group restrictions

Matching models prohibit pairs through NA in `group_exogamy` (D135): all
exogamy in a dimension (`NA`) or selected symmetric pairs (`"B:C" = NA`) are
ineligible, alongside age and kin restrictions. Compiled
specifications carry either a same-group flag or blocked-pair matrix; the
exact draw skips restricted cells. Partner/father searches use their applicable
model's restrictions. See [group matching](group-matching.md).

## Milestone 2 construction

`matching_model()` accepts a formula with age and group predictors,
named utility coefficients, and optional symmetric additive pair effects (`A:B` in both directions).
Group predictors are named by dimension (`female_race`, `same_race`) with the
dimension's categories as levels (D134); see
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
Rate mode adds a search process. Rate mode selects the current status fertility schedule; an overall schedule
applies to both statuses with one warning per simulation. On-demand requires
overall fertility. Birth spacing persists across status transitions.
Failed searches leave the woman in her search cell; age-band gaps have zero hazard.

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

Omitting `matching` uses `age_gap_matching_model(0, 5)`. The helper
converts a preferred age gap and half-weight distance in years to the existing
quadratic utility model. `diagnostics = TRUE` adds weighted choice's applied
`.log_weight_<group>` columns to new union records in grouped populations,
without consuming random draws. Diagnostics default off; see the
[union specification](union-formation-models.md#optional-union-diagnostics).


## Typed partnerships (D086)

Optional cohab_prob on union formation enables cohab/marriage histories in rate
mode. `cohab_prob` optionally includes flat half-open calendar bands (D122);
active probability schedules update before same-time searches. Probability-only
changes refresh marriage/cohabitation contributions without touching clocks,
while rate changes follow the D118/D128 rules. Cohab matching and
conversion are part of that process; dissolution
optionally stratifies by type. See [cohabitation](cohabitation.md) for complete
inputs, chronology, linked episodes, legal-history recovery and observations.
Fertility input strata are now union_type (the old input union_status errors).
Calibration has marriage-search, father and cohab-search channels;
preference tables may be shared, but role weights use each role's intensity.


## Persistent fertility weights (D102)

Hetfert is internal (D176): users cannot enable it, and fertility always uses
the uniform path unless the internal `enable_het_fert()` is called.
`het_fert(cv, alpha)` describes mean-one lognormal founder draws and maternal
inheritance on the log scale. Without it, fertility keeps the uniform path. Enabled simulations initialize missing weight columns after seed
setup, preserve supplied positive female weights, and return weights for all
women in `persons$hetfert`. Sons receive NA. Zero-year runs initialize weights
without entering the engine. Existing weights survive enabled restarts. Disabled output omits the weight
column, including when it was supplied, so re-enabling draws fresh weights.

`FertilityWeightTree` stores eligible weights by swap-remove slot with geometric
capacity growth and logarithmic updates/draws. Cell intensity uses member count,
not weights. Matching calibration aggregates the same individual hazards using
per-bin compensated weight sums, and internal scan checks verify their totals.
No tree is allocated on the uniform path. Daughter weights use independent
normal innovations, except exact-copy or zero-innovation endpoints. No
per-generation marginal normalization is imposed. See [the specification](fertility-heterogeneity-notes.md).

## Weight models (D129)

`weight_model(formula, coefficients)` or
`weight_model(fun = f)` builds a weight model, validated by
`validate_weight_model()`; `weights` arguments on `set_mortality()`,
`set_fertility()`, `set_union_dissolution()`, `set_group_transition()` and
`set_union_formation()` (with `weights_conversion`) accept one model or a
`time_start`/`time_end`/`model` schedule. `prepare_weight_model()` decodes a
formula into variable codes and coefficients for the engine, and simulation
warns about terms that cancel within every cell. The engine stores formula
weights at each refresh and selects by rejection sampling against a per-cell
weight bound; custom functions are called for the chosen cell at selection,
except for fertility and search under matching calibration, where stored
weights are needed. See [weight models](binned-group-hazards.md#weight-models-d129).

## Individual group transitions

`group_transition(rates, group)` and `set_group_transition(model, rates, group)`
configure optional competing person-level moves in one group dimension, one
process per dimension (D134). Required direction columns are `origin` and
`destination`; one age or duration band type applies across the table,
optionally stratified by `sex`, `union_type` and other declared dimensions. The transition
registry includes initially empty groups and all other grouped specifications
must cover it. Native events update pools, calibration, fertility cells and
union cells affected by either partner's identity. Duration uses latest entry,
or age without history; ancestry is retained. No household moves or target
adjustments are implemented.

Grouped populations always include a `transitions` tibble, including empty
histories. Restart retains the history and reconstructs duration. Observation
kernels use private indexed history adapters to split group exposure at moves
and attribute events at their occurrence. See [the contract](group-transitions.md).
