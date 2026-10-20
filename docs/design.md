# Design and implemented interface

Updated October 19, 2026. Milestone 2 was accepted as complete by the user
on this date (D050). This document describes implemented architecture.
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
| Union process construction and formula matching | `R/union_formation.R` |
| Shared kinship traversal and extraction | `src/kinship.h`, `src/kinship.cpp`, `R/kinship.R` |
| HMD/HFD downloads and preparation | `R/construct_demographic_rates.R` |
| Simulation validation, RNG isolation, result assembly | `R/simulate.R` |
| Continuous scheduling and events | `src/simulate.cpp` |
| Demographic extraction and measures | `src/demography.cpp`, `R/calc_demographic_measures.R` |
| Software and analytical validation | `tests/testthat/` |
| Replication assessment | `validation/milestone1.qmd`, `validation/milestone2.qmd` |

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
Fertility also stores `birth_interval`; source attributes are copied into the
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
conditional-logit on-demand and rate-based union processes. Other processes are rejected.

## Execution and output

`simulate(model, years, seed = NULL)` requires population, mortality,
fertility, and an explicit union-formation process. It uses continuous years, with no `dt`, execution step, monthly
compensation grid, or birth-spacing argument on the execution function.

C++ constructs midpoint-adjusted fertility bands, samples cumulative-hazard
thresholds, and processes the earliest event in a sparse ordered queue.
R validates inputs and assembles the returned `hammel_population`, retaining
histories and existing extra person columns. Newborns inherit the mother's optional `group`; other custom attributes are
missing. Group-specific schedules and baseline group sampling are documented in
[groups](groups.md). The `simulation` component records `start_time`, `end_time`,
`birth_interval`, `births`, and `deaths` for the latest run.

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
resamples candidates and is not guaranteed to match an uninterrupted fixed-seed
run. Fertility adjustment after changing rates ignores the previous schedule;
that transition issue is [documented and deferred](scientific-contract.md#known-limitations-and-deferred-issues).

## Data and observation interfaces

`get_mortality_hmd(country, year, ...)` and `get_fertility_hfd(country, year, ...)`
return rate tibbles, without configuring a model. Users can save/load them with
RDS and reuse them without downloads. HMDHFDplus is optional for these helpers.
There is no automatic disk cache or rate refresh.

`calculate_lexis_counts()`, `calculate_tfr()`, `calculate_life_table()`, and
`calculate_population_size()` operate on inspectable outputs. See
[demographic analysis](demographic-analysis.md) for exposure, missing cells,
replications, and age alignment. `calculate_lexis_counts(include_unions = TRUE)`
returns overall/partnered/unpartnered rows with exposure, births, deaths, and
rates. Union counts and marriage rates appear only on unpartnered rows,
using the recorded union histories (D046/D047).

`calculate_lexis_counts(by_group = TRUE)` and
`calculate_population_size(by_group = TRUE)` optionally split observations by
recorded group. The native count kernels retain full histories and add a group
axis to their output arrays, preserving cross-group parent and spouse links.
TFR and life-table calculations automatically retain input group strata. See
[groups](groups.md) for attribution, ordering, and empty-output rules (D053).

## Extension questions

These remain proposals or unresolved requirements:

- Custom process execution and permitted state mutations/invalidation.
- Non-death union dissolution, poor-pool waiting, previous-father continuity,
  and transition-aware fertility compensation; see [union formation](union-formation-models.md)
  and the [roadmap](roadmap.md).
- Heterogeneity, parity-specific fertility, and user-defined newborn attributes.
- Observation retention controls and exact stochastic checkpoints.
- Calendar-varying rate representations and additional database adapters.

Do not infer future demographic rules from SOCSIM or implement arbitrary
callbacks that can bypass consistency checks.

## Milestone 2 construction

`create_matching_model()` accepts an age-difference formula with named utility
coefficients. `union_formation()` and `set_union_formation()` specify initiation,
candidate sampling, and a matching specification. The current C++ on-demand
increment creates unions at births using conditional-logit selection from
sampled living men after minimum-male-age screening (default 15 continuous
years) and known-kin screening (default two generations on each side of a
common ancestor). Only living, unpartnered men enter the sampling pool. Only age difference
(male minus female, continuous years) is supported initially. See the
[union specification](union-formation-models.md) for worked examples and
deferred extensions. Rate-based initiation schedules searches for unpartnered women.

The engine always builds a call-local current-partner lookup for both sexes once
from ongoing union histories. Per-person union links allow death to close all
affected unions. Formation/death update these caches immediately; restart
reconstructs them from histories. This gives constant-time partner lookup at
birth without adding redundant public `persons` columns. Multiple ongoing
unions for either sex are rejected. Formation removes the selected man from
the pool; widowhood returns a surviving husband. Pool updates use constant-time
swap removal with per-person slot tracking.

There is no union-enable flag or no-union execution path. Both initiation modes share the same partner and union state and search helper.
Rate mode adds a search event type and per-person event versions to invalidate
queued births after widowhood. Rate mode selects the current status fertility schedule; an overall schedule
applies to both statuses with one warning per simulation. On-demand requires
overall fertility. Birth spacing persists across status transitions.
Failed searches reschedule at attained age; age-band gaps have zero hazard.

Candidate sampling uses a reusable array of male-pool positions and a partial
Fisher–Yates shuffle, undoing swaps after each sample. This avoids repeated
hash allocations while retaining distinct sampled candidates and the existing
RNG draw sequence. Workspace is local to each simulation; it tracks positions
in the current male vector rather than person identity.

Matching uses the entered formula and named coefficients. R prepares an exact
native linear/raw-quadratic evaluator where applicable; other formulas retain
`model.matrix()` evaluation through an R callback receiving age differences
only. C++ computes stable softmax probabilities and samples one candidate. No
formula callback receives mutable population state, and no candidate sets are
retained. The native path avoids an R model-matrix call at every search for
the common linear-plus-squared specification.

`min_male_age` is stored in the union process and passed to native execution.
A screening helper filters the temporary sample at the exact search time before
matching. The male vector still contains all living, unpartnered men; age screening does not
change the sampling frame or trigger refilling. Existing partner caches bypass
new-search restrictions.

The shared `Kinship` module resolves parent IDs to positions, rejects cycles,
and lazily caches distinct ancestors to the configured generation bound.
Screening marks one woman's ancestors and checks candidate membership directly
in C++, skipping these queries at zero generations; R's `extract_kin()`
uses the same kernel and an inverted ancestor index for selected relatives.
Shortest-distance lineage and within-depth completeness accompany extraction
rows. No all-pairs kin matrix or full ancestral-path expansion is constructed.
Newborn append leaves existing ancestry caches valid; restart rebuilds them.


Nonmarital rate-mode births and union searches share a `find_candidate` helper
covering sampling, screening and matching. Only `form_union` updates partner
caches and removes the man from the available pool. Parentage alone leaves
both parties unpartnered; the newborn's father index is appended to kinship
independently of current-partner state. No candidate sets or additional
persistence state are retained. `father_matching = NULL` resolves to the general model at execution; an explicit
model supplies a separately prepared scorer only for nonmarital father selection.
Both scorers use the native evaluator or R formula fallback. Previous-father
preference is deferred (D042).

Omitting `matching` uses `create_age_gap_matching_model(0, 5)`. The helper
converts a preferred age gap and half-weight distance in years to the existing
quadratic utility model. `diagnostics = TRUE` adds `.n_eligible` and
`.utility_rank` to new union records without retaining candidate sets or
consuming additional random draws. Diagnostics default off; see the
[union specification](union-formation-models.md#optional-union-diagnostics).
