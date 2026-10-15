# Design and implemented interface

Updated October 13, 2026. This document describes implemented architecture.
Future extension ideas are listed separately and do not authorize new behavior.
The [scientific contract](scientific-contract.md) defines demographic semantics.

## Responsibilities and source locations

| Responsibility | Implementation |
| --- | --- |
| Model construction and population attachment | `R/hammel_model.R` |
| Population schema, validation, baseline generation | `R/hammel_population.R`, `R/checks.R` |
| Mortality/fertility specification and rate validation | `R/set_process.R` |
| Union process construction and formula matching | `R/union_formation.R` |
| HMD/HFD downloads and preparation | `R/construct_demographic_rates.R` |
| Simulation validation, RNG isolation, result assembly | `R/simulate.R` |
| Continuous scheduling and events | `src/simulate.cpp` |
| Demographic extraction and measures | `src/demography.cpp`, `R/calc_demographic_measures.R` |
| Software and analytical validation | `tests/testthat/` |
| Replication assessment | `inst/examples/milestone1_bootstrap.R`, `validation/milestone1.qmd` |

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

Setters validate built-in specifications again, including edited objects.
Mortality constructors permit partial nonoverlapping schedules, but simulation
requires complete coverage. Fertility bands must be finite; gaps have zero risk
at execution. Time-varying rate columns are rejected. Arbitrary named process
specifications can be stored; execution supports mortality/fertility and the
conditional-logit on-demand union process. Other processes are rejected.

## Execution and output

`simulate(model, years, seed = NULL)` requires population, mortality,
fertility, and an explicit union-formation process. It uses continuous years, with no `dt`, execution step, monthly
compensation grid, or birth-spacing argument on the execution function.

C++ constructs midpoint-adjusted fertility bands, samples cumulative-hazard
thresholds, and processes the earliest event in a sparse ordered queue.
R validates inputs and assembles the returned `hammel_population`, retaining
histories and existing extra person columns. Newborn custom attributes are
missing. The `simulation` component records `start_time`, `end_time`,
`birth_interval`, `births`, and `deaths` for the latest run.

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
replications, and age alignment.

## Extension questions

These remain proposals or unresolved requirements:

- Custom process execution and permitted state mutations/invalidation.
- Union initiation, partner eligibility, matching, dissolution, and fertility
  denominators; see [union-formation notes](union-formation-models.md).
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
years). Other pair restrictions are not yet applied. Only age difference
(male minus female, continuous years) is supported initially. See the
[union specification](union-formation-models.md) for worked examples and
unresolved execution rules. `simulate()` rejects rate-based initiation.

The engine always builds a call-local current-partner lookup for women once
from ongoing union histories. Per-person union links allow death to close all
affected unions. Formation/death update these caches immediately; restart
reconstructs them from histories. This gives constant-time partner lookup at
birth without adding redundant public `persons` columns. Multiple ongoing
unions for a woman are unsupported; male exclusivity is deliberately not
enforced until eligibility screening is implemented.

There is no union-enable flag or no-union execution path. Rate-based initiation
remains rejected until implemented; its eventual implementation will share the
same partner and union state.

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
matching. The male vector still contains all living men; age screening does not
change the sampling frame or trigger refilling. Existing partner caches bypass
new-search restrictions.
