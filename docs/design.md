> Current timing update (D024): event execution is continuous, with no truncation or execution step. The monthly SOCSIM-style fertility compensation is retained separately. Earlier monthly-engine descriptions below are superseded where they concern execution. See the continuous-time scientific contract.

# Design specification

Status: initial proposal. The project direction in the brief is accepted; the architecture and API below are candidates requiring further specification.

## Candidate architecture

Separate four responsibilities:

1. **R model specification:** explicit objects for populations, relationships, processes, parameters, and observation requests. Constructors compose without executing a simulation.
2. **Simulation engine:** owns population state, relationships, scheduling, and random-number generation independently for each run. The first mortality/fertility engine uses C++17 and Rcpp with R's RNG.
3. **Process extension contract:** built-in and custom mechanisms declare eligibility, rates or schedules, permitted state changes, dependencies, and invalidation or rescheduling behavior.
4. **Observation system:** selects event histories, snapshots, and summaries to retain so large simulations need not record everything.

Persons, demographic time, and kinship are candidate core concepts. Demographic mechanisms should be replaceable. Evaluate vectorized R hooks and compiled extensions against correctness and performance requirements before choosing either.

## Open decisions

| Question | What the specification must establish |
| --- | --- |
| Time model | Monthly time accepted for the first milestone; calendar representation, age calculation, and boundary rules remain open; configurable units are a longer-term goal |
| Scheduling | Competing events, simultaneous events, ordering, and invalidation after changes |
| Rate representation | Hazards, probabilities, central rates, units, intervals, conversions, interpolation |
| Population state | Required person fields, custom attributes, identity, and dead or departed persons |
| Relationships | Biological and social parenthood, partnership histories, validity constraints |
| Processes | Eligibility, heterogeneity, mutation operations, dependencies, custom event types |
| Random numbers | Seed behavior, independent runs, R RNG interaction, replication strategy |
| Observations | Snapshot timing, event schemas, retention, and result structure |
| Initialization | Supplied populations, history requirements, missing information, burn-in |
| Implementation | Language, R/native boundary, dependencies, performance targets |

Do not infer answers from SOCSIM's implementation: document reference behavior and make an explicit Hammel choice.

## Interface direction

The accepted initial table schema is in [population schema](population-schema.md) (D012). Person and union dates are exposed in numeric years with fractional years; the first engine retains integer monthly scheduling internally. `persons` has no marital-status column. Relationship state is derived from union histories and may be cached internally for performance. Population metadata uses `end_time` in years; `set_population()` accepts only a population object and retains its reference time.

Simulation outputs must also be usable as inputs to further simulation (D010). The pipeable workflow should support initializing a subsequent model or run from a `hammel_population`, retaining person and union identities and histories. The result must provide enough endpoint information to establish the new starting time and current population. Subsequent runs can configure demographic mechanisms; the syntax and configuration carryover rules remain open.

Distinguish reuse of demographic output from exact checkpoint continuation. Reusing output must not inherently depend on retaining a live compiled engine. Whether pending events and RNG state are saved, and whether split runs reproduce an uninterrupted run, requires a separate decision.

Simulation results use a named-list object of class `hammel_population`. The core tables are `persons` (one row per person, including deceased persons) and `unions` (one row per union episode, linked by person IDs). Access them as `result$persons` and `result$unions`. Column schemas, data-frame implementation, optional observations and exposure summaries, metadata, and methods remain to be specified. See D009.

Model construction should support R pipes, with simulation execution as a separate operation. D011 establishes the basic configuration API below. Simulation execution and observation function names remain proposals.

Before implementing the first simulator, specify a minimal end-to-end example with input schemas, time and rate units, observation requests, returned objects, and reproducibility behavior.

## First milestone scope

The accepted first milestone combines calendar-time-constant fertility and mortality inputs with parent assignment. For the first executable version (D019), sample a living male at each birth without forming unions; include minimum birth spacing with SOCSIM-style rate compensation. Build it incrementally, but evaluate milestone completion against the combined demographic validation targets in the scientific contract. Partnership formation and richer matching are deferred.

The minimum observation design must support age-specific births, deaths, and the corresponding person-time exposures, as well as accounting and relationship checks. Specify these observations before deciding the result schema. Basic partnership rules must be compatible with the intended fertility denominators and target rates.

## Time unit

Use a monthly increment for the first milestone, represented by an explicit step duration `dt = 1/12` years. Record public times and durations in years and exposure in person-years. Do not hard-code month assumptions: demographic algorithms use configured durations, and any integer step indexing and conversion belong in centralized time-handling operations. Supporting other increments can come later, but the initial implementation must preserve this separation. The first engine schedules events ahead in a sparse queue and randomly executes pending events within each step (D019).

## Rate input interface

Accept age-specific fertility and mortality rates expressed as events per person-year. Keep rates per person-year and multiply by the relevant duration in years when computing cumulative hazard or interval probabilities. For the initial monthly setting this is equivalent to division by twelve; do not scatter that conversion through the implementation. Users should not need to construct monthly rate files. Output rate estimates used for validation should be available in the corresponding person-year units.

Keep exposure-unit conversion separate from the mapping of rates to event probabilities or schedules. The latter depends on the rate interpretation and the accepted demographic algorithm.

## Rate-to-probability direction

Use `1 - exp(-rate * dt)` for a single-process event probability under a constant hazard within the time span, with `dt` obtained from time configuration (`1/12` years for the initial monthly setting). This resolves the previously open conversion direction; rate boundaries, process competition, and fertility event generation remain to be specified. See D006 and the scientific contract for the mortality interpretation and validation assumptions.

## Candidate waiting-time sampling

Sample directly on the hazard scale: draw `E = -log(U)` and consume that threshold through successive age bands using `rate * duration`. In the band containing the event, divide the remaining threshold by its rate to obtain the remaining waiting time. Known age-specific rate changes are thus included in the original sample. See D007 and the scientific contract for units, a worked example, and validation requirements. D008 resolves monthly event assignment and candidate ties: truncate waiting times to integer months and prioritize death, birth, formation, then dissolution for equal-month candidates. D019 specifies random sequential execution and post-birth resampling for the first engine.

## Basic configuration API

Implemented under D015:

```r
population <- hammel_population(persons, unions, end_time = 0)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)
```

`hammel_population()` normalizes absent unions and validates the tables. `set_population(model, population)` accepts only a `hammel_population`, validates it, and stores it intact as `model$population`. Its `end_time` supplies the next run's starting time. No separate model starting-time field is stored. Additional components are preserved. Unknown reference times remain `NULL`.

The common operation is `set_process(model, process)`, using the name carried by a `hammel_process` object. `mortality(rates)` and `fertility(rates)` validate age intervals and construct those objects. `set_mortality()` and `set_fertility()` are constructor-plus-setter shortcuts. Setting the same process name replaces its specification while retaining others. Rates remain per person-year. Mortality requires separate `m` and `f` schedules. Calendar-time intervals and richer process behavior remain open; see the scientific contract for implemented rate checks and coverage limits.

## First execution interface

`simulate(model, years, seed = NULL, birth_interval = 0.75)` returns an updated `hammel_population`. Required model components are population, mortality, and fertility. The population's `end_time` supplies the starting time. C++ handles cumulative-hazard sampling, adjusted fertility schedules, a sparse event queue, and constant-time random selection from currently living males. R handles configuration, optional seed isolation, custom-column preservation, and result assembly. Engine state is local to each call.

Public units remain years and person-year rates; only explicit defaults express the initial monthly increment and nine-month birth spacing. See the scientific contract for boundary, competition, coverage, and adjustment details. The example in `inst/examples/milestone1.R` runs 50,000 initial persons for 100 years with downloaded constant USA/2000 rates. Scientific accuracy beyond the initial analytical checks remains to be assessed.
