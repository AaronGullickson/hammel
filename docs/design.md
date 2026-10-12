# Design specification

Status: initial proposal. The project direction in the brief is accepted; the architecture and API below are candidates requiring further specification.

## Candidate architecture

Separate four responsibilities:

1. **R model specification:** explicit objects for populations, relationships, processes, parameters, and observation requests. Constructors compose without executing a simulation.
2. **Simulation engine:** owns population state, relationships, scheduling, and random-number generation independently for each run. A compiled engine is proposed; its language is undecided.
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

Simulation outputs must also be usable as inputs to further simulation (D010). The pipeable workflow should support initializing a subsequent model or run from a `hammel_results`, retaining person and union identities and histories. The result must provide enough endpoint information to establish the new starting time and current population. Subsequent runs can configure demographic mechanisms; the syntax and configuration carryover rules remain open.

Distinguish reuse of demographic output from exact checkpoint continuation. Reusing output must not inherently depend on retaining a live compiled engine. Whether pending events and RNG state are saved, and whether split runs reproduce an uninterrupted run, requires a separate decision.

Simulation results use a named-list object of class `hammel_results`. The core tables are `persons` (one row per person, including deceased persons) and `unions` (one row per union episode, linked by person IDs). Access them as `result$persons` and `result$unions`. Column schemas, data-frame implementation, optional observations and exposure summaries, metadata, and methods remain to be specified. See D009.

Model construction should support R pipes, with simulation execution as a separate operation. D011 establishes the basic configuration API below. Simulation execution and observation function names remain proposals.

Before implementing the first simulator, specify a minimal end-to-end example with input schemas, time and rate units, observation requests, returned objects, and reproducibility behavior.

## First milestone scope

The accepted first milestone combines calendar-time-constant fertility and mortality inputs with basic partnerships and parent assignment. Build it incrementally, but evaluate completion against the combined demographic validation targets in the scientific contract. Rich partnership matching is deferred.

The minimum observation design must support age-specific births, deaths, and the corresponding person-time exposures, as well as accounting and relationship checks. Specify these observations before deciding the result schema. Basic partnership rules must be compatible with the intended fertility denominators and target rates.

## Time unit

Use a fixed monthly time unit for the first milestone. Keep unit assumptions explicit and conversions in defined operations to support future configurable units. This does not decide whether the engine scans the population monthly or schedules events ahead, or how events within a month are ordered.

## Rate input interface

Accept age-specific fertility and mortality rates expressed as events per person-year. Convert to events per person-month internally by dividing by twelve, keeping the input units and denominator explicit. Users should not need to construct monthly rate files. Output rate estimates used for validation should be available in the corresponding person-year units.

Keep exposure-unit conversion separate from the mapping of rates to event probabilities or schedules. The latter depends on the rate interpretation and the accepted demographic algorithm.

## Rate-to-probability direction

Use `1 - exp(-rate * dt)` for a single-process event probability under a constant hazard within the time span, with `dt = 1/12` for person-year rates and monthly time. This resolves the previously open conversion direction; rate boundaries, process competition, and fertility event generation remain to be specified. See D006 and the scientific contract for the mortality interpretation and validation assumptions.

## Candidate waiting-time sampling

Sample directly on the hazard scale: draw `E = -log(U)` and consume that threshold through successive age bands using `rate * duration`. In the band containing the event, divide the remaining threshold by its rate to obtain the remaining waiting time. Known age-specific rate changes are thus included in the original sample. See D007 and the scientific contract for units, a worked example, and validation requirements. D008 resolves monthly event assignment and candidate ties: truncate waiting times to integer months and prioritize death, birth, formation, then dissolution for equal-month candidates. Execution order across people and resampling remain open.

## Basic configuration API

For a person-table input, `set_population(model, persons)` and explicit
`unions = NULL` both initialize an unmarried population. The model stores an
empty union data frame. This does not infer or alter person attributes;
demographic schema consistency checks remain to be specified.

Implemented under D011:

```r
model <- hammel_model() |>
  set_population(persons, unions, start_month = 0) |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)

restarted_model <- hammel_model() |>
  set_population(results) |>
  set_mortality(mortality_rates)
```

`results` is a `hammel_results` object. Its `end_month` becomes the model's `start_month`. These optional whole-number indices do not yet establish calendar or exposure boundary conventions. Unknown times remain `NULL`; setters do not infer time from the latest event.

The common operation is `set_process(model, name, process)`. The demographic shortcuts use names `mortality` and `fertility`. Setting the same name replaces its specification and retains other processes. Specifications are stored unchanged, including their person-year input units. Rate schemas and the executable process contract are deferred.
