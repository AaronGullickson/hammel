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

Model construction should support R pipes, with simulation execution as a separate operation. Concrete function names and object classes have not been chosen. The illustrative API in the conversation is not a public API commitment.

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
