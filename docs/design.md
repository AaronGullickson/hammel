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
| Time model | Monthly or continuous time; calendar representation; age calculation; boundary rules |
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
