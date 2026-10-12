# Decision log

Keep accepted choices separate from proposals. Only record user agreement as acceptance of a demographic or major architectural decision. Each new entry should include status, date, rationale, alternatives, implications, and supporting evidence where available. Update superseded entries with a pointer to their replacement.

## D001: Standalone rewrite

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** build a standalone R package named Hammel, informed by SOCSIM; do not build a SOCSIM wrapper or rsocsim companion package.
- **Rationale:** modernize the user interface and enable custom demographic mechanisms and methodological innovation beyond the legacy file-oriented workflow.
- **Implications:** investigate SOCSIM algorithms and assumptions; legacy input compatibility is optional.
- **Evidence:** project conversation handoff, summarized in `project-brief.md`.

## D002: Scientific ownership

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** the user guides demographic requirements and scientific interpretation; agents support software investigation, implementation, documentation, and validation.
- **Implications:** expose unresolved demographic assumptions and obtain decisions before implementing dependent behavior.

## D003: First demographic validation milestone

- **Status:** accepted October 12, 2026.
- **Decision:** the first major milestone simulates a population with demographic rates constant over calendar time, including fertility, mortality, and a basic partnership mechanism to assign parents.
- **Rationale:** establish demographic accuracy before investing in richer mechanisms.
- **Validation targets:** output fertility, including TFR, and life-table measures should agree with expected input quantities under explicitly defined assumptions and sampling uncertainty.
- **Scope:** a simple partnership model is sufficient. Random marriage after a minimum age was offered as an example; the pairing algorithm and eligibility rules remain open.
- **Implications:** specify rate denominators, partnership effects on fertility eligibility, observation exposure, initialization, and statistical acceptance criteria before implementing dependent behavior.

## Pending choices

The design specification lists open choices. In particular, C++, any example function names, and the four-part candidate architecture are not accepted decisions.

## D004: Monthly time for the first milestone

- **Status:** accepted October 12, 2026.
- **Decision:** use a fixed month as the simulation time unit for the first milestone, following SOCSIM's monthly convention.
- **Long-term direction:** make time units configurable; supported units and the interface remain unspecified. This flexibility is not required for the first milestone.
- **Implications:** make unit assumptions explicit and centralize conversions, including age, event dates, input rates, and person-time exposure.
- **Still open:** within-month ordering, competing events, boundary conventions, and rate-to-probability or scheduling conversions. Monthly time does not settle the scheduling algorithm.

## D005: Person-year rate inputs

- **Status:** accepted October 12, 2026.
- **Decision:** users can supply age-specific fertility and mortality rates as events per person-year. For the monthly engine, convert these internally to events per person-month by dividing by twelve; users need not prepare monthly rate files.
- **Rationale:** support familiar demographic inputs and remove unnecessary input preparation.
- **Implications:** retain explicit rate units and denominators; report demographic estimates in person-year units when comparing them with these inputs. The factor of twelve treats a simulation year as twelve simulation months.
- **Boundary:** division by twelve converts exposure units, not rates into probabilities. The mapping from fertility intensities or mortality central rates/hazards to simulation events still requires specification. Annual death probabilities are a different input type and cannot use this conversion as probabilities.

## D006: Rate-based event probability conversion
- **Status:** accepted October 12, 2026.
- **Decision:** use rates as the input convention and derive an event probability over a time span `dt` using `1 - exp(-rate * dt)`. For person-year rates and a month, `dt = 1/12`.
- **Assumption:** this expression is exact for a hazard constant over the relevant span. Rates may vary between age bands; handling spans that cross rate boundaries still requires specification.
- **Mortality interpretation:** the intended model treats the supplied mortality rate as a piecewise constant hazard. An empirical life-table central rate `nmx` is not universally identical to that hazard; using it directly is a modeling assumption that must be documented and reflected in validation references.
- **Boundary:** this defines a single-process probability, not the treatment of competing events, fertility eligibility, repeat births, or within-month ordering.
- **Implementation guidance:** use the numerically stable equivalent `-expm1(-rate * dt)` when implementing this conversion.

## D007: Sample waiting times directly on the hazard scale

- **Status:** accepted October 12, 2026.
- **Decision:** use the mathematical equivalence between log survival and cumulative hazard to sample candidate waiting times directly from rates, without converting to probabilities and back to logarithms.
- **Method:** draw `E = -log(U)` for `U` uniform on `(0, 1)`. Accumulate `rate * duration` across successive age bands until the cumulative hazard reaches `E`. Carry the remaining threshold into each subsequent band; within the final positive-rate band, the remaining waiting time is `remaining_threshold / rate`.
- **Units:** rates per person-year use durations in years; one simulation month contributes `rate / 12` to cumulative hazard.
- **Relationship to D006:** this is the intended candidate waiting-time method. D006's probability formula remains the equivalent single-process interval probability and may be used for interpretation and validation.
- **Resolved subsequently:** D008 adopts SOCSIM-style integer truncation and candidate tie priority. Rescheduling after state changes, terminal age handling, and fertility eligibility and repeat-event semantics remain open.

## D008: Integer-month candidate dates and competing-event ties

- **Status:** accepted October 12, 2026.
- **Decision:** adopt SOCSIM-style integer truncation of nonnegative candidate waiting times to months and its priority for competing candidate events in the same month.
- **Monthly assignment:** convert the sampled waiting time to months and truncate its fractional part (equivalently, take the floor). Zero means the current simulation month. Exact age-band boundary conventions still need worked specifications.
- **Tie priority:** for a person's equal-month candidates, select death before birth before marriage/partnership formation before divorce/partnership dissolution. Group transitions follow these; specify their ordering when introduced. Compare integer-month dates, not hidden fractional waiting times.
- **Scope:** retain one selected next event per person. This priority does not establish a global order among different people's events.
- **Implications:** an earlier fractional birth time does not beat a later fractional death time if both truncate to the same month. Validate and document this deliberate discretization effect.
- **Still open:** execution order across people, resampling and same-month follow-up events, checkpoint behavior, and terminal age handling. This decision does not automatically adopt all SOCSIM scheduling behavior.
- **Reference:** the SOCSIM rate and competing-event audits at rsocsim commit `125450cff968bb3ff762fa99d4c926a8136062a2`.
