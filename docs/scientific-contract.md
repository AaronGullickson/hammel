# Scientific contract

## Baseline population generator

`create_baseline_population(size, max_age = 70, end_time = 0, female_probability = 0.5)` creates a simple starter population inspired by the user's SOCSIM initialization example. All persons are alive at `end_time`, have IDs `1` through `size`, unknown parents, and no unions. Sex is independently sampled with probability `female_probability` of code 1 (female), otherwise code 0 (male); this does not enforce an exact sex balance.

Ages are independently uniform on `[0, max_age]` years and birth times are `end_time - age`. Unlike the example's integer-month birth dates, the generator uses continuous year-valued ages, with no hard-coded increment or rounding. `max_age = 0` produces persons born at the reference time. An empty population is permitted with `size = 0`. Size must be a nonnegative whole number representable as an R integer, maximum age must be finite and nonnegative, the reference time must be finite, and the sex probability must lie in `[0, 1]`.

The generator uses R's current RNG state. `set.seed()` makes generation reproducible; generation advances that state normally. It does not set or restore the user's seed itself.

For example, `size = 3`, `max_age = 0`, `end_time = 2000`, and `female_probability = 1` gives IDs 1–3, sex code 1, birth times 2000, missing death and parent fields, and an empty union table. Validation includes this known-answer case, reproducibility, empty populations, bounds, schema consistency, and distribution checks on a fixed-seed sample.

This population is not claimed to be stable, stationary, or representative of a real population. Burn-in and realistic initial age or kinship structures require separate specifications. Parameter defaults follow the user's example (maximum age around 70 and equal-probability sexes); the continuous-year sampling is an explicit simplification for this helper.

## Output reused as input

D010 requires that a result can initialize further simulation. The eventual contract must specify endpoint time, the living population at that endpoint, retention of dead persons and past unions, unique ID allocation, and continuity of dates and relationships. Validate that subsequent births and unions refer to valid persons and do not reuse existing IDs, and that deceased persons do not become eligible for future events.

Exact fixed-seed equality between split and uninterrupted runs is not yet a requirement. Define restart RNG and pending-event behavior before making such a claim. Existing demographic histories must remain available for initialization even if a particular output observation window is restricted.

Status: requirements for specifying and validating behavior. No mortality, fertility, partnership, or scheduling algorithm has yet been accepted.

## General requirements

Every consequential mechanism needs a plain-language definition and worked examples before implementation. Explicitly distinguish hazards, probabilities, and central rates. State units, interval definitions, conversions, and interpolation assumptions. Do not apply a conversion merely because it is conventional: identify the assumptions that make it valid.

Specifications must define eligibility, event timing, competing-event semantics, state transitions, and consequences for other pending events. Record intentional differences from SOCSIM and their scientific implications.

## Mechanism specification template

Use this structure for each proposed mechanism:

1. **Status and scope:** proposed or accepted; behavior covered and excluded.
2. **Inputs:** schema, units, valid ranges, missing-value policy, and provenance.
3. **Eligibility:** population at risk and dependence on age, history, relationships, or attributes.
4. **Timing:** rate or probability interpretation, scheduling, interval boundaries, ties, and competition.
5. **Transition:** person and relationship changes; dependent events invalidated or rescheduled.
6. **Observations:** emitted events, recorded attributes, and snapshot consequences.
7. **Worked examples:** ordinary cases and boundary cases with expected outcomes.
8. **Validation:** analytical expectations, accounting and relationship invariants, and reference comparisons.
9. **Open questions:** decisions needed from the user before dependent implementation.

## Foundation requiring specification

The accepted first milestone is a population simulation with fertility and mortality rates constant over calendar time, basic partnerships for assigning parents, and demographic validation of the outputs. Rates may vary with age; their exact representation remains open. Implementing mortality first is a possible incremental step, not the full milestone.

Before implementation, resolve input requirements, time semantics, rate interpretation, competing events, observation timing, and RNG behavior. Fertility requires explicit eligibility, newborn initialization, and parenthood rules. Specify a minimal partnership mechanism, including age eligibility, pairing constraints, what happens when no partner is available, and consequences of a partner's death. Random pairing after a minimum age is a candidate, not a settled rule.

## First milestone demographic checks

These are validation requirements; the estimators, reference cases, run sizes, and numerical acceptance thresholds still need specification.

- **Fertility:** record births and exposure by the age of the person giving birth. Estimate age-specific fertility rates using the same denominator and age intervals as the target inputs. Calculate TFR by summing age-specific fertility rates times their interval widths, with consistent units.
- **Fertility eligibility:** distinguish population-based fertility rates from rates conditional on partnership or another eligibility criterion. Applying population-based rates only to partnered people generally changes the resulting population TFR. Resolve this explicitly before implementing the partnership–fertility interaction.
- **Mortality:** record deaths and person-time exposure by age, estimate the corresponding age-specific death rates, and construct life tables using a documented method consistent with the input interpretation. Compare survival and life expectancy with reference quantities derived from those inputs. Do not equate a central death rate, hazard, and interval death probability without the required assumptions.
- **Population structure:** constant input rates do not imply an immediately stationary or stable population. Distinguish period measures from cohort measures and specify initialization, observation windows, and any burn-in. Raw birth or death counts depend on population composition and are not sufficient rate checks.
- **Sampling uncertainty:** predefine replication strategy and statistical criteria. Simulated TFR and life-table estimates fluctuate; deterministic equality to inputs is not an appropriate general acceptance rule. Specify treatment of sparse ages, zero exposure, incomplete cohorts, and the oldest age interval.
- **Mechanism checks:** verify birth and death accounting, parent identities, partnership consistency, and eligibility at event time. Include known-answer mortality and fertility cases alongside combined population simulations so failures can be localized.

## Validation principles

- Derive checks from the scientific specification, rather than reproducing implementation logic.
- Use analytical cases with known answers, including deterministic boundary cases where applicable.
- Check population accounting: for a closed population, living population changes equal births minus deaths. Extend the identity explicitly when migration or other entry and exit mechanisms are introduced.
- Check valid event dates, person identities, and relationship consistency under the chosen model.
- Verify reproducibility under the declared seed contract and independence between simulations.
- Compare demographic outcomes against a pinned SOCSIM reference over multiple replications when relevant. Define comparison measures and acceptance criteria before interpreting results.
- Require exact fixed-seed comparisons only when algorithms and random-number consumption are intentionally preserved. For changed algorithms, use statistical comparisons with justified uncertainty and tolerances.

Passing software checks does not establish research readiness. Before research use, seek targeted review of scheduling, RNG handling, and the R/native boundary by an experienced simulation developer, alongside demographic review.

## Monthly time decision

The first milestone uses a monthly simulation increment, following SOCSIM's convention, with its duration explicitly configured as `dt = 1/12` years. Record times and durations in years and exposure in person-years. Do not hard-code month assumptions in demographic algorithms or output schemas. Other increments are a longer-term feature, but time handling must be isolated from the outset. Age and event-date conventions, exposure accounting, within-month ordering, and conversion of input rates to monthly event probabilities or schedules still require specification. A monthly time unit does not imply simply dividing annual probabilities by twelve.

## Person-year input convention

Accepted inputs include age-specific fertility and mortality rates expressed as events per person-year. Convert their exposure units internally: a rate of 0.12 events per person-year is 0.01 events per person-month. This is a unit conversion, not a statement that the monthly event probability is 0.01.

For example, if the input represents a constant mortality hazard of 0.12 per year, the one-month death probability is `1 - exp(-0.12 / 12)`, approximately 0.00995. A central death rate requires a stated relationship to the underlying mortality model; an annual death probability requires a different conversion. These interpretations must be specified before implementation.

For fertility, define whether the denominator is all women or an eligible subset, and how birth events are generated from the converted rate. Validate rates by dividing event counts by exposure in the corresponding units. Include a unit-consistency check: annualized monthly rates equal the original person-year rates, and equivalent exposure units produce equivalent estimates.

## Accepted rate-to-probability direction

Derive a single-process event probability as `1 - exp(-rate * dt)`, with time units matching the rate denominator. Person-year rates use `dt = 1/12` for a month. This assumes a constant hazard over the span; the numerically stable implementation is `-expm1(-rate * dt)`.

The intended mortality model interprets supplied rates as piecewise constant hazards. When users supply empirical central death rates (`nmx`), interpreting them as hazards is an explicit modeling assumption, not a universal life-table identity. Derive validation expectations consistently with this assumption and identify any differences from conventional life-table conversions. Boundary-crossing rules remain open.

This decision resolves the earlier open single-process conversion direction. It does not specify competition between processes, fertility eligibility or birth spacing, or within-month event ordering. For repeatable events, the probability of at least one event and the expected event count are distinct quantities; fertility event generation still needs its own specification.

## Accepted hazard-scale waiting-time direction

Sample candidate event times directly using cumulative hazard under piecewise constant hazards:

1. Draw `U` uniformly on `(0, 1)` and set the positive threshold `E = -log(U)`.
2. For the remaining duration `d` in the current age band, calculate `h = rate * d`, using matching rate and time units.
3. If `E > h`, subtract `h` from `E` and continue into the next band with its rate.
4. Otherwise, for a positive rate, the event lies `E / rate` time units into the current band.

A zero-rate band contributes no hazard and is crossed without dividing by zero. If the specified future schedule never accumulates enough hazard, no event occurs within that schedule; do not invent a terminal event. Exact boundary assignment and terminal-age policy still require specification.

For person-year rates, durations are in years and a month contributes `rate / 12`. Survival across bands is `exp(-sum(rate * duration))`. This is equivalent to accumulating `log(1 - p)` when `p = 1 - exp(-rate * duration)`, without a probability conversion round trip. The known future age schedule is incorporated at sampling time, so crossing an age band alone does not require a new random draw.

Worked example: the current rate is 0.12 per year for two remaining months, then 0.24 per year. With threshold `E = 0.05`, the first band contributes `0.12 * 2/12 = 0.02`. The remaining threshold is 0.03, reached after `0.03 / 0.24 = 0.125` years (1.5 months) in the next band. The candidate continuous waiting time is 3.5 months. Under D008, truncate 3.5 months to an offset of 3 simulation months.

Validate constant-hazard survival analytically and test piecewise schedules against `exp(-sum(rate * duration))`, including zero-rate bands and events on either side of an age boundary. D008 specifies monthly truncation and candidate tie priority; their boundary cases need worked examples and validation. Other competing-event behavior remains open.

## Accepted monthly truncation and candidate ties

Under D008, truncate nonnegative waiting times in months to integer offsets, with offset zero indicating the current month. Retain one next event per person. Equal-month candidates are prioritized as death, birth, partnership formation, then partnership dissolution. Compare the truncated dates, not fractional times. This is a per-person selection rule; execution order across people and resampling after events remain open.

Worked example: birth at offset 2.1 months and death at offset 2.9 months both map to offset 2; death wins. Birth at 1.9 months and death at 2.1 months map to offsets 1 and 2; birth wins. Candidate offsets below one month map to the current month. Test exact integer boundaries and age-band boundaries explicitly, along with these tie cases. Document differences from continuous-time event competition in demographic validation.

## Time-unit independence requirement

Canonical public time and duration units are years; exposure is recorded in person-years and rate inputs are per person-year. The first milestone's step duration is `dt = 1/12` years. Algorithms must consume the configured step duration rather than assume twelve steps per year. Integer scheduling indices, if used, are generic step indices.

D008's truncation is expressed generally as `step_offset = floor(waiting_time_years / dt)`; its corresponding public duration is `step_offset * dt` years. This preserves the accepted monthly behavior when `dt = 1/12`, without encoding months throughout the engine. Equal-step candidates follow the accepted tie priority. Numerical boundary policy remains to be specified.

Age boundaries and exposure accumulation must use elapsed years. Observation frequency must be explicitly expressed or converted in the time-handling layer. Validate unit consistency and review for scattered month constants; support for additional increments need not be implemented in the first milestone.
