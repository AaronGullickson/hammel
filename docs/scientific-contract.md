# Scientific contract

## Database rate helpers

Download one population and explicitly selected year for a constant schedule; do not silently select the latest year or average rates. HMD uses `Mx_1x1`, retaining male and female central rates per person-year under the accepted hazard interpretation. Request raw age labels using `fixup = FALSE`, while respecting HMD OpenInterval flags because HMDHFDplus 2.0.8 still parses its web output; closed single-year ages map to `[age, age + 1)` and an upper open HMD age maps to `[age, Inf)`. Do not drop missing rates or fill gaps. Retain database, population code, year, product, retrieval timestamp, and reader version, without credentials.

HMDHFDplus is optional and required only for download helpers. Verify conversions with synthetic data and the package's local readers; authenticated web requests require user credentials and are not part of ordinary tests. Source interfaces were checked against HMDHFDplus 2.0.8. HFD uses `asfrRR`, retaining total ASFR in births per woman-year without scaling. At the user's direction, omit open-age rows (`12-`, `55+`) and map closed ages to one-year bands. The retained-age TFR excludes endpoint contributions; do not claim full HFD TFR agreement. No fertility rates are inferred outside the retained age range. An optional HFD revision date is passed to the reader and recorded as requested; NULL means the current revision was requested, not that a specific revision has been pinned.

References: [HMDHFDplus source](https://github.com/timriffe/HMDHFDplus), [HFD FAQ](https://www.humanfertility.org/Project/FAQ), and [HFD basic indicator methods](https://www.demogr.mpg.de/papers/technicalreports/tr-2011-001.pdf).

## Baseline population generator

`create_baseline_population(size, max_age = 70, end_time = 0, female_probability = 0.5)` creates a simple starter population inspired by the user's SOCSIM initialization example. All persons are alive at `end_time`, have IDs `1` through `size`, unknown parents, and no unions. Sex is independently sampled with probability `female_probability` of label `"f"` (female), otherwise `"m"` (male); this does not enforce an exact sex balance.

Ages are independently uniform on `[0, max_age]` years and birth times are `end_time - age`. Unlike the example's integer-month birth dates, the generator uses continuous year-valued ages, with no hard-coded increment or rounding. `max_age = 0` produces persons born at the reference time. An empty population is permitted with `size = 0`. Size must be a nonnegative whole number representable as an R integer, maximum age must be finite and nonnegative, the reference time must be finite, and the sex probability must lie in `[0, 1]`.

The generator uses R's current RNG state. `set.seed()` makes generation reproducible; generation advances that state normally. It does not set or restore the user's seed itself.

For example, `size = 3`, `max_age = 0`, `end_time = 2000`, and `female_probability = 1` gives IDs 1–3, sex label `"f"`, birth times 2000, missing death and parent fields, and an empty union table. Validation includes this known-answer case, reproducibility, empty populations, bounds, schema consistency, and distribution checks on a fixed-seed sample.

This population is not claimed to be stable, stationary, or representative of a real population. Burn-in and realistic initial age or kinship structures require separate specifications. Parameter defaults follow the user's example (maximum age around 70 and equal-probability sexes); the continuous-year sampling is an explicit simplification for this helper.

## Output reused as input

D010 requires that a result can initialize further simulation. The eventual contract must specify endpoint time, the living population at that endpoint, retention of dead persons and past unions, unique ID allocation, and continuity of dates and relationships. Validate that subsequent births and unions refer to valid persons and do not reuse existing IDs, and that deceased persons do not become eligible for future events.

Exact fixed-seed equality between split and uninterrupted runs is not yet a requirement. Define restart RNG and pending-event behavior before making such a claim. Existing demographic histories must remain available for initialization even if a particular output observation window is restricted.

Status: the first mortality/fertility engine is specified below; richer mechanisms and full demographic validation remain to be developed.

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

## Initial process specification objects

`mortality(rates)` and `fertility(rates)` construct inspectable `hammel_process` specifications without executing events. Both require nonempty data frames with `age_start`, `age_end`, and `rate`. Ages use years and half-open intervals `[age_start, age_end)`; starts must be finite and nonnegative, ends must be greater than starts and may be positive infinity. Rates must be finite and nonnegative, per person-year. Mortality also requires character `sex` labels and separate schedules for both `m` and `f`.

Reject overlapping intervals within each schedule, including duplicate and nested bands; adjacent bands are valid. Preserve input order and custom columns. Do not infer rates outside supplied intervals or fill gaps; coverage must be resolved before simulation. Reject `time_start` or `time_end` columns in this initial implementation to avoid pretending to support calendar-varying rates. The schedule interface is intended to gain that dimension later.

`set_process(model, process)` stores a process under its own name, replacing an existing process with that name and preserving others. `set_mortality(model, rates)` and `set_fertility(model, rates)` construct the appropriate object and delegate to that setter. Objects contain `name`, `rates`, and `rate_unit = "person_year"`, with a process-specific subclass and `hammel_process` as their common class. Validate built-in schedules again when attaching objects so edited invalid specifications cannot bypass checks.

Known structural cases include adjacent age bands (accepted), overlapping bands (rejected), identical age bands across different sexes (accepted), and rates of zero or greater than one (accepted: these are rates, not probabilities). These constructors do not resolve fertility eligibility, population versus eligible-person denominators, competing-event execution, or state transitions.

## First executable simulator contract

Current user-authorized scope: mortality and fertility with all living females at risk under the age-specific population fertility schedule, random male parentage at each birth, and no new unions. Use C++ for scheduling and event execution with R setup and inspectable population output. The initial setup check requires a population plus mortality and fertility processes. Additional custom process execution is not supported yet.

`simulate(model, years, seed = NULL, step_duration = 1/12, birth_interval = 0.75)` advances from the population's finite `end_time`. Public units are years. Duration and minimum birth interval must be whole multiples of the explicitly configured step duration. Process a half-open run interval `[start, start + years)`; events are assigned to step starts by flooring candidate waiting times in step units. Equal-step death and birth candidates select death. Within a step, select pending persons randomly and execute sequentially, resampling a mother and newborn after birth. Minimum spacing prevents repeated maternal births within the same step when positive. Known age-specific changes are included at scheduling time.

Mortality covers `[0, Inf)` separately for each sex, without gaps. Fertility outside supplied bands is zero. No maximum age forces death. A zero hazard can yield no future event. Histories supply the last recorded maternal birth for spacing on restart. New IDs exceed all supplied IDs; dead persons remain dead. Every birth chooses uniformly from all currently living males, with no minimum age or kinship exclusions in this deliberately basic version; if none exist, father ID is NA and birth still occurs. Newborn sex has equal probabilities. No unions are formed; existing ongoing unions close when a partner dies.

Birth-spacing adjustment follows the audited SOCSIM method, generalized to configured steps. For age step `k`, average the supplied annual rate over that age interval to obtain `f_k`, then form `q_k = 1 - exp(-f_k * dt)`. With spacing of `B` steps, estimate eligibility as `e_k = 1 - sum(q_j)` over the preceding `B` age steps. The eligible annual hazard is `f_k / e_k`. Reject nonpositive eligibility estimates instead of clamping them. Zero spacing leaves rates unchanged. Construct this schedule once to cover ages reachable during the run, not using a fixed maximum age or a month constant.

This is SOCSIM's approximate population-rate compensation, not an exact guarantee of observed ASFRs. Averaging into step-sized age bands is explicit, while mortality retains its supplied bands. Initially unknown maternal histories mean everyone with no recorded birth is treated as eligible; initialization transients and age-structured validation require burn-in. Integer truncation, tie priority, birth spacing, and approximate adjustment can cause departures from input rates, which must be quantified before research use. Fertility event dates can precede the underlying continuous candidate by less than one step, including near age boundaries.

An explicit seed controls R's RNG for this call and restores the caller's RNG state afterward; without a seed, use and advance the current RNG. Each run owns its state; restarting resamples pending candidates and is not an exact checkpoint. Preserve original custom columns for existing persons; newborn custom attributes are missing, pending a user-defined initialization interface.

Validation: analytical constant and piecewise mortality survival, known threshold waiting times, death-before-birth ties, spacing between every mother's recorded births, father availability, unique IDs, population accounting, restart and input immutability, fixed-seed reproducibility, zero-risk and empty-population runs. Test constant fertility after burn-in against the renewal expectation induced by adjusted hazards and integer truncation, separately from eventual full ASFR/TFR and life-table accuracy checks.

## Demographic extraction contract

Accept a `hammel_population` or a persons data frame for one independent history. IDs must be unique within that history; analyze stacked simulations separately. No user-supplied grouping interface is implemented. For populations, require a known endpoint and do not allow observation beyond it. For bare data frames, the caller must ensure complete history through the requested finite observation endpoint; missing death is treated as censored there. Do not pool replicated person IDs implicitly.

`calculate_lexis_counts()` uses explicit increasing `time_breaks` and `age_breaks`, returning person-years, deaths, births, ASMR, and ASFR for each period, sex, and age band. Births and ASFR are NA for males. Calendar intervals and age intervals are half-open. Exposures are exact overlaps of `[birth_time, death_time)` with each period and each person's age band, censored at the observation endpoint. No monthly assumptions, integer-year loops, or birthday rounding are used. Zero-duration lives contribute zero exposure, although recorded deaths remain events. Retain zero-count/zero-exposure cells; rates with zero exposure are NA, not zero or infinity. Exposures use recorded event times; the engine's step-start truncation can therefore affect estimated central rates.

Mortality counts deaths by period, age at death, and sex, then divides by person-years. Fertility counts each child's birth once using the mother's age at that birth, then divides by all female person-years in that band, including postpartum exposure. Unknown maternal IDs cannot yield age-specific birth counts; report their counts in an `unassigned_births` attribute by period. Do not substitute father's age or discard female exposure during birth spacing.

`calculate_tfr()` takes the fertility-rate table and sums `asfr * (age_end - age_start)` within each period, automatically selecting female rows when sex is present. An optional `age_range` selects complete bands from a summary that also contains an open mortality interval. Require finite contiguous age bands. Missing rates yield missing TFR rather than omitting bands. The selected age range is explicit; it is not necessarily full HFD TFR. No factor of five is assumed.

`calculate_life_table()` takes a mortality-rate table (or a standard rate schedule with `rate`) and constructs a period life table independently by sex and period. Require contiguous age bands starting at zero and a final explicitly open interval (`age_end = Inf`); do not silently turn a finite final age into an open band. Interpret central rates as piecewise constant hazards, consistent with the engine; this differs from conventional HMD infant adjustments. For finite width n and rate m, q = 1 - exp(-m*n), survival propagates from the selected radix, and L = l*(1-exp(-m*n))/m. At m=0 use L=l*n. Compute T as reverse cumulative L, and e=T/l, not simply T. Include ax derived consistently from the same hazard model.

For a positive terminal rate, q=1, L=l/m, and ax=1/m. A zero terminal rate represents immortality: q=0, and L and remaining life expectancy are infinite if survival to that band is positive. Missing rates propagate uncertainty rather than being filled with zero; zero survivorship has undefined e (NA). Use stable small-rate calculations. These are period synthetic-cohort tables, not empirical cohort follow-up estimates.

Validation uses known fractional overlaps, birthdays inside periods, event/boundary inclusion, unknown mothers, zero exposure and empty populations, rejection of duplicate IDs, unequal-width TFR, and constant-hazard life tables with e=1/m at every age. Test zero and missing rates and explicit open-interval requirements. The user's SOCSIM analysis code informs the workflow, but interval arithmetic and life-table computation are implemented for Hammel's year-based schema.

Age-break clarification: `calculate_lexis_counts()` appends Inf when age breaks end at a finite value, creating an open band starting at that last breakpoint. Existing explicit Inf breaks are unchanged. This includes exposures and deaths above that age rather than discarding them. Life tables from manually supplied rates must still identify an actual open interval; no terminal rate is fabricated.

TFR clarification: `calculate_tfr()` no longer takes `age_range`. Select female rows when sex is present, omit rows with infinite `age_end`, and sum ASFR times width across the remaining finite contiguous age bands for each period. Missing ASFRs in retained bands still produce NA. Report the actual finite range used.

TFR missing-value clarification: at the user's direction, sum ASFR times interval width with `na.rm = TRUE`, omitting missing ASFR contributions. An all-missing set of retained bands therefore sums to zero, following R's sum semantics. Missing mortality-rate handling in life tables is unchanged.

## Continuous-time execution (D024)

Supersedes step-start truncation and random within-step execution above. Public times remain years. Sample continuous mortality and fertility candidates by cumulative-hazard inversion and execute the earliest candidate chronologically, without flooring or rounding. Keep one next event per person; exact equal-time candidates for a person select death. Across persons, exact ties execute in insertion order (initial population row order, then newborn before mother when rescheduling). This deterministic fallback consumes no additional random draws. No execution step or annual batching exists.

Run boundaries remain half-open `[start, start + years)`; any finite nonnegative duration is allowed. After a birth at time `t`, maternal fertility sampling starts at `max(now, t + birth_interval)`, while mortality remains active throughout spacing. Newborns enter risk immediately. Events update father availability immediately. Restart behavior remains resampling rather than exact checkpoint continuation.

Retain the SOCSIM-style fertility compensation on its existing monthly age grid (`1/12` years), with nine preceding intervals for the default spacing of 0.75 years. This grid describes the adjusted hazard, not execution times. Birth intervals still must be whole multiples of that compensation interval; zero spacing uses the unadjusted schedule. The adjustment remains approximate and is not a continuous-time calibration guarantee.

Examples: a birth candidate at 0.02 years precedes a death candidate at 0.06 years, even though both lie in the first month. A birth at 1.03 with spacing 0.75 makes fertility eligible at 1.78, without rounding; mortality may occur before then. A candidate exactly at the endpoint is excluded. Constant adjusted hazard `h` with spacing `b` has long-run renewal rate `1 / (b + 1/h)`, unlike the previous truncated-wait expectation.

Validation: exact sampled mortality dates against R exponential draws, chronological competing events within a former monthly bin, continuous birth spacing and restart spacing, half-open boundaries, constant-hazard survival, renewal expectation, accounting, reproducibility, and father availability at recorded birth times. Full age-specific demographic validation remains outstanding.
