# Decision log

Keep accepted choices separate from proposals. Only record user agreement as acceptance of a demographic or major architectural decision. Each new entry should include status, date, rationale, alternatives, implications, and supporting evidence where available. Update superseded entries with a pointer to their replacement.

## D023: Tibble table outputs

- **Status:** accepted October 12, 2026.
- **Decision:** return all public tables as tibbles, including person/union histories, process schedules, demographic summaries, and empty or auxiliary tables. Continue accepting base data frames as inputs.
- **Implications:** add tibble as a package dependency and normalize table classes while preserving column contents and order. This supersedes earlier promises to preserve arbitrary input data-frame subclasses unchanged.

## Demographic extraction implementation

At the user's request, add person-year exposure, mortality/fertility rate, TFR, and period life-table functions based on the supplied SOCSIM analysis workflow. Support stacked person histories with explicit group columns. Implement exact recorded-lifetime overlaps in C++, with R tables and summaries. Life tables use the already accepted piecewise constant hazard assumption, with consistent exponential person-years and `ex = Tx/lx`; this corrects the supplied code's omission of division by survivorship. See `demographic-analysis.md` and the scientific contract for interfaces, missing-cell behavior, and validation scope.

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

The design specification lists open choices. The four-part candidate architecture is not fully accepted. D019 selects C++ for the first engine; the scientific contract specifies its initial execution interface and rules.

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

## D009: Simulation result class and core tables

- **Status:** accepted October 12, 2026.
- **Decision:** expose simulation results as a named list with class `hammel_results`. Its two core datasets are `persons` and `unions`, accessed as `result$persons` and `result$unions`.
- **Structure:** `persons` has one row per person who entered the simulation, including those who died, with unique person IDs. `unions` has one row per union episode, with references to the participating person IDs; repeated unions are separate episodes.
- **Rationale:** retain distinct person and relationship histories in inspectable, related tables without duplicating person records or nesting histories inside each person.
- **Still open:** exact columns, base data frames versus tibbles, observation and exposure-summary components, metadata, and result methods. No requirement to add a tibble dependency follows from this decision.

## D010: Reuse simulation output as input

- **Status:** accepted October 12, 2026.
- **Decision:** simulation output must be reusable as input for further simulation, preserving the useful SOCSIM workflow. A `hammel_results` supplies the person and union histories and population state needed to initialize a subsequent run.
- **Implications:** preserve unique person and union IDs, deceased persons and completed union histories, and enough endpoint information to identify the starting time and living population. Further runs must allocate new IDs without collisions and maintain consistent event dates and relationship references.
- **Interface direction:** support a pipeable route from a result into further simulation, with demographic processes configurable for the subsequent run. Function names and whether this uses a model constructor or direct result dispatch remain open.
- **Boundary:** reusable demographic output does not by itself promise an exact stochastic checkpoint. Retention of pending event times, RNG state, or compiled engine state, and equality between split and uninterrupted runs remain separate open decisions. Do not require a live engine merely to reuse demographic output.
- **Still open:** endpoint metadata schema, how histories are retained in subsequent results, process configuration carryover, restart RNG semantics, and serialization requirements.

## D011: Basic model configuration API and result rename

- **Status:** accepted October 12, 2026.
- **Decision:** rename the result class and constructor to `hammel_results`. Use `hammel_model()` for an empty specification and pipeable `set_population()`, `set_process()`, `set_mortality()`, and `set_fertility()` functions.
- **Population inputs:** accept person and union tables or a `hammel_results` object. Results supply their tables and endpoint time; existing model process settings are retained.
- **Unmarried initialization:** accepted clarification: with a person-table input, omitted or `NULL` unions represent no initial unions and are normalized to an empty data frame.
- **Process interface:** mortality and fertility setters are shortcuts to the same `set_process()` operation. These constructors and setters do not execute simulations.
- **Current implementation choices:** `set_process(model, name, process)` stores a non-NULL specification unchanged and replaces an existing process of the same name. Results carry optional `end_month`; models carry optional `start_month`, both whole-number indices with an unspecified user-supplied calendar origin. Unknown times remain `NULL`.
- **Scope:** this updates D009's class name and resolves D010's population configuration route. Table schemas, executable process contracts, calendar conventions, and simulation execution remain open. No demographic default or rate conversion is implemented by these setters.

## D012: Initial table schema and year-based dates

- **Status:** accepted October 12, 2026.
- **Decision:** use the initial person and union history schema in `population-schema.md`, leaving marital status off `persons`. Expose times as numeric years with fractional years rather than month indices.
- **Rationale:** year-based dates are convenient for exposure estimates and support future simulation increments other than months. Relationship state can be derived from union histories; caching it internally remains possible.
- **Scope:** table columns use `birth_time`, `death_time`, `start_time`, and `end_time`; these names are implementation choices reflecting the year-based convention. Group and engine linkage pointers are omitted for now.
- **Relationship to earlier decisions:** D004's monthly engine and D008's integer-month truncation remain in effect. D011's current month metadata API requires a future explicit migration to year-based endpoint and starting times.
- **Still open:** calendar origin convention, input time alignment, floating-point tolerances, exposure boundary semantics, and detailed schema validation. The agreed direction is documented; the existing constructors have not yet been migrated or given column validation.

## D013: Canonical year units and explicit step duration

- **Status:** accepted October 12, 2026.
- **Decision:** record public times and durations in years, exposure in person-years, and inputs as rates per person-year. Use monthly increments for the initial milestone but do not hard-code month assumptions throughout the implementation.
- **Implementation direction:** obtain step duration `dt` in years from explicit time configuration; the initial setting is `1/12`. Work directly with rates times durations. Centralize conversions to generic internal step indices if needed.
- **Clarification of D004–D008:** months are the first milestone's configured increment, not the canonical public time unit. D005's division by twelve is the monthly instance of multiplication by interval duration. D008's truncation is `floor(waiting_time_years / dt)` and candidate ties compare step indices; its priority is unchanged.
- **Scope:** supporting other step durations remains future work, but avoiding embedded month assumptions is required from the initial implementation. Existing month-based constructor arguments still need explicit migration; this documentation change does not reinterpret them.

## Implementation update: population validation and year metadata

Both `set_population()` and `hammel_results()` now apply the shared structural checks documented in `population-schema.md`. Public metadata uses `start_time` and `end_time` in years; the earlier month arguments are removed, not reinterpreted. D011–D013 statements describing migration or validation as pending are historical and superseded by this implementation update. Demographic eligibility and time-boundary policies remain open.

## D014: Store the population as a result object

- **Status:** accepted October 12, 2026.
- **Decision:** `hammel_model` stores one `population` component of class `hammel_results`, alongside `processes`. Do not duplicate `persons`, `unions`, or the starting time at the model's top level.
- **Behavior:** table inputs construct and retain a result object; existing result inputs are validated and retained intact. The population's `end_time` is the next run's starting time. `set_population(..., start_time = ...)` remains the table-input interface and sets that reference time.
- **Rationale:** use the same population representation for initial setup and restarts, preserving additional result components and avoiding duplicate state.

## D015: Population class and single population setter input

- **Status:** accepted October 12, 2026.
- **Decision:** rename the class and constructor from `hammel_results` to `hammel_population`. Simplify the setter to `set_population(model, population)`, accepting only that class.
- **Behavior:** create and normalize tables through `hammel_population(persons, unions = NULL, end_time = NULL)`; the setter validates and stores the object unchanged. Initial and simulated populations share this representation. `end_time` remains the population reference time in years and supplies the next run's starting time.
- **Supersedes:** earlier class names and table-input setter signatures in D009–D014. Direct table inputs, `unions`, and `start_time` arguments are removed from `set_population()`.

## D016: Character sex labels

- **Status:** accepted October 12, 2026.
- **Decision:** record sex as character values `"m"` (male) and `"f"` (female), replacing SOCSIM's numeric codes in Hammel's public interface. Population validation requires these nonmissing labels and baseline generation produces them.
- **Process direction:** mortality schedules must use the same labels and provide separate age-specific schedules for males and females from the first milestone. Rate table validation remains to be implemented.
- **Supersedes:** D012's initial numeric sex encoding and the baseline generator's numeric labels. Numeric codes are not silently converted.

## D017: Mortality and fertility process constructors

- **Status:** accepted October 12, 2026.
- **Decision:** add `mortality(rates)` and `fertility(rates)` constructors inheriting from `hammel_process`; use `set_process(model, process)` and make demographic setters constructor-plus-setter shortcuts.
- **Rate schema:** `age_start`, `age_end`, and `rate`, with years and person-year units. Mortality also requires `sex` with separate `m` and `f` schedules. Intervals are half-open and must not overlap within a schedule.
- **Implementation scope:** validate structure and rates; retain supplied schedules. No simulation behavior, interpolation, missing-coverage default, or time-varying rates are implemented. Calendar interval columns are rejected until supported.
- **Supersedes:** D011's explicit-name setter signature and raw process payloads for mortality/fertility.

## D018: HMD and HFD download helpers

- **Status:** accepted October 12, 2026.
- **Decision:** use HMDHFDplus download functions to configure mortality and fertility from a population code and one explicitly selected year, producing constant-over-time schedules for the first milestone.
- **Inputs:** HMD `Mx_1x1` for both sexes; HFD total `asfrRR`. Preserve person-year units and interpret mortality central rates under the accepted hazard model.
- **Fertility endpoints:** at the user's direction, omit HFD open-age rows and retain only closed single-year bands. Compare TFR for that retained age range, not full HFD TFR.
- **Implementation:** optional HMDHFDplus dependency; preserve source metadata and never credentials. Missing retained rates are rejected, not imputed. Authenticated USA/2000 downloads passed the opt-in online test on October 12, 2026; see `testing.md`.

## D019: First executable simulator scope

- **Status:** accepted October 12, 2026.
- **Decision:** implement the initial simulator in C++ with R interfaces. All living females face population age-specific fertility risk; assign a uniformly sampled currently living male parent at each birth, without forming unions. If no males are living, retain the birth with unknown father.
- **Birth spacing:** include SOCSIM-style minimum spacing and its fertility-rate adjustment. The initial default follows SOCSIM's nine months, expressed as 0.75 years. Details and approximation limits are recorded in the scientific contract.
- **Coverage:** fertility outside supplied bands is zero; mortality schedules must cover all ages and end in an open interval.
- **Milestone update:** basic parent assignment replaces partnership formation for this first executable version. Richer matching and union formation are deferred.
- **Implementation choices:** generic step duration, random sequential execution within steps, birth rescheduling, half-open run boundaries, call-local seed behavior, preserved histories, and missing newborn custom attributes are documented in the scientific contract. Full model diagnostics and scientific accuracy assessment are subsequent work.

## D020: Unified demographic summary

- **Status:** accepted October 12, 2026.
- **Decision:** expose one verb-named `summarize_demography()` function returning person-years, counts, ASMR, and ASFR by period, age, and sex. Male births and ASFR are NA. Remove the three separate exposure/mortality/fertility extraction functions and the user-supplied `groups` interface.
- **Use:** feed that same summary into `calculate_tfr()` and `calculate_life_table()`. TFR selects female rows and accepts an optional finite age range matching complete summary bands; life tables separate built-in periods and sexes.
- **Scope:** summarize one independent population history per call; custom strata and simulation-group support are deferred. This supersedes the grouping and separate-extractor portions of the earlier demographic extraction implementation note.

## D021: Automatic open-age handling in demographic summaries

- **Status:** accepted October 12, 2026.
- **Decision:** append Inf to summary age breaks when absent. For TFR, omit open-ended age rows and use all remaining finite age bands; remove the age_range argument.
- **Implications:** report the finite age range used and preserve missing-rate behavior. Users can subset the summary explicitly when desired. This supersedes D020's optional TFR age-range argument.

## D022: Lexis counts function name

- **Status:** accepted October 12, 2026.
- **Decision:** rename `summarize_demography()` to `calculate_lexis_counts()`.
- **Scope:** retain existing exposure, counts, ASMR/ASFR, and interval behavior. Life-table missing-rate handling remains unchanged. This supersedes the function name in D020.

## D024: Continuous-time event execution

- **Status:** accepted October 12, 2026.
- **Decision:** on a branch from master, remove event-time truncation and the execution step; execute continuous candidates chronologically. Retain existing monthly SOCSIM-style fertility compensation separately.
- **Supersedes:** D004/D008/D013/D019 requirements for discrete event execution, truncated candidate comparison, and random within-step order. Year units, cumulative-hazard sampling, and demographic scope remain in effect.
- **Implementation rules:** scientific contract's continuous-time section specifies exact ties, boundaries, spacing, and analytical validation. Monthly compensation and its approximation remain; this does not establish demographic accuracy.

## D025: Reusable database rate downloads

- **Status:** accepted October 12, 2026.
- **Decision:** separate HMD/HFD downloading and preparation from model configuration. Replace database process constructors and setters with `get_mortality_hmd()` and `get_fertility_hfd()`, returning validated rate tibbles accepted by ordinary process constructors and setters.
- **Reuse:** users explicitly save/load tables with RDS; model configuration and repeated simulations need no database requests. Preserve source metadata as a `source` attribute on the table and carry it into constructed processes. Credentials are never retained.
- **Scope:** retain D018's single-year selection, products, units, and age handling. Automatic caching and calendar-varying schedules are outside this change.

## D026: Birth interval belongs to the fertility process

- **Status:** accepted October 13, 2026.
- **Decision:** move `birth_interval` from `simulate()` to `fertility()` and
  `set_fertility()`. Store it in the fertility process and read it at execution.
- **Scope:** retain the 0.75-year default, monthly compensation restrictions,
  continuous execution, and existing spacing and rate-adjustment semantics.
  Run metadata continues to record the configured interval.
- **Supersedes:** earlier execution signatures exposing `birth_interval`.

## D027: Continuous birth-spacing compensation

- **Status:** accepted October 13, 2026.
- **Decision:** use `h(a) = f(a) / (1 - integral_{a-b}^{a} f(u) du)` with
  fertility zero outside supplied bands, replacing monthly SOCSIM compensation.
- **Scope:** fixed age schedule per run; retain current initialization and
  recorded-history spacing, with no transient correction. Remove whole-month
  spacing restrictions; retain the 0.75-year default and reject nonpositive
  eligibility.
- **Implementation:** piecewise linear eligibility and analytic cumulative
  hazard integration/inversion; no numerical age grid or major dependencies.
- **Supersedes:** monthly compensation portions of D019/D024/D026.

## D028: Midpoint approximation using the shared hazard sampler

- **Status:** accepted October 13, 2026.
- **Decision:** calculate one adjusted rate at each original fertility band's
  midpoint, using exact overlaps of the preceding birth interval with the
  original input schedule. Use ordinary Band records and waiting_time() for
  fertility and mortality.
- **Rationale:** prefer a smaller implementation over exact within-band
  adjustment; quantify resulting ASFR/TFR approximation separately.
- **Scope:** retain continuous execution, arbitrary year-valued spacing,
  the 0.75-year default, and existing initialization. Reject nonpositive
  eligibility at evaluated midpoints. At the user's direction, reject all
  open-ended fertility input bands, including zero rates/zero spacing.
- **Supersedes:** D027's age-varying hazard representation and exact sampler.

## Deferred concern: fertility compensation across rate changes

- **Status:** noted and deferred October 13, 2026, at the user's request.
- **Issue:** sequential runs retain actual maternal birth histories, while
  midpoint compensation uses only the new fertility schedule. Recently
  excluded women may reflect earlier rates, producing a transition mismatch
  with the new ASFR target. Its magnitude and duration remain unverified.
- **Current action:** document the limitation; do not implement a correction.
- **Proposal only:** retaining recent rate schedules and effective dates could
  support a history-aware compensation calculation. No algorithm is accepted.
- **Boundary:** this deferral does not adopt built-in calendar-varying rates
  or establish a new scientific guarantee for incremental runs.
