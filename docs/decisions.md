# Decision log

Updated October 19, 2026. Entries retain the choices made at their original
dates; each numbered entry now states whether and how it remains current.
Read [the scientific contract](scientific-contract.md) for implemented behavior,
not a superseded API or timing statement buried in a historical entry.

Only user agreement establishes acceptance of demographic or major architectural
choices. New entries should record status, date, rationale, implications,
alternatives, and evidence. Proposals and deferred concerns are not authorization.

## Current decision map

| Topic | Governing decisions |
| --- | --- |
| Standalone project and scientific ownership | D001–D002 |
| Initial engine scope and milestone completion | D003 as narrowed by D019; D029 |
| Milestone 2 union-formation scope and completion | D030–D043; scope narrowed by D044/D045; completion D050 |
| Union-status summaries and rate-column names | D046–D048; October 17 `overall` label clarification |
| Default age-gap matching | October 18 default age-gap matching entry below |
| Shared kinship and generation-distance exclusions | D036 |
| Person-year inputs and hazard sampling | D005–D007, D013, D024 |
| Population representation, sex labels, public tibbles | D012, D015, D016, D023 |
| Process constructors/setters | D017, D026 |
| Database rate products and reusable downloads | D018, D025 |
| Demographic extraction and function names | D020–D022 and missing-ASFR clarification |
| Continuous execution | D024 |
| Fertility spacing/compensation and finite age bands | D026, D028 (superseding D027's exact sampler) |

## Recorded decisions

## D001: Standalone rewrite

**Current status:** Current project direction.

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** build a standalone R package named Hammel, informed by SOCSIM; do not build a SOCSIM wrapper or rsocsim companion package.
- **Rationale:** modernize the user interface and enable custom demographic mechanisms and methodological innovation beyond the legacy file-oriented workflow.
- **Implications:** investigate SOCSIM algorithms and assumptions; legacy input compatibility is optional.
- **Evidence:** project conversation handoff, summarized in `project-brief.md`.

## D002: Scientific ownership

**Current status:** Current scientific ownership rule.

- **Status:** accepted in the October 11, 2026 conversation.
- **Decision:** the user guides demographic requirements and scientific interpretation; agents support software investigation, implementation, documentation, and validation.
- **Implications:** expose unresolved demographic assumptions and obtain decisions before implementing dependent behavior.

## D003: First demographic validation milestone

**Current status:** Milestone 1 accepted as complete October 15, 2026 (D029); D019
replaced initial partnership formation with random male parent assignment.

- **Status:** accepted October 12, 2026.
- **Decision:** the first major milestone simulates a population with demographic rates constant over calendar time, including fertility, mortality, and a basic partnership mechanism to assign parents.
- **Rationale:** establish demographic accuracy before investing in richer mechanisms.
- **Validation targets:** output fertility, including TFR, and life-table measures should agree with expected input quantities under explicitly defined assumptions and sampling uncertainty.
- **Scope:** a simple partnership model is sufficient. Random marriage after a minimum age was offered as an example; the pairing algorithm and eligibility rules remain open.
- **Implications:** specify rate denominators, partnership effects on fertility eligibility, observation exposure, initialization, and statistical acceptance criteria before implementing dependent behavior.

## D004: Monthly time for the first milestone

**Current status:** Superseded for execution by D024. There is no monthly engine or
configurable execution increment.

- **Status:** accepted October 12, 2026.
- **Decision:** use a fixed month as the simulation time unit for the first milestone, following SOCSIM's monthly convention.
- **Long-term direction:** make time units configurable; supported units and the interface remain unspecified. This flexibility is not required for the first milestone.
- **Implications:** make unit assumptions explicit and centralize conversions, including age, event dates, input rates, and person-time exposure.
- **Still open:** within-month ordering, competing events, boundary conventions, and rate-to-probability or scheduling conversions. Monthly time does not settle the scheduling algorithm.

## D005: Person-year rate inputs

**Current status:** Person-year input units remain current; D024/D028 remove internal
monthly execution/compensation conversions.

- **Status:** accepted October 12, 2026.
- **Decision:** users can supply age-specific fertility and mortality rates as events per person-year. For the monthly engine, convert these internally to events per person-month by dividing by twelve; users need not prepare monthly rate files.
- **Rationale:** support familiar demographic inputs and remove unnecessary input preparation.
- **Implications:** retain explicit rate units and denominators; report demographic estimates in person-year units when comparing them with these inputs. The factor of twelve treats a simulation year as twelve simulation months.
- **Boundary:** division by twelve converts exposure units, not rates into probabilities. The mapping from fertility intensities or mortality central rates/hazards to simulation events still requires specification. Annual death probabilities are a different input type and cannot use this conversion as probabilities.

## D006: Rate-based event probability conversion

**Current status:** The constant-hazard identity and mortality interpretation remain
current; current execution rules are in the scientific contract.
- **Status:** accepted October 12, 2026.
- **Decision:** use rates as the input convention and derive an event probability over a time span `dt` using `1 - exp(-rate * dt)`. For person-year rates and a month, `dt = 1/12`.
- **Assumption:** this expression is exact for a hazard constant over the relevant span. Rates may vary between age bands; handling spans that cross rate boundaries still requires specification.
- **Mortality interpretation:** the intended model treats the supplied mortality rate as a piecewise constant hazard. An empirical life-table central rate `nmx` is not universally identical to that hazard; using it directly is a modeling assumption that must be documented and reflected in validation references.
- **Boundary:** this defines a single-process probability, not the treatment of competing events, fertility eligibility, repeat births, or within-month ordering.
- **Implementation guidance:** use the numerically stable equivalent `-expm1(-rate * dt)` when implementing this conversion.

## D007: Sample waiting times directly on the hazard scale

**Current status:** Cumulative-hazard sampling remains current. D024 resolves continuous
dates and ties; the implementation draws Exp(1) thresholds.

- **Status:** accepted October 12, 2026.
- **Decision:** use the mathematical equivalence between log survival and cumulative hazard to sample candidate waiting times directly from rates, without converting to probabilities and back to logarithms.
- **Method:** draw `E = -log(U)` for `U` uniform on `(0, 1)`. Accumulate `rate * duration` across successive age bands until the cumulative hazard reaches `E`. Carry the remaining threshold into each subsequent band; within the final positive-rate band, the remaining waiting time is `remaining_threshold / rate`.
- **Units:** rates per person-year use durations in years; one simulation month contributes `rate / 12` to cumulative hazard.
- **Relationship to D006:** this is the intended candidate waiting-time method. D006's probability formula remains the equivalent single-process interval probability and may be used for interpretation and validation.
- **Resolved subsequently:** D008 adopts SOCSIM-style integer truncation and candidate tie priority. Rescheduling after state changes, terminal age handling, and fertility eligibility and repeat-event semantics remain open.

## D008: Integer-month candidate dates and competing-event ties

**Current status:** Integer truncation and equal-month competition are superseded by
D024's continuous chronology and exact ties.

- **Status:** accepted October 12, 2026.
- **Decision:** adopt SOCSIM-style integer truncation of nonnegative candidate waiting times to months and its priority for competing candidate events in the same month.
- **Monthly assignment:** convert the sampled waiting time to months and truncate its fractional part (equivalently, take the floor). Zero means the current simulation month. Exact age-band boundary conventions still need worked specifications.
- **Tie priority:** for a person's equal-month candidates, select death before birth before marriage/partnership formation before divorce/partnership dissolution. Group transitions follow these; specify their ordering when introduced. Compare integer-month dates, not hidden fractional waiting times.
- **Scope:** retain one selected next event per person. This priority does not establish a global order among different people's events.
- **Implications:** an earlier fractional birth time does not beat a later fractional death time if both truncate to the same month. Validate and document this deliberate discretization effect.
- **Still open:** execution order across people, resampling and same-month follow-up events, checkpoint behavior, and terminal age handling. This decision does not automatically adopt all SOCSIM scheduling behavior.
- **Reference:** the SOCSIM rate and competing-event audits at rsocsim commit `125450cff968bb3ff762fa99d4c926a8136062a2`.

## D009: Simulation result class and core tables

**Current status:** Core person/union tables remain current; the class is now
hammel_population (D015), with tibble outputs (D023).

- **Status:** accepted October 12, 2026.
- **Decision:** expose simulation results as a named list with class `hammel_results`. Its two core datasets are `persons` and `unions`, accessed as `result$persons` and `result$unions`.
- **Structure:** `persons` has one row per person who entered the simulation, including those who died, with unique person IDs. `unions` has one row per union episode, with references to the participating person IDs; repeated unions are separate episodes.
- **Rationale:** retain distinct person and relationship histories in inspectable, related tables without duplicating person records or nesting histories inside each person.
- **Still open:** exact columns, base data frames versus tibbles, observation and exposure-summary components, metadata, and result methods. No requirement to add a tibble dependency follows from this decision.

## D010: Reuse simulation output as input

**Current status:** Output reuse is implemented. D019/D024 specify history retention,
resampling, and seed behavior; exact checkpointing remains unsupported.

- **Status:** accepted October 12, 2026.
- **Decision:** simulation output must be reusable as input for further simulation, preserving the useful SOCSIM workflow. A `hammel_results` supplies the person and union histories and population state needed to initialize a subsequent run.
- **Implications:** preserve unique person and union IDs, deceased persons and completed union histories, and enough endpoint information to identify the starting time and living population. Further runs must allocate new IDs without collisions and maintain consistent event dates and relationship references.
- **Interface direction:** support a pipeable route from a result into further simulation, with demographic processes configurable for the subsequent run. Function names and whether this uses a model constructor or direct result dispatch remain open.
- **Boundary:** reusable demographic output does not by itself promise an exact stochastic checkpoint. Retention of pending event times, RNG state, or compiled engine state, and equality between split and uninterrupted runs remain separate open decisions. Do not require a live engine merely to reuse demographic output.
- **Still open:** endpoint metadata schema, how histories are retained in subsequent results, process configuration carryover, restart RNG semantics, and serialization requirements.

## D011: Basic model configuration API and result rename

**Current status:** Historical API stage. D015 replaces population inputs/class; D017
replaces process setter signatures; public metadata is in years.

- **Status:** accepted October 12, 2026.
- **Decision:** rename the result class and constructor to `hammel_results`. Use `hammel_model()` for an empty specification and pipeable `set_population()`, `set_process()`, `set_mortality()`, and `set_fertility()` functions.
- **Population inputs:** accept person and union tables or a `hammel_results` object. Results supply their tables and endpoint time; existing model process settings are retained.
- **Unmarried initialization:** accepted clarification: with a person-table input, omitted or `NULL` unions represent no initial unions and are normalized to an empty data frame.
- **Process interface:** mortality and fertility setters are shortcuts to the same `set_process()` operation. These constructors and setters do not execute simulations.
- **Current implementation choices:** `set_process(model, name, process)` stores a non-NULL specification unchanged and replaces an existing process of the same name. Results carry optional `end_month`; models carry optional `start_month`, both whole-number indices with an unspecified user-supplied calendar origin. Unknown times remain `NULL`.
- **Scope:** this updates D009's class name and resolves D010's population configuration route. Table schemas, executable process contracts, calendar conventions, and simulation execution remain open. No demographic default or rate conversion is implemented by these setters.

## D012: Initial table schema and year-based dates

**Current status:** Table structure and year-valued dates remain current.
Migration/validation are implemented; monthly timing references are superseded by D024.

- **Status:** accepted October 12, 2026.
- **Decision:** use the initial person and union history schema in `population-schema.md`, leaving marital status off `persons`. Expose times as numeric years with fractional years rather than month indices.
- **Rationale:** year-based dates are convenient for exposure estimates and support future simulation increments other than months. Relationship state can be derived from union histories; caching it internally remains possible.
- **Scope:** table columns use `birth_time`, `death_time`, `start_time`, and `end_time`; these names are implementation choices reflecting the year-based convention. Group and engine linkage pointers are omitted for now.
- **Relationship to earlier decisions:** D004's monthly engine and D008's integer-month truncation remain in effect. D011's current month metadata API requires a future explicit migration to year-based endpoint and starting times.
- **Still open:** calendar origin convention, input time alignment, floating-point tolerances, exposure boundary semantics, and detailed schema validation. The agreed direction is documented; the existing constructors have not yet been migrated or given column validation.

## D013: Canonical year units and explicit step duration

**Current status:** Canonical year units remain current. Explicit monthly steps and
integer execution indices are superseded by D024.

- **Status:** accepted October 12, 2026.
- **Decision:** record public times and durations in years, exposure in person-years, and inputs as rates per person-year. Use monthly increments for the initial milestone but do not hard-code month assumptions throughout the implementation.
- **Implementation direction:** obtain step duration `dt` in years from explicit time configuration; the initial setting is `1/12`. Work directly with rates times durations. Centralize conversions to generic internal step indices if needed.
- **Clarification of D004–D008:** months are the first milestone's configured increment, not the canonical public time unit. D005's division by twelve is the monthly instance of multiplication by interval duration. D008's truncation is `floor(waiting_time_years / dt)` and candidate ties compare step indices; its priority is unchanged.
- **Scope:** supporting other step durations remains future work, but avoiding embedded month assumptions is required from the initial implementation. Existing month-based constructor arguments still need explicit migration; this documentation change does not reinterpret them.

## D014: Store the population as a result object

**Current status:** Single stored population remains current; its class/setter are
replaced by D015.

- **Status:** accepted October 12, 2026.
- **Decision:** `hammel_model` stores one `population` component of class `hammel_results`, alongside `processes`. Do not duplicate `persons`, `unions`, or the starting time at the model's top level.
- **Behavior:** table inputs construct and retain a result object; existing result inputs are validated and retained intact. The population's `end_time` is the next run's starting time. `set_population(..., start_time = ...)` remains the table-input interface and sets that reference time.
- **Rationale:** use the same population representation for initial setup and restarts, preserving additional result components and avoiding duplicate state.

## D015: Population class and single population setter input

**Current status:** Current population class and setter. D023 normalizes tables to
tibbles while preserving contents and extra components.

- **Status:** accepted October 12, 2026.
- **Decision:** rename the class and constructor from `hammel_results` to `hammel_population`. Simplify the setter to `set_population(model, population)`, accepting only that class.
- **Behavior:** create and normalize tables through `hammel_population(persons, unions = NULL, end_time = NULL)`; the setter validates and stores the object unchanged. Initial and simulated populations share this representation. `end_time` remains the population reference time in years and supplies the next run's starting time.
- **Supersedes:** earlier class names and table-input setter signatures in D009–D014. Direct table inputs, `unions`, and `start_time` arguments are removed from `set_population()`.

## D016: Character sex labels

**Current status:** Current sex labels.

- **Status:** accepted October 12, 2026.
- **Decision:** record sex as character values `"m"` (male) and `"f"` (female), replacing SOCSIM's numeric codes in Hammel's public interface. Population validation requires these nonmissing labels and baseline generation produces them.
- **Process direction:** mortality schedules must use the same labels and provide separate age-specific schedules for males and females from the first milestone. Rate table validation remains to be implemented.
- **Supersedes:** D012's initial numeric sex encoding and the baseline generator's numeric labels. Numeric codes are not silently converted.

## D017: Mortality and fertility process constructors

**Current status:** Current constructors and setter names. D028 subsequently prohibits
open-ended fertility inputs; execution resolves coverage and eligibility.

- **Status:** accepted October 12, 2026.
- **Decision:** add `mortality(rates)` and `fertility(rates)` constructors inheriting from `hammel_process`; use `set_process(model, process)` and make demographic setters constructor-plus-setter shortcuts.
- **Rate schema:** `age_start`, `age_end`, and `rate`, with years and person-year units. Mortality also requires `sex` with separate `m` and `f` schedules. Intervals are half-open and must not overlap within a schedule.
- **Implementation scope:** validate structure and rates; retain supplied schedules. No simulation behavior, interpolation, missing-coverage default, or time-varying rates are implemented. Calendar interval columns are rejected until supported.
- **Supersedes:** D011's explicit-name setter signature and raw process payloads for mortality/fertility.

## D018: HMD and HFD download helpers

**Current status:** Data products and endpoint omission remain current. D025 replaces
database-specific process/setter interfaces with reusable rate downloads.

- **Status:** accepted October 12, 2026.
- **Decision:** use HMDHFDplus download functions to configure mortality and fertility from a population code and one explicitly selected year, producing constant-over-time schedules for the first milestone.
- **Inputs:** HMD `Mx_1x1` for both sexes; HFD total `asfrRR`. Preserve person-year units and interpret mortality central rates under the accepted hazard model.
- **Fertility endpoints:** at the user's direction, omit HFD open-age rows and retain only closed single-year bands. Compare TFR for that retained age range, not full HFD TFR.
- **Implementation:** optional HMDHFDplus dependency; preserve source metadata and never credentials. Missing retained rates are rejected, not imputed. Authenticated USA/2000 downloads passed the opt-in online test on October 12, 2026; see `testing.md`.

## D019: First executable simulator scope

**Current status:** C++/R implementation, fertility scope, and parent assignment remain
current. D024 replaces discrete scheduling; D026 moves spacing to the process; D028
replaces compensation.

- **Status:** accepted October 12, 2026.
- **Decision:** implement the initial simulator in C++ with R interfaces. All living females face population age-specific fertility risk; assign a uniformly sampled currently living male parent at each birth, without forming unions. If no males are living, retain the birth with unknown father.
- **Birth spacing:** include SOCSIM-style minimum spacing and its fertility-rate adjustment. The initial default follows SOCSIM's nine months, expressed as 0.75 years. Details and approximation limits are recorded in the scientific contract.
- **Coverage:** fertility outside supplied bands is zero; mortality schedules must cover all ages and end in an open interval.
- **Milestone update:** basic parent assignment replaces partnership formation for this first executable version. Richer matching and union formation are deferred.
- **Implementation choices:** generic step duration, random sequential execution within steps, birth rescheduling, half-open run boundaries, call-local seed behavior, preserved histories, and missing newborn custom attributes are documented in the scientific contract. Full model diagnostics and scientific accuracy assessment are subsequent work.

## D020: Unified demographic summary

**Current status:** Shared analysis and built-in period/sex grouping remain current.
D021 changes age-range handling; D022 renames the function calculate_lexis_counts().
Missing-ASFR handling was clarified subsequently.

- **Status:** accepted October 12, 2026.
- **Decision:** expose one verb-named `summarize_demography()` function returning person-years, counts, ASMR, and ASFR by period, age, and sex. Male births and ASFR are NA. Remove the three separate exposure/mortality/fertility extraction functions and the user-supplied `groups` interface.
- **Use:** feed that same summary into `calculate_tfr()` and `calculate_life_table()`. TFR selects female rows and accepts an optional finite age range matching complete summary bands; life tables separate built-in periods and sexes.
- **Scope:** summarize one independent population history per call; custom strata and simulation-group support are deferred. This supersedes the grouping and separate-extractor portions of the earlier demographic extraction implementation note.

## D021: Automatic open-age handling in demographic summaries

**Current status:** Open-age summary/TFR behavior remains current. D059
supersedes the later clarification that omitted missing ASFR contributions.

- **Status:** accepted October 12, 2026.
- **Decision:** append Inf to summary age breaks when absent. For TFR, omit open-ended age rows and use all remaining finite age bands; remove the age_range argument.
- **Implications:** report the finite age range used and preserve missing-rate behavior. Users can subset the summary explicitly when desired. This supersedes D020's optional TFR age-range argument.

## D022: Lexis counts function name

**Current status:** Current summary function name.

- **Status:** accepted October 12, 2026.
- **Decision:** rename `summarize_demography()` to `calculate_lexis_counts()`.
- **Scope:** retain existing exposure, counts, ASMR/ASFR, and interval behavior. Life-table missing-rate handling remains unchanged. This supersedes the function name in D020.

## D023: Tibble table outputs

**Current status:** Current public table convention and dependency.

- **Status:** accepted October 12, 2026.
- **Decision:** return all public tables as tibbles, including person/union histories, process schedules, demographic summaries, and empty or auxiliary tables. Continue accepting base data frames as inputs.
- **Implications:** add tibble as a package dependency and normalize table classes while preserving column contents and order. This supersedes earlier promises to preserve arbitrary input data-frame subclasses unchanged.

## D024: Continuous-time event execution

**Current status:** Continuous event execution remains current. Its retained monthly
compensation is superseded by D027/D028.

- **Status:** accepted October 12, 2026.
- **Decision:** on a branch from master, remove event-time truncation and the execution step; execute continuous candidates chronologically. Retain existing monthly SOCSIM-style fertility compensation separately.
- **Supersedes:** D004/D008/D013/D019 requirements for discrete event execution, truncated candidate comparison, and random within-step order. Year units, cumulative-hazard sampling, and demographic scope remain in effect.
- **Implementation rules:** scientific contract's continuous-time section specifies exact ties, boundaries, spacing, and analytical validation. Monthly compensation and its approximation remain; this does not establish demographic accuracy.

## D025: Reusable database rate downloads

**Current status:** Current reusable database download interface.

- **Status:** accepted October 12, 2026.
- **Decision:** separate HMD/HFD downloading and preparation from model configuration. Replace database process constructors and setters with `get_mortality_hmd()` and `get_fertility_hfd()`, returning validated rate tibbles accepted by ordinary process constructors and setters.
- **Reuse:** users explicitly save/load tables with RDS; model configuration and repeated simulations need no database requests. Preserve source metadata as a `source` attribute on the table and carry it into constructed processes. Credentials are never retained.
- **Scope:** retain D018's single-year selection, products, units, and age handling. Automatic caching and calendar-varying schedules are outside this change.

## D026: Birth interval belongs to the fertility process

**Current status:** Spacing belongs to fertility. D027/D028 remove its earlier whole-
month restriction.

- **Status:** accepted October 13, 2026.
- **Decision:** move `birth_interval` from `simulate()` to `fertility()` and
  `set_fertility()`. Store it in the fertility process and read it at execution.
- **Scope:** retain the 0.75-year default, monthly compensation restrictions,
  continuous execution, and existing spacing and rate-adjustment semantics.
  Run metadata continues to record the configured interval.
- **Supersedes:** earlier execution signatures exposing `birth_interval`.

## D027: Continuous birth-spacing compensation

**Current status:** Exact age-varying compensation/specialized sampler superseded by
D028. Year-valued fractional spacing remains current.

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

**Current status:** Current midpoint compensation, shared sampler, and finite fertility
inputs.

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

## D029: Milestone 1 completion

**Current status:** Milestone 1 complete; milestone 2 scope subsequently selected in
D030.

- **Status:** accepted October 15, 2026 by the user, after rendering the report.
- **Decision:** close the initial mortality/fertility demographic-validation
  milestone on the basis of Aaron Gullickson's completed assessment in
  `validation/milestone1.qmd` and its rendered `validation/milestone1.html`.
- **Evidence:** replicated demographic comparisons, age-specific coverage and
  deviations, two-sex Leslie projection and stable age structure comparisons,
  extended runs, fertility-rate restart cases, and performance/test results
  presented in the report. The agent did not independently rerun the report
  when recording acceptance.
- **Implications:** completion reflects the user's scientific assessment of
  the reported cases; it does not establish universal accuracy or add numerical
  thresholds. Existing engine assumptions remain unchanged.
- **Remaining limits:** initialization transients, sparse oldest ages, and
  fertility compensation after rate changes. The latter remains deferred;
  the report proposes a GitHub issue, but no issue is created by this decision.
  No controlled executable SOCSIM comparison is claimed.
- **Scope:** no new feature, union-formation implementation, or next milestone
  is authorized by this completion decision.

## D030: Milestone 2 union-formation scope

**Current status:** Implemented under D031–D043 with waiting and transition-aware
compensation deferred by D044/D045. The report now includes demographic analysis,
parameter recovery and performance comparisons; D050 records user acceptance
of milestone completion on October 19, 2026.
D032 supersedes the original requirement to keep no-union execution available.

- **Status:** accepted October 15, 2026 through the user's milestone 2 request.
- **Decision:** separate female search initiation from partner matching. Support
  birth-triggered on-demand search and independently scheduled union-search
  events. Use conditional-logit partner choice, initially based on age difference,
  with an extensible matching specification and eligibility/kinship screening.
  Scheduled mode must support partnered and unpartnered fertility. End unions by
  widowhood; other dissolution is outside this milestone.
- **Terminology:** use union/union formation and partnered/unpartnered in code
  and documentation, agnostic about legal marriage versus cohabitation. The
  user's conversational use of “marriage” does not introduce a legal-status model.
- **Matching rationale:** probabilistic choice allows rare eligible matches
  through model coefficients and later extensions such as group preferences.
  The user also wants to validate matching by fitting the same choice model to
  simulated decisions. Candidate sets will not be retained. Approximate model recovery will use
  eligible men reconstructed in time slices; slice definitions and recovery
  criteria remain unresolved.
- **Deferred:** parity-specific fertility remains for a later milestone;
  milestone 2 does not add it.
- **Open choices:** full versus sampled candidate pools, optional waiting,
  failed-search behavior, exact eligibility/kin prohibitions, age-gap model,
  search-rate interpretation/retries, conditional-fertility compensation, and
  birth-time classification/parent assignment. No defaults are accepted for these.
- **Implications:** write specifications and worked examples and obtain dependent
  demographic decisions before implementation. Keep milestone 1 behavior
  available. SOCSIM matching and kinship behavior require a source audit before
  claiming equivalence; no legacy algorithm is adopted by this decision.
- **Evidence and next specification:** [milestone 2 notes](union-formation-models.md).

## D031: Temporary on-demand union execution

**Current status:** Temporary uniform/unrestricted matching is superseded by
D034–D037; D043 replaces default unmatched-birth preservation. Both initiation
modes and conditional fertility are now implemented (D038–D041).

- **Status:** accepted October 15, 2026 through the user's implementation request.
- **Decision:** at birth use a woman's ongoing partner, otherwise sample living
  men without replacement at the configured size and choose one uniformly to
  form a union. Skip pair screening and conditional-logit choice for now.
- **Empty pool:** the user explicitly selected preserving the birth with an
  unknown father and no union. A future toggle for missing father versus birth
  cancellation is requested but deferred.
- **Implications:** the temporary pool includes partnered men, relatives, and
  children. Maintain current partners, union histories, widowhood, restarts,
  and unique IDs. No union-status fertility dependence is introduced, so
  existing birth/death candidates of surviving partners remain valid.
- **Boundary:** scheduled initiation, pair eligibility, utility-based matching,
  and conditional fertility remain pending; the final eligibility policy is
  not established by this intentionally unrestricted increment.
- **Specification:** [on-demand increment](union-formation-models.md#on-demand-conditional-logit-execution-specification).

## D032: Union formation required for execution

**Current status:** Accepted and implemented; supersedes optional execution in
D031. D038 subsequently adds rate-based initiation.

- **Status:** accepted October 15, 2026 at the user's direction.
- **Decision:** every simulation requires an explicit union-formation process.
  Initialize partner caches and union histories unconditionally; do not add an
  enable flag. Initiation mode determines search timing, not state tracking.
- **Implications:** remove the no-union father-assignment path. Missing union
  processes fail even for zero-duration runs. On-demand is currently the only
  executable initiation mode; scheduled initiation remains rejected.
- **Compatibility:** existing workflows need `set_union_formation()`. Preserve
  the accepted milestone 1 report as historical evidence at its recorded commit;
  do not rewrite its analysis for the changed milestone 2 model.

## D033: Sex-specific union partner columns

**Current status:** Accepted and implemented; supersedes unordered union roles in D012.

- **Status:** accepted October 15, 2026 at the user's request.
- **Decision:** replace `partner1_id` and `partner2_id` with `wife_id` and
  `husband_id`. Wife references a female person and husband a male person,
  for both ongoing and ended histories. Continue using “union” for the
  relationship; these role names do not establish legal marriage status.
- **Implications:** population validation enforces sex roles. C++ histories,
  partner lookup, death handling, and output use explicit wife/husband roles.
  Empty union tables have the same renamed columns. Existing unordered
  tables must be mapped by sex, not just renamed by position.
- **Validation:** valid noncontiguous IDs, swapped-role rejection for ongoing
  and ended unions, empty schema, and existing birth/death/restart cases.

## D034: Execute conditional-logit partner selection

**Current status:** Current conditional-logit implementation; eligibility is now
implemented under D035–D037. The milestone 2 report now includes parameter recovery.

- **Status:** accepted October 15, 2026 through the user's implementation request.
- **Decision:** use the entered matching formula and coefficients to generate
  candidate utilities on continuous male-minus-female age differences. Convert
  utilities to softmax probabilities and sample one partner. Keep existing
  sampling-before-screening order; eligibility screening remains deferred.
- **Worked test model:** linear coefficient 0.072 and raw squared coefficient
  -0.014 have a utility maximum at 2.5714286 years. These are requested test
  inputs, not defaults or a target mean gap for every candidate pool.
- **Implementation:** exact C++ evaluation for linear/raw quadratic terms,
  existing R model-matrix semantics for other formulas, stable softmax and one
  uniform selection draw. No candidate histories or new dependencies. Reject
  nonfinite formula values rather than silently dropping candidates.
- **Validation:** analytical utilities/probabilities, choice frequencies,
  formula fallback, numerical extremes, birth-time parentage, and existing
  reproducibility/widowhood/restart checks. Scientific model recovery from
  reconstructed time slices is still future validation work.

## D035: Minimum male age for partner searches

**Current status:** Implemented; D043 supersedes the original failed-screen birth-
preservation default.

- **Status:** accepted October 15, 2026 through the user's eligibility request
  and explicit choices of default age and failed-screen behavior.
- **Decision:** add user-selectable `min_male_age`, default 15 continuous years.
  Finite nonnegative scalar; fractional ages supported, zero disables exclusion.
  Include men at or above the threshold at search time. No additional female
  age restriction; fertility/search rates govern female initiation.
- **Ordering:** sample first, exclude underage sampled men, then score and
  normalize matching probabilities. No refilling or search outside the sample.
- **Failed screen:** preserve the birth with unknown father and no new union
  when no sampled men pass, as explicitly selected by the user.
- **Scope:** applies to new searches; retain existing ongoing unions. Other
  age/kinship/availability policies remain outside this increment.
- **Validation:** constructor/edited-process checks, exact/fractional age
  boundaries, aging into eligibility, no added female restriction, retained
  existing partners, and sampling-before-screening failures.

## D036: Shared kinship module and generation-distance exclusions

**Current status:** Implemented; D037 adds exclusivity and D043 supersedes the original
failed-screen birth-preservation default.

- **Status:** accepted October 15, 2026 through the user's implementation request
  and explicit incomplete-ancestry decision.
- **Decision:** use one separate native `kinship` module for union screening and
  R extraction. Exclude candidates sharing a recorded ancestor within
  `kinship_generations` parent steps of EACH person, including each person at
  distance zero. Default two; zero imposes no cross-person exclusion.
- **Incomplete histories:** exclude demonstrated relationships only. Unknown
  links stop that branch and are not evidence of unrelatedness; do not skip
  known paternal relationships because a mother is missing.
- **Consequences:** two includes known half/full siblings and first cousins,
  aunts/uncles, parents and grandparents. Removed/second cousins and direct
  ancestors beyond two generations require a larger bound. No unlimited
  direct-ancestor exclusion or new male-exclusivity rule is inferred.
- **Engineering:** iterative cycle-safe traversal, lazily cached distinct
  ancestors and minimum distances, shortest-route maternal/paternal evidence,
  newborn append and restart reconstruction. No dense all-pairs kin matrix or
  enumeration of all ancestry paths. No major new dependency or copied SOCSIM
  implementation.
- **R interface:** `extract_kin()` returns compact common-ancestor evidence for
  selected egos/relatives as a tibble, using the same distance kernel. Output
  retains closest tied witnesses rather than assigning ambiguous relationship
  labels; lineage and within-depth completeness are explicit.
- **Scope:** screening follows sampling and male-age checks, before utilities,
  without refilling. Failed searches preserve births with unknown father.
  Existing unions are retained. Social kin and genetic coefficients are not
  implemented. See [kinship](kinship.md) for worked examples and limitations.

## Additional clarifications and deferred concerns

### Demographic extraction implementation

The original request adapted the user's SOCSIM analysis workflow to exact
recorded-lifetime exposure and R-native summaries. Life tables follow the
accepted hazard interpretation and use ex=Tx/lx. The earlier stacked-history
interface was removed by D020; current functions analyze one independent
history per call. See [demographic analysis](demographic-analysis.md).

### Deferred concern: fertility compensation across rate changes

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

### Missing-ASFR clarification

**Superseded by D059:** missing ASFRs now propagate into TFR.

At the user's direction, `calculate_tfr()` sums with `na.rm = TRUE`, omitting
missing contributions; an all-missing retained age range therefore sums to zero.
This supersedes the earlier missing-TFR statements in D020/D021 and the
earlier scientific contract (available in Git history).
It does not change life-table missing-rate propagation. This restates an
existing accepted clarification; it is not a new decision from the doc audit.

### Milestone 2 construction clarification (October 15, 2026)

At the user's direction, candidate sampling precedes pair-specific screening.
Age difference is male age minus female age in continuous years. Process
construction now accepts explicit initiation, sample size (default 1000, zero for full pool),
and a matching formula with named coefficients. R model-matrix expansion and
stable softmax provide the construction/scoring infrastructure without a new
major dependency. At that construction-only stage, execution rules were unresolved and the engine
rejected union processes. D031–D043 subsequently implemented both initiation
modes and their eligibility, matching, failure and fertility-transition rules;
D044/D045 defer waiting and transition-aware compensation.

The user subsequently selected `sample_size = 1000` as the default for both
`union_formation()` and `set_union_formation()` on October 15, 2026.

### Sex-specific union death reasons (October 15, 2026)

At the user's direction, union `end_reason` now uses `death_husband` and
`death_wife` in place of `death`. Existing `dissolution` and ongoing `NA`
remain valid. The executed death event records the terminating partner;
subsequent deaths do not overwrite an ended union, including equal-time
events under existing queue ordering. Supplied ended histories and restart
reasons are retained. This supersedes generic death labels in earlier entries.


## D037: Exclusive ongoing unions and pre-sampling male availability

**Current status:** Implemented; D043 supersedes the referenced birth-preservation
default.

- **Status:** accepted by the user October 15, 2026.
- **Decision:** require at most one ongoing union per person. Sample from living,
  unpartnered men before applying age and kin restrictions. Remove a man when
  selected; return a surviving husband when widowed. Never return dead men.
  Reject conflicting ongoing input unions for either sex; retain ended history.
- **Consequences:** scarcity can leave births with unknown fathers under the
  existing failed-search policy. Fertility/mortality schedules remain unchanged.
  This supersedes the temporary concurrent-male-union behavior in D031/D036.
  Polygamy/polyamory remains a possible later explicit model, not current behavior.
- **Validation:** pre-sampling exclusion, scarcity, widow re-entry, death removal,
  conflicting input unions, restart consistency, and no concurrent unions.


## D038: Female rate-based search events and partnered-only births

**Current status:** Search scheduling remains current; D039/D040 replace partnered-only
births with status-specific fertility and nonmarital father selection.

- **Status:** accepted by the user October 15, 2026.
- **Decision:** use age-specific search hazards as competing events for
  unpartnered women only, with no differentiation by previous/current legal
  marital status. Reuse on-demand matching. Failed attempts draw another waiting
  time from the same schedule, starting at failure. Rate mode permits births
  only while partnered; union formation enables risk subject to birth spacing,
  and widowhood cancels pending births and restores search risk immediately.
- **Boundary:** retain the existing adjusted fertility schedule, without
  claiming recovery of population ASFRs under this gate. Separate conditional
  fertility schedules and the wait alternative remain unimplemented.
- **Validation:** analytical competing search/death hazards, age support/gaps,
  retries, shared partner choice, no unpartnered births, widowhood invalidation,
  restart histories and reproducibility.

## D039: Stacked union-status fertility schedules and overall fallback

**Current status:** Implemented; D040/D041 replace unconditional unknown-father
assignment and D043 governs failed matches.

- **Status:** accepted by the user October 15, 2026.
- **Decision:** keep the fertility API, adding an optional `union_status` column
  with both `partnered` and `unpartnered`. Without it, use one overall schedule.
  Reject partial/mixed status tables; explicit zero schedules express suppression.
  On-demand rejects conditional tables. Rate mode accepts overall rates for both
  statuses with one warning per call, superseding D038's partnered-only gate.
- **Transitions:** reschedule on formation/widowhood, preserving maternal spacing;
  births to unpartnered women in rate mode have unknown father and do not search.
  Existing midpoint adjustment operates separately on each schedule; transition
  compensation remains deferred. Death wins ties, then birth, then search.

## D040: Nonmarital father selection shares union matching

**Current status:** Implemented with D041 separate father matching and D043 default
failed-match cancellation.

- **Status:** accepted by the user October 15, 2026.
- **Decision:** rate-mode unpartnered births select fathers from currently living,
  unpartnered men using the same sampling, eligibility and matching model as
  partner selection. Assigning a father does not form a union or reserve the man.
  Empty screened pools preserve births with unknown father. This supersedes
  D039's unconditional unknown-father rule for unpartnered births.
- **Deferred:** separate father-selection and partner-selection models are an
  explicitly deferred extension. Preferential reuse of the previous child's
  father at later births or union formation is under discussion, not accepted;
  current selections have no special persistence preference.
- **Validation:** father assignment without unions, repeated availability,
  partnered/dead/young/kin exclusions, matching probabilities, known paternal
  links in kinship, retained births when no candidate qualifies, and restarts.

## D041: Optional separate nonmarital father matching model

**Current status:** Implemented; previous-father continuity is deferred under D042.

- **Status:** accepted by the user October 15, 2026; supersedes D040's deferral.
- **Decision:** add `father_matching = NULL` to union constructors/setters. NULL
  uses the current general `matching` model; an explicit model applies only to
  unpartnered rate-mode births. All union searches retain the general model.
  Sampling, eligibility and parentage/state semantics remain shared.
- **Validation:** constructor/edited-object checks, exact fallback equivalence,
  opposite preference models separating births from unions, and native versus
  R formula evaluation. Preferential reuse of previous fathers remains unresolved.


## D042: Revised milestone 2 completion scope

**Current status:** Completion scope subsequently narrowed by D044/D045; cancellation is
implemented (D043), and the report includes performance comparisons. Milestone 2
was accepted as complete on October 19, 2026 (D050).

- **Status:** accepted by the user October 15, 2026.
- **Required work:** implement a failed-match birth cancellation toggle and
  address fertility compensation across status/rate transitions in milestone 2.
  This supersedes their earlier deferral; exact cancellation semantics and the
  compensation algorithm still require specification and user decisions.
- **Under consideration:** the poor-pool wait alternative may remain in this
  milestone or be postponed; do not treat it as a settled completion requirement.
- **Later milestone:** preferential reuse of previous nonmarital fathers at
  subsequent births or union formation is explicitly deferred. Current draws
  retain no special previous-father preference.
- **Validation:** benchmark R formula-fallback performance in the final milestone
  validation report rather than as a separate implementation task. No benchmark
  result is claimed. Parity conditioning, non-death dissolution and multiple
  ongoing unions remain outside current scope.

## D043: Failed-match birth cancellation defaults on

- **Status:** toggle and default accepted by the user October 15, 2026.
- **Decision:** add `cancel_birth_if_no_match = TRUE` to the union process. Apply
  to failed on-demand partner matching and failed rate-mode nonmarital father
  matching. FALSE explicitly permits new births with unknown father; input
  histories are preserved. This supersedes default birth preservation in D031/D040.
- **Retry:** accepted by the user: redraw from the attempt time without starting
  a new birth-spacing interval or changing the last actual birth. No correction
  should erase the intended fertility effect of failed matches.


## D044: Defer poor-pool waiting to a later milestone

**Current status:** Waiting remains deferred; D045 also defers compensation. The report
includes performance comparisons; milestone 2 is complete (D050).

- **Status:** accepted by the user October 15, 2026; supersedes D042's open scope.
- **Decision:** postpone the option to decline a nonempty eligible pool and wait
  to a later milestone. It is not a milestone 2 completion requirement.
- **Current behavior:** a nonempty screened pool always yields a selection via
  conditional-logit matching. Empty searches retain existing retry behavior;
  unmatched births follow `cancel_birth_if_no_match` (default TRUE).
- **Subsequent disposition:** D045 deferred compensation; the report includes
  demographic validation and R formula-fallback timing. D050 records completion.


## D045: Defer transition-aware fertility compensation

**Current status:** Compensation remains deferred; the report now assesses fertility
discrepancies and performance. Milestone 2 was accepted as complete (D050);
this deferral remains in force.

- **Status:** accepted by the user October 15, 2026; supersedes D042/D044's
  remaining implementation requirement.
- **Decision:** postpone correction of fertility compensation across union-status
  transitions and changes of rates between simulation runs to later work. Retain
  per-schedule midpoint compensation and spacing based on actual recorded births.
  No individual rate-history or dynamic exposure correction is implemented.
- **Evidence boundary:** the report assesses discrepancies for its input case;
  it does not establish generally negligible bias or calibrated conditional
  ASFR recovery. Retain this limitation after acceptance.
- **Milestone 2:** the report includes demographic assessment and R formula-fallback
  timing; D050 records user acceptance. Poor-pool waiting and previous-father
  continuity remain deferred.

## D046: Optional union summaries in Lexis counts

- **Status:** union exposure and formation summaries requested by the user
  October 16, 2026; classification at formation explicitly accepted that day.
- **Decision:** expose partnered/unpartnered person-years and union formation
  counts for both sexes through an optional Lexis summary extension. Include
  maternal births by status to support conditional fertility denominators.
  Union status is current partnership, not first-marriage history or legal
  marital status. Births at the union start are partnered, including those
  triggering on-demand formation; use half-open union episodes.
- **Implementation:** superseded by D047 below. Initially, `include_unions = TRUE` added columns to the existing
  period/age/sex grid. Count repeat union starts for both spouses at their
  respective ages. Integrate coverage without duplicating overlapping
  person-time; existing overall measures remain unchanged. Require a population
  object with union histories. Missing earlier records do not prove historical
  singlehood; users must choose windows with complete follow-up.
- **Validation:** worked exact exposure, boundary births, spouse-specific ages,
  widowhood/remarriage, overlapping and zero-duration episodes, empty histories,
  conservation, and unchanged default summaries. No demographic calibration or
  milestone acceptance is implied.

## D047: Long union-status Lexis summaries and marriage rates

**Current status:** Implemented with D048 rate-column names and the October 17 pooled
label `overall` replacing `total`.

- **Status:** accepted by the user October 16, 2026; supersedes D046's wide layout.
- **Decision:** `include_unions = TRUE` gives `total`, `partnered`, and
  `unpartnered` rows identified by `union_status`. Split exposure, births,
  and deaths, with ASFR/ASMR using each row's exposure. `unions` is NA on
  total/partnered rows; on unpartnered rows count formations for each spouse.
  Add `marriage_rate = unions / person_years` on unpartnered rows, NA elsewhere
  and whenever exposure is zero. All rates retain person-year units.
- **Observation convention:** deaths use status immediately before death,
  including partnered status for a death ending a union. Births retain D046's
  half-open status. See the scientific contract and worked analysis example.
- **Compatibility:** default FALSE stays unchanged. Existing TFR status grouping
  handles long results without changes to its API. Life tables preserve status
  grouping too. The total duplicates status components and must not be added
  to them. Complete histories are still required.

## D048: Consistent rate column names

- **Status:** requested by the user October 16, 2026.
- **Decision:** rename Lexis rate columns to `rate_fertility`, `rate_mortality`,
  and `rate_marriage` (formerly `asfr`, `asmr`, and `marriage_rate`). This applies
  to both default and union-status summaries. TFR and life-table extraction
  consume the new names; process schedules retain `rate` and life tables
  retain `mx`. Units and calculation rules are unchanged.

## D049: Direct-rate approximation for USA 2001 marriage inputs

- **Status:** requested by the user October 16, 2026.
- **Decision:** provide separate female/male marriage schedules in
  `validation/data`, with `age_start`, `age_end`, and `rate`, using published
  transitions per 1,000 divided by 1,000 directly as rates per person-year.
  Retain source bands, including 70+, and the original combined source table.
- **Limits:** the source denominator is the unmarried population at the start
  of the period, not person-years. This approximation does not apply a
  probability-to-hazard transformation. The female schedule can configure
  current union-search initiation; male rates are reference data. Search
  failures can make realized formation rates differ from these inputs.

The user renamed D047's pooled `union_status` label from `total` to `overall`
on October 17, 2026, to avoid confusion with total fertility rate terminology.
This changes the label in Lexis counts and downstream summaries, without changing
counts, denominators, or demographic behavior.

## Default age-gap matching preferences

Accepted by the user October 18, 2026: `union_formation()` and
`set_union_formation()` default to `create_age_gap_matching_model(0, 5)` when
matching is omitted. The preferred male-minus-female age gap is zero years;
selection weight halves at a distance of five years on either side. Explicit
matching specifications override this default. Initiation mode remains explicit.


## D050: Milestone 2 accepted as complete

- **Status:** accepted by the user October 19, 2026, following review of
  `validation/milestone2.qmd` and its rendered HTML and user revisions.
- **Decision:** mark milestone 2 union formation complete. The user accepts
  the demographic assessment of both initiation modes, status-specific fertility,
  formation rates, union diagnostics, matching-parameter recovery and performance.
  The report and status documentation are to be committed together.
- **Limits:** acceptance retains deferred transition-aware fertility compensation
  (D045), poor-pool waiting (D044), and previous-father continuity (D042). It does
  not establish universal accuracy, independent replication, or completion of
  every review suggestion. No new scientific acceptance thresholds are imposed.
- **Next scope:** groups remain the selected milestone 3 theme. Detailed rules
  require user decisions; milestone completion does not authorize implementation
  of unspecified mechanisms or other roadmap proposals.

## D051: Distribute group support and attach inheritance to fertility

- **Status:** accepted by the user October 20, 2026; design direction only,
  not implemented functionality.
- **Decision:** store group data in the population and configure group effects
  in the relevant demographic processes. Group-specific schedules and
  proportional rate effects belong with their rate processes; group effects on
  partner selection belong with matching. Do not introduce a separate general
  group process for milestone 3.
- **Inheritance:** attach a distinct inheritance model object to fertility to
  specify newborn group assignment and, if supported, ancestry inheritance.
  Allow simple rules through character-string modes rather than a boolean,
  with following the mother, following the father, and random assignment as
  examples. Richer model specifications may be supported through the same
  object, analogous to matching models.
- **Open details:** exact constructors and mode names, random assignment rules
  and probabilities, richer model inputs, missing-parent handling, ancestry
  semantics, and interactions among group effects and existing schedules
  require specification before implementation. No default inheritance rule or
  probability is selected by this decision.


## D052: Initial group schedules, maternal inheritance, and baseline sampling

- **Status:** requested and accepted by the user October 20, 2026.
- **Decision:** permit optional `group` columns in population persons and
  mortality, fertility, and female union-search rate tables. Ungrouped rate
  tables apply to all groups. Grouped tables must cover every population group;
  the user selected an error for missing group schedules.
- **Fallback:** without population groups, select each grouped table's first
  group alphabetically and issue a simulation warning. Do not create population
  group identities from that fallback. Group values accept nonmissing character
  labels or numeric codes; finite numeric codes and nonempty labels are required.
  Alphabetical lookup uses radix order of character representations.
- **Births and matching:** initially inherit each newborn's group from the
  mother automatically. Group introduces no matching barriers or utility terms.
  This is an explicitly requested initial rule, not the deferred configurable
  inheritance model described in D051.
- **Initialization:** add optional named numeric `group_weights` to
  `create_baseline_population()`. Names give groups; normalize finite nonnegative
  weights with at least one positive weight and independently sample every
  person's identity. A named vector is the chosen interface; omitted weights
  preserve the existing generator and its RNG draws.
- **Scope:** group/status-specific fertility compensation uses the existing
  midpoint method. No group transitions, ancestry, proportional multipliers,
  inheritance model API, or group summary extension is included. See
  [groups](groups.md) for rules and worked validation cases.

## D053: Group-specific demographic observations

- **Status:** requested by the user October 20, 2026.
- **Decision:** add optional group counts to `calculate_lexis_counts()` and
  `calculate_population_size()`. Automatically calculate TFR and life tables
  separately by group whenever the rate input includes a `group` column.
- **Interface:** `by_group = FALSE` preserves pooled counting; TRUE requires
  group data and retains a group column in the output. Group labels retain
  their type and first-occurrence order. No extra pooled group rows are added.
  Lexis group counts combine with existing union-status summaries.
- **Existing attribution rules:** births follow the mother (including her
  group), exposure and deaths follow the individual, and union formation counts
  follow each spouse. Unknown-mother births remain in the single period-level
  attribute; no maternal group is inferred from the child. Full histories must
  remain available for cross-group parent and partner references.
- **Validation:** exact exposure and event boundaries, group and union-status
  accounting identities, analytical TFR and exponential life tables, input
  validation, label types, empty tibbles, and unchanged pooled defaults. See
  [groups](groups.md) and [demographic analysis](demographic-analysis.md).

## D054: Group rate multipliers, including conditional fertility

- **Status:** requested and specified by the user October 20, 2026.
- **Decision:** add optional `group_prop_hazard` to mortality, fertility, and
  rate-based union-formation constructors and setters. Named vector entries
  multiply the selected shared or group-specific rates; missing groups use 1.0.
  Factors are dimensionless and not normalized. Require finite nonnegative
  factors with unique, nonmissing, nonempty group names; zero suppresses risk.
- **Fertility order:** the user selected multiplication of input ASFRs before
  existing midpoint birth-spacing compensation. This is equivalent to explicitly
  scaled input group schedules, rather than scaling adjusted eligible hazards.
- **Conditional fertility:** the user requested either a shared named vector
  or a named list of partnered/unpartnered group vectors. The user selected
  1.0 for an omitted status, as for an omitted group. Status lists require
  conditional fertility input rates; a shared vector applies to both statuses
  or to overall fertility. Mortality and search processes accept vectors only.
- **No population groups:** the user selected a simulation warning and ignoring
  multipliers. Existing alphabetical fallback for grouped rate tables remains
  separate. Factors do not create groups or supply missing schedules.
- **Execution:** scale an execution copy of input schedules before native
  construction; retain original tables and provenance in processes. Reject
  nonfinite effective rates or infeasible effective fertility compensation.
  Search factors apply to women and do not modify partner utilities; on-demand
  formation has no search hazard and rejects non-NULL multipliers.
- **Validation and scope:** require fixed-seed equivalence to explicit schedules,
  omitted/unit/zero factors, conditional-status defaults, no-group warnings,
  edited-object validation, RNG preservation, restarts, and empty populations.
  Input-rate summaries still describe stored tables; population-dependent
  factors are applied by simulation, not those summaries. See [groups](groups.md).

## D055: Track ancestry shares at birth

- **Status:** requested and specified by the user October 20, 2026.
- **Decision:** track expected group ancestry separately from assigned identity
  in an optional `ancestry` tibble on `hammel_population`. Without supplied
  proportions, every starting person is 100% of their recorded group, even when
  existing parent histories are present. Otherwise retain supplied proportions.
  No ancestry table is created automatically without population groups.
- **Births:** each component is half the mother's share plus half the father's.
  Missing parental contributions enter an explicit unknown share; parental
  unknown shares transmit by the same averaging rule. Group identity still
  follows the mother. No stochastic genetic transmission is modeled.
- **Restarts and state:** preserve existing proportions, including deceased
  people. Store call-local ancestry state in the engine and append at birth,
  with no genealogy traversal or additional RNG draws. Ancestry does not yet
  affect rates or matching.
- **Table representation:** `person_id`, known shares named `group_<label>`,
  and `unknown`. Require exactly one row per person, finite shares in `[0,1]`,
  and row sums within `1e-8` of one; do not renormalize. Reference groups may
  differ from current identity groups. Align supplied rows to person order;
  retain supplied ancestry even without current group identities.
- **Validation:** worked pure/mixed parents and unknown fathers, multiple
  generations, supplied values and shuffled IDs, restarts, deaths, empty
  populations, numeric labels, run independence, and unchanged demographic RNG
  consumption. See [ancestry](ancestry.md).

## D056: Add group effects to matching scores

- **Status:** requested and specified by the user October 20, 2026.
- **Decision:** attach optional `group_exogamy` effects to matching models, adding their
  log-scale value after all formula terms are evaluated. Named-vector syntax
  `A:B` means female A/male B; `A::B` supplies the same effect in both directions.
  The user clarified that omitted pairs have odds ratio 1.0, hence zero log
  effect. Within-group effects are allowed. Reject overlapping directed entries
  after expansion rather than choosing precedence. A directional table supports
  group labels containing colons.
- **Formula interface:** expose `female_group`, `male_group`, and numeric 0/1
  `same_group` alongside continuous male-minus-female age difference. Categorical
  formulas use explicitly supplied `group_levels` and fixed treatment contrasts;
  group columns and coefficients must not vary with sampled availability.
- **Execution:** pair effects apply after either native age scoring or custom
  formula evaluation, before softmax and diagnostic utility ranking. Both
  initiation modes and nonmarital father choice use the applicable full matching
  model. An explicit `father_matching` replaces general preferences for father
  choice; NULL reuses them. Group preferences require population identities.
- **Scope and validation:** finite utility effects change preference, not
  eligibility or search hazards. Validate analytical probabilities, direction,
  symmetry, factor stability and interactions, native/R agreement, father
  overrides, ranking, errors, empty populations and unchanged RNG consumption
  for zero effects. Multiple identities and exact group prohibitions are not
  introduced. The user requested later diagnostic tools to flag potential
  omissions; these remain future work. See [group matching](group-matching.md).


## D057: Hard group exclusions after candidate sampling

- **Status:** selected by the user October 20, 2026.
- **Decision:** add `prevent_exogamy` to matching models: FALSE by default,
  TRUE to prohibit every different-group pairing, or a character vector of
  cross-group pairs with female-first `:` and symmetric `::` syntax. A
  directional table also supports labels containing colons.
- **Sampling:** draw from the unchanged male pool, then exclude prohibited
  pairs alongside age/kin eligibility without refilling. The user selected this
  to preserve encounter scarcity for minorities and the same encounter frame
  used by finite penalties. Users may increase candidate sample size.
- **Failure and scope:** empty eligible samples use existing retry and
  cancellation rules. Restrictions override finite scores, but do not alter
  existing partners. The applicable matching/father model owns its restrictions;
  NULL father matching inherits the general specification. Diagnostics count
  and rank the final eligible set. No NA or infinite utility effects are added.
- **Validation:** direction/symmetry, no refilling, compatible partners outside
  an unsuccessful sample, complete exclusion, finite-penalty contrast, father
  overrides, retained partners, diagnostics and unchanged RNG when disabled.
  See [group matching](group-matching.md).


## D058: Configurable newborn identity inheritance on fertility

- **Status:** specified and authorized by the user October 20, 2026; maternal
  assignment remains the explicit default when users supply no rule.
- **Decision:** attach a `create_group_inheritance()` object to fertility.
  Parameter-free method names may also be supplied directly as strings to
  fertility constructors/setters, which normalize them to the same model object
  (user requested October 20, 2026).
  Built-ins follow mother, father, random parent (0.5/0.5), same-sex parent,
  opposite-sex parent, known child ancestry, or a strict hierarchy of parental
  identities. Missing required information falls back to the mother's group,
  including missing fathers and unranked parents. Ancestry sampling conditions
  out unknown shares and falls back to mother if entirely unknown.
- **Hierarchy:** highest first in `order`; optional `passing` is the per-birth
  probability of reversing the hierarchy (default zero). Equal identities or
  missing information retain the maternal fallback. Complex hierarchies/ties
  use custom functions.
- **Custom rules:** receive parental groups and ancestry, child ancestry, and
  child sex and logical `parents_in_union` at birth (including on-demand
  formation; FALSE for nonmarital/unknown-father births); return one valid
  identity. The union-status input was requested October 20, 2026; callbacks
  without the new argument or `...` remain supported. Ordinary R RNG draws share the native
  stream. Only copied inputs are exposed; functions must not mutate simulation
  state. Errors and invalid/undeclared returns fail explicitly.
- **Possible groups:** infer starting identities plus reference groups with
  positive ancestry; optional `groups` adds possible outputs. Validate all
  grouped schedules and categorical matching levels before execution, and keep
  a fixed registry during a run. All-zero ancestry columns add no identity.
- **Timing:** calculate sex and ancestry, assign identity, then schedule newborn
  hazards. Assignment does not change ancestry, existing people's identities,
  or default maternal random-number consumption. Rules are reapplied only to
  new births, with supplied proportions retained on restart.
- **Mixed labels:** the user explicitly noted that custom rules can assign AB
  from A/B parents and configure independent AB rates/matching. AB identity
  need not introduce an AB ancestry component. Multiple simultaneous identities
  remain deferred.
- **Validation:** all rules/fallbacks, sampling probabilities, passing endpoints,
  custom/native RNG agreement and error restoration, future identity rate and
  matching coverage, AB neonatal mortality and within-run matching, numeric
  labels, empty populations, and restarts. See [group inheritance](group-inheritance.md).


## D059: Propagate missing ASFRs into TFR

- **Status:** accepted by the user October 21, 2026.
- **Decision:** `calculate_tfr()` returns `NA` if any ASFR in the included
  finite age bands is missing, independently within each period/group/status.
  All-missing bands also yield `NA`; known zero rates still contribute zero.
  Open age bands remain excluded. This supersedes omission of missing ASFR
  contributions from TFR sums.
- **Rationale:** sparse partnered exposure for minority-group teenagers can
  leave age-specific rates unobserved. Summing remaining bands produced
  incomplete TFRs that looked complete and depressed replication means.
  Returning `NA` exposes missing information; it does not remove variability
  from tiny positive exposures. Summaries omitting missing TFRs describe only
  complete estimates and should report their nonmissing count.
- **Validation:** partial/all-missing schedules, zero rates, excluded open
  bands, and isolation across period/group/status strata.
