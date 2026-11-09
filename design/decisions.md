# Decision log

Updated November 6, 2026. Entries retain the choices made at their original
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
| Person-year inputs and hazard sampling | D005–D007, D013, D024; cell scheduling D128 |
| Population representation, sex labels, public tibbles | D012, D015, D016, D023 |
| Process constructors/setters | D017, D026 |
| Noun names for specification constructors | D147 |
| Database rate products and reusable downloads | D018, D025 |
| Demographic extraction and function names | D020–D022, D081 and missing-ASFR clarification |
| Optional Poisson intervals for observed rates | D124; argument name D148 |
| Continuous execution | D024 |
| Fertility spacing (cell accounting) | D026, D101 (superseding D027/D028/D099 compensation) |
| Two-sex union models and union API restructuring (deferred) | D073 |
| Duration-dependent dissolution and symmetric group pairs | D078–D079 |
| Age-specific searches by previous union; separate marriage observations | D080 |
| Actual-hazard calibration contributions | D082, mechanism replaced by cells under D128 |
| Age/duration remarriage in one rate table | D083 |
| Optional duration bands for marriage observations | D084 |
| Cohabitation design and implementation | D085–D086 |
| Union/marriage entry, dissolution and conversion observations | D087–D091 |
| Unpartnered birth extraction | D092 |
| Milestone 4 completion and retained validation evidence | D093 |
| Parity schedule bands and HFD preparation | D094 |
| Open fertility age bands, HFD endpoints and TFR warnings | D095 |
| Missing conditional HFD rates replaced with zero | D096 |
| PATFR helper and automatic TFR dispatch for parity tables | D097 |
| History-authoritative parity initialization | D098 |
| Parity schedule validation and history parity in execution | D099 (correction superseded by D101) |
| Cell-accounted fertility, the only spacing method | D101; unrealized-birth counts D128 |
| Persistent lognormal fertility weights and maternal inheritance | D102 |
| Milestone order and groups-based characteristics scope | D103 |
| Disabled fertility heterogeneity removes output weights | D104 |
| Milestone 5 acceptance | D105 |
| Group-transition inputs are annual hazards | D106 |
| Directed transition schedules, clocks and optional strata | D107 |
| No transition-count target adjustment | D108 |
| Individual transitions core; household moves deferred to regions | D109 as superseded by D114 |
| Transition rate-table columns and common band type | D110 |
| Grouped output always includes transition history tibble | D111 |
| Individual transition implementation, duration and ancestry | D112 |
| Directed transition Lexis counter | D113 |
| Household moves deferred; regional applications under multiple-group milestone 9 | D114, D121, D131 |
| Detailed union partner groups refer to formation | D115 |
| Milestone 6 acceptance; time-varying inputs next | D116 |
| Two-sex marriage model as tentative milestone 12 | D117, D121 |
| Calendar-varying rate inputs | D118 |
| Multi-year HMD/HFD preparation and native HMD shapes | D119 |
| Optional HMD open-age pooling from deaths and exposure | D120 |
| Milestone 8 direction and milestones 9–12 revision | D121, reordered by D131 |
| Calendar-varying flat cohabitation probabilities | D122 |
| Complete matching/inheritance model-list calendar schedules | D123 |
| Detailed and pooled previous-union inputs and summaries | D126 |
| Milestone 7 acceptance and retained validation evidence | D127 |
| Cell scheduling of all rate processes, one clock per process | D128 |
| Weight models (formulas, custom functions, schedules) | D129 |
| Milestone 8 acceptance | D130 |
| Milestones 9–11 reordered: group columns, immigration, then summaries | D131 |
| `stratify_rates()` replaces `group_prop_hazard` | D132 |
| One `by` argument for rate calculators; shared observation engine | D133 |
| Multiple group dimensions: declaration, names, staging | D134 |
| Matching prohibitions written as NA in `group_exogamy` | D135 |
| First public release after milestone 9; immigration and two-sex models post-release | D136 |
| GitHub milestones and Projects; multiple issues and PRs per milestone | D137 |
| Matching calibration within separate markets | D138 |
| `group_exogamy` table form for cross-dimension pair effects (`_x`/`_y` columns) | D139 |
| Household moves for group transitions | D140 |
| Milestone 9 accepted complete | D141 |

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

**Current status:** Superseded by D101 on October 29, 2026: midpoint
compensation was removed and births use cell-accounted fertility. The shared
waiting-time sampler remains for other events. D095 had superseded the
prohibition on upper-open fertility inputs.

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
default, and D072 replaces candidate sampling (exclusive unions are unchanged).

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

**Current status:** Resolved by D101 (October 29, 2026): cell-accounted
fertility uses the current population composition, including women who changed
union type during spacing. Milestone 2 was accepted as complete (D050) under
the earlier deferral.

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

- **Status:** selected by the user October 20, 2026; sampling semantics superseded by D072.
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


## D060: Configurable female fraction at birth

- **Status:** requested and accepted by the user October 22, 2026.
- **Decision:** add `ffab` to `set_fertility()` and `fertility()`, defaulting
  to the existing 0.4886. It is a finite scalar probability in `[0, 1]`,
  applying to every actual birth under that process. Sex precedes group
  inheritance and newborn scheduling; existing people are unchanged.
- **Software contract:** append the argument to preserve positional calls;
  validate at construction, attachment, and execution. Missing fields in older
  process objects use the original default. Preserve one uniform draw per
  birth at all probabilities, including endpoints, and record the effective
  setting in simulation metadata. Default seeded histories remain unchanged.
- **Scope:** no sex-specific fertility input or group/time-varying sex
  probabilities are introduced. This does not establish a cause for the
  inheritance trajectories under investigation.


## D061: Prototype sequential two-sided temporary matching markets

- **Status:** experiment authorized by the user October 22, 2026; superseded
  by D062 and removed with the two-sided method in D067.
- **Decision:** add optional `matching_method = "two_sided"` to union formation,
  preserving `"female_choice"` as the default. Sample M - 1 alternate living,
  unpartnered women with the same age cutoff as ego's M eligible sampled men.
  Uniformly reduce men if too few women exist. Do not weight alternate women
  by demographic hazards or apply kinship to their hypothetical pairs.
- **Procedure:** repeatedly flip a fair sex coin, uniformly draw an initiator,
  and select their partner with the existing conditional-logit pair utilities.
  Remove hypothetical pairs locally and stop when ego matches. Enforce hard
  group exclusions; an alternate woman with no allowed man leaves unmatched.
  All scores retain female-first/male-second orientation.
- **Scope:** shared selection for unions and nonmarital fathers, preserving
  explicit father preferences, actual eligibility, and existing event timing.
  No hypothetical pair changes population state. The user requested runtime
  comparisons and will assess input benchmarks using milestone 3 analyses.
- **Specification:** `docs/two-sided-matching.md`, removed in D067; see Git
  history.


## D062: Fixed-pool reciprocal probabilities for the two-sided experiment

- **Status:** implementation authorized by the user October 22, 2026;
  removed by the user October 23, 2026 (D067).
- **Decision:** supersede D061's sequential temporary depletion. Retain its
  sampling, eligibility, shortage, hard exclusion, and father-selection rules.
  Score the fixed market, sum ego's row-normalized partner probabilities and
  the men's column-normalized probabilities of choosing ego, normalize, and
  draw once. This integrates the agreed fair-coin, uniform-initiator procedure
  with replacement until ego participates.
- **Custom formulas:** evaluate each woman's allowed male row once, then use
  those fixed pair utilities for both normalizations. No column-specific
  re-evaluation of data-dependent transforms occurs.
- **Validation:** analytical probabilities, independent rejection simulation,
  extreme utilities, exclusions, native/R parity, actual-state invariants,
  default compatibility, and measured timing against the previous prototype.
- **Specification:** `docs/two-sided-matching.md`, removed in D067; see Git
  history.


## D063: Geometric rescaling experiment removed

- **Status:** tested and removed at the user's request October 22, 2026.
- **Experiment:** replace D062's additive probabilities with exponentiated
  pair scores divided by the square root of their row sum times column sum,
  then normalize across ego's male candidates. Sampling and eligibility were
  unchanged. Fixed-matrix cross-product odds were preserved analytically.
- **User-reported assessment:** aggregate log odds -5.79 (SE 0.0165), against
  -6, and directional ratio 0.943 (SE 0.0140), against 1. The fixed-matrix
  property did not establish aggregate recovery or directional balance in
  the simulation. The user requested removal and further consideration.
- **Current behavior:** restore D062's additive fixed-pool selection. Its
  documented limitations remain; no replacement mechanism is accepted.

## D064: Two-sided matching is the default for new specifications

- **Status:** implemented at the user's request October 22, 2026; superseded
  by D067.
- **Decision:** default `union_formation()` and `set_union_formation()` to
  `matching_method = "two_sided"`, using D062's existing additive fixed-pool
  probabilities. Retain explicit `"female_choice"` and its inheritance warning.
- **Rationale:** preliminary user-reported conditional-logit parameter recovery
  and improved directional pairing balance support using reciprocal competition
  by default. Aggregate log cross-products need not equal twice symmetric input
  coefficients. This default change does not establish general scientific
  validation; further parameter-recovery and inheritance assessment remains open.
- **Compatibility:** new calls omitting the argument use two-sided matching and
  may change seeded histories and runtime. Existing specifications explicitly
  storing a method retain it; older objects without the field retain female
  choice. Initiation hazards, birth timing, and existing unions are unchanged.
- **Validation:** omitted and explicit two-sided matching agree in both
  initiation modes for constructors and setters; explicit female choice agrees
  with legacy objects. `TZ=UTC Rscript -e 'devtools::test(reporter = "summary")'`
  passed with no failures or warnings; the opt-in live database test was skipped.

## D065: Candidate sample size belongs to the matching model

- **Status:** implemented at the user's request October 22, 2026; the
  default size of 250 is superseded by 1000 in D067; `sample_size` removed in D072.
- **Decision:** both matching-model constructors accept `sample_size = 250`.
  Remove the argument from union formation and its setter, without a deprecated
  alias. The user prefers direct API changes while the package is unpublished.
- **Semantics:** sample up to this many men without replacement before screening;
  zero uses the full pool. Both matching methods use the selected model's size.
  A separate nonmarital father model supplies its own sample size; NULL inherits
  the general model. This supersedes the earlier shared process-level size and
  default 1000. Initiation hazards and score formulas are unchanged.
- **Examples and limits:** sizes 250 for general matching and 50 for father
  matching apply to their respective searches. Screening does not refill samples;
  female shortage can further reduce two-sided markets. Reducing the sample
  changes encounters and seeded outcomes, not just runtime. No scientific
  convergence at 250 is claimed. Rebuild older specifications using the new API.
- **Validation:** constructor/edited-object validation, screening and full-pool
  behavior, independent union/father sample sizes, NULL inheritance, and both
  initiation/matching modes pass. After a clean optimized native build,
  `TZ=UTC Rscript -e 'devtools::test(reporter = "summary")'` passed with no
  failures or warnings; the opt-in live database test was skipped. The report's
  three explicit size-250 settings were moved into its matching models without
  executing the report or downloading rates.

## D066: Experimental weighted-choice matching

- **Status:** specified with and authorized by the user October 23, 2026;
  made the default in D067, which also supersedes its retention of two-sided
  matching.
- **Decision:** add `matching_method = "weighted_choice"`: female choice among
  screened candidates, with each man's selection weight multiplied by a weight
  for his group. Weights are recalibrated at every search, in continuous time,
  so expected male union flows by group meet a target: each unpartnered man
  receives the per-capita intensity of an unpartnered woman of his group
  (crude counts at or above `min_male_age`), interpolated toward female group
  shares by the matching model's `imbalance_absorption` (default 1).
- **Intensity:** search hazards for rate-mode searches; unpartnered fertility,
  zero during birth spacing, for on-demand searches. Rate-mode nonmarital
  fathers use separate fertility-driven weights with the father model.
- **Rationale:** female choice leaves male group flows unspecified, producing
  directional imbalance; additive two-sided probabilities (D062) attenuate
  preferences; per-market Sinkhorn balancing attenuated age and group effects
  in an exploratory static analysis. Group-only weights leave the preference
  association unchanged. Age-standardized targets and male input marriage
  rates were not selected: they would impose a male age pattern or targets
  women's timing cannot meet.
- **Retained:** additive two-sided matching remains available and the default;
  explicit female choice is unchanged. With one group the new method equals
  female choice exactly.
- **Consequences:** in two-group populations directed effects change weights,
  not expected group flows; directional counts are sensitive to the sexes'
  group composition among unpartnered people. Approximate weights are used,
  with one warning per simulation, when no exact solution exists.
- **Specification:** [weighted-choice matching](weighted-choice-matching.md).
- **Validation:** see that document; the full suite passed after
  implementation. This also records the user's October 23 clarification that
  candidate `sample_size` is a speed shortcut, not an encounter model, which
  qualifies the encounter-scarcity rationale recorded under D057.

## D067: Weighted choice is the default; two-sided matching removed

- **Status:** decided and requested by the user October 23, 2026, after
  reviewing the replicated results recorded in
  [weighted-choice matching](weighted-choice-matching.md).
- **Decision:** `union_formation()` and `set_union_formation()` default to
  `matching_method = "weighted_choice"`, and both matching-model constructors
  default to `sample_size = 1000`. Remove the additive two-sided method
  (D061/D062/D064) entirely: its engine code, inspection functions, tests,
  specification, and benchmark scripts and results. Explicit
  `"female_choice"` remains available; older objects without a method field
  retain it. Objects storing `"two_sided"` are rejected with an instruction to
  rebuild, following the user's preference for direct API changes while the
  package is unpublished.
- **Rationale:** the user judged two-sided matching a dead end: its additive
  probabilities attenuate the input association. Weighted choice balanced
  directional pairing and recovered preferences in the user's runs. A
  250-man sample attenuated pooled associations under both remaining methods;
  1000 removed most of this at modest cost for weighted choice.
- **Consequences:** default seeded histories change. With hard exclusions that
  leave a male group unreachable, the default now warns that male-group
  targets could not be met. The milestone 3 report's chunk using
  `matching_method = "two_sided"` must be changed by the user.
- **Validation:** full suite passed with no failures or warnings; the opt-in
  live database test was skipped. See [testing](testing.md).

## D068: Matching method belongs to the matching model

- **Status:** requested by the user October 23, 2026.
- **Decision:** move `matching_method` from `union_formation()` and
  `set_union_formation()` to `create_matching_model()` and
  `create_age_gap_matching_model()`, defaulting to `"weighted_choice"`. A
  separate `father_matching` model applies its own method; NULL inherits the
  general model's method, as for `sample_size` (D065) and
  `imbalance_absorption`.
- **Rationale:** the method describes how a matching model turns utilities
  into choices, not how unions are initiated. It sits with the other
  matching-model settings and allows, for example, weighted unions with
  female-choice nonmarital fathers.
- **No legacy handling:** at the user's direction, there are no compatibility
  checks while the package is unpublished. Matching models lacking
  `matching_method` or `imbalance_absorption` fail ordinary validation; a
  `matching_method` field left on an old union process is ignored. The earlier
  rule that process objects without a method retain female choice is removed.
- **Engine:** the calibrator is built when either model is weighted in a
  grouped population; each search applies weights only if its own model is
  weighted.
- **Validation:** constructor defaults and validation, removal of the union
  argument, NULL father inheritance of the method, and mixed union/father
  methods; the full suite passed with no failures or warnings, the live
  database test skipped. The milestone 3 report's chunks passing
  `matching_method` to `set_union_formation()` must be updated by the user.

## D069: Advisory warnings for matching-method interactions

- **Status:** requested by the user October 23, 2026.
- **Weighted choice and directed effects:** constructing a weighted-choice
  matching model warns once when `group_exogamy` contains directed pairs whose
  two directions differ (including one-sided entries) or the formula uses
  `male_group`. Weighted choice pins male group flows, so with two groups only
  each pair's sum affects pairing, directed effects are partly reflected with
  more groups, and `male_group` main effects are absorbed entirely. Symmetric,
  equal-valued, within-group, and `same_group` terms do not warn.
- **Female choice and inheritance:** a model with fertility, union formation,
  and a population of at least two groups warns when a female-choice union or
  father model is combined with the `mother`, `father`, `same_sex`, or
  `opposite_sex` inheritance rule. The user proposed warning for every built-in
  rule except `random_parent`; `ancestry` and `hierarchy` were excluded because
  they treat mother and father symmetrically, so directional imbalance does
  not change their outcomes. Custom rules are not checked. The check runs in
  `set_process()` and `set_population()`, so it fires while building a model;
  `simulate()` re-validates without repeating it.
- **Validation:** warning and silent cases for both checks, each setter path,
  father-only female choice, and no repetition at simulation; tests that
  intentionally combine these settings silence only the relevant warning. The
  full suite passed with no failures or warnings, the live database test
  skipped.

## D070: Weighted-choice refinements: strict endogamy and weight diagnostics

- **Status:** requested by the user October 23, 2026.
- **Strict endogamy:** a matching model with `prevent_exogamy = TRUE` skips
  weighted-choice calibration, because a group weight cannot change a choice
  among same-group candidates. It reproduces female choice exactly and no
  longer raises target warnings for unreachable groups. Partial restrictions
  still calibrate.
- **Diagnostics:** with `diagnostics = TRUE`, a weighted general matching
  model in a grouped population records `.log_weight_<group>` on new unions:
  the applied log group weights (reference group 0; unweighted searches 0;
  earlier rows NA). These serve as exact conditional-logit offsets.
- **Coverage:** weighted choice with identities created during a run (custom
  AB inheritance), previously tested only under female choice after D069, is
  now tested directly, including a male group with no counted women.
- **Validation:** full suite passed with no failures or warnings; the live
  database test skipped.

## D071: Weighted-choice calibration accuracy fixes

- **Status:** accepted by the user October 23, 2026, after an external
  review (Codex) reported two issues, both reproduced here.
- **Unmet targets reported as success:** women whose only allowed men belonged
  to zero-target groups (log weight -1000, underflowing to zero) dropped out
  of the solve and final flows, so an impossible target was reported as met.
  Final flows are now computed in log space including those women, and any
  miss above 1e-8 marks the search approximate, triggering the warning.
  Actual partner selection was already correct.
- **Narrow rate bands:** bin intensities used the hazard at each bin's centre
  age, so a band narrower than a quarter-year bin could be missed, silently
  returning unit weights. Intensities now average the hazard over each bin's
  age span. The whole-bin age cutoff and bin-centre pair utilities remain
  documented approximations; measured sensitivity is in the specification.
- **Validation:** both reported cases added as tests (the impossible target is
  flagged; a 0.05-year band at a woman's age gives the analytical weight), plus
  one-year bands against an exact person-level reference. Full suite passed
  with no failures or warnings; the live database test skipped.

## D072: Exact full-pool partner draws replace candidate sampling

- **Status:** decided by the user October 23, 2026, after a prototype
  benchmark and accuracy checks.
- **Decision:** remove candidate sampling. Every partner and nonmarital father
  search draws exactly from conditional-logit choice over all eligible men, for
  both matching methods (female choice uses zero group weights). Remove
  `sample_size` and its sampling code and test hooks, and remove the
  `.n_eligible` and `.utility_rank` diagnostics; `.log_weight_<group>`
  remains. Nonmarital fathers use the same draw with the father model's
  settings.
- **Method:** unpartnered men are kept in cells by group and quarter-year birth
  bin. A cell is chosen in proportion to its size times the largest weight any
  member could have for the woman, a member uniformly, and he is accepted with
  probability exp(U - that maximum); only he is screened for age and kinship.
  Rejection sampling makes the result exact; after 200 rejections the eligible
  pool is enumerated. The native age model's within-bin maximum is exact.
  Custom R formulas are bounded on a 1/32-year grid (accepted by the user for
  now); bound violations are counted and warned about.
- **Rationale:** the user's runs showed a 250-man sample attenuated pooled
  associations; sampling was only a speed shortcut, so exact draws give the
  intended model without a size to tune. A prototype matched exact full-pool
  probabilities in frozen-pool tests and full-pool sampling over 40 seeds.
- **Consequences:** a search fails only when no eligible man exists in the
  whole pool; this supersedes the sample-then-screen semantics of D037 and
  D057, the encounter-scarcity rationale, and D065/D067's sample sizes.
  Seeded histories change. Random-number use depends on how steeply utility
  varies within a bin. The user noted possible future options for a minimum
  match threshold or a woman cancelling her search; none is implemented.
- **Performance:** optimized builds, one seed: female choice 0.10 s and
  weighted choice 0.25 s for 20,000 persons over 50 years (0.30 s and 0.47 s
  with 1,000-man samples); 0.39 s and 1.12 s for 50,000 persons over 100 years
  (1.39 s and 2.18 s before).
- **Validation:** full suite passed with no failures or warnings; the live
  database test skipped.

## D073: Keep current union formation; defer two-sex models and API restructuring

- **Status:** decided by the user October 24, 2026.
- **Decision:** keep weighted-choice matching for the current milestone. Do
  not add user-specified directional (gender-asymmetric) intermarriage ratios,
  male-group union propensities or directional ratio targets to weighted
  choice; users who need gender-asymmetric intermarriage should be served by a
  future, fuller harmonic-mean approach. Leave the `union_formation()` and
  `set_union_formation()` API unchanged until a single working implementation
  incorporates all union-formation components.
- **Rationale:** with balanced sexes, weighted choice forces equal directional
  counts (AB - BA equals the difference between male and female group shares),
  so asymmetry would require changing male-group targets, and female timing
  still could not respond to male availability (no marriage squeeze). The user
  distinguished models whose union event is independent of matching, which
  preserve input rates, from full two-sex models, which allow composition
  feedback without necessarily preserving rates. Cohabitation and nonmarital
  father selection must be designed consistently with any two-sex model, so
  restructuring the API now would be premature.
- **Future work, not authorized:** a participation function giving the
  probability that a woman forms a union under current market conditions; a
  Schoen harmonic-mean model, possibly modified along Dagsvik's lines; possibly
  Choo–Siow; a loose translation from log odds ratios to harmonic-mean
  parameters; a union-model abstraction for `set_union_formation()`.
- **Evidence:** discussion notes and calculations in
  [two-sex union models](two-sex-union-models.md). No code changed.


## D074: Replace imbalance absorption with male target adjustments

- **Status:** accepted by the user October 26, 2026.
- **Decision:** remove `imbalance_absorption` from matching constructors and
  execution. Always use the mirrored male baseline (formerly kappa = 1), then
  normalize baseline[g] * exp(delta[g]). Expose `delta` as a named finite
  numeric vector, omitted groups zero, on each matching model, including
  separate father models. Zero adjustments retain the previous default.
- **Rationale:** kappa cannot tune direction when female and mirrored male
  shares match. Delta can represent systematic directional differences while
  preserving static group-only cross-products. It does not guarantee a
  pair-specific ratio or change female event timing.
- **Supersession:** replaces D066's target interpolation and D068's absorption
  argument, and supersedes D073's deferral of male-group propensities only.
  Full two-sex models and participation feedback remain deferred.
- **Specification:** [weighted-choice matching](weighted-choice-matching.md).


## D075: One matching mechanism and symmetric group-pair effects

- **Status:** accepted by the user October 26, 2026.
- **Decision:** remove `matching_method` and explicit female-choice execution.
  All partner and father matching uses calibrated male-group weights with
  `delta` target adjustments (D074). Preserve unit-weight shortcuts for one/no
  groups and strict endogamy, where group weights cannot affect selection.
- **Pair effects:** `group_exogamy = c("A:B" = beta)` applies beta in both
  directions. Reject `::` syntax and repeated unordered pairs, even if their
  values agree. A pair table retains `female_group`, `male_group`, `log_odds`
  columns, but each row now specifies an unordered pair and expands both ways;
  within-group entries expand once. Directed pair effects are removed.
- **Scope:** formula predictors and hard eligibility restrictions retain their
  existing contracts; directional hard exclusions are not preference effects.
  Separate father models retain their own preferences, restrictions and delta.
  No event timing, rates, inheritance, or demographic acceptance changes.
- **Validation:** symmetric vector/table expansion, reverse-pair rejection,
  native/R score agreement, unchanged delta analytical cases, both initiation
  modes, father overrides, strict endogamy and RNG/restart invariants.
- **Supersession:** replaces D068's method selection and D069's directed-pair
  warning; symmetric pair effects replace D056's directed utility interface.


## D076: Default to random-parent group inheritance

- **Status:** accepted by the user October 27, 2026.
- **Decision:** default newborn identity inheritance to `random_parent`, selecting
  either parent's group with equal probability. Retain maternal fallback when
  the father is unknown. Equal parental identities and fallback cases consume
  no identity draw. Explicit `mother` remains available.
- **Rationale:** the user selected the default after assessing inheritance
  scenarios and sensitivity to shared starting populations. This choice does
  not establish that random-parent inheritance eliminates baseline sensitivity.
- **Compatibility:** missing inheritance fields use the new default; stored
  explicit maternal rules remain maternal. Ungrouped simulations remain valid
  and consume no additional identity draws. Grouped default histories can change.
- **Specification:** [group inheritance](group-inheritance.md).


## D077: Milestone 3 accepted as complete

- **Status:** accepted by the user October 27, 2026, following report review,
  user revisions and a closing request while the user rerenders the report.
- **Decision:** close milestone 3 (groups), based on the user's assessment in
  [milestone3.qmd](../validation/milestone3.qmd). Scope includes group-specific
  schedules and multipliers, group observations, ancestry, weighted matching
  with symmetric pair preferences and delta, exact full-pool draws, and
  built-in/custom newborn inheritance, including emergent mixed identities.
- **Evidence:** the report assesses rates, matching preferences and directional
  ratios, inheritance and baseline sensitivity, custom probabilities, performance,
  and the test suite. The user's conclusion is that measured behavior meets
  expectations while computation remains fast. This documentation update is
  not an independent rerun or a claim that every review suggestion was applied.
- **Retained limitations:** small or imbalanced available pools can produce
  unexpected directional ratios; delta does not fix pairwise ratios over time.
  Calibration uses age bins and omits kinship; pooled cross-products need not
  exactly equal group-only benchmarks. No marriage-squeeze feedback is present.
  Fertility-transition bias (D045), sparse exposure and baseline sensitivity
  remain documented. No new scientific acceptance thresholds are adopted.
- **Next scope:** proposed dissolution/remarriage and full two-sex models remain
  future work, requiring explicit specifications and user decisions. Fixed newborn
  assignment is deferred to group transitions; mixed-parent assignment to a
  chosen category remains possible through custom inheritance. No fixed method
  is implemented or selected. Completion does not authorize future features.
- **Provenance:** revised conclusion at source commit
  `7514916bb8cf0c06faef24dad7fa013304093779`; final rendering is handled by the
  user. No final HTML hash or combined acceptance commit is assigned here.


## D078: Duration-dependent dissolution carried by women

- **Status:** accepted by the user October 27, 2026; dissolution precedes work
  on separate remarriage rates.
- **Decision:** use a continuous annual dissolution hazard, indexed by elapsed
  union duration in years. One clock per union is carried by its female partner;
  there is no second male risk. Shared schedules and optional schedules and
  proportional hazards were initially by the woman's group; D079 supersedes
  this with symmetric pair stratification.
- **Ties:** death first, then birth, then dissolution for the same woman.
  A tied birth occurs while partnered; the pending dissolution follows if still
  scheduled. Births preserve the union's pending dissolution date. Births
  exactly at a recorded dissolution end count as partnered in summaries;
  exposure remains half-open and death-ending birth boundaries are unchanged.
- **Rationale:** duration follows the SOCSIM baseline idea while retaining
  Hammel's person-year hazard units and female-driven organization. Under
  identical couple hazards carrier sex does not change risk; group-specific
  schedules now follow the unordered pair under D079.
- **Supersession:** replaces the October 16 female-age-specific dissolution
  direction in the roadmap. Input duration bands use `dur_start` and `dur_end`
  (matching the existing `age_start`/`age_end` naming pattern).
- **Scope:** optional separate dissolution process. Existing unpartnered search
  rules continue after dissolution; separate remarriage schedules, cohabitation,
  couple covariates and general transition-aware fertility compensation remain
  deferred. No demographic rates are defaulted.
- **Specification:** [union dissolution](union-dissolution.md).


## D079: Symmetric group-pair dissolution schedules and observations

- **Status:** accepted by the user October 27, 2026, replacing the initial
  request for separate partner axes with gender-symmetric pairings.
- **Decision:** dissolution rate strata are unordered pairs such as A:A, A:B,
  and B:B. A:B applies in both directions. Use `group_pair` on rate tables and
  unordered pair names on `group_prop_hazard`. Omitted factors use 1; there is
  still one annual duration-dependent hazard per union, carried by the woman.
- **Observations:** `calculate_union_counts(by_group = TRUE)` returns the same
  symmetric `group_pair` strata, counting each episode once. No partner-specific
  strata or multiplier-list interface is introduced.
- **Validation/coverage:** canonicalize reversed labels, reject conflicting
  reverse multiplier keys and overlapping pair bands, and require schedules
  for all possible unordered identity pairs. No missing-pair risk is invented.
- **Supersession:** replaces D078's woman-only group schedule and multiplier
  interpretation, and woman-only union-count grouping. Other process `group`
  columns and their multiplier conventions are unchanged.
- **Specification:** [union dissolution](union-dissolution.md).


## D080: Age-specific union searches by previous union and separate observations

- **Status:** accepted by the user October 28, 2026. Duration-dependent
  remarriage remains deferred.
- **Decision:** add optional `previous_union` to age-specific female search
  rate tables, using `none`, `dissolved`, and `widowed`. Require all three
  categories within each scheduled group. Without the column, apply the same
  rates to all previous types. Group multipliers retain their current form.
- **State:** the most recent recorded episode determines the unpartnered
  category. Cache state independently per simulation, recover it from complete
  restart histories, and update it on endings. No public person column or male
  search hazard is introduced. No recorded episode means `none`, conditional
  on supplied history; missing histories cannot establish never-partnered status.
- **Calibration:** maintain category counts within female group/age bins and
  sum their category-specific intensities before the existing weighted-choice
  solver. Matching preferences and mirrored male-group targets keep their
  current definition; no separate male previous-status targets are added.
- **Observation pooling clarification (October 28):** default `by_previous_union = FALSE`
  pools counts and exposure and omits the category column. TRUE returns the three
  categories, independently of `by_group`; input rate stratification does not
  choose the observation setting.
- **Observations:** add `calculate_marriage_counts()` for unpartnered
  person-years, starts, and `rate_marriage` by sex and previous category, with
  optional person group. The user explicitly requested removing `unions` and
  `rate_marriage` from `calculate_lexis_counts()` rather than keeping compatibility
  columns. Its status-specific birth, death and exposure summaries remain.
- **Validation convention:** exclusive histories are required for previous
  categories; reject overlaps rather than invent a preceding status. Touching
  episodes are valid; tied starts follow union ID order. Historical absence
  retains its documented interpretation, rather than imputing past unions.
- **Supersession:** updates D046/D047's location of formation counts and
  marriage rates and D078/D079's deferral of separate remarriage schedules.
  No duration-remarriage clocks, cohabitation, cooldowns, former-partner exclusion,
  or fertility-transition correction is added.
- **Specification:** [age-specific remarriage](remarriage.md).


## D081: Name demographic observation helpers for their returned rates

- **Status:** accepted by the user October 28, 2026.
- **Decision:** rename `calculate_lexis_counts()` to `calculate_demog_rates()`,
  `calculate_marriage_counts()` to `calculate_marriage_rates()`, and
  `calculate_union_counts()` to `calculate_dissolution_rates()`.
- **Scope:** names change; arguments, exposure, event counts, rate definitions
  and returned columns retain their existing contracts. Counts and exposure
  remain available alongside rates. Old function names are removed, without
  compatibility aliases. Historical decision entries and accepted validation
  reports retain their original terminology; current docs and examples use
  the new names.
- **Supersession:** updates public function names in D022, D046–D048,
  D079 and D080 without changing their demographic rules.


## D082: Aggregate actual female hazards using contribution and boundary updates

- **Status:** implementation/performance trial authorized October 28, 2026.
  The user explicitly selected actual current hazards for calibration; no
  milestone 4 scientific acceptance is inferred. *Current (November 5, 2026):*
  actual hazards remain the calibration rule, but D128 replaced the per-woman
  tracker and its boundary queue with contributions taken from the search and
  fertility cells.
- **Decision:** maintain individual current female search intensity and sum it
  into group/birth-cohort calibration cells. Replace the previous category
  count arrays and age-band averaging with contributions updated at age-band
  boundaries, birth-spacing release and demographic state changes. Maintain
  separate marriage-search and fertility-driven totals with one pending
  bookkeeping update per available woman.
- **Interpretation:** mirrored male targets still use actual group-level
  female intensities. No reference male schedule or male previous-status strata
  are introduced. Matching utility/cutoff bins and the exact partner draw remain.
  Female hazard evaluation now uses actual age instead of bin averages.
- **Scheduling:** bookkeeping consumes no RNG and does not resample individual
  demographic events. New weights can alter partner selections and subsequent
  grouped histories. Ungrouped/strict-endogamy paths skip the tracker. Pool exits
  erase pending updates; restarts reconstruct contributions from histories.
- **Scope:** apply to existing age schedules only and benchmark optimized old
  and new builds. Duration-based remarriage, parity, region, custom state
  mechanisms and transition-aware fertility compensation remain deferred.
- **Supersession:** replaces D080's previous-category count representation and
  the rate-bin averaging introduced under D071, while preserving their public
  schedule inputs and all other calibration conventions.
- **Specification:** [current-hazard calibration](calibration-intensity.md).


## D083: Age or duration remarriage schedules in one rate data frame

- **Status:** accepted by the user October 28, 2026.
- **Decision:** retain a single rate data frame; optional `dur_start`/`dur_end`
  accompany required age columns. Each row uses one complete clock pair with
  both unused boundaries NA. `none` must use age; `dissolved` and `widowed`
  each independently choose age or duration, consistently for all bands and
  groups within that category. No simultaneous age/duration risk surface and
  no list-based schedule API are introduced. All-age tables may omit duration
  columns. Existing three-category coverage and group multipliers remain.
- **Origin/timing:** duration starts at the most recent union ending, whether
  dissolution or partner death, and resets at every new ending. Reconstruct
  origin and category from complete histories on restart. Search waiting times
  integrate annual hazards on the chosen clock. Gaps have zero risk; zero-rate
  bands can express delays, but no default delay or age constraint is invented.
- **Calibration:** D082's per-woman contribution tracker evaluates the search
  schedule using its own age/duration origin, independently of birth-cohort
  matching cells and fertility ages. The same bounded queue handles duration
  boundaries. Male targets, utilities, exact draws and event priority remain.
- **Observations/scope:** current marriage summaries remain age-specific with
  optional previous category. Duration-band observation denominators need a
  separate interface decision; no age/duration surface or new observation
  default is silently added. No cohabitation, former-spouse exclusion, parity,
  calendar-varying rates or fertility-transition correction is introduced.
- **Supersession:** implements the duration-remarriage portion deferred in
  D080/D082 while preserving their other rules. The duration clock is based
  on the latest ending, not first marriage or the beginning of the simulation.
- **Specification:** [remarriage](remarriage.md).


## D084: Duration-binned observations for previously partnered people

- **Status:** accepted by the user October 28, 2026, including automatic
  previous-category stratification when duration breaks are supplied.
- **Decision:** add optional `dur_breaks = NULL` to `calculate_marriage_rates()`.
  When supplied, `none` uses `age_breaks`, while `dissolved` and `widowed` use
  duration since the latest ending. Automatically include/separate
  `previous_union`, even when its flag is FALSE. Return both age and duration
  boundary pairs with unused columns NA; do not pool unlike clocks.
- **Semantics:** split unpartnered person-years and formation events across
  half-open period and selected-clock bands, resetting duration at each ending
  and using pre-window history. Duration rows include all ages. Preserve sex
  and optional individual-group strata, censor at death/observation end,
  count zero-duration starts, and return NA rates at zero exposure. Append Inf
  after finite final duration breaks, matching existing observation conventions.
- **Compatibility:** NULL preserves existing age-only output, pooling settings
  and calculations. No automatic inference from the original model is needed.
- **Supersession:** implements the duration-observation interface deferred in
  D083. No age-by-duration surface or simultaneous clock stratification is added.
- **Specification:** [remarriage](remarriage.md).

## D085: Cohabitation classification and accepted stretch-goal direction

- **Status:** cohabitation design decisions accepted by the user October 28,
  2026; not implemented. Remaining interface/history choices are unresolved.
- **Entry:** retain overall partnership-search waiting times. At each attempt,
  an age/group/previous-union probability table selects intended cohabitation
  versus direct marriage before matching. Support matching roles for marriage,
  cohabitation and unpartnered father selection, sharing a common default with
  optional overrides. Type-specific calibration uses the corresponding female
  intensities. A failed match leaves the woman unpartnered and reschedules an
  attempt; its intended type is drawn anew when it occurs.
- **Conversion:** an existing cohabiting couple may convert to marriage using
  an age- or cohabitation-duration-specific annual hazard, optionally by group.
  Keep the current partners; conversion does not invoke matching. The user
  chose to omit previous-union stratification from conversion inputs.
- **Previous category:** default to SOCSIM-style recovery of earlier legal
  marital history when a cohabitation ends. Someone with no earlier marriage
  returns to `none` after either breakup or cohabitant death. Earlier legal
  divorce/dissolution restores `dissolved`; earlier legal spouse death restores
  `widowed`.
- **Toggle:** expose a boolean on `union_formation()`, working API name
  `cohabitation_counts_as_previous_union = FALSE`. With TRUE, cohabitation breakup
  assigns `dissolved` and cohabitant death assigns `widowed`, irrespective of
  earlier legal history. This is classification policy, not a conversion-risk
  stratum or an implicit choice of duration-clock origin. Conversion followed
  by a marriage ending uses the marriage history under either setting.
- **Bookkeeping:** preserve enough history to recover legal category and the
  selected policy at initialization/restart; do not discard earlier marriage
  history merely because conversion rates do not use it. Observation policy
  must be made consistent with the chosen simulation classification.
- **Still open:** exact type-probability and conversion APIs/coverage, shared
  versus type-specific fertility/dissolution defaults, representation of type
  changes and durations, and observation interfaces. SOCSIM's table fallback,
  rate inheritance, birth-triggered legalization, counters and apparent queue
  omission are not adopted. No cohabitation implementation or milestone closure
  is implied by these partial decisions.
- **Reference:** [static SOCSIM cohabitation audit](socsim-cohabitation-audit.md).


## D086: Implement optional typed partnerships in existing union processes

- **Status:** authorized October 28, 2026. Rate-mode only; on-demand remains the
  simplified kinship mechanism without cohabitation. No milestone acceptance
  or new empirical default is inferred.
- **API:** union formation accepts age/group/previous-category `cohab_prob`
  (prob, not rate), optional `cohab_matching`, `union_convert_rate` and the D085
  classification toggle. Probability gaps mean zero, as explicitly selected
  by the user. NULL probability preserves untyped execution/RNG. Conversion
  uses age or current-cohab duration, optionally woman's group; missing schedule
  warns and gives no conversions. The conversion argument corrects the prompt's
  `union_covert_rate` typo to `union_convert_rate`.
- **State/history:** typed histories use cohab/marriage and conversion creates
  same-partner linked episodes, resetting dissolution duration. The user selected
  dissolution before conversion on exact ties: death, birth, dissolution,
  conversion, search. Cohabiting people remain unavailable for other matching;
  genuine endings release both partners. Conversion preserves partnership,
  parentage and birth spacing, without a new matching draw.
- **Calibration:** split female search intensity by p/(1-p), with separate role
  weights and probability-boundary updates in the common queue. Father intensity
  stays fertility-driven. Models may share preference tables without sharing
  role-specific intensity/weight solutions.
- **Rates:** fertility input `union_status` is replaced by `union_type`, with
  complete unpartnered/partnered or unpartnered/cohab/marriage coverage per group.
  The user explicitly confirmed matching list names for fertility group hazards;
  reject list keys not present in rates and retain omitted-entry multiplier 1.
  Dissolution can supply both type schedules or one shared schedule, retaining
  symmetric group-pair strata and multipliers. Generic partnered observation
  labels remain union_status; typed fertility inputs/TFR retain union_type.
- **Observation/restart:** exclude conversion from unpartnered entry counts,
  preserve continuous partnered exposure and episode-duration clocks, and recover
  legal category from typed histories. Record the classification policy with
  simulation metadata for observation; changing policy reclassifies supplied
  history rather than reconstructing old policy changes. Typed cohab histories
  require enabled cohabitation for continuation.
- **Specification:** [cohabitation](cohabitation.md). D085's implementation
  deferral is superseded; unresolved general mechanisms remain out of scope.


## D087 — Union entry rate naming and type-specific observation (October 28, 2026)

- **Status:** accepted by the user.
- **Decision:** rename `calculate_marriage_rates()` to `calculate_union_rates()`
  and `rate_marriage` to `rate_union`. Add `by_union_type = FALSE`; TRUE separates
  cohabitation and direct-marriage entries over the same unpartnered exposure.
  Conversion is excluded. Type-specific rates sum to pooled rates at positive
  exposure; exposure must not be summed over the repeated type rows. Untyped
  histories support pooling; requesting a type split requires typed histories.
- **Scope:** retain existing age/duration, previous-union and group options.
  No compatibility alias is retained during pre-release development.


## D088 — Separate marriage entry measure (October 28, 2026)

- **Status:** accepted by the user.
- **Decision:** add `calculate_marriage_rates()` alongside
  `calculate_union_rates()`. Marriage starts include direct marriage and
  cohabitation conversion; exposure includes unpartnered and cohabiting time,
  excludes married time, and returns `marriages` and `rate_marriage`.
  Initial interface supports age, period, sex and optional individual group.
  Untyped histories represent marriages. Partnership-entry semantics are unchanged.


D088 extension: the user requested `by_previous_union` on the marriage helper
as well. Default FALSE pools; TRUE separates exposure and starts by the
pre-entry category. During cohabitation and at conversion this remains the
category held before the cohabitation began. Real endings follow the recorded
cohabitation policy. Duration stratification is unchanged.


## D089 — Dissolution observation by union type (October 28, 2026)

- **Status:** accepted by the user.
- **Decision:** `calculate_dissolution_rates()` gains `by_union_type = FALSE`.
  TRUE separates cohabitation and marriage episodes, crossed with group pairs
  when requested. Each type uses its own union-years and dissolution events.
  Conversion stops cohabitation exposure, is not dissolution, and resets
  duration at the linked marriage start. Require typed histories for splitting;
  retain pooled behavior by default and NA rates at zero exposure.


## D090 — Demographic observations by discovered union type (October 28, 2026)

- **Status:** accepted by the user.
- **Decision:** replace `calculate_demog_rates(include_unions = TRUE)` with
  `by_union_type = TRUE`; return `union_type` strata for unpartnered and types
  discovered in history, without an overall row. Default FALSE remains pooled.
  Legacy untyped histories use partnered/unpartnered. Typed histories require
  exclusive episodes; untyped overlapping coverage remains merged.
- **Timing:** exposure is half-open, death uses the left-limit type, and a birth
  tied to conversion or dissolution retains the ending type. Preserve existing
  formation-boundary birth classification. Unknown mothers retain the single
  period-level unassigned attribute. Derived measures keep type strata separate.


## D091 — Conversion observation helper (October 28, 2026)

- **Status:** requested by the user.
- **Decision:** add `calculate_cohab_convert_rate()` with mutually exclusive
  `dur_breaks` or `age_breaks`, using typed population histories. Conversion
  events use cohabitation union-years as exposure. Duration is since cohabitation
  entry; age uses the woman, matching input conversion hazards. Optional
  `by_group` uses her group. Report `union_years`, `conversions`, and
  `rate_cohab_convert`, keeping empty cells and NA at zero exposure. Dissolution
  and death end exposure without conversion events; marriage is excluded.
- **Implementation:** share the episode exposure accumulator with dissolution
  observation, supplying the appropriate axis origin and event reason. No
  simulation scheduling or demographic defaults change.


## D092 — Unpartnered birth extraction (October 28, 2026)

- **Status:** requested by the user.
- **Decision:** add `get_unpartnered_births()` using the previously discussed
  known-parent birth filter, enriched with parent ages and optional groups.
  Exclude marriage and cohabitation coverage, including births tied to
  dissolution or conversion; ordinary coverage is half-open. Death-ending
  births are unpartnered unless another episode covers the instant. Preserve
  person order, custom columns and repeated parental pairs. Omit unknown-parent
  records and founders. No time window is imposed; full histories determine
  the interpretation. No simulation behavior changes.


## D093 — Milestone 4 acceptance and closeout (October 29, 2026)

- **Status:** accepted by the user. After correcting the report and adding a
  duration-specific conversion assessment, the user designated it the final
  report and requested milestone closeout while rendering it.
- **Scope:** duration-based dissolution with symmetric group-pair effects;
  age/duration search schedules by previous union; actual-hazard calibration
  contributions; optional rate-mode cohabitation, separate matching roles,
  typed fertility/dissolution, linked same-partner conversion and legal-history
  recovery; demographic/entry/dissolution/conversion rate observations and
  unpartnered-birth extraction (D078–D092).
- **Evidence:** [final report source](../validation/milestone4.qmd) and
  [rendered assessment](../validation/milestone4.html). It uses independent
  seeded baselines and replication-level MCSEs to assess rates, group-pair
  dissolution, age/duration remarriage and conversion, typed fertility,
  separate matching contexts, and runtime. The report also executes the offline
  test suite. The final HTML was still rendering when closeout began; no new
  render, simulation or database download is implied by this documentation edit.
- **Limits:** acceptance is the user's assessment, not a guarantee of universal
  empirical accuracy or exact pooled cross-tab recovery. Transition-aware
  fertility compensation (D045), thin-pool behavior, initialization transients,
  and approximate calibration ages remain limitations. Richer couple covariates,
  marriage order, former-partner/father continuity, parity and calendar-time
  schedules are not implemented or authorized by closeout. No executable
  SOCSIM comparison or source copying is claimed.
- **Cleanup:** keep reports 1–4 and their rate/provenance inputs. The author
  removed the exploratory `milestone4_tests.R`; the report supersedes it.
  Move the 14 cited benchmark/recovery scripts and result files into
  `validation/milestone4/`, updating links and reproduction commands. Keep
  rendered-report execution intermediates intact while the author's render is
  active; ignore them in version control. See the validation directory guide.


## D094 — Half-open parity bands and conditional HFD import (October 29, 2026)

- **Status:** the user selected parity_start/parity_end and authorized extending
  the HFD helper after confirming the milestone 4 PR was merged.
- **Representation:** current maternal parity uses half-open intervals and Inf
  for the upper open category. HFD birth orders 1,2,3,4,5+ map to [0,1), [1,2),
  [2,3), [3,4), [4,Inf). Rates are per woman-year in that current parity.
- **Importer:** `get_fertility_hfd(..., by_parity = TRUE)` requests conditional
  period product mi, retaining the default total asfrRR import when FALSE.
  Select one year/revision, omit open age endpoints as before, preserve missing
  conditional cells as NA with a warning, and record product/mapping provenance.
  No interpolation, zero imputation, probability conversion or spacing adjustment
  occurs in this data-preparation operation. Unavailable products/years do not
  fall back to unconditional rates.
- **Boundary:** parity-dependent execution is not yet implemented. Fertility
  constructors reject parity columns rather than silently using them as pooled
  age rates. Initialization, state transitions, rate coverage, compensation,
  observations and heterogeneity remain to be specified before engine changes.
- **Specification:** [parity fertility design/preparation](parity-fertility.md).
  Offline synthetic fixtures verify mapping and metadata; no live download or
  demographic validation is implied by those tests.


## D095 — Open fertility ages, HFD endpoints and truncated TFR (October 29, 2026)

- **Status:** accepted by the user. Upper-open bands remain allowed; their
  interior-limit compensation was superseded by D101 cell accounting.
- **Decision:** allow upper-open fertility age bands. Their birth-spacing
  compensation is the interior limit h=f/(1-b*f), with zero rate staying zero,
  zero spacing unchanged and nonpositive eligibility rejected. Finite bands
  retain D028's original midpoint calculation. The first b years in an open band
  approximate the effect of preceding rates; no artificial midpoint is assigned.
- **TFR:** omit upper-open age bands as before, warning when a known omitted
  female rate is nonzero that the resulting TFR underestimates full-age fertility.
  Zero omitted rates do not warn. Missing omitted rates remain unknown, not zero.
- **HFD:** the user additionally selected preserving endpoint rows: map 12- to
  [12,13) using its supplied rate, and 55+ to [55,Inf). This applies to total
  asfrRR and conditional mi preparation. Record the changed endpoint convention
  in source metadata. The lower-open value is represented at its labeled age,
  not recomputed into a new empirical estimate. Missing conditional cells still
  remain NA under D094. No parity engine execution is added.
- **Supersession:** replaces D028's blanket finite-age restriction and D094's
  earlier inherited omission of HFD open endpoints. Other timing/spacing,
  missing-data and parity-preparation decisions remain in force.


## D096 — Zero-fill unavailable HFD conditional rates (October 29, 2026)

- **Status:** accepted by the user after inspecting the conditional-rate import.
- **Decision:** replace unavailable conditional HFD age/parity rates with zero
  on import. Preserve explicit zero/nonmissing values and all selected cells;
  no interpolation. Record the zero-fill policy and replacement count in source
  metadata. Do not warn merely because this selected policy replaced missing
  cells. Total HFD rates and arbitrary input validation retain their existing
  missing-data rules. No parity-dependent engine or PATFR is implemented here.
- **Supersedes:** D094's preserve-NA-with-warning rule for HFD conditional cells.
  Zero fill is the user's modeling choice, not an empirical reclassification
  of missing HFD observations as measured zero hazards.


## D097 — PATFR and parity-aware TFR dispatch (October 29, 2026)

- **Status:** the user authorized a separate calculate_patfr(), automatic
  calculate_tfr() delegation for parity schedules, and the same stratification
  as calculate_tfr(). They explicitly selected exclusion of upper-open age bands
  with nonzero-rate warning, and NA results for missing included rates.
- **Measure:** deterministic synthetic-cohort parity progression using the HFD
  Methods Protocol (2015), sections 8.3–8.4, midpoint fertility-table convention.
  Output is patfr, including when called through calculate_tfr(), to distinguish
  it from ordinary TFR. Start at parity zero at the reported finite lower age;
  count repeated births within the terminal open parity category. Each supplied
  period/group/union-type stratum is calculated independently, selecting female
  rows if sex is present. Hold other strata fixed; no mortality or explicit
  birth-spacing state is included.
- **Numerics:** use annual age boundaries plus input boundaries, with midpoint
  entry/exit exposure and q=d*m/(1+d*m/2); fractional steps generalize the annual
  convention. Reject invalid progression probabilities rather than clipping.
  This is not the exact continuous-time hazard solution or a new engine time
  grid. All age/parity cells require unique contiguous coverage; broader finite
  parity bands apply their rate to each contained state. See the worked formula
  in [parity fertility](parity-fertility.md).
- **Limits:** distinguish raw missing rates (NA result) from earlier explicit
  HFD zero filling (D096). Upper-age omission, source endpoints and missing-cell
  assumptions can prevent exact equality to HFD's full published PATFR. No
  parity-dependent simulation or heterogeneity is authorized by this helper.


## D098 — Initialize parity from recorded maternal birth histories (October 29, 2026)

- **Status:** accepted by the user during parity-execution design.
- **Decision:** on detecting parity-specific fertility, initialize the internal
  parity vector from recorded births by mother_id. Women with no recorded births
  start at zero, including random baselines with no birth histories. An optional
  parity column must not override observed history. Count children regardless
  of survival or known father; retain actual counts beyond an open band's lower
  bound. Successful births increment parity, cancelled attempts do not.
- **Architecture:** per-simulation integer storage indexed by stable engine
  person position, avoiding repeated joins. Restart reconstructs counts from
  retained recorded histories; no additional authoritative public state is added.
- **Limit:** incomplete histories understate parity under this explicit rule.
  Birth-spacing compensation and remaining execution contracts are unresolved;
  this decision does not itself implement the engine or fertility heterogeneity.


## D099 — Execute parity fertility with SOCSIM-style compensation (October 29, 2026)

- **Current status:** the preceding-parity correction was superseded by D101 on
  October 29, 2026. Its parity partition validation, history-derived parity
  counts and schedule switching at births remain current.
- **Status:** accepted by the user after reviewing the pinned SOCSIM source trace.
- **Decision:** implement age/current-parity fertility with D098 history-derived
  integer counts. Parity bands partition [0,Inf) within each group/type; age gaps
  retain existing zero-risk meaning. Increment only on successful births and
  select the new schedule, retaining actual-birth spacing and refreshing any
  unpartnered fertility contribution used for father matching. No new public
  authoritative parity column or parity axis in matching bins is added.
- **Correction:** no adjustment at actual parity zero. For higher parity, divide
  its current input band rate by one minus the preceding actual parity's
  unadjusted rate integral over the birth interval, evaluated at the current
  age midpoint. Upper-open ages use the preceding schedule's tail rate times
  spacing. Zero spacing and zero input hazards remain unchanged. Positive
  rates require positive eligibility and representable adjustment. Group factors
  apply before correction; prior group/union-type history is not reconstructed.
- **Highest parity:** retain actual counts; entry to the final open category
  corrects from the prior category, while subsequent births correct from its
  own input schedule. Broad finite parity bands also use actual k and k-1.
- **Intentional departure:** use annual inputs, continuous execution and existing
  midpoint/interior-limit approximations rather than SOCSIM's monthly
  probability-sum grid. This mathematical implementation copies no SOCSIM source.
- **Limits:** the correction is not the exact current-parity eligible fraction
  and does not guarantee input conditional-rate recovery. Synthetic checks must
  report observed differences. No fertility heterogeneity, automatic prehistory
  synthesis or parity-stratified public observation API is added in this change.
- **Specification:** [parity fertility](parity-fertility.md);
  [static SOCSIM trace](socsim-parity-spacing-audit.md). Milestone 5 scientific
  acceptance remains with the user.


## D100 — Separate mortality and fertility observation APIs (October 29, 2026)

Accepted by the user: replace the public combined `calculate_demog_rates()`
with `calculate_mortality_rates()` and `calculate_fertility_rates()`. Both retain
period/age splitting and optional group/union-type strata. Fertility additionally
accepts optional `parity_breaks` defining half-open observed current-parity bands.
The shared combined implementation remains internal; there is no exported
compatibility alias. Fertility returns women only; mortality returns both sexes.
Observed parity uses maternal birth history, with pre-birth event assignment
and post-birth exposure progression; it does not change simulation calibration.
See [demographic analysis](demographic-analysis.md#fertility-rates-by-observed-parity).

## D101 — Cell-accounted fertility under birth spacing (October 29, 2026)

- **Current (November 5, 2026):** cell accounting and spacing are unchanged;
  D128 moved births to one fertility clock, pools groups when schedules are
  not group-specific, and reports unrealized births as counts.
- **Status:** accepted by the user as a specification; implemented October 29,
  2026, first as the opt-in `spacing_method = "cell"`. After reviewing the
  validation, the user made it the only method on October 29, 2026: midpoint
  compensation (D028, D095's interior limit, D099's correction) and the
  `spacing_method` argument were removed. Midpoint comparisons are reproducible
  at commit `89d451b`. The user revised decision (5) after validation showed
  that input-band cells recover only band averages, and later agreed to keep
  low-rate upper-open bands as one cell (`birth_interval * rate < 0.001`).
- **Decision:** generate births per fertility cell (age sub-band of at most one
  year x actual parity x union type x group, after group multipliers) at intensity
  `m_c * N_c`, where `N_c` counts all living members including those inside
  their birth interval. Assign each birth to an eligible member with
  probability proportional to an individual weight (all 1 until heterogeneity
  is accepted). Equivalent individual hazard: `m_c * N_c * w_i / W_c^e`.
- **Agreed choices:** (1) accept population-dependent individual hazards as a
  calibration device; (2) drop intensity in cells with members but no eligible
  weight and report unrealized expected births; (3) add as an option alongside
  D028/D099 first, then make it the only rule once validated; (4) weights
  determine selection probability within a cell, not a direct rate multiplier
  (option A), so cell rates stay pinned to inputs; (5) revised after
  validation: split age bands into sub-cells of at most one year at whole ages
  and band edges, and split pooled parity bands by actual parity, so the rate
  is recovered at every age and contained parity (input-band cells recovered
  only band averages); (6) report
  unrealized births per period and cell in run metadata as a tibble.
- **Rationale:** input rates count all woman-years, including spacing-excluded
  time. Precomputed corrections assume a mix of eligible and ineligible women
  and fail when the simulated mix differs: young high-parity cells (milestone
  5; also reproduced in SOCSIM), union-type transitions (milestone 2, D045),
  restarts and start-up (milestone 1). Accounting of the simulated population
  recovers inputs in expectation for any mix while eligible women exist.
- **Alternatives:** a deterministic composition-based parity correction from
  the target synthetic cohort (exact only when the simulated composition
  matches it); weight-scaled rate multipliers (option B), which double-count
  parity selection and move observed rates away from inputs.
- **Implementation:** a `BoundaryQueue` class extracted from the calibration
  tracker (no behavior change) serves both the tracker and the cells; each keeps
  its own instance, and under cell mode the tracker no longer tracks fertility
  boundaries. Calibration fertility totals come from cells at search time. Cell
  clocks are kept in a separate indexed heap compared with the event queue, and
  cell boundaries are processed in time order before events. Uniform selection;
  D102 adds weighted selection for persistent lognormal heterogeneity. Unrealized births are reported
  in `simulation$unrealized_births` by one-year period from the run start.
- **Evidence:** toy diagnostics in `validation/milestone5/`
  (`check-composition-spacing.R`, `check-runtime-spacing.R`) and the SOCSIM
  reference experiment. Engine validation follows the specification's plan.
- **Specification:** [cell-accounted fertility](fertility-cell-accounting.md);
  [investigation](parity-spacing-investigation.md). D102 subsequently selects lognormal weights and maternal inheritance.


## D102 — Persistent lognormal fertility weights and maternal inheritance (October 29, 2026)

The user selected the API and lognormal inheritance construction:
`set_fertility(..., hetfert = het_fert(cv, alpha))` enables
mean-one lognormal weights with requested CV and maternal inheritance on the
latent normal scale. `hetfert = NULL` retains uniform sampling. The population's
existing female `persons$hetfert` values are retained; if the column is absent,
simulation draws values for women. Daughters use maternal alpha plus independent
innovation scaled by sqrt(1-alpha^2); no generational normalization is applied.
Weights determine eligible-mother selection within fixed-intensity cells, not
direct rate multipliers. Output always includes the column when enabled.

The [heterogeneity specification](fertility-heterogeneity-notes.md) records
units, strict-positive supplied female weights, seed handling, zero CV, endpoint
inheritance, restarts and the reproductive-selection caveat. Distribution
preservation is conditional on the reproducing mothers' distribution, not a
universal guarantee. The source audit of SOCSIM does not select its cubic
distribution or betaT transformation for Hammel.


## D103 — Milestone order and groups-based characteristics scope (October 29, 2026)

The user confirmed milestone 4 is complete and milestone 5 is current. Move
group transitions to milestone 6 and time-varying inputs to milestone 7. Use
"group transitions" rather than "group conversions"; this does not rename
cohabitation-to-marriage conversion. Remove the separate individual-
characteristics milestone: most categorical applications should use groups,
with multiple types of groups as an option for more complex applications.

This selects planning order and design direction only. Transition semantics,
multiple-group representation and public APIs remain to be specified; it does
not authorize arbitrary state mutation or mark milestone 5 complete. See the
[current roadmap](roadmap.md). Milestones 8–10 retain their existing numbering.


## D104 — Disabled fertility heterogeneity removes output weights (October 29, 2026)

The user selected that `hetfert = NULL` produces no `persons$hetfert` column,
even if the starting population includes one. This applies to zero-year and
empty runs as well. The input population is not mutated. A subsequent enabled
run therefore draws fresh weights for all recorded women when using that
output; it does not silently assign unit weights to daughters born while the
feature was disabled. Weight persistence is guaranteed across enabled restarts,
not across an explicit disabled interval. This replaces D102's original
ordinary-custom-column treatment for disabled runs.


## D105 — Milestone 5 accepted complete (October 29, 2026)

- **Status:** accepted complete by the user after reviewing and revising the
  milestone 5 report. Its final HTML render was running when the user requested
  documentation closeout; preserve the source, rendered artifact and active
  intermediates. Closeout does not independently rerun or certify that render.
- **Scope:** parity-specific conditional fertility inputs and history-derived
  parity (D094–D099), HFD preparation and PATFR/TFR dispatch, separate mortality
  and fertility observations with parity bands (D100), cell-accounted fertility
  replacing precomputed spacing compensation (D101), persistent lognormal
  fertility weights with maternal inheritance (D102), and the disabled-weight
  lifecycle (D104). The report revisits baseline, restart and union-type
  fertility, and assesses parity rates, heterogeneity and performance.
- **Evidence:** the user-authored [report](../validation/milestone5.qmd), its
  [rendered assessment](../validation/milestone5.html), and supporting
  [diagnostics/provenance](../validation/milestone5/README.md). Acceptance is the
  user's scientific assessment, not a new agent-defined numerical threshold.
- **Retained limitations:** cells with members but no eligible women lose
  intensity, reported as unrealized expected births. Adding it back is an
  accounting diagnostic, not realized births. Sparse or structurally infeasible
  young/high-parity cells need not recover their target. Missing birth
  histories cannot be reconstructed. Inherited weight distributions may change
  under reproductive selection; alpha and raw mother–daughter correlations
  are not generally identical. Restarts are not exact checkpoints, and failed
  partner searches retain their existing cancellation behavior.
- **Next:** milestone 6 is group transitions, using that term for life-course
  group changes. Multiple types of groups remain an optional complexity
  direction; transition semantics, histories, tracker/event invalidation,
  ancestry effects and public API require specification before implementation.
  Milestone 7 remains time-varying inputs under D103. Completion does not
  authorize unspecified future demographic behavior.

## D106 — Group-transition inputs use annual hazards (October 29, 2026)

**Status:** accepted by the user; transition execution remains unimplemented.

Group-transition inputs continue Hammel's rate convention: hazards per
person-year, rather than SOCSIM's monthly probabilities. Users with an annual
probability p can approximate a rate with `-log(1 - p)` (numerically,
`-log1p(-p)` in R). This conversion is exact for a constant hazard over one
year without competing exits; a cause-specific probability observed with
competing events does not generally identify that hazard by this formula.

The user explicitly selected rates following the SOCSIM transition audit.
This settles input units and interpretation only. Eligibility, clocks,
histories, linked moves, ancestry effects, event invalidation and the public
API still require specification before implementation.

## D107 — Directed transition schedules and optional strata (October 29, 2026)

**Status:** accepted by the user; transition execution remains unimplemented.

For n groups, up to n(n-1) directed origin-to-destination schedules are
possible, excluding self-transitions. For a person in an origin group, the
applicable outgoing schedules compete with each other and other applicable
events. For example, A-to-B and A-to-C compete for a person currently in A;
rates are annual hazards under D106 and need not sum to one.

Each schedule can be age-specific or duration-specific. Schedules may
optionally stratify by sex and marital status. An omitted sex dimension
applies to all sexes; an omitted marital-status dimension applies to all
marital statuses. Either dimension can be supplied independently.

The duration origin and initialization, marital-status categories, treatment
of missing directed schedules, and validation of overlapping universal and
stratified inputs remain to be specified. This decision does not select a
public API, joint age-by-duration surfaces, or household moves.

## D108 — No transition-count target adjustment (October 29, 2026)

**Status:** accepted by the user; transition execution remains unimplemented.

Use supplied transition hazards directly, without scaling them to attain
requested origin-to-destination event counts. Do not add SOCSIM-style
`transittarget` adjustment. The user rejected this mechanism after reviewing
its target/expected-count scaling in the reference audit. Realized counts
follow the supplied hazards, changing exposure, competing events and random
variation; this decision adds no count-calibration API.

## D109 — Household transitions are a stretch goal (October 29, 2026)

**Status:** original scope direction superseded by D114 on October 30:
household moves are deferred to milestone 10, internal migration and regions.

Individual group transitions form milestone 6's core scope. Linked household
moves are a stretch goal, potentially useful when groups represent regions.
Core completion does not depend on delivering household moves.

Household membership, follower eligibility, destination assignment and
interactions with partners' pending events require a separate specification
before implementation. This decision does not adopt SOCSIM's household rules
or select a residence/custody schema.

## D110 — Transition rate-table schema and common band type (October 29, 2026)

**Status:** accepted by the user; transition execution remains unimplemented.

Transition rate tables require `group_origin`, `group_dest`, and `rate`, plus
either age bands or duration bands. All origin-to-destination pairs in the
table must use the same band type: age-based or duration-based. Do not mix
age-based and duration-based directed schedules in one transition rate table.
Rates are hazards per person-year (D106).

The optional stratification columns are `sex` and `union_type`. If a column
is present, schedules apply by its categories; if absent, schedules apply
universally across that dimension. Either column may be supplied independently.
`union_type` replaces D107's preliminary marital-status wording; it does not
introduce a separate marital-status stratification column. The common-band
rule narrows D107's earlier per-schedule age/duration choice.

Duration origin and initialization, category coverage and labels (including
unpartnered people), absent directed-pair treatment, and detailed band
validation remain to be specified before implementation. No within-table
universal/category override scheme is selected.

## D111 — Grouped output includes a transitions tibble (October 29, 2026)

**Status:** accepted requirement; transition execution and output remain
unimplemented.

When groups are in use, simulation output must include a `transitions` tibble
recording individual group changes, including an empty tibble when no changes
occurred. Presence does not depend on configuring a transition process or
observing a transition. Each change must retain the person, transition time,
origin and destination so group membership can be reconstructed over time.
An empty table must retain the defined column types.

Exact history column names, initialization/restart validation and observation
boundary semantics remain to be specified. This requirement does not select
whether ungrouped output includes the component.

## D112 — Individual group-transition execution (October 30, 2026)

**Status:** implementation authorized by the user and implemented; milestone 6
subsequently accepted complete under D116.

The user requested `group_transition(rates)` and pipeable
`set_group_transition(model, rates)`, with person-level competing waiting times
and updates to all state affected by a move. D106–D111 govern annual hazards,
directed schedules, one common age/duration band type, optional sex and
union-type strata, no count targeting, individual-only core scope and retained
transition histories. Household switching is explicitly excluded from this
implementation.

The user confirmed duration since latest group entry, falling back to age
without recorded transitions; omitted pairs and gaps mean zero risk. The user
also confirmed ancestry proportions stay unchanged. Initial ancestry may be
inferred from starting identities when shares are absent; newborn shares
continue to follow parental ancestry independently of identity.

The [implemented contract](group-transitions.md) specifies the history schema,
half-open bands, state/event invalidation, restarts and observations. Columns
are `person_id`, `time`, `group_origin`, and `group_dest`; grouped populations
retain a typed empty table when no moves occur. Software validation does not
establish scientific acceptance or authorize household moves.

## D113 — Directed transition Lexis counter (October 30, 2026)

**Status:** user requested implementation; implemented under the
[transition observation contract](group-transitions.md#observed-transition-rates).

`calculate_transition_rates()` reports person-years, transition events and rates
by directed origin/destination, calendar period and either age or group
duration. The user requested optional sex and union-type disaggregation;
`by_sex` and `by_union_type` follow the existing sex field and union-type
conventions. Defaults pool both dimensions. Duration follows D112's latest
entry/initial-age rule. Rates divide each directed event count by exposure
in its origin, so outgoing destinations share denominators. This is an
observation helper and changes neither event execution nor model assumptions.

## D114 — Defer household moves to the regions milestone (October 30, 2026)

**Status:** accepted by the user; supersedes D109's milestone 6 stretch goal.

Household moves were deferred to the then-planned milestone 10, internal
migration and regions. D121 subsequently removes that standalone milestone and
retains these unresolved applications under tentative milestone 10, multiple
group columns; no household behavior is selected by that reordering.
The user identified ambiguity between individual initiation hazards and total
move hazards induced by household members, and questioned conditional
rescheduling after externally triggered moves. Discussion also identified
membership/custody, partial moves and mixed origins, changing households,
atomic updates, shared dependence, and the distinction between household
initiations and individual observed moves.

The [regions milestone notes](roadmap.md#household-moves-deferred-design-questions-d114)
record these issues and suggested validation. Correct conditional resampling
is possible while retaining attained ages, durations and persistent state;
it does not resolve the input-rate or household-definition questions. No
household model, rate carrier, schema or implementation is selected.

## D115 — Detailed union groups refer to formation (October 30, 2026)

**Status:** correction requested by the user and implemented.

The user identified that `get_detailed_unions()` joined final person groups,
which relabeled earlier pairings after transitions and distorted formation
cross-tabs. `group_w` and `group_h` now reconstruct partner identities at the
union's start from retained histories. A same-time transition follows
formation under the existing event priority, so its origin applies. Without
transition history, use the recorded group. Dates, row order, custom columns
and label types retain their existing behavior. The native indexed-history
lookup avoids repeated per-union scans in R. Regression tests verify return
moves, repeated episodes, exact-time moves, unchanged earlier pairings after
later transitions, empty outputs and character/numeric labels. Correcting
attribution does not independently establish matching-parameter recovery.

## D116 — Milestone 6 accepted complete (October 30, 2026)

**Status:** accepted complete by the user after report review, observation
corrections and external code-review follow-up. The user explicitly requested
completion documentation while compiling the final report.

- **Scope:** individual directed group transitions using annual hazards,
  common age or duration bands, optional sex/union-type strata, latest-entry
  duration (age without earlier history), zero risk in omitted pairs/gaps,
  unchanged ancestry, consistent event/pool/cell/union updates, retained
  transition histories and restart reconstruction (D106–D112). Directed Lexis
  rates and formation-time partner group extraction are included (D113/D115).
- **Assessment:** the user's [report](../validation/milestone6.qmd) and
  [rendered assessment](../validation/milestone6.html) examine transition-rate
  recovery, demographic summaries, matching cross-tabs at formation,
  population shares and performance. Scientific acceptance belongs to the
  user; no numerical threshold is added by this closeout. The final render
  was active at acceptance. Preserve all report artifacts and intermediates;
  this entry does not certify render completion or invent final results.
- **Retained limitations:** finite-sample/sparse-cell effects, complete
  histories for historical observations, initial duration falling back to age
  without earlier transitions, non-checkpoint restarts, and the documented
  life-table age-aggregation caveat. An attached zero-rate process does not
  promise fixed-seed equality with omitting it. Existing fertility, matching
  and history limitations remain documented in their contracts.
- **Deferred:** household moves belong to milestone 10, internal migration
  and regions (D114). Multiple group dimensions and fixed newborn assignment
  remain unselected future work; custom mixed-parent assignment remains
  available. Completion does not authorize those mechanisms.
- **Next:** milestone 7, time-varying inputs (D103). Its scientific contract,
  schedule boundaries, event invalidation and public API require specification
  before consequential implementation. Earlier proposals are not accepted
  design decisions or automatic implementation authorization.

## D117 — Tentative milestone 11: two-sex marriage model (October 30, 2026)

**Status:** user-selected tentative planning milestone; renumbered to milestone
12 by D121 on October 31, 2026. The model and API remain unselected.

Add a two-sex marriage model as tentative milestone 11, following the regions
milestone. Existing [two-sex research notes](two-sex-union-models.md) and D073
provide the starting discussion. No particular model, parameterization, input
semantics, cohabitation integration or public API is selected by assigning this
milestone. Current union formation remains unchanged; milestone 7 remains next.


## D118 — Calendar-varying rates within one run (October 30, 2026)

**Status:** user-approved semantics and implementation authorization; implemented.
Milestone 7 was subsequently accepted complete under D127.

The user approved flat rate tables with `time_start`/`time_end`, half-open
calendar intervals, piecewise constant rates without interpolation, coverage
checks, boundary changes before same-time events, and updates to pending risks,
fertility cells and matching calibration while preserving demographic histories.
The user narrowed this implementation step to the six existing annual-rate
inputs: mortality, fertility, union search, dissolution, conversion and individual
group transitions. Matching and inheritance model schedules are deferred to a
second step. Other arguments, including `cohab_prob` and rate multipliers,
remain static. Fractional calendar years and simulation durations are supported.

For the deferred model-schedule step, the user also approved selecting matching
at partner selection and inheritance at birth, without rewriting existing unions
or identities. Complete-model schedule tibbles are the selected representation;
constructor names and public arguments still need design. This step does not
implement those schedules.

See [the rate contract](time-varying-rates.md) for validation, worked examples,
clock origins, observation semantics and scope. Conditional exponential
resampling at calendar boundaries is the engine implementation choice; histories
and duration origins do not restart. Static inputs remain supported. Neither
implementation authorization nor software checks establish demographic acceptance.


## D119 — Multi-year database preparation and native HMD shapes (October 30, 2026)

**Status:** user-requested implementation, extending D118's rate preparation.

The user requested both getters accept one year (static output without calendar
columns) or a vector of years (output with `time_start`/`time_end`), and requested
additional source shapes while retaining `1x1` as default. All requested years
are selected locally from one download of each required product.

HMD supports its six native age-by-time shapes through `shape`; grouped `year`
requests select source period starts and preserve complete source intervals.
The getter retains source pooled rates and shortened endpoint periods, without
estimating or interpolating rates. HMDHFDplus 2.0.8's loss of end labels requires
an auxiliary annual coverage download for parsed grouped-time products; this is
an implementation workaround recorded in provenance. HFD's existing products
remain native annual/single-age, with no equivalent shape option or invented
local aggregation. The existing HFD endpoint and conditional zero-fill policies
are retained.

See [database preparation](database-rates.md) for interval mapping, examples,
validation and scope. This does not establish milestone 7 scientific acceptance.


## D120 — Optional pooled HMD open age (October 30, 2026)

**Status:** user-selected pooling rule and implementation authorization.

After Swedish historical Mx downloads exposed missing oldest-age rates, the
user requested optional age truncation and explicitly selected pooling deaths
and exposure instead of extending the cutoff-age rate. The getter accepts
`open_age = NULL` by default. A supplied cutoff retains source rates below it
and replaces its tail with `[open_age, Inf)` at total deaths / total person-years,
separately by sex and selected source period. Download matching native deaths
and exposure products once each. Missing Mx in the replaced tail are allowed;
missing/invalid retained rates or tail counts/exposure remain errors. Pooled
exposure must be positive; no imputation or forced maximum lifespan is added.

The cutoff must align with native source bands; no splitting counts or silently
moving the cutoff. Validate matching complete age tails and source periods.
Provenance retains cutoff, method, auxiliary products and a tibble of pooled
counts, exposure and rates. Singleton versus vector-year output schemas are
unchanged. See [database preparation](database-rates.md) for worked examples
and validation. Stationary-year investigation and historical simulation setup
were explicitly deferred while fixing the getter.

## D121 — Binned-group hazard direction and roadmap revision (October 31, 2026)

**Status:** user-selected engine direction and milestone priorities. *Current
(November 5, 2026):* implemented under D128/D129 and accepted under D130. The
single-year subdivision below was revised during implementation: only
fertility splits to single years (D128). The milestone numbering below was
reordered by D131.

The user will finish milestone 7 with the current engine, then make milestone 8
use binned-group hazards for all six existing rate processes and allow custom
formula weights on all rates. Input rates determine aggregate cell intensity;
weights redistribute event risk within cells. Use single-year age or duration
subdivisions even for wider input bands, retaining fractional boundaries and
all other input strata. Fertility retains its actual-parity cells. Calendar
changes update active rates rather than add historical-period membership bins.

The user distinguished separate rate targets from relative formula stratification
and recognized that weighting an already stratified group has no effect when
the resulting factor is constant within a cell. Redundancy warnings should be
considered; their exact diagnostic policy is still a proposal. Formula weights
provide a route to future user-defined characteristics, but do not specify their
initialization, inheritance, longitudinal updates or missing-data rules.
See [the discussion record](binned-group-hazards.md) for equations, examples,
composition effects, process-specific exposure/eligibility distinctions,
unresolved contracts and proposed validation. No new scheduler is implemented.

The user revised the roadmap as follows:

| Milestone | Revised scope | Status |
| --- | --- | --- |
| 7 | Finish existing time-varying-input work with the current engine | In progress; scientific acceptance remains pending |
| 8 | Binned-group hazards for all rate processes and custom formula weights | User-selected direction; specification/API and implementation pending |
| 9 | Former milestone 8: summaries, diagnostics and helpers | Proposed details retained |
| 10 | Possibility of more than one group column | Tentative; schema and mechanisms unselected |
| 11 | Former milestone 9: external immigration and generation tracking | Tentative; no emigration mechanism implied |
| 12 | Former milestone 11: two-sex marriage model | Tentative; model and API unselected |

The former standalone milestone 10 for internal migration and regions is removed
because those applications may be handled through multiple group columns in new
milestone 10. This supersedes D114's dedicated milestone placement, not its
unresolved household-move questions, which remain in the roadmap. A separate
mandatory `region` schema is no longer the roadmap direction. Household
membership, regional matching restrictions and linked-move semantics remain
unselected. D117's two-sex priority moves from milestone 11 to 12; its scientific
model remains unselected. No existing validation report, result or acceptance
status is changed by this reordering.


## D122 — Calendar-varying cohabitation probabilities (October 31, 2026)

**Status:** user-approved semantics and explicit implementation authorization.

The user selected `cohab_prob` to work like flat rate tables, rather than whole-
model schedules. Optional `time_start`/`time_end` supply half-open absolute-year
intervals, including fractional years. Retain age/group/previous-union strata,
probabilities in [0,1], static-table compatibility and age-gap probability zero.
Validate every active slice, reject overlapping risks and require calendar
coverage during execution. The active probability at search execution selects
intended type before matching, with changes preceding exact same-time searches.

Refresh marriage/cohabitation calibration contributions as lambda*(1-p) and
lambda*p. Probability-only boundaries preserve pending person/union/fertility
clocks and consume no randomness; a probability change does not change the
overall search hazard. Shared rate boundaries apply both changes before D118's
existing rescheduling. Existing unions keep their type; histories, ancestry,
spacing, parity, duration origins and prior-union classification remain unchanged.
The configured table retains typed mode even during probability-zero intervals.

This extends D118's original rate-only scope. Matching and inheritance model
schedules remain future milestone 7 work, with APIs still unselected; no
milestone 8 engine change is introduced. See [cohabitation](cohabitation.md#calendar-varying-type-probabilities-d122)
for worked examples and validation. Implementation checks are distinct from
milestone 7 scientific acceptance.


## D123 — Complete model schedules for all matching and inheritance arguments (October 31, 2026)

**Status:** user-selected API and explicit implementation authorization; implemented.
Milestone 7 was subsequently accepted complete under D127.

The user selected either a direct static model or a tibble with `time_start`,
`time_end` and list-column `model` containing complete objects for `matching`,
`father_matching`, `cohab_matching` and fertility's `inheritance`. Base data
frames with list columns are accepted and normalized to tibbles. Preserve full
objects, formula environments, callbacks, extra columns and provenance.

D118's event-time selection applies: active matching at partner/father choice,
active inheritance at each successful birth, with half-open absolute-year bands,
no interpolation, overlap rejection and coverage across a run. Validate every
row, including unused/future intervals. NULL entries and nested schedules are
not complete models; existing top-level optional matching NULLs dynamically
follow general matching. Static inheritance method aliases remain supported.

Use one registry covering identities from all inheritance models, including
future custom outputs; validate grouped schedules and matching levels against
that registry before execution. Keep past identities, unions, ancestry,
spacing, durations, weights and pending event dates at model-only boundaries.
Refresh active role-specific preference/draw structures and weighted-calibration
sums from the current pool, without demographic RNG draws. Custom inheritance
executes only at births; deterministic matching formulas may be evaluated for
cache preparation under the existing purity contract.

See [the model-schedule contract](time-varying-models.md) for worked examples,
validation, retained defaults and limitations. Direct static execution remains
supported; caches belong to the simulation call. Model schedules complete this
implementation component of milestone 7, not its scientific acceptance. Other
non-rate arguments remain static, and the milestone 8 hazard rewrite remains
unimplemented.

## D124 — Optional exact Poisson intervals for observed event rates (November 1, 2026)

**Status:** User-requested interface and Poisson method; implemented. This does
not establish milestone 7 scientific acceptance.

**Decision:** Add `conf.level = NULL` to mortality, fertility, union, marriage,
dissolution, cohabitation conversion and transition rate calculators. A supplied
confidence level requests `conf_low` and `conf_high` using two-sided
`stats::poisson.test()` intervals for each row's observed event count and its
existing person-year or union-year denominator. NULL preserves existing output.

**Implementation choices:** Require a finite numeric scalar strictly between
zero and one. Zero exposure yields NA bounds; zero events with positive exposure
retain the Poisson upper limit. Empty tables retain numeric confidence columns.
Fertility bounds refer to observed births, including with `unrealized = TRUE`;
fractional expected unrealized intensity is not treated as an observed count.

**Rationale and limits:** Supply within-run counting uncertainty without requiring
replications. Normal SE approximations are unsuitable for very small counts.
Exactness refers to the Poisson model, not guaranteed coverage for all simulator
feedback or shared exposure. These are marginal intervals, not Monte Carlo SEs
or simultaneous confidence bands. See [the analysis contract](demographic-analysis.md#optional-poisson-rate-intervals).

## D125 — Widowed search-rate fallback to dissolved (November 2, 2026)

**Status:** Superseded by D126, which replaces the ambiguous detailed
`dissolved` label with `separated` and reserves `dissolved` for pooled inputs.

**Decision:** Union-search rate tables stratified by `previous_union` may supply
`none` and `dissolved` without `widowed`. For each group and active calendar
slice, omitted widowed rates use the complete dissolved schedule, including
its age/duration clock. Explicit widowed schedules take precedence as a whole;
their gaps remain zero. Effective clocks retain the existing consistency rules.
Actual widowhood histories and observation categories do not change.
This supersedes D080's mandatory third category for search-rate inputs only;
cohabitation-probability tables retain their existing three-category contract.
Preserve the public input and expand only a private execution copy. See
[remarriage schedules](remarriage.md) for examples and validation.

## D126 — Distinct detailed and pooled previous-union categories (November 2, 2026)

**Status:** User-requested, implemented and software-tested; replaces D125.

**Decision:** Detailed previous-union categories are `none`, `separated`, and
`widowed`; pooled categories are `none` and `dissolved`, where dissolved
includes separation and widowhood. Use the same two complete alternatives for
union-search rate and cohabitation-probability inputs, independently per group
and active calendar slice. Never mix the two alternatives within one slice and
group. Expand pooled inputs privately to the two detailed engine states; retain
actual histories, duration origins, eligibility and event timing. Effective
categories must retain consistent clocks across groups and calendar periods.
No compatibility aliases, migration messages or automatic reinterpretation of
old three-category tables are required for this unreleased package.

Union and marriage calculators retain `by_previous_union = TRUE` for detailed
output and add `pool_dissolved = TRUE` for the bipartite output. The user
selected this additional flag rather than overloading `by_previous_union`;
pooling requires previous-union stratification (explicitly enabled or implied
by duration breaks). Sum separated/widowed events and exposure in
each otherwise identical cell before calculating rates and exact Poisson
confidence intervals. With duration breaks, the pooled dissolved duration is
still elapsed time since the most recent ending, irrespective of cause. Empty
cells, grouping, union types, sex and calendar/age boundaries retain their rules.
For example, 2 separated events/10 years plus 3 widowed events/20 years gives
5 dissolved events/30 years and a confidence interval from that count/exposure,
not an average of rates or confidence limits. Validate equivalence of pooled
inputs to explicitly equal detailed inputs and accounting/Poisson identities
for pooled outputs, including calendars, groups, durations and conversions.

## D127 — Milestone 7 accepted complete (November 2, 2026)

**Status:** Accepted complete by the user's explicit request to record milestone
7 completion. Date uses America/Los_Angeles; recorded November 3 UTC.

**Decision:** Calendar-varying inputs are complete with the current engine:
six annual-rate processes under D118, calendar-banded cohabitation probabilities
under D122, and complete matching/inheritance model-list schedules under D123.
Retain supporting database preparation (D119/D120), optional Poisson rate
intervals (D124), and consistent detailed/pooled previous-union categories and
confidence intervals (D126). This does not select or implement milestone 8's
binned-group hazard and custom formula-weight mechanisms (D121).

**Evidence:** The user's [milestone 7 report](../validation/milestone7.qmd),
[rendered assessment](../validation/milestone7.html), its conclusion that the
expectations were met, and the user's explicit completion instruction establish
acceptance. The code/report revision at acceptance is
`3b65d2202b006ef6f3e232c98595814b1544f9fe`. Software checks and the calendar-review
regressions are separately recorded in [testing](testing.md). This documentation
update does not rerun simulations, render the report, or independently certify
all numerical findings or a remote CI rerun. Preserve all seven accepted report
sources, rendered artifacts, source comments, provenance and active render
intermediates.

**Retained limits:** Acceptance applies to the assessed scenarios and documented
contracts. Finite samples and sparse cells, complete-history requirements,
failed searches and matching calibration, non-checkpoint restarts, and finite
resolution of unrealized-birth reporting remain relevant. The historical-data
demonstrations use documented approximations; completion does not establish
universal rate recovery, historical reconstruction, SOCSIM equivalence, or
release readiness. Non-rate settings outside D122/D123 remain static. Current
calendar coverage and strata rules are unchanged. No additional acceptance
thresholds or demographic defaults are inferred.

**Next priority:** Milestone 8 remains the selected future direction under D121;
its scientific specification, formula API and unresolved mechanisms require
separate work and user decisions.

## D128 — Cell scheduling for every rate process (November 4, 2026)

**Status:** User-directed during milestone 8 (step 1 and follow-up directions,
November 3–4, 2026); accepted with milestone 8 under D130.

**Decision:** Mortality, fertility, rate-mode search, dissolution, conversion
and each group-transition slot are scheduled by cells of at-risk units
(persons, or unions for dissolution and conversion). A cell produces events at
its input rate times its members. Each process keeps one Exp(1) clock spent at
its total intensity; a fired clock picks a cell in proportion to intensity and
a unit within it. Band-edge moves are bookkeeping without random draws. The
engine keeps no individual waiting-time clocks: one approach for every process,
chosen by the user to simplify the code base rather than keep a hybrid.
User directions within the decision:

- Only fertility splits bands to single years and actual parity, because
  spacing is an eligibility rule (D101); other processes use input bands.
- Cells stratify only as finely as the input rates: by group only when a
  process's effective schedules are group-specific.
- Dissolution units are unions, not their female partners.
- Unrealized births are counted when the clock draws a cell with no eligible
  woman, by one-year period, group, union type, parity band and age band;
  `years_without_eligible` and expected-value records are removed (revising
  D101's reporting; realized births are unchanged in distribution).
- Calendar rate changes update cell rates in place when the cell layout is
  unchanged; otherwise cells are rebuilt.
- Matching calibration takes search and fertility contributions from the cells,
  replacing D082's per-woman tracker; a probability tracker remains for
  cohabitation probabilities.

**Rationale and evidence:** Equivalent in distribution to independent
individual hazards; validated against the previous engine over replications,
the milestone 5 recovery checks and analytical mortality. Run times were
measured against the previous engine and SOCSIM 1.9.18: slower for simple
mortality/fertility models with many birthday moves, faster for models with
calendar rate changes. See [cell scheduling](binned-group-hazards.md#cell-scheduling-d128).

**Implications:** Fixed-seed histories differ from earlier versions and between
differently stratified but equivalent models. Same-time ties between processes
follow death, birth, dissolution, conversion, search, then transitions.

## D129 — Weight models (November 4, 2026)

**Status:** Specified with the user November 4, 2026 and implemented; accepted
with milestone 8 under D130.

**Decision:** `create_weight_model()` defines log-linear formula weights
(coefficients, treatment-coded categorical variables with fixed levels,
`group_levels` for group variables) or a custom function of a cell's members
and its description. Models attach through `weights` on mortality, fertility,
rate-mode search, dissolution and transitions, and `weights_conversion`, as one
model or a calendar schedule using the D123 tibble structure. Weights are
normalized within cells, so cell intensities are unchanged; exp(beta) is a
within-cell hazard ratio. Variables are fixed person and union characteristics
that change only at events; user-supplied columns are excluded. Missing
formula variables and NA function results contribute zero with a warning;
non-finite coefficients and infinite results are errors. Fertility weights
combine with `hetfert` on the log scale. Formula weights are stored and
selected by rejection sampling; search weights under matching calibration
enter the calibration sums. Terms that cancel within every cell warn.

**Alternatives deferred:** unnormalized weights (frailty selection),
eligibility rules beyond spacing, `NULL` schedule entries (only in parallel for
matching, inheritance and weight schedules), weights varying within bands and
user-defined characteristics. See
[weight models](binned-group-hazards.md#weight-models-d129).

## D130 — Milestone 8 accepted complete (November 5, 2026)

**Status:** Accepted complete by the user's explicit request to close the
milestone as accepted.

**Decision:** Milestone 8 is complete: cell scheduling for every rate process
(D128) and weight models (D129), carrying out D121's direction.

**Evidence:** The user's [milestone 8 report](../validation/milestone8.qmd) and
[rendered assessment](../validation/milestone8.html), and the user's explicit
completion instruction. The code/report revision at acceptance is the
`milestone8-binned-intensity` commit that adds this entry, following `8ae2209`
(review fixes in `9176466`, coverage tests in `8ae2209`); the closeout commit
changes only documentation. Software checks are
recorded in [testing](testing.md). This documentation update does not rerun
or rerender the report or independently certify its numerical findings.

**Retained limits:** Weights are normalized within cells, so population-level
effects depend on cell size, homogeneity and composition change and cannot
produce frailty selection. Performance depends on model structure; simple
models with many band-edge moves are slower than the previous engine. Fixed-seed
histories are not comparable with earlier versions. Fertility spacing remains
the only eligibility rule. No acceptance thresholds or demographic defaults are
inferred.

**Next priority:** Milestone 9, multiple group columns, as reordered by D131.

## D131 — Milestones 9–11 reordered (November 5, 2026)

**Status:** User-selected planning priority, at milestone 8 closeout. No scope,
schema or API is selected by the reordering.

**Decision:** Move multiple group columns from tentative milestone 10 to
milestone 9, the next milestone, and external immigration from tentative
milestone 11 to 10. Summaries, diagnostics and post-simulation helpers move
from milestone 9 to 11. The two-sex marriage model remains tentative milestone
12. This revises D121's numbering only; D114's household questions remain
attached to multiple group columns.

| Milestone | Scope | Status |
| --- | --- | --- |
| 9 | More than one group column, potentially including regions and internal migration | Next; schema and mechanisms unselected |
| 10 | External immigration and generation tracking | Tentative |
| 11 | Summaries, diagnostics and helpers | Proposed details retained |
| 12 | Two-sex marriage model | Tentative; model and API unselected |

## D132 — `stratify_rates()` replaces `group_prop_hazard` (November 5, 2026)

**Status:** User-selected API change at the start of milestone 9, specified in
discussion November 5, 2026. Supersedes the `group_prop_hazard` arguments of
the October 20 multiplier request, D079's pair multipliers and the fertility
per-status multiplier lists.

**Decision:** Remove `group_prop_hazard` from every rate process. A new
function `stratify_rates(rates, by)` adds proportional-hazard strata to a rate
table: `by` is a named list of named factor vectors, for example
`list(group = c(A = 1.2, B = 1))`. Each input row is copied per category with
its rate multiplied by the factor. User choices:

- A column already present in `rates` is an error.
- Several columns in `by` are crossed, multiplying their factors.
- The function cannot know the process, so each setter checks the table's
  actual stratum columns, however they were added, against those its process
  supports. An attribute recording added columns was rejected, because tables
  can be pre-stratified or entered directly. `by` may only name stratum
  columns.
- Interactions with existing strata, such as different group effects by union
  status, are built by stratifying subsets and binding them.

**Implications:** Stratified tables are ordinary group-specific tables. Every
needed category must be listed; an omitted group no longer defaults to factor
one. Factors of one now split cells by group (D128), so weight terms on that
variable cancel with a warning. Without population groups, the standard
alphabetical-first fallback applies instead of ignoring multipliers. Results
equal the previous multipliers in distribution, with equal rates. Conversion
rates, which had no multiplier argument, can now be stratified the same way.
The accepted milestone 3 report uses `group_prop_hazard` and reproduces at its
recorded commit, as the milestone 1 report already does.
See [groups](groups.md#proportional-stratification-with-stratify_rates-d132).

## D133 — One `by` argument for rate calculators (November 5, 2026)

**Status:** User-selected API change in milestone 9, specified in discussion
November 5, 2026.

**Decision:** Replace the `by_group`, `by_union_type`, `by_previous_union` and
`by_sex` arguments of every rate calculator and `calculate_population_size()`
with `by`, a character vector of categorical variables. Variables may go beyond
what a process accepts as rate strata when they are easily computed from the
histories. User choices:

- Band arguments (`age_breaks`, `dur_breaks`, `parity_breaks`) stay separate;
  `by` holds categories only. `parity_breaks` stays in fertility only.
- Sex stays always split where it is now (mortality, fertility, union and
  marriage rates) and is a `by` value for transitions and population size.
- Dissolution and conversion rates use the wife's and husband's groups,
  `group_w` and `group_h`, instead of the unordered `group_pair`, plus
  `union_type` for dissolution.
- `pool_dissolved` stays a separate option, now available wherever
  `previous_union` is in `by`; `dur_breaks` still adds `previous_union` in
  union rates.
- One row-order rule for all calculators.
- New variables now: `previous_union` for mortality, fertility, transitions and
  population size; `union_type` for population size.

**Implementation:** A shared native observation engine builds each person's
timeline of states (group, union type, previous union, parity) with same-time
changes ordered as in the simulation, and one counter splits exposure and
assigns events. It replaces the per-calculator counting functions and is the
basis for milestone 9's additional group columns. Results equal the previous
calculators' output except for row order, the fixed `union_type` level order
and the new variables.

**Implications:** The old arguments are removed, not deprecated. Accepted
reports for milestones 1 and 3–8 use them and reproduce at their recorded
commits. See
[demographic analysis](demographic-analysis.md#stratifying-results-with-by-d133).

## D134 — Multiple group dimensions (November 5, 2026)

**Status:** User-selected design for milestone 9 (D131), decided in discussion
November 5, 2026. Step 1 implemented in stages 1–5 (November 5–6, 2026), each
committed by the user, and stage 6 (tests, articles, documents). Step 2
(matching across several dimensions, calibrated over their combinations)
implemented November 6, 2026, following the user's choices recorded in
`multiple-groups.md`: the combination covers every dimension any matching
model of the run uses, effects add, and custom cross-dimension functions are
left for later. Accepted complete November 6, 2026 (D141).

**Decision:** Populations may carry several categorical group dimensions,
always declared with `set_groups(model, race = group(...), ...)`; there is no
implicit single `group` column. `group()` holds a dimension's optional ordered
`categories` (default inferred; replaces inheritance's `groups` and the
`group_levels` of matching and weight models, with the first category as
reference), `ancestry` (default `FALSE`) and `inheritance` (moved from
`set_fertility()`). The declaration is stored with the population and carried
by results. Dimensions use their own names, with reserved names and derived
names as listed in the design. Rate tables may be stratified by any declared
dimension; transitions are set per dimension with
`set_group_transition(model, rates, group = ...)` and use `origin` and
`destination`. Matching settings become per-dimension lists, calibrated
jointly across the dimensions a matching model uses. Baseline populations take
`groups` as independent weights or a joint table and carry a default
declaration. Ancestry is a named list of tibbles.

**Staging:** step 1 delivers everything except joint matching calibration,
with group effects in matching limited to one dimension; step 2 adds joint
calibration.

**Later user decisions (stage 2, November 5, 2026):** dissolution and
conversion rates are directional, with `<dim>_w` and `<dim>_h` columns (the
wife's and husband's categories); the symmetric `group_pair` is removed, and a
symmetric rate lists both directions. Search rates and cohabitation
probabilities use the dimension's own name (the woman's category).

**Implementation choices** (detailed in [multiple group dimensions](multiple-groups.md)):
one transition process per dimension, stored by dimension, with histories
holding `group`, `origin` and `destination` as character labels; uncovered
transition strata have zero risk; analysis splits by several dimensions
through one combined internal category; matching models are translated at
simulation to an internal one-dimension form; a bare matching setting applies
to the only declared dimension; summaries keep every non-measure column as a
stratum.

**Implications:** Single-group code, tests, articles and reports change
(`set_groups(model, group = group())` reproduces a single dimension named
`group`). Accepted reports reproduce at their recorded commits. See
[multiple group dimensions](multiple-groups.md).

## D135 — Prohibitions as NA in `group_exogamy` (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** Remove `prevent_exogamy`. `group_exogamy` carries both pair
effects and prohibitions: a finite value is a symmetric log-odds effect, NA
prohibits that unordered pair in both directions, and NA for a whole dimension
(`list(region = NA)`, or a bare `NA` with one dimension) prohibits every union
between different categories. Unlisted pairs have no effect. A within-category
prohibition is an error.

**Rationale:** One argument, and one meaning of `:`. `prevent_exogamy` was
directional (`"A:B"` female-first, `"A::B"` both), unlike the symmetric pair
effects of D075, which was inconsistent.

**Supersedes:** the directional restriction syntax and the separate
`prevent_exogamy` argument (October 20, 2026 restrictions; D072's pool-wide
eligibility rule is unchanged). Directional prohibitions are no longer
available.


## D136 — First public release after the current milestone (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** Finish the current milestone (9, multiple group dimensions),
then work toward the first version bump, coinciding with the first public
release. External immigration and two-sex marriage models move to post-release
projects; they are not prerequisites for this release.

**Rationale:** The user considers the current feature set sufficiently close
to the intended SOCSIM scope, with added flexibility, to prioritize a public
release over these additional mechanisms. This is a scope judgment, not a
verified claim of complete SOCSIM equivalence.

**Implications:** Supersedes D131's immediate post-milestone-9 sequence and
D117/D121's placement of two-sex marriage models before the first release.
Milestone 9 still requires completion and user acceptance. The version number,
release channel, release criteria, and the scope/timing of proposed summaries
and helpers remain unselected. This decision does not authorize publication.
Existing immigration, family-population construction and generation discussions
remain future design material; no immigration mechanism or API is implemented
or authorized by this decision.


## D137 — Align development with GitHub milestones and Projects (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** Use GitHub milestones and Projects to organize development. A
milestone may be delivered through multiple issues and pull requests; it does
not require a single PR or a branch spanning the entire milestone.

**Working approach:** Milestones group work toward a feature or release goal.
Issues describe bounded tasks, their dependencies and completion criteria. PRs
deliver reviewable changes and link the relevant issues. Projects track the
status and priority of that work. Keep scientific contracts, worked examples,
accepted decisions and validation evidence in the repository, with links from
GitHub items rather than duplicate specifications. Merging a component PR does
not establish scientific acceptance of the milestone.

**Current mapping:** GitHub has open milestones for multiple groups (9),
external immigration (10), cleanup/documentation/helpers/polish (11), and
two-sex marriage (12). Under D136, immigration and two-sex marriage are
post-release work; an open milestone does not make them release prerequisites.
The cleanup milestone is a starting point for release planning, whose concrete
scope and version number remain to be selected. This entry records the
workflow; it does not claim that GitHub items or a Project have been updated.

## D138 — Matching calibration within separate markets (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** When matching prohibitions split the matching categories into
blocks that cannot match across (such as regions under `region = NA`),
calibrate each block separately, with mirrored male targets normalized
within the block. Each block is a separate marriage market.

**Rationale:** Women search only within their block, so the block is where
male selection should be balanced. The previous population-wide targets were
unreachable whenever a block's men and women differed in their shares, which
made the solver fail and slowed simulation about tenfold. Sex symmetry of
intermarriage is preserved within blocks and so also when summed over them;
what is given up is only the unreachable population-wide target.

**Implications:** Results change only for models whose prohibitions split the
categories into several blocks. A block with more men than women of marrying
age leaves more men unmatched.

## D139 — Table form of `group_exogamy` for cross-dimension effects (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** `group_exogamy` takes either the list form (per-dimension named
vectors or NA, adding across dimensions) or a table form: a data frame with
`<dim>_x`/`<dim>_y` columns and `log_odds` that fully specifies symmetric
pair effects over combinations of categories, NA prohibiting. The forms are
not combined. Each row is one unordered pair; unlisted pairs have no effect;
rows have no wildcards. The earlier per-dimension pair table
(`female_group`/`male_group`) is removed.

**Rationale:** One place per model for pair effects. Combining a table with
per-dimension vectors (each adding to the other) was judged easy to misread;
the list form remains the shortcut when there are no cross-dimension effects.
Neutral `_x`/`_y` names keep the symmetry visible. The table is translated
into the same native effect table, so it is faster than formula group terms,
which are evaluated in R during draws.

## D140 — Household moves (November 6, 2026)

**Status:** User decision, November 6, 2026. Specified in
[group transitions](group-transitions.md#household-moves-d140); implemented
November 6, 2026, with two user additions after testing (the party is built
recursively, so a dependent's own dependents follow; a dependent re-attached
to a parent elsewhere joins that parent at once, recorded as the parent's
move), with the rule as its own object (`household()`, per
the user) and the cutoff named `age_independence`.

**Decision:** A household rule (`household = household(age_independence = 18)` on
`set_group_transition()`, `NULL` by default for individual moves; each process
takes its own) makes a move carry the mover's partner and dependents. Dependents are living,
unpartnered people under `age_independence` (user-settable, default 18) attached
to their mother while she lives, otherwise to their father; people with
neither parent alive are independent (no grandparents or stepparents). One
clock per household: the woman in a union carries it, a partnered man and
dependents have none, everyone else has their own. The moving party is the
carrier, their partner and the dependents of either. The transitions history
gains `move_id` and `initiator_id`.

**Rationale:** Household moves were deferred by D114 with open questions.
SOCSIM's assembly (audited November 6, 2026) has no household entity, lets
every member initiate (so input rates are initiation rates and moves are
double counted relative to observed individual rates), uses different rules
for male and female movers, and has a defective child-age test. One clock per
household lets observed individual rates be used as input; one attachment
rule (mother, then father) treats children consistently.

**Supersedes:** the deferral in D114 for transitions; household moves remain
unavailable for other processes.

## D141 — Milestone 9 accepted complete (November 6, 2026)

**Status:** User decision, November 6, 2026.

**Decision:** The user accepted milestone 9, multiple group dimensions, as
complete, following the user-authored assessment in
[`validation/milestone9.qmd`](../validation/milestone9.qmd), which the user
maintains and renders. Scope accepted: declared group dimensions (D134, steps
1 and 2), directional union strata, `stratify_rates()` (D132) and one `by`
argument (D133), NA prohibitions in `group_exogamy` (D135), separate-market
calibration (D138), the table form of `group_exogamy` (D139), household moves
(D140), and baseline populations declaring their categories. Code at
acceptance: `fc4b9dbf29dba9782cec9265704f9e87fef4022c` plus the closeout
documentation.

**Limitations retained:**

- Matching calibration cost grows steeply with the number of matching
  category combinations (about 4 to 4.5 times per doubling); separate markets
  (D138) make prohibited dimensions cheap.
- Household moves (D140) require age bands, and followers move less often
  than the input rate when households are split across regions.
- Uncovered transition strata have zero risk, unlike other rate tables.
- Matching formula group terms are evaluated in R during draws (slower than
  `group_exogamy`) and can express directional effects; there are no custom
  cross-dimension matching functions (the table form covers pair effects).
- Transition history labels are always character.

**Also recorded (user preference):** the user writes all vignettes and
articles and makes all edits to the top-level `README.md`; agents do not edit
them.

## D142 — Fixed newborn group category (November 7, 2026)

**Status:** User decision, November 7, 2026.

**Decision:** Add a built-in rule assigning every newborn to a user-specified
category of a group dimension. The motivating case is education: all children
start in the bottom education category regardless of parental education.
Existing identities and ancestry averaging are unchanged; later transitions
remain available. Assignment requires no R callback or identity random draw.

**Interface:** `create_group_inheritance("fixed", fixed_value = "none")`, attached
through `group(inheritance = ...)`. `fixed_value` names a category in the dimension's
registry, rather than implicitly selecting the first category. Calendar schedules
use the ordinary birth-time selection contract. See [the specification](group-inheritance.md).

**Supersedes:** earlier deferral of fixed newborn assignment; random-parent
inheritance remains the default for dimensions without an explicit rule.


## D143 — R source organization (November 9, 2026)

**Status:** User decision and implementation authorization, November 9, 2026.

**Decision:** Organize R source files by cohesive responsibility, keeping related public functions and their specific helpers together. Use the exact lead function name for a filename when there is a clear lead, preferring the constructor over its convenience setter; otherwise use a descriptive responsibility name. The agreed map is implemented in [R source organization](r-source-organization.md).

Observed event-rate calculators and rate-specific helpers share `calculate_rates.R`. Shared history-analysis preparation lives in `analysis_utils.R`; input-rate validation and engine encoding live in `rates_util.R`. HMD/HFD import functions and all their helpers share `get_database_rates.R`. TFR/PATFR share `calculate_tfr.R`. Matching is separate from union-formation process configuration; transition histories are separate from transition-process configuration. Keep `model_schedules.R` and `set_population()` in `hammel_model.R` as agreed. The `simulate()` body remains intact.

**Future additions:** When a new function has no logical home in an existing R file, explain the proposed responsibility and filename and ask the user before creating a new R script. Do not create one automatically. This requirement is also in `AGENTS.md`.

**Scope:** Mechanical organization, preserving function definitions, constants, source comments, exports, and help topics. Help-topic redesign, pkgdown reference grouping, new demographic behavior, and simulation-body refactoring are separate work. Explicit roxygen ordering preserves the existing shared-help display order after files are split.


**Within-file organization (user instruction, November 9, 2026):** Start each hand-maintained R script with a brief purpose comment above the section headings. Put exported functions before internal constants and helper functions, with `# PUBLIC FUNCTIONS ------`, `# INTERNAL CONSTANTS ------`, and `# HELPER FUNCTIONS ------` headings where applicable. Preserve function bodies and documentation blocks. Keep generated wrappers untouched; package documentation uses its own section.


## D144 — Shared help for process constructors and setters (November 9, 2026)

**Status:** User decision and implementation authorization, November 9, 2026.

**Decision:** Each demographic process constructor shares a help topic with its corresponding convenience setter. The five pairs are mortality, fertility, union formation, union dissolution, and group transitions. Each page explains the specification once, shows both interfaces, and distinguishes the returned process specification from the updated model. Neither interface executes a simulation.

Keep generic `set_process()` on its own topic for attaching an existing specification, with links to the process-specific topics. Separate the formerly combined mortality/fertility constructor page and remove their setters from the generic attachment page. The other three process pairs already shared topics; clarify their constructor/setter relationship. This is documentation organization, with no function, argument, default, or demographic behavior changes.


## D145 — Shared help for observed event-rate calculators (November 9, 2026)

**Status:** User decision and implementation authorization, November 9, 2026.

**Decision:** Combine the seven observed event-rate calculators into the `calculate_rates` help topic: mortality, fertility, union formation, marriage, dissolution, cohabitation conversion, and group transitions. Preserve each function's own help alias. Document common arguments once, provide an event/exposure comparison, and retain separate sections describing each function's event definitions, eligibility, timing, strata, exposure denominator, and output. Include examples covering all seven calculators without running a simulation or downloading rates.

Keep population size, TFR/PATFR, and life tables on separate help topics. This changes documentation organization only; function names, arguments, defaults, exports, and executable behavior remain unchanged.


## D146 — Shared fertility-summary and group-declaration help (November 9, 2026)

**Status:** User decision and implementation authorization, November 9, 2026.

**Decision:** Combine `calculate_tfr()` and `calculate_patfr()` on the `calculate_tfr` help topic, explaining the age-only sum, the conditional age/parity fertility table, automatic dispatch when parity columns are present, and the differing output column names. Preserve the distinct scientific conventions of both calculations.

Combine `group()` and `set_groups()` on the `group` help topic, distinguishing a single dimension's settings from declaring named dimensions on a model. Retain each function's help alias, arguments, and relevant examples. Leave `create_matching_model()` and `create_age_gap_matching_model()` on separate help topics, as explicitly requested. No executable behavior, defaults, or exports change.


## D147 — Noun names for specification constructors (November 9, 2026)

**Status:** User decision and implementation authorization, November 9, 2026.

**Decision:** Replace the four `create_*` specification constructors with noun-based names:

| Former name | Current name |
| --- | --- |
| `create_matching_model()` | `matching_model()` |
| `create_age_gap_matching_model()` | `age_gap_matching_model()` |
| `create_weight_model()` | `weight_model()` |
| `create_group_inheritance()` | `group_inheritance()` |

Keep `create_baseline_population()` unchanged. Verbs remain appropriate for pipeable actions such as `set_population()`; standalone specification constructors use nouns. Rename the three corresponding source files to match their lead functions. Matching constructors retain separate help topics, as requested under D146.

**Compatibility and scope:** This is a direct pre-release API rename, without aliases under the former names. Defaults, classes, validation rules, stochastic behavior, and model semantics remain unchanged. Update package code, tests, generated help, and current design specifications. Preserve user-maintained README/vignettes, historical validation reports and supporting scripts, and dated verification/handoff records. Those sources can retain historical calls; running them against the renamed API requires updating the calls or using their recorded software version. Earlier decision entries retain their original constructor spellings; this entry governs current names.


## D148 — Snake-case confidence argument and function-style review (November 9, 2026)

**Status:** User instruction and implementation authorization, November 9, 2026.

**Decision:** Rename Hammel's `conf.level` argument to `conf_level` on all seven observed event-rate calculators and their internal helpers. Update tests, generated help, diagnostics, and current specifications. Keep the external `stats::poisson.test(conf.level = conf_level)` call intact. This is a direct pre-release argument rename without an alias for the former spelling; default output, interval calculation, column names, and exports are unchanged. This supersedes the argument spelling in D124, not its statistical contract.

**Review:** Check the public API and hand-maintained R code against tidyverse naming, function-writing, and syntax conventions, retaining the user's noun constructors, purpose comments, and uppercase section headings. Record the scope, findings, and limits in [R source organization](r-source-organization.md). Broader formatting and unrelated API changes are not part of the confidence-argument rename. User-owned articles and historical validation evidence remain intact.


## D149 — Eighty-column R source convention (November 9, 2026)

**Status:** User instruction and implementation authorization, November 9, 2026.

**Decision:** Wrap hand-maintained `R/` source to at most 80 characters per line, including comments and roxygen documentation, and follow this convention in future R code. For long function signatures use one argument per line with a two-space indent. Wrap long calls and expressions at syntax boundaries; preserve explanatory comments, public-first sectioning, documentation content, names, defaults, and demographic behavior. Generated wrappers and user-owned articles/README/historical validation evidence remain untouched. Markdown prose need not be hard-wrapped.

Long diagnostic literals may use `paste0()` of literal pieces so their complete text is unchanged. Small braced branches are permitted where needed to lay out long conditional expressions. Verify the parsed expressions modulo constant string concatenation and single-expression conditional braces, compare generated help content and example expressions, and run the offline suite. The tooling may use a temporary formatting library/cache without adding a package dependency. Other control-flow style findings from D148 remain separate work.
