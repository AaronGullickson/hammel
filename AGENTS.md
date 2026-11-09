# Working on Hammel

Start with the [documentation index](design/README.md). Read [the project
brief](design/project-brief.md), [the design specification](design/specification.md), and [the
scientific contract](design/scientific-contract.md) before consequential implementation.
Consult [the decision log](design/decisions.md) for the status and rationale of choices
and [the audit plan](design/socsim-audit.md) when investigating SOCSIM.

## Purpose and scope

Hammel is a new, standalone demographic microsimulator with a modern R package
interface, informed by SOCSIM's approach and features, but offering more flexibility to
users to create customizable model parameters. It is a substantial rewrite, not a SOCSIM
wrapper, rsocsim companion package, or upstream contribution. Legacy file compatibility
must not determine its architecture.

The user guides demographic requirements and scientific interpretation. Agents
investigate, implement, document, and validate software. Do not silently resolve
questions that change the demographic model.

## Before implementing a feature

1. Identify the requested behavior and the relevant accepted decisions.
2. For consequential simulation behavior, write or update its plain-language scientific specification and worked examples first. Identify units, eligibility, event timing, edge cases, and expected outcomes.
3. Identify unresolved modeling choices. Explain their consequences to the user and obtain a decision before implementing behavior that depends on them. Continue independent investigation or infrastructure work where possible.
4. Define validation from the specification, including analytical cases or invariants where applicable.
5. Implement a narrow, reviewable change and keep the documents consistent with the resulting behavior.

Suggestions, illustrative APIs, and proposed milestones are not accepted decisions or
authorization to implement them. Routine software choices within an authorized task can
be made autonomously; record choices with substantial architectural consequences.

## Engineering expectations

- Organize R functions by cohesive responsibility using
  [the source organization guide](design/r-source-organization.md). Place new
  functions in an existing file when its responsibility fits. If no existing
  file is a logical home, explain the proposed responsibility and filename and
  ask the user before creating a new R script; do not create one automatically.
- Start each hand-maintained R file with a brief purpose comment. Place
  exported functions first and internal
  constants/helpers below, using `# SECTION NAME ------` headings. Preserve
  documentation blocks and leave generated wrappers untouched.
- Use noun-based names for specification constructors, such as `matching_model()`,
  `age_gap_matching_model()`, `weight_model()`, and `group_inheritance()`.
  Prefer verbs for pipeable actions such as `set_population()`. Keep
  `create_baseline_population()` as the explicitly retained exception.
- Keep hand-maintained `R/` source lines at most 80 characters, including
  comments and roxygen documentation. Follow the same limit for new or edited
  R code elsewhere. Wrap long calls and use one argument per line with a
  two-space indent for long function signatures. Split long diagnostic strings
  without changing their text. Leave generated wrappers untouched.
- Use braces for all loop bodies and statement-style or multiline `if`/`else`
  branches. Keep `else` on the same line as the preceding closing brace.
  Compact, single-line conditional value expressions may remain unbraced when
  they have no side effects or control-flow modifiers. Put `return()`, `stop()`,
  `break`, and `next` on their own lines inside the relevant block; brace
  control-flow/error-handler callback bodies as well.
- Use snake_case for package-owned function and argument names (for example,
  `conf_level`). Preserve external APIs when calling them, such as
  `stats::poisson.test(conf.level = conf_level)`. Follow tidyverse function-writing
  conventions subject to the explicit project naming and section-heading rules;
  see the source organization guide for the current style review.
- Return all public tables as tibbles, including empty tables, population histories,
  process schedules, and demographic summaries. Accept base data frames as inputs.
- Keep model construction separate from simulation execution. Aim for composable R workflows and R-native, inspectable results at the user-facing boundary. This does not require an R-only implementation: performance-critical computation and internal state may use compiled code, such as C++, with results exposed as ordinary R objects for inspection and analysis.
- Each simulation must own independent state. Avoid hidden process-wide population state and undocumented effects on later runs.
- Make time, rate units, random-number behavior, and observation semantics explicit.
  Public times/durations are years, exposure is person-years, and input rates are
  per person-year. The current engine executes continuously (D024), with no
  monthly increment or execution-step parameter. Every rate process is
  scheduled by cells with one exponential clock per process (D128); there are
  no individual waiting-time clocks. Births use cell-accounted
  fertility (D101): each one-year age x actual parity x union type (x group,
  when schedules are group-specific) cell produces births at its input rate
  times all members, assigned to a random eligible member, weighted by
  persistent hetfert and any weight model when enabled (D102, D129). Draws of
  cells with no eligible woman are counted as unrealized births. See
  design/binned-group-hazards.md and design/fertility-cell-accounting.md. Midpoint compensation
  (D028/D095/D099) was removed. Do not reintroduce superseded monthly timing,
  midpoint compensation or the specialized exact fertility sampler without an
  explicit modeling decision.
- Design built-in and custom demographic mechanisms around a documented contract; arbitrary callbacks must not bypass state consistency or event invalidation rules.
- Preserve scientifically relevant provenance. Before copying or adapting source, establish its version, license, and applicable obligations.
- C++17/Rcpp is selected for the first engine (D019); tibble is the public table
  convention (D023). Do not replace the engine language or add major dependencies
  merely because future notes mention a candidate.
- Do not invent demographic defaults, fabricate audit findings, or claim scientific validity from plausible-looking output.
- Do not start sub-agents solely because the project notes describe possible agent roles. Use delegation when the current task authorizes it.

## Validation and reporting

Use tests appropriate to the behavior: analytical cases, accounting identities, valid
dates and relationships, reproducibility, and independence between runs. Compare with
SOCSIM over multiple replications when a reference comparison is appropriate. Exact
fixed-seed agreement is meaningful only when algorithms and random-number consumption
are preserved; otherwise use justified statistical comparisons.

Report what changed, what was checked, and remaining limitations. Distinguish
implemented behavior, verified findings, proposals, and unresolved assumptions. Document
intentional departures from SOCSIM. Never mark a scientific decision accepted without
the user's agreement.

## Documentation and current priorities

- Current specifications describe implemented behavior. Use Git history for
  superseded specifications and verification records; do not maintain snapshot
  archives in `design/`. Keep rationale and supersession in `design/decisions.md`.
  Do not make readers reconcile contradictory current rules with appended
  corrections.
- Milestone 1 was accepted as complete on October 15, 2026 (D029). Milestone 2
  union formation was accepted as complete on October 19, 2026 (D050). Milestone 3
  groups was accepted as complete on October 27, 2026 (D077). Milestone 4
  dissolution/remarriage/cohabitation was accepted October 29, 2026 (D093).
  Milestone 5 parity-specific fertility, cell accounting and heterogeneity was
  accepted October 29, 2026 (D105). Individual group transitions (milestone 6)
  were accepted October 30, 2026 (D116). Milestone 7 calendar-varying rates are
  implemented under the user-approved D118 contract in
  `design/time-varying-rates.md`; the user accepted milestone 7 as complete on
  November 2, 2026 (D127), based on `validation/milestone7.qmd` and its retained
  rendered assessment. Matching
  and inheritance model schedules are implemented under D123; see
  `design/time-varying-models.md` for the shared model-list tibble interface,
  event-time selection, dynamic NULL fallbacks and future identity registry.
  D122 adds calendar
  bands to `cohab_prob`; probability-only changes update type selection and
  matching calibration while preserving pending event dates. Other arguments
  remain static. Milestone 7 was completed with the previous engine. Milestone
  8, cell scheduling for every rate process (D128) and weight models (D129),
  was accepted complete on November 5, 2026 (D130), based on
  `validation/milestone8.qmd`; see `design/binned-group-hazards.md`.
  Milestone 9, multiple group dimensions (D134; see
  `design/multiple-groups.md`), including matching across several dimensions,
  separate-market calibration (D138) and household moves (D140), was accepted
  complete on November 6, 2026 (D141), based on `validation/milestone9.qmd`.
  Under D136 the next focus is the first public release; external immigration
  and the two-sex model are post-release projects.
  Distinguish current
  specifications from deferred mechanisms and future proposals.
- `validation/milestone1.qmd` is the accepted historical assessment; reproduce it
  at its recorded commit because current execution requires a union process.
  `validation/milestone2.qmd` contains the user's accepted assessment, including
  matching-parameter recovery and performance comparisons.
  `validation/milestone3.qmd` contains the accepted group-rate, matching, ancestry,
  inheritance, and performance assessment. `validation/milestone4.qmd` is the
  accepted dissolution/remarriage/cohabitation assessment.
  `validation/milestone5.qmd` is the accepted parity/spacing/heterogeneity
  assessment. `validation/milestone6.qmd` is the accepted individual-transition
  assessment (D116), with its final render active at acceptance.
  `validation/milestone7.qmd` and `validation/milestone8.qmd` are the accepted
  calendar-input and cell-scheduling/weight-model assessments, and
  `validation/milestone9.qmd` the accepted multiple-dimension assessment.
  Preserve all nine reports, rendered artifacts and any active render
  intermediates; do not
  fabricate results or set scientific acceptance thresholds. Milestone
  completion retains documented limitations and does not authorize future features.
- Milestone 4 includes duration-dependent dissolution, age/duration searches by
  previous union, and optional rate-mode cohabitation with linked conversion
  episodes and type-specific matching/rates (D078–D092). Unions are the units
  in dissolution and conversion cells (D128); exact ties prioritize death,
  birth, dissolution, conversion, search, then group transition.
  See `design/union-dissolution.md`, `design/remarriage.md`, and `design/cohabitation.md`.
  Supporting benchmarks/recovery evidence lives in `validation/milestone4/`;
  top-level reports and `validation/data/` remain in place. Preserve active render
  intermediates until rendering finishes, including milestone 5's active final
  render at closeout and milestone 6's active final render. Parity, fertility
  heterogeneity and individual group transitions are implemented and accepted.
  See `design/group-transitions.md`; household moves, deferred by D114, are
  implemented in milestone 9 as `household()` on a transition process
  (D140), with regional applications under multiple group columns (D121,
  D131). Calendar-time rates follow `design/time-varying-rates.md`
  (D118); matching and inheritance model schedules are implemented under D123. Most categorical characteristics should
  use groups; multiple group dimensions follow D134.
- Use `devtools::test()` for the testthat suite. Live database tests require
  explicit opt-in; ordinary verification must not download rates or expose
  credentials. See `design/testing.md` for commands and evidence.
- Treat session handoffs as dated status records, not standing pause instructions.
  The current user request controls work. The user writes all vignettes and
  articles (`vignettes/`) and makes all edits to the top-level `README.md`
  (D141); do not edit them unless asked. Preserve the user's edits, including
  explanatory source comments and the accepted reports. Fixed newborn assignment
  is available through `group_inheritance("fixed", fixed_value = ...)` (D142); custom inheritance already supports assigning
  mixed-parent births to a specified category.
