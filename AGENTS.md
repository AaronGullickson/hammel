# Working on Hammel

Start with the [documentation index](docs/README.md). Read [the project
brief](docs/project-brief.md), [the design specification](docs/design.md), and [the
scientific contract](docs/scientific-contract.md) before consequential implementation.
Consult [the decision log](docs/decisions.md) for the status and rationale of choices
and [the audit plan](docs/socsim-audit.md) when investigating SOCSIM.

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

- Prefer verb-based names for user-facing action functions, such as `create_baseline_population()` and `set_population()`.
- Return all public tables as tibbles, including empty tables, population histories,
  process schedules, and demographic summaries. Accept base data frames as inputs.
- Keep model construction separate from simulation execution. Aim for composable R workflows and R-native, inspectable results at the user-facing boundary. This does not require an R-only implementation: performance-critical computation and internal state may use compiled code, such as C++, with results exposed as ordinary R objects for inspection and analysis.
- Each simulation must own independent state. Avoid hidden process-wide population state and undocumented effects on later runs.
- Make time, rate units, random-number behavior, and observation semantics explicit.
  Public times/durations are years, exposure is person-years, and input rates are
  per person-year. The current engine executes continuously (D024), with no
  monthly increment or execution-step parameter. Births use cell-accounted
  fertility (D101): each one-year age x actual parity x union type x group cell
  produces births at its input rate times all members, assigned to a random
  eligible member, weighted by persistent hetfert when enabled (D102). Births
  in cells with no eligible woman are reported as
  unrealized. See docs/fertility-cell-accounting.md. Midpoint compensation
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
  archives in `docs/`. Keep rationale and supersession in `docs/decisions.md`.
  Do not make readers reconcile contradictory current rules with appended
  corrections.
- Milestone 1 was accepted as complete on October 15, 2026 (D029). Milestone 2
  union formation was accepted as complete on October 19, 2026 (D050). Milestone 3
  groups was accepted as complete on October 27, 2026 (D077). Milestone 4
  dissolution/remarriage/cohabitation was accepted October 29, 2026 (D093).
  Milestone 5 parity-specific fertility, cell accounting and heterogeneity was
  accepted October 29, 2026 (D105). Group transitions (milestone 6) are next;
  detailed semantics and API need specification before implementation.
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
  assessment. Preserve all five reports and rendered
  artifacts; do not fabricate results or set scientific acceptance thresholds. Milestone
  completion retains documented limitations and does not authorize future features.
- Milestone 4 includes duration-dependent dissolution, age/duration searches by
  previous union, and optional rate-mode cohabitation with linked conversion
  episodes and type-specific matching/rates (D078–D092). One woman carries each
  union clock; exact ties prioritize death, birth, dissolution, conversion, search.
  See `docs/union-dissolution.md`, `docs/remarriage.md`, and `docs/cohabitation.md`.
  Supporting benchmarks/recovery evidence lives in `validation/milestone4/`;
  top-level reports and `validation/data/` remain in place. Preserve active render
  intermediates until rendering finishes, including milestone 5's active final
  render at closeout. Parity and fertility heterogeneity are implemented; group
  transitions (milestone 6) and calendar-time schedules (milestone 7) require
  their own scientific specifications. Most categorical characteristics should
  use groups, with multiple types of groups an optional complexity direction.
- Use `devtools::test()` for the testthat suite. Live database tests require
  explicit opt-in; ordinary verification must not download rates or expose
  credentials. See `docs/testing.md` for commands and evidence.
- Treat session handoffs as dated status records, not standing pause instructions.
  The current user request controls work. Preserve the user's edits, including
  explanatory source comments and the accepted reports. Fixed newborn assignment
  is deferred to group transitions; custom inheritance already supports assigning
  mixed-parent births to a specified category.
