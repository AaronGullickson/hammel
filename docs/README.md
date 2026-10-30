# Hammel documentation

Updated October 30, 2026. Milestones 1–5 are complete by user acceptance
(D029, D050, D077, D093, D105), including parity-specific fertility,
cell-accounted birth spacing and persistent lognormal fertility weights.
Current behavior uses continuous event times with no monthly grid. Acceptance
is based on the user-authored validation reports and retains their limitations.
**Milestone 6 is in progress:** individual group transitions are implemented;
scientific milestone acceptance remains pending. See the
[group-transition contract](group-transitions.md). Household moves are deferred to
the regions milestone (D114); milestone 7 is time-varying inputs (D103).

## Reading order

| Document | Purpose |
| --- | --- |
| [Project brief](project-brief.md) | Purpose, current scope, and milestone status |
| [Design](design.md) | Implemented architecture, API, and source locations |
| [Scientific contract](scientific-contract.md) | Current demographic assumptions and event rules |
| [Population schema](population-schema.md) | Person/union tables and structural validation |
| [Group transitions](group-transitions.md) | Individual transitions, histories, state updates and age/duration Lexis rates |
| [Groups](groups.md) | Group-specific schedules, multipliers, newborn identity, and baseline sampling |
| [Ancestry](ancestry.md) | Starting proportions, parental averaging, unknown shares, and restarts |
| [Group inheritance](group-inheritance.md) | Built-in/custom identity rules, passing, and possible groups |
| [Calibration intensity updates](calibration-intensity.md) | Actual female hazards, aggregate cells, and bookkeeping boundary queue (D082) |
| [Weighted-choice matching](weighted-choice-matching.md) | Calibrated male-group weights and delta target adjustments (D074/D075) |
| [Group matching](group-matching.md) | Symmetric pair effects and fixed group formula predictors |
| [Demographic analysis](demographic-analysis.md) | Separate mortality/fertility rates, parity ASFRs, union rates, TFR, and life tables |
| [Cell-accounted fertility](fertility-cell-accounting.md) | Birth-spacing mechanism (D101): cells, unrealized births, validation |
| [Parity spacing investigation](parity-spacing-investigation.md) | Evidence that led to D101, including the SOCSIM reference |
| [Testing and validation](testing.md) | Test commands, evidence, and remaining validation |
| [Current handoff](session-handoff.md) | Recent findings, active work, and deferred issues |

Read the brief, design, and scientific contract before consequential changes.
[AGENTS.md](../AGENTS.md) contains the working rules. Model choices require the
user's agreement; proposed APIs and research notes do not authorize implementation.

## Explanations, proposals, and provenance

| Document | Status |
| --- | --- |
| [Waiting-time sampling](waiting-time-sampling.md) | Explanation of the current shared hazard sampler |
| [Union-formation models](union-formation-models.md) | Implemented milestone 2 rules, matching defaults, diagnostics, and deferred work |
| [Cohabitation](cohabitation.md) | Rate-mode type probabilities, type-specific matching/rates, conversion and histories (D085/D086) |
| [Remarriage schedules](remarriage.md) | Age/duration previous-union search strata, calibration, and marriage observations (D080/D083) |
| [Fertility heterogeneity](fertility-heterogeneity-notes.md) | D102: persistent lognormal weights, maternal inheritance, validation and selection caveat |
| [Parity fertility and PATFR](parity-fertility.md) | History-derived parity, schedules, HFD recovery and PATFR; milestone 5 accepted (D105) |
| [Union dissolution](union-dissolution.md) | Milestone 4 duration hazards, female carrier, state changes, and worked examples (D078) |
| [Kinship infrastructure](kinship.md) | Shared native distance screening and R extraction |
| [Two-sex union models](two-sex-union-models.md) | Proposal notes: directional asymmetry, marriage squeeze, Schoen/Choo–Siow/Dagsvik models, API considerations (D073) |
| [Future roadmap](roadmap.md) | Accepted milestones and future priorities with unresolved specifications |
| [Decision log](decisions.md) | Accepted choices and supersession history |
| [SOCSIM audit plan](socsim-audit.md) | Scope and remaining reference investigation |
| [SOCSIM fertility heterogeneity](socsim-fertility-heterogeneity-audit.md) | Pinned fmult distribution, inheritance, defaults, switch checks and loader caveat |
| [SOCSIM rates](socsim-rate-audit.md) | Pinned static source findings |
| [SOCSIM competing events](socsim-competing-events-audit.md) | Pinned static scheduling findings |
| [SOCSIM union formation](socsim-union-formation-audit.md) | Pinned static search, matching, and eligibility findings |
| [SOCSIM dissolution and remarriage](socsim-dissolution-audit.md) | Pinned static divorce clocks, remarriage schedules, state updates, and cohabitation findings |
| [SOCSIM cohabitation](socsim-cohabitation-audit.md) | Pinned static formation/type assignment, conversion, breakup, legal history and output trace |
| [SOCSIM group transitions](socsim-group-transitions-audit.md) | Pinned static transition clocks, household moves, histories, rescheduling and source caveats |

Git retains earlier specifications and verification records; do not keep
snapshot copies in the docs tree. The decision log records rationale and
supersession. Current contracts take precedence over earlier monthly timing,
class names, and fertility compensation descriptions.

To inspect prior versions, use `git log -- docs` and
`git show <commit>:docs/scientific-contract.md`.

## Validation artifacts and examples

- [Milestone 1 validation report](../validation/milestone1.qmd): completed
  assessment by Aaron Gullickson; the user rendered `validation/milestone1.html`
  and accepted milestone 1. Reproduce at the recorded acceptance commit; current
  execution requires a union process absent from that historical assessment.
- [Milestone 2 validation report](../validation/milestone2.qmd): user-authored
  assessment accepted October 19, 2026 (D050), using USA circa 2000–2001 inputs; includes both initiation
  modes, matching-parameter recovery, and performance comparisons. Rendering
  executes simulations and an HMD download. Deferred limitations remain documented.
- [Union-formation example](../inst/examples/milestone2.R): current construction
  and execution workflow.
- [Milestone 3 validation report](../validation/milestone3.qmd): user-authored
  assessment accepted October 27, 2026 (D077), covering group-specific rates,
  weighted matching and delta, conditional-logit recovery, ancestry, built-in
  and custom inheritance, baseline sensitivity, and performance. Rendering
  downloads HMD inputs and runs substantial simulations and the test suite.
- [Milestone 4 validation report](../validation/milestone4.qmd) and
  [rendered assessment](../validation/milestone4.html): user-authored assessment
  accepted October 29, 2026 (D093), covering dissolution, age/duration remarriage,
  cohabitation, age/duration conversion, fertility by type, separate matching
  roles, and performance. Rendering downloads HMD mortality and runs 100
  independent-baseline replications per scenario plus the offline test suite.
- [Milestone 5 validation report](../validation/milestone5.qmd) and
  [rendered assessment](../validation/milestone5.html): accepted October 29,
  2026 (D105), covering parity-specific rates, cell-accounted spacing, baseline,
  restart and union-type fertility checks, heterogeneous weights, maternal
  inheritance and performance. The final render was active at closeout;
  preserve all intermediates. Rendering installs the package, downloads HFD/HMD
  inputs and runs substantial replications and offline tests.
- [Validation directory guide](../validation/README.md): retained reports,
  rate provenance, and supporting benchmark/conditional-choice evidence.
The former `milestone1.R` and `milestone1_bootstrap.R` examples have been
removed from the current tree; use Git history for those historical workflows.

`inst/examples/` ships with the installed R package. `validation/` and `docs/`
are development material excluded from package builds. Preserve the completed
report sources and rendered artifacts. Some local inputs and generated files
are Git-ignored; reproducible evidence needs its input provenance and software
environment as well as a documented archive or sharing location.
