# Validation reports and supporting evidence

Milestones 1–6 are accepted by the user (D029, D050, D077, D093, D105, D116). Keep the
source reports, final rendered HTML, explanatory comments, and input provenance.
Acceptance retains the limitations documented in each report and the scientific
contract; software checks and empirical assessment are different evidence.

- `milestone1.qmd` / `.html`: accepted October 15, 2026. Historical execution
  requires its recorded commit; the current engine requires a union process.
- `milestone2.qmd` / `.html`: accepted October 19, 2026; union formation.
- `milestone3.qmd` / `.html`: accepted October 27, 2026; groups, matching,
  ancestry and inheritance.
- `milestone4.qmd` / `.html`: accepted October 29, 2026; dissolution,
  age/duration remarriage, cohabitation, conversion, matching roles and rate
  observation.
- `data/`: local rate inputs, source documents, metadata and preparation scripts
  used by the reports and shipped examples. Retain these for reproducibility.
  See the [data directory guide](data/README.md) for country/period folders and
  the path mapping for existing reports, which were not edited during relocation.
- [England and Wales, 1992 onward](data/england_wales_1992_onward/README.md):
  fertility by three union types, annual marriage evidence and population
  denominators, cohabitation composition, and duration-specific divorce
  probabilities for a proposed second-transition demonstration. Observations
  and British supplementary evidence are retained separately from unselected
  scenario assumptions; includes an offline rebuild and source manifest.
- `milestone5.qmd` / `.html`: accepted October 29, 2026 (D105); parity-specific
  fertility, cell-accounted birth spacing, lognormal fertility weights with
  maternal inheritance, and baseline/restart/union-type recovery checks. The
  final render was active at closeout; preserve all active intermediates.
- `milestone6.qmd` / `.html`: accepted October 30, 2026 (D116); individual
  age/duration group transitions, sex/union-type strata, demographic rate and
  formation-group recovery, population shares and performance. The final render
  is active at acceptance; preserve all intermediates. Actual replication and
  population settings are in the report, not inferred from earlier drafts.
- `milestone5/`: supporting method investigations and diagnostics, including
  superseded-correction evidence, runtime cell-accounting checks, reference
  SOCSIM experiments and version-specific reproduction instructions.
- `milestone4/`: supplementary calibration benchmarks and conditional-choice
  recovery evidence cited by the current methods and testing documentation.

Render the final milestone 4 report from this directory with
`quarto render milestone4.qmd`. It downloads HMD mortality using configured
access, runs 100 independent starting populations per scenario, and executes
the offline test suite. It is a substantial analysis, not an ordinary package
check. Do not run it merely to verify documentation edits. Milestone 5 renders
with `quarto render milestone5.qmd`; it installs the package, downloads HFD/HMD
inputs, runs 100 independently initialized replications per scenario and the
offline tests. Its assessment uses 150-year runs and a 50-year burn-in, with
additional restart follow-up. Reproduction settings remain in the report.

Quarto intermediates (`*.rmarkdown`, `*.knit.md`, `*_files/`, `.quarto/` and
`.cache/`) are disposable after rendering finishes. Never remove them during
an active render. `milestone4_tests.R` was an exploratory script superseded by
the final report and removed by the author; Git history retains it.
