# Validation reports and supporting evidence

Milestones 1–4 are accepted by the user (D029, D050, D077, D093). Keep the
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
  observation. The author is producing the final render at closeout.
- `data/`: local rate inputs, source documents, metadata and preparation scripts
  used by the reports and shipped examples. Retain these for reproducibility.
- `milestone4/`: supplementary calibration benchmarks and conditional-choice
  recovery evidence cited by the current methods and testing documentation.

Render the final milestone 4 report from this directory with
`quarto render milestone4.qmd`. It downloads HMD mortality using configured
access, runs 100 independent starting populations per scenario, and executes
the offline test suite. It is a substantial analysis, not an ordinary package
check. Do not run it merely to verify documentation edits.

Quarto intermediates (`*.rmarkdown`, `*.knit.md`, `*_files/`, `.quarto/` and
`.cache/`) are disposable after rendering finishes. Never remove them during
an active render. `milestone4_tests.R` was an exploratory script superseded by
the final report and removed by the author; Git history retains it.
