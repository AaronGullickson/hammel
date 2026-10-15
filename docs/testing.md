# Testing and milestone validation

Updated October 15, 2026. The user accepted milestone 1 as complete after
rendering the validation report. Software test history and the user's demographic
assessment are distinct evidence; recording acceptance did not rerun either.
Earlier implementation-specific runs and counts remain available in Git history.

## Run the suite

The package uses testthat edition 3 with testthat >= 3.2.0 for scoped mocks.
Tests live in `tests/testthat/`, with the standard runner `tests/testthat.R`.
From the repository root:

```r
devtools::test()
devtools::test(filter = "simulation")
```

Package checks use the standard runner. Tests isolate their case fixtures;
optional HMDHFDplus reader checks skip if it is unavailable, and mocked reader
bindings restore on exit. The live database test skips unless explicitly enabled.
The old top-level base-R test scripts have been converted and removed.

| Test file | Main evidence |
| --- | --- |
| `test-model-api.R` | Construction, replacement, schemas, history and metadata retention |
| `test-baseline-population.R` | Reproducibility, bounds, known cases, distributions, empty populations |
| `test-processes.R` | Rate validation, typed specifications, units, spacing, edited objects |
| `test-hazard-sampler.R` | Deterministic hazard inversion and analytical survival |
| `test-fertility-adjustment.R` | Annual/five-year midpoint overlaps, gaps, shared sampler, feasibility, open-age rejection, independent mortality |
| `test-simulation.R` | Chronology, spacing, renewal rate, accounting, parent availability, seed isolation, restarts, unions closing on death |
| `test-demography.R` | Exact exposure, event boundaries, missing/zero cells, TFR widths, analytical life tables, population counts, tibbles |
| `test-database-processes.R` | Synthetic/local conversions, metadata, RDS reuse without downloads |
| `test-database-online.R` | Opt-in authenticated download contracts |

Minimum-spacing checks group births by mother, sort dates, and check successive
intervals against the configured spacing with a small numerical tolerance.
They include default 0.75 years and nondefault 1.13 years across restarts.

Most recent recorded implementation verification: the development suite,
optimized package installation, installed-package runner, and documentation
examples passed for D028. The installed suite recorded 1,832 expectations,
zero failures/warnings, and one opt-in live-test skip. Those expectation counts
are an execution record, not a scientific acceptance criterion. R's timezone
service was unavailable in the sandbox, so verification set
`TZ=America/Los_Angeles` explicitly.

The development suite was rerun during this documentation audit and passed,
with the opt-in live test skipped. Synthetic construction, simulation, analysis,
and restart examples also ran successfully; no data downloads were made.

## Optional live database checks

Install HMDHFDplus and configure `HMD_USERNAME`, `HMD_PASSWORD`, `HFD_USERNAME`,
and `HFD_PASSWORD` in the R environment or `~/.Renviron`. Then run:

```sh
HAMMEL_RUN_ONLINE_TESTS=true Rscript -e 'devtools::test(filter = "database-online")'
```

Without opt-in, the test skips before inspecting credentials. With opt-in,
missing credentials/dependency, authentication errors, or invalid data fail it.
It checks USA/2000 rate preparation and provenance, not demographic simulation
accuracy. It does not print credentials or raw downloaded tables.

Authenticated USA/2000 access was verified October 12, 2026 against HMDHFDplus
2.0.8. Later midpoint verification used offline fixtures; that historical live
check is not a guarantee of present database availability. Retained missing
rates are rejected; neither tests nor helpers impute them.

## Reproducible demographic assessment

The [replication example](../inst/examples/milestone1_bootstrap.R) currently uses
100,000 initial persons, one shared baseline, and 100 simulation seeds for
100 years with USA/2000 schedules. It examines TFR, ASFR, life expectancy,
ASMR with 20-year periods, and population growth. This conditions on the same
starting population; it does not resample initial populations or observed data.
Its filename uses “bootstrap,” but the runs are simulation replications.

The completed [Quarto report](../validation/milestone1.qmd) contains the user's
analysis and interpretation. Execution is enabled; rendering reruns simulations
and input downloads and can require substantial memory. From the repository
root:

```sh
quarto render validation/milestone1.qmd
```

Rate downloads should be saved explicitly with RDS and loaded for later renders;
record source attributes, checksums, seeds/settings, R/package versions, and
commit plus any uncommitted code changes. CSV does not preserve provenance
attributes. Local validation data/rendered output are Git-ignored, and the
validation directory is excluded from package builds. Preserve completed inputs
and evidence in a documented archive when sharing the report.

## Accepted milestone evidence and remaining limits

The user rendered `validation/milestone1.html` and accepted milestone 1 on
October 15, 2026 (D029). The source report includes replicated ASFR/ASMR coverage
and deviations, TFR and life-expectancy comparisons, two-sex Leslie projection,
stable growth and age structure, extended runs, fertility-rate restart cases,
and performance/test results. Scientific conclusions belong to the report's
author; no new acceptance thresholds are imposed here.

At acceptance, the clean working tree was at commit
`dd02ac681ca2cb3c14f98b8b73a5b0f625637445`. The local HTML artifact had SHA-256
`66d3cb3ba35f86f44031fc1044abc594b76f967c259f1e7162ad948c46622974`.
This identifies the inspected artifact, not an independent replication. The
HTML and local validation data remain Git-ignored; they have not been deposited
in an external archive as part of completion.

Limitations remain: initial fertility transients, sparse oldest-age cells and
normal-interval reliability, and fertility compensation after rate changes.
The report defers correction of the restart issue beyond milestone 1. It does
not establish accuracy for all rate schedules or supply a controlled executable
SOCSIM comparison. Further studies and archive/publication work require their
own scope; these are not outstanding conditions on the user's acceptance.

## Milestone 2 construction checks (October 15, 2026)

After adding union-process construction and formula-based matching utilities,
`devtools::test()` passed with 1,870 expectations, zero failures/warnings,
and the opt-in live database test skipped. No rates were downloaded. The new
checks cover analytical softmax probabilities, quadratic utilities, numerical
stability, formula/coefficient validation, process attachment/replacement,
provenance, edited-specification rejection, and the engine's explicit rejection
of union processes. Documentation generation and `git diff --check` passed.
These checks validate construction/scoring infrastructure, not executed union
formation or milestone 2 scientific completion.
