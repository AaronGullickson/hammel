# Testing and milestone validation

Updated October 13, 2026. Software/analytical tests pass for midpoint fertility
compensation. Broader demographic validation is ongoing; acceptance criteria
and a completed milestone report are still needed. Earlier implementation-specific
runs and counts are available in Git history; they do not verify a later engine.

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

The [Quarto template](../validation/milestone1.qmd) is the development report.
Execution is disabled by default while the user drafts it. From the repository
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

## Remaining milestone checks

- Quantify replicated ASFR/ASMR differences and TFR/life-table agreement,
  distinguishing bias from variability and uncertainty in the mean.
- Align input/output age bands and report sparse cells, especially oldest ages;
  keep descriptive 100+ life-table comparisons distinct from aligned rates.
- Choose observation/burn-in windows and numerical acceptance criteria with the
  user. Initial excess fertility does not require an automatic correction.
- Broaden beyond USA/2000 to higher fertility, five-year fertility bands, and
  alternative spacing. An African historical case and new source integration
  are currently postponed, not selected.

User-reported visual agreement and analytical unit checks are useful evidence,
but do not establish full demographic accuracy or milestone completion. A
controlled executable SOCSIM comparison has not been performed; use a pinned
reference and statistical comparisons if that becomes appropriate.
