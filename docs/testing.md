# Testing and milestone validation

Updated October 19, 2026. The user accepted milestone 1 as complete after
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
| `test-groups.R` | Group rate lookup, analytical survival, coverage errors, fallback warnings, maternal inheritance, cross-group unions, restarts, RNG equivalence, and weighted initialization |
| `test-group-demography.R` | Exact group exposure/counts, cross-group parent and spouse attribution, conservation, analytical group TFR/life tables, population sizes, validation, and typed empty tables |
| `test-model-api.R` | Construction, replacement, schemas, history and metadata retention |
| `test-baseline-population.R` | Reproducibility, bounds, known cases, distributions, empty populations |
| `test-processes.R` | Rate validation, typed specifications, units, spacing, edited objects |
| `test-hazard-sampler.R` | Deterministic hazard inversion and analytical survival |
| `test-fertility-adjustment.R` | Annual/five-year midpoint overlaps, gaps, shared sampler, feasibility, open-age rejection, independent mortality |
| `test-simulation.R` | Chronology, spacing, renewal rate, accounting, parent availability, seed isolation, restarts, unions closing on death |
| `test-lexis-unions.R` | Exact long status exposure/counts, pre-death status, boundary births/deaths, formation rates, conservation, empty histories, TFR/life-table grouping, unchanged default |
| `test-demography.R` | Exact exposure, event boundaries, missing/zero cells, TFR widths, analytical life tables, population counts, tibbles |
| `test-database-processes.R` | Synthetic/local conversions, metadata, RDS reuse without downloads |
| `test-database-online.R` | Opt-in authenticated download contracts |

Minimum-spacing checks group births by mother, sort dates, and check successive
intervals against the configured spacing with a small numerical tolerance.
They include default 0.75 years and nondefault 1.13 years across restarts.

Historical suite counts and implementation-specific timings are available in
Git history. They are execution records, not scientific acceptance criteria.
Set `TZ=UTC` (or the local timezone) when sandbox timezone discovery is unavailable.

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

The current executable example is [milestone2.R](../inst/examples/milestone2.R).
The former milestone 1 single-run and replication scripts have been removed;
Git history retains them. Analyze simulation replications separately before
combining summaries; replication uncertainty does not resample the observed
input data or imply uncertainty over initial populations.

The accepted [milestone 1 report](../validation/milestone1.qmd) contains the
user's analysis and interpretation. Reproduce it from the acceptance commit
below, in a separate checkout: the current engine requires a union process,
while this historical report used the former no-union path. Rendering executes
simulations and input downloads and can require substantial memory. From that
checkout's repository root:

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

## Current milestone 2 software coverage

| Test files | Main evidence |
| --- | --- |
| `test-union-formation.R`, `test-conditional-logit.R`, `test-age-gap-matching.R` | Constructors, defaults, named coefficients, analytical utilities/probabilities, native/R formula agreement and age-gap helper |
| `test-on-demand-unions.R`, `test-union-eligibility.R`, `test-union-kinship.R` | Sampling, birth-time unions, age/kin screening, exclusivity, widowhood and restarts |
| `test-kinship.R` | Shared recorded-parent traversal, distance witnesses, incomplete ancestry and cycles |
| `test-rate-unions.R`, `test-status-fertility.R` | Search hazards, retries, status-specific fertility, spacing and event invalidation |
| `test-nonmarital-fathers.R`, `test-cancel-births.R` | Separate father matching, parentage without unions, failed-match cancellation and retry semantics |
| `test-union-diagnostics.R` | Eligible counts, utility ranks, tied utilities, historical retention and unchanged RNG/results |

These checks establish software behavior and analytical cases, not universal
demographic accuracy. Use `devtools::test()` for ordinary verification, with
live database checks disabled.

## Accepted milestone 2 demographic assessment

The user-authored [report](../validation/milestone2.qmd) uses USA circa
2000–2001 HMD mortality, NCHS overall/status-specific fertility, and SIPP-based
female marriage inputs with male rates as a reference. Input preparation and
limitations are recorded in [data provenance](../validation/data/usa_2000_rates.README.md).
England and Wales 2015 prepared data are an additional example, not the case
used in the current report.

The report now includes both initiation modes, conditional fertility and
formation rates, ever-married summaries, eligible-candidate diagnostics,
age-gap sensitivity, and conditional-logit recovery using 5,000 on-demand
unions and 50 reconstructed alternatives per union. The user reports recovery
of the generating linear/quadratic coefficients 0.072 and -0.014. Reconstruction
uses each union's formation time, age, survival and partnership eligibility;
it does not recover the actual sampled choice set or apply pair-specific kin
screening. See [matching validation](union-formation-models.md#rationale-and-matching-model-validation).

### Matching performance comparison

The report includes timing comparisons for 50,000 initial persons over 100
years, varying initiation, kin-screening depth and a custom formula
requiring R evaluation. The author's interpretation reports custom matching
as the largest timing cost. This comparison is no longer an unstarted task;
these are report findings, not an independently rerun benchmark or performance
guarantee. Formula and resulting population trajectories can affect timings.

The user accepted milestone 2 as complete on October 19, 2026 (D050), after
report review and revisions. Acceptance records the user's demographic assessment;
this documentation update does not independently rerun simulations or tests,
set acceptance thresholds, or claim every review suggestion was implemented.
The report and documentation are to be committed together; no separate final
artifact hash or acceptance commit is recorded here. Rendering executes an HMD
download and substantial simulations. Transition-aware compensation remains
deferred under D045; the report examines discrepancies, but no general
correction or guarantee of conditional ASFR recovery exists.
