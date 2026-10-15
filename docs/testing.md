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

## Temporary on-demand execution checks (October 15, 2026)

After implementing D031, `devtools::test()` passed with 1,918 expectations,
zero test failures/warnings, and the opt-in live database test skipped. New
cases check sample sizes, uniqueness and inclusion frequencies, current-partner
reuse, birth-time union creation/parentage, empty pools, historical IDs/custom
columns, widowhood/repartnering, restarts, unsupported active histories,
explicitly skipped screening and utilities, RNG isolation, and independence
from later milestone 1 runs. Rcpp bindings and Rd documentation were regenerated;
`git diff --check` passed. No rates were downloaded. The debug C++ build emitted
the toolchain's `_FORTIFY_SOURCE` warning about compilation without optimization;
compilation and installation succeeded. These are software checks of the
temporary uniform implementation, not milestone 2 scientific acceptance.

After D032 made union formation mandatory, the suite passed with 1,920
expectations, zero test failures/warnings, and one opt-in live-test skip. The
missing-process checks cover positive and zero duration. Simulation fixtures
now supply an explicit union process; restart cases retain it. Current examples
were updated, while the accepted milestone 1 report remains unchanged and
should be reproduced at its recorded commit. The native enable flag and
legacy no-union father-assignment branch were removed.

## Candidate-sampling performance (October 15, 2026)

A synthetic 50,000-person, 100-year run isolated candidate sampling overhead:
with constant mortality 0.02 for both sexes and fertility 0.06 at ages 15–50,
seed 42 for baseline and simulation, the hash-map sampler took 6.026 seconds
with sample size 1000 versus 0.152 seconds with sample size 1. Replacing its
per-candidate hash nodes with a reusable identity-position array and reversible
partial Fisher–Yates swaps reduced sample-size-1000 elapsed time to 0.962 seconds
(sample size 1: 0.163 seconds). Both measurements used the development/debug
build on the same machine; they are single-run diagnostics, not general timing
guarantees or measurements of the user's USA case. Sampling size, sampling
without replacement, uniform choice, and RNG draw order are preserved. Scratch
storage grows with the largest partial-search male pool and remains call-local.

## Sex-specific union schema checks (October 15, 2026)

After D033 renamed union roles to `wife_id`/`husband_id`, `devtools::test()`
passed 1,920 expectations with no test failures/warnings and one live-test skip.
An added role-validation test was then checked with
`devtools::test(filter = "model-api")`: all 106 expectations passed, including
five new checks for valid roles, swapped roles in ongoing/ended histories, and
the empty-table schema. Roxygen documentation and `git diff --check` passed.
The user's edits to `inst/examples/milestone2.R` were preserved.

After adding sex-specific union death reasons, `devtools::test()` passed all
1,941 expectations with zero test failures/warnings and one opt-in live-test
skip. Cases cover wife/husband deaths in both supplied and newly formed
unions, exact end times, and preservation of ended records on restart.
`git diff --check` passed.

## Conditional-logit execution checks (October 15, 2026)

`devtools::test()` passed 1,969 expectations with zero test failures/warnings
and one opt-in live database skip after D034. New tests cover the requested
0.072/-0.014 quadratic peak, analytical utilities/probabilities, 50,000 choice
draws within six marginal standard errors, native/R formula agreement,
coefficient ordering, custom transformations, extreme/nonfinite utilities,
empty/singleton pools, and birth-time selection in both age-gap directions.
Existing partner retention, history, widowhood, restart, and RNG checks pass.
This is algorithm validation; model recovery on demographic simulation data
and milestone 2 scientific acceptance remain pending.

The same synthetic 50,000-person/100-year case used in the sampling performance
record took 2.374 seconds with a sample of 1000 and the requested quadratic
coefficients (19,777 unions). This is a single development/debug-build timing,
not the user's USA case or a performance guarantee. No rates were downloaded.
Rcpp bindings and Rd documentation were regenerated.

## Deferred matching performance comparison

At the user's request on October 15, 2026, defer benchmarking the R formula
fallback against the native linear/raw-quadratic matching path. A future check
should include a cubic term such as `I(age_difference^3)`, controlling inputs,
candidate sample size, run duration, and build settings. Distinguish R
model-matrix overhead from changes in model complexity and simulated population
trajectory. No fallback slowdown has been measured; this note does not
authorize running the benchmark now.

## Minimum male age checks (October 15, 2026)

After D035, `devtools::test()` passed 2,007 expectations with zero test
failures/warnings and one opt-in live-test skip. Checks cover the default 15,
fractional/zero thresholds, invalid and edited specifications, exact threshold
entry, aging during a run, sampling before screening without refilling, births
with missing fathers when all sampled men fail, retaining supplied underage
partners, and no added minimum age for female initiation. Existing
conditional-logit, widowhood, restart, and demographic tests pass.
No rates were downloaded or deferred formula benchmarks run. Rcpp/Rd
generation and `git diff --check` passed.

## Shared kinship checks (October 15, 2026)

The D036 development suite initially passed 2,152 expectations with zero test
failures/warnings and one opt-in live database skip. Subsequent extraction
checks add base-data-frame input and half-first-cousin evidence; the targeted
`devtools::test(filter = "kinship")` run passed all 149 expectations with no
failures or warnings. Hand-built
families cover direct ancestry, half/full siblings, aunts/uncles, first and
removed/second cousins, missing maternal/paternal links, deceased ancestors,
unrelated co-parents, shuffled/noncontiguous IDs, repeated ancestry, shortest
lineage, empty schemas, invalid settings/IDs and ancestry cycles. A 5,000-person
chain verifies iterative deep traversal. Pair-by-pair integration checks agree
between R extraction and native partner screening. Newborn/restart cases
verify sons cannot be selected by their mothers and recorded parent links
are rebuilt consistently. Existing union histories are retained.

A synthetic 50,000-person/100-year run (the existing constant-rate case,
quadratic matching, sample 1000, minimum male age 15, seed 42) took 3.626 seconds
with generation bound zero and 4.197 seconds with bound two. These are
single-run debug-build diagnostics with different resulting stochastic
histories (33,656 versus 33,554 births), not isolated profiling or the user's
USA case. No rates were downloaded and the deferred R formula-fallback
benchmark was not run. Rcpp/Rd regeneration and `git diff --check` passed.

### Kinship screening optimization (October 15, 2026)

After the user reported 8.92 seconds for their case, temporary native timing
in the synthetic benchmark found roughly 0.035 seconds building ancestry
caches versus 1.88 seconds in candidate screening. The repeated-comparison
path was replaced with one marked ancestor set per woman and direct candidate
membership checks. An additional zero-bound short circuit skips kinship
queries/cache building when no cross-person relationship can be excluded.
Cycle validation still runs. Temporary timing code/output was removed.

The final synthetic benchmark took 2.608 seconds at zero and 3.481 at two
generations, versus 3.859 and 4.196 in the instrumented original. These are
single-run diagnostics; instrumentation and different simulated histories limit
exact timing comparisons. The user subsequently reported approximately 4.4
seconds at zero after reloading, close to their earlier 4.2-second baseline.

Before the final zero-bound short circuit, the full `devtools::test()` suite
passed all 2,160 expectations with no test failures/warnings and one live-test
skip, including a regression for switching between women without stale
ancestor markers. A saved 3,000-person, 100-year fixed-seed simulation was
identical before/after the membership optimization. No scientific rule or
random-number consumption was changed. The user's example edits were retained.
After the zero-bound short circuit, all 174 targeted kinship/eligibility
expectations passed, and the saved fixed-seed result remained identical.


### Exclusive unions (D037; October 15, 2026)

`devtools::test()` passed 2195 expectations with no failures. One live database
test was skipped as intended; one local-reader dependency warning came from
`lubridate` attempting unavailable sandbox `timedatectl` timezone discovery.
New checks cover partnered men excluded before singleton sampling, male scarcity
with preserved unknown-father births, widower re-entry, duplicate ongoing male
union rejection, and restart consistency. Existing death and sampling tests pass.
