# Running tests

The suite uses testthat edition 3. Tests live in `tests/testthat/test-*.R`,
with shared helpers in `tests/testthat/helper-*.R` and the standard package
runner in `tests/testthat.R`. Run all tests from the project root with:

```r
devtools::test()
```

Run a subset with, for example, `devtools::test(filter = "simulation")`.
Package checks run the suite through the standard testthat runner. Analytical
thresholds, seeds, and demographic expectations are retained from the original
base-R scripts. Each named test case has local fixtures; reader replacements
use testthat's scoped mocks and restore the original bindings on exit.
Optional HMDHFDplus reader tests explicitly skip when it is unavailable.

## Simulation foundation

`test-demography.R` checks exact fractional-year exposure, event boundary rules, maternal-age birth assignment, rejection of duplicated IDs in stacked histories, unequal-width TFR, empty cells, and analytical constant-hazard life tables, including missing and zero terminal rates. It also checks stable small- and high-rate calculations. These test extraction and estimation algorithms, separately from whether the simulator reproduces its configured demographic rates.

Verified October 12, 2026: all offline test scripts and documentation examples passed after adding extraction. A synthetic 50,000-person, 50-year simulation also completed both rate-extraction paths, TFR, and period- and sex-specific life tables in approximately 0.14 seconds of analysis. Period death totals matched the records in the observation window. This is a workflow and accounting check, not a formal demographic validation over replications.

`test-hazard-sampler.R` checks known cumulative-hazard thresholds and analytical survival across age bands. `test-simulation.R` checks constant-hazard survival, minimum spacing, adjusted fertility renewal behavior after burn-in, exact continuous mortality dates, within-month chronological competition, half-open endpoints, accounting, ID and history preservation, unknown fathers, union closure, reproducibility, and restart behavior. These are initial foundation checks, not a completed assessment of age-specific rates, TFR, or full life tables.

The live example `inst/examples/milestone1.R` downloads USA/2000 schedules and runs a population of 50,000 for 100 years. It prints aggregate counts and elapsed simulation time; it does not write credentials or downloaded tables into the repository.

Verified October 12, 2026: the live example completed from 2000 to 2100 with 75,698 births, 67,291 deaths, and 58,407 living persons at the endpoint. Simulation elapsed time was approximately 0.23 seconds on this machine, excluding downloads. Counts satisfy `50000 + 75698 - 67291 = 58407`; this is not a completed ASFR/TFR or life-table validation. Hazard sampler, simulation, existing model/process/population tests, documentation examples, and Rd checks passed.

The live run also revealed that HMDHFDplus 2.0.8 parses HMD web ages despite `fixup = FALSE`. The adapter now respects `OpenInterval` as well as raw `+` labels; local-reader regression tests and the opt-in online checks cover the open terminal mortality interval.

## Live database integration

`tests/testthat/test-database-online.R` is opt-in. It downloads USA mortality and fertility rates for 2000, passes the prepared rate tibbles to the ordinary setters, and checks the imported process structure, units, sex schedules, rates, and source metadata. It does not simulate a population or establish demographic accuracy. It never prints credentials or raw database tables.

Install HMDHFDplus and supply `HMD_USERNAME`, `HMD_PASSWORD`, `HFD_USERNAME`, and `HFD_PASSWORD` through your environment or `~/.Renviron`. R reads the user environment file at startup (unless overridden by its startup configuration). Then run:

```sh
HAMMEL_RUN_ONLINE_TESTS=true Rscript -e 'devtools::test(filter = "database-online")'
```

With no opt-in flag, the test skips before inspecting credentials or loading the optional dependency. With the flag enabled, missing credentials, authentication failures, unavailable data, and conversion errors fail the test. Live data can change; the test checks contracts rather than exact rate values. Offline tests use synthetic fixtures and temporary reader replacements to check exact conversions independently of database availability.

Verified October 12, 2026: authenticated USA/2000 downloads through both setters passed with HMDHFDplus 2.0.8. This verifies the current access and conversion paths, not future database availability or simulation correctness.

Continuous-time update verified October 12, 2026 on `continuous-execution`: all offline scripts passed after installation in a temporary library. The online script skipped without opt-in; live downloads were not rerun. New checks cover sampled continuous mortality times, chronological birth-before-death competition inside a former monthly bin, exact endpoint exclusion, fractional run durations, father availability, and the continuous renewal expectation with retained monthly compensation. Historical live-run counts above describe the former discrete engine.

D025 verified October 12, 2026: all offline scripts passed with the new download helpers. Synthetic reader tests check RDS preservation of rates and source metadata for both databases and reuse in multiple models with readers replaced by functions that fail on any download. Live requests were not rerun.

Download once and persist outside the repository as appropriate:

```r
saveRDS(get_mortality_hmd("USA", 2000), "mortality-USA-2000.rds")
saveRDS(get_fertility_hfd("USA", 2000), "fertility-USA-2000.rds")

model <- hammel_model() |>
  set_mortality(readRDS("mortality-USA-2000.rds")) |>
  set_fertility(readRDS("fertility-USA-2000.rds"))
```

RDS retains provenance; CSV does not preserve the `source` attribute. The former `mortality_hmd()`, `fertility_hfd()`, and database setters have been replaced by the download helpers.

D026 verified October 13, 2026: package installation, all seven offline test
scripts, and generated documentation examples passed. The online test skipped
without opt-in. Added checks cover fertility process spacing validation and
18-month spacing through simulation and restart. Default renewal and zero-spacing
chronological checks still pass. No live database requests were made.

Testthat conversion verified October 13, 2026: all eight legacy test scripts
were migrated to named test cases and testthat expectations, and removed from
the top-level tests directory. Live tests retain explicit opt-in and fail on
missing credentials or dependency when requested. No live downloads were run.

Both `devtools::test()` and the installed-package `tests/testthat.R` runner
passed: 1,791 expectations, zero failures or warnings, and one opt-in live-test
skip across 52 named cases. Verification set `TZ=America/Los_Angeles` to avoid
the sandbox's unavailable system timezone service. Testthat >= 3.2.0 is required
for the scoped reader-mocking API.
