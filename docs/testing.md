# Running tests

The scripts in `tests/` use base R and an installed Hammel package. Run them with `Rscript tests/<name>.R` after installing the package.

## Simulation foundation

`demography.R` checks exact fractional-year exposure, event boundary rules, maternal-age birth assignment, rejection of duplicated IDs in stacked histories, unequal-width TFR, empty cells, and analytical constant-hazard life tables, including missing and zero terminal rates. It also checks stable small- and high-rate calculations. These test extraction and estimation algorithms, separately from whether the simulator reproduces its configured demographic rates.

Verified October 12, 2026: all offline test scripts and documentation examples passed after adding extraction. A synthetic 50,000-person, 50-year simulation also completed both rate-extraction paths, TFR, and period- and sex-specific life tables in approximately 0.14 seconds of analysis. Period death totals matched the records in the observation window. This is a workflow and accounting check, not a formal demographic validation over replications.

`hazard-sampler.R` checks known cumulative-hazard thresholds and analytical survival across age bands. `simulation.R` checks constant-hazard survival, minimum spacing, adjusted fertility renewal behavior after burn-in, exact continuous mortality dates, within-month chronological competition, half-open endpoints, accounting, ID and history preservation, unknown fathers, union closure, reproducibility, and restart behavior. These are initial foundation checks, not a completed assessment of age-specific rates, TFR, or full life tables.

The live example `inst/examples/milestone1.R` downloads USA/2000 schedules and runs a population of 50,000 for 100 years. It prints aggregate counts and elapsed simulation time; it does not write credentials or downloaded tables into the repository.

Verified October 12, 2026: the live example completed from 2000 to 2100 with 75,698 births, 67,291 deaths, and 58,407 living persons at the endpoint. Simulation elapsed time was approximately 0.23 seconds on this machine, excluding downloads. Counts satisfy `50000 + 75698 - 67291 = 58407`; this is not a completed ASFR/TFR or life-table validation. Hazard sampler, simulation, existing model/process/population tests, documentation examples, and Rd checks passed.

The live run also revealed that HMDHFDplus 2.0.8 parses HMD web ages despite `fixup = FALSE`. The adapter now respects `OpenInterval` as well as raw `+` labels; local-reader regression tests and the opt-in online checks cover the open terminal mortality interval.

## Live database integration

`tests/database-online.R` is opt-in. It downloads USA mortality and fertility rates for 2000, exercises both pipeable database setters, and checks the imported process structure, units, sex schedules, rates, and source metadata. It does not simulate a population or establish demographic accuracy. It never prints credentials or raw database tables.

Install HMDHFDplus and supply `HMD_USERNAME`, `HMD_PASSWORD`, `HFD_USERNAME`, and `HFD_PASSWORD` through your environment or `~/.Renviron`. R reads the user environment file at startup (unless overridden by its startup configuration). Then run:

```sh
HAMMEL_RUN_ONLINE_TESTS=true Rscript tests/database-online.R
```

With no opt-in flag, the test skips before inspecting credentials or loading the optional dependency. With the flag enabled, missing credentials, authentication failures, unavailable data, and conversion errors fail the test. Live data can change; the test checks contracts rather than exact rate values. Offline tests use synthetic fixtures and temporary reader replacements to check exact conversions independently of database availability.

Verified October 12, 2026: authenticated USA/2000 downloads through both setters passed with HMDHFDplus 2.0.8. This verifies the current access and conversion paths, not future database availability or simulation correctness.

Continuous-time update verified October 12, 2026 on `continuous-execution`: all offline scripts passed after installation in a temporary library. The online script skipped without opt-in; live downloads were not rerun. New checks cover sampled continuous mortality times, chronological birth-before-death competition inside a former monthly bin, exact endpoint exclusion, fractional run durations, father availability, and the continuous renewal expectation with retained monthly compensation. Historical live-run counts above describe the former discrete engine.
