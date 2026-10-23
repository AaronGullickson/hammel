# Fixed-pool versus sequential two-sided matching: timings

These measurements precede the native bulk optimization. See
[native performance](native_two_sided_performance.md) for current measurements.

Measured October 22, 2026 on aarch64 Linux, R 4.6.1, GCC 16.1.1. This is
software performance evidence; recovery of demographic input parameters
remains under investigation. The milestone 3 report was not executed or edited.

The new fixed-pool method was slower than sequential depletion in these
workloads. Every pair must contribute to column probabilities; the sequential
method scores only selected remaining rows/columns and stops when ego matches.
Removing the partner-draw loop therefore does not remove the dominant score
work. Both are substantially more expensive than default female choice.

## Complete simulations

Three trials per method/sample size, medians in seconds. Each starts with the
same seeded 2,000-person population, ages 0–80, 75% A / 25% B. Mortality is
0.01/year, fertility 0.08/year at ages 15–45 with default spacing, duration
20 years, on-demand formation, preferred gap zero, half-weight distance five,
and symmetric cross-group effect −2.5. Sample size precedes eligibility
screening. Histories/event counts differ because algorithms consume different
random streams; raw CSVs retain births, deaths, unions, and final history sizes.

| Sample size | Female choice | Sequential two-sided | Fixed-pool two-sided | Fixed / sequential |
| --- | ---: | ---: | ---: | ---: |
| 50 | 0.016 | 0.019 | 0.024 | 1.26 |
| 200 | 0.017 | 0.053 | 0.112 | 2.11 |
| 1000 | 0.025 | 0.662 | 1.685 | 2.55 |

These are synthetic timings, not predictions for the user's 10,000-person
rate-mode workload. Baseline durations of a few milliseconds have coarse
resolution. Timing comparisons do not assess scientific performance.

## Fixed markets

Milliseconds per fresh market, medians of three trials. M women and M men,
ego age 25, other ages uniform 15–65, independent 75% A / 25% B groups, same
preferences as above. Native production paths are compared. The fixed-pool
implementation streams female rows and retains column sums plus ego's row,
using linear scratch space rather than allocating the complete matrix.

| M | Sequential | Fixed pool |
| --- | ---: | ---: |
| 50 | 0.015 | 0.030 |
| 100 | 0.050 | 0.110 |
| 300 | 0.380 | 0.920 |
| 1000 | 3.850 | 10.350 |

The new script also compares streamed scores with explicitly materialized
matrices and verifies identical choices and RNG states. Raw measurements
include M=10, whose timings approach timer resolution.

## R formulas

Five-year simulations with the same population and three trials. A pointwise
equivalent `I(abs(age_difference)^2)` formula forces R evaluation. Each fixed
female row is evaluated once; the sequential method evaluates only each chosen
remaining row/column. Median seconds:

| Sample size | Sequential native | Fixed native | Sequential R | Fixed R |
| --- | ---: | ---: | ---: | ---: |
| 50 | 0.015 | 0.016 | 0.445 | 0.901 |
| 200 | 0.028 | 0.047 | 1.692 | 3.629 |

## Reproduction and provenance

Run `TZ=America/Los_Angeles Rscript validation/two_sided_matching_benchmark.R`
from the package root. It uses synthetic inputs only. The script now removes
existing compiled objects before `compile_dll(debug = FALSE, force = TRUE)`
and loads with `compile = FALSE`. Requesting `force = TRUE` alone can reuse
debug objects; earlier recorded timings did not establish a clean optimized
build and should not be compared directly with this assessment.

The sequential comparison used a separate temporary package tree with the
source from commit `03eba9b` (the previous experimental implementation), no
pre-existing compiled objects, and the same benchmark workload. The current
fixed-pool implementation is an uncommitted change on that commit. Three
trials alternate method order within each build; builds were measured
sequentially on the same machine, so machine-load differences remain possible.

Raw CSVs and R session/build metadata are under
[fixed_pool](two_sided_matching_benchmark/fixed_pool/simulations.csv) and
[sequential](two_sided_matching_benchmark/sequential/simulations.csv), with
`fixed_markets.csv`, `simulations.csv`, `callbacks.csv`, and `metadata.rds` in
each directory. Earlier raw files at the parent directory retain the original
prototype measurements; their optimized-build status is unverified.

Source MD5 hashes at measurement time:

- Fixed pool `src/simulate.cpp`: `443af8c24c73495159a78801b9edb921`
- Fixed pool `R/union_formation.R`: `1a9c5ef5c1a97272f8edafbee7e8cc8e`
- Fixed pool `R/simulate.R`: `cc513bfc4e0f3e98215a57672b5479fa`
- Sequential `src/simulate.cpp`: `ebccacbf934fe2d468bffb9e7b1fde85`

## Verification

The full `devtools::test()` suite passed; live database tests remained skipped.
Checks include analytical probabilities, an independent literal rejection
simulation with replacement, fixed-matrix normalization, extreme utilities,
exclusions, pointwise native/R agreement, row-context custom formulas, one
selection draw, state isolation, eligibility, shortages, father overrides,
restart/pool maintenance, and diagnostic invariance. Scientific odds-ratio
recovery is a separate milestone 3 assessment.
