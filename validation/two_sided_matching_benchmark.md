# Sequential two-sided matching: synthetic timing record

Measured October 22, 2026 on aarch64 Linux with R 4.6.1. This is a software
performance experiment, not scientific validation or acceptance of input
parameter benchmarks. It did not execute or alter milestone3.qmd.

Run [the benchmark script](two_sided_matching_benchmark.R) from the package
root. It uses only synthetic/local inputs. Raw measurements are in
[fixed_markets.csv](two_sided_matching_benchmark/fixed_markets.csv),
[simulations.csv](two_sided_matching_benchmark/simulations.csv), and
[callbacks.csv](two_sided_matching_benchmark/callbacks.csv).

## Build and provenance

Base commit: `4d03f55d814344fccff258a0f3f3912a88711db0`. Experimental changes were uncommitted.
Native code was rebuilt with `pkgbuild::compile_dll(debug = FALSE, force = TRUE)`
using the installed R compiler flags. Methods in each comparison used the
same build. The compiler was GCC 16.1.1. Reversing method order between trials
reduces simple order effects; three trials do not provide a broad hardware
or workload assessment. Millisecond-scale baseline timings are coarse.

MD5 hashes at measurement time (later documentation-only edits may differ):

- `src/simulate.cpp`: `ebccacbf934fe2d468bffb9e7b1fde85`
- `R/union_formation.R`: `d64fb2dc269c800a6e5c7c07e8e96291`
- `R/simulate.R`: `10f6f14670ee97b3910e97be60cac168`

## Native on-demand scoring versus a full precomputed matrix

The market has M women and M men, with ego aged 25, other ages uniform on
15–65, independently sampled 75% A/25% B identities, preferred gap zero,
half-weight distance five, and symmetric cross-group effect -2.5. For each
size and trial the same RNG seed is used for both scoring strategies. The
full matrix is rebuilt per fresh market, not amortized over repeated draws.
All recorded choices, choosing sexes, steps, and RNG-dependent results agree
exactly between strategies. The simulation uses the on-demand strategy.

Median milliseconds per fresh market:

| M | On demand | Precomputed | Precomputed / on demand |
| --- | ---: | ---: | ---: |
| 10 | 0.005 | 0.015 | 3.00 |
| 50 | 0.080 | 0.180 | 2.25 |
| 100 | 0.270 | 0.675 | 2.50 |
| 300 | 2.160 | 5.640 | 2.61 |
| 1000 | 23.100 | 62.200 | 2.69 |

## Complete native-scored simulations

Each run starts with the same seeded 2,000-person population, ages 0–80,
75% A/25% B. Mortality is 0.01/year for both sexes; fertility is 0.08/year
at ages 15–45 with default spacing. Runs last 20 years and use on-demand
formation with the same age/group preference model as above. Three seeds
are used for each method and sample size. Matching methods consume different
random streams and produce slightly different event counts; timings are not
comparisons of identical histories. Candidate sample size precedes screening,
so actual temporary markets can be smaller than the requested sample.

Median wall-clock seconds and event counts:

| Sample size | Female choice s | Two-sided s | Ratio | Female choice unions | Two-sided unions |
| --- | ---: | ---: | ---: | ---: | ---: |
| 50 | 0.020 | 0.037 | 1.85 | 346 | 344 |
| 200 | 0.029 | 0.232 | 8.00 | 336 | 342 |
| 1000 | 0.060 | 3.708 | 61.80 | 349 | 342 |

Final history sizes and births/deaths are retained in the raw CSV. These
measurements demonstrate substantial sample-size-dependent overhead; they
are not runtime predictions for the user's larger validation workloads.

## R formula overhead

A pointwise-equivalent R expression `I(abs(age_difference)^2)` forces the R
fallback while preserving the native curve. These two-sided runs use the
same population and synthetic inputs but last five years, with three seeds.
Both paths score only the chosen remaining row/column. Median seconds:

| Sample size | Native | R formula |
| --- | ---: | ---: |
| 50 | 0.023 | 0.465 |
| 200 | 0.095 | 1.861 |

## Software verification

The complete `devtools::test()` suite passed after implementation; live database
tests remained skipped. Additional checks established exact default seeded
history agreement with pre-experiment runs in both initiation modes. Tests
cover analytical two-by-two probabilities, an independent R reference, ego
eligibility, shortage reduction, group exclusions, selected-vector R callbacks,
actual-state isolation, father overrides, ageing-in, deaths, widowhood, restart
pool reconstruction, and diagnostics/RNG invariance. No inference about
scientific parameter recovery follows from these software checks.
