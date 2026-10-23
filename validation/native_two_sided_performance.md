# Native two-sided matching: bulk calculation performance

Measured October 22, 2026, aarch64 Linux, R 4.6.1, GCC 16.1.1. The additive
fixed-pool D062 selection rule is unchanged. These are software checks and
synthetic timing measurements, not scientific parameter-recovery validation.

## Implementation

Native linear/quadratic models now calculate whole matrix columns in reusable
plain C++ arrays. This removes per-row R allocations and per-cell Rcpp object
access, and uses a maximum followed by exponentiation/summation for each
column. Raw age-gap arithmetic and directed group offsets are retained. Only
linear scratch storage is needed; a dense matrix is not allocated. The final
additive directional probability calculation and single random draw are shared
with the reference path. Custom R formula evaluation remains row-based.

This still evaluates every pair: the work is quadratic. R matrix operations
also visit every entry internally. Their potential advantage over interpreted
R loops does not imply an advantage over a compiled bulk calculation.

## Full simulations

Three-trial median seconds. Both builds use the same seeded 2,000-person
population (uniform ages 0–80, 75% A / 25% B), mortality 0.01/year,
fertility 0.08/year at ages 15–45 with default spacing, 20-year on-demand
formation, preferred gap zero, half-weight distance five, cross-group effect
-2.5. Sampling precedes eligibility screening. Per-trial event counts and
history sizes are included in the raw data.

| Sample size | Previous native | Bulk native | Speedup |
| --- | ---: | ---: | ---: |
| 50 | 0.024 | 0.019 | 1.26× |
| 200 | 0.114 | 0.063 | 1.81× |
| 1000 | 1.697 | 0.883 | 1.92× |

These workloads do not predict the user's larger rate-mode timings. Small
runtimes have coarse resolution; only three trials were run. Before/after
builds were measured sequentially, so changing machine load remains possible.

## Explicit R matrix calculation

[The R comparison script](native_two_sided_matrix_benchmark.R) constructs age
gaps with `outer`, applies the quadratic coefficients and group effects to the
matrix, normalizes its columns and ego row, combines probabilities, and draws
once. It uses stable column shifts. Markets have uniform ages 15–65, ego age
25, independently sampled 75/25 groups, and the preferences above. The native
and R methods recompute the market each draw, alternate order over three
trials, and produce identical sampled choices using the same random seeds.

Median milliseconds per fresh market:

| Market size | Bulk native | Vectorized R matrix |
| --- | ---: | ---: |
| 50 | 0.01 | 0.16 |
| 250 | 0.34 | 2.48 |
| 1000 | 5.15 | 55.70 |

The native method outperformed this R implementation. This does not establish
that every possible dense matrix/BLAS implementation would be slower.

## Reproduction and provenance

Run `TZ=America/Los_Angeles Rscript validation/two_sided_matching_benchmark.R`
for the full current synthetic benchmark, and
`TZ=America/Los_Angeles Rscript validation/native_two_sided_matrix_benchmark.R`
for the explicit R comparison. Both clean old objects, rebuild optimized native
code, and load with compilation disabled. No rate downloads or report execution
are involved. The prior implementation is in commit `bc18a8f`; the before
benchmark used a separate temporary source tree with that native implementation
and the same benchmark script. Source hashes verify the comparison despite
that temporary tree not containing Git metadata.

Raw CSVs and `metadata.rds` files with session information and source hashes
are under [before](native_two_sided_benchmark/before/simulations.csv),
[after](native_two_sided_benchmark/after/simulations.csv), and
[matrix](native_two_sided_benchmark/matrix/matrix_comparison.csv).

Native source MD5:

- Before: `443af8c24c73495159a78801b9edb921`
- After: `a49ff7fbcebe64a96ac7c31ccc287ec5`

## Verification

The reference materialized matrix path remains available internally. Tests
compare it with bulk columns across directed group offsets, hard restrictions,
linear/quadratic scores, extreme utilities, rectangular diagnostic markets,
and RNG states. Native and equivalent custom formulas agree. Six complete
fixed-seed histories (both initiation modes, sample sizes 50, 200, 1000) match
exactly before and after optimization. Column sums can differ at floating-point
roundoff, so universal bitwise probability equality is not claimed.
