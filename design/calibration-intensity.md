# Current-hazard contributions for matching calibration

Updated November 5, 2026. D082 (October 28, 2026) replaced age-band-averaged
female intensities with actual current hazards; D083 extended it to
duration-specific remarriage. Since milestone 8 (D128) the contributions come
from the search and fertility cells themselves, rather than from a separate
per-woman tracker with its own boundary queue (the former
`FemaleIntensityTracker`, removed November 4, 2026).

## Scientific behavior

Each living unpartnered woman contributes her current female hazards to her
group/birth-cohort calibration cell, in two channels kept separately:

- **Marriage searches (rate mode).** Her search hazard is
  `rate_c * N_c * w_i / W_c` in her search cell c (`rate_c` with weights of
  one), split by her current cohabitation probability p into a cohabitation
  share `p` and a direct-marriage share `1 - p`. The search origin, age or
  duration since the latest union ended, is independent of fertility age.
- **Fertility-driven partner and father searches.** Each unpartnered eligible
  woman's cell-derived fertility hazard (D101: cell intensity times her weight
  over the eligible weight sum), zero during birth spacing.

Sum contributions into the existing group/cohort cells. Derive male mirrored
targets and expected flows from these totals. Male-group target
normalization, delta, matching utilities, kinship conventions, and the solver
are unchanged. Quarter-year cohort bins still approximate matching utilities
and cutoff eligibility for calibration; their center-based inclusion rule is
retained. Calibrated weights can change partner identities and downstream
histories. One/no-group paths and paths with all exogamy prohibited do not
calibrate or consume additional draws.

## Updates and event timing

- **Search cells** keep, per calibration bin, their member count and the sums
  of `w (1 - p)` and `w p` (`SearchCalibration`). A cell whose contents
  changed is marked and recomputed before the next search solve, so the
  totals change only through cell membership, weights and rates.
- **Fertility cells** keep the weight of unpartnered eligible members per
  calibration bin and contribute `Lambda_c / W_c^e` per unit weight.
- **Cohabitation probabilities** are the one quantity that can change between
  cell moves, because their age bands can differ from the search clock. A small
  `ProbabilityTracker` keeps each unpartnered woman's current p and a
  bookkeeping boundary at her next probability band edge; it exists only with
  rate-mode search and cohabitation. Calendar probability changes refresh p
  for everyone before same-time searches.

Cell moves at band edges, spacing releases and probability boundaries are
bookkeeping: no threshold is drawn and no demographic event is rescheduled.
At an exact boundary the new half-open band applies. Restarts rebuild cells and
spacing from supplied histories.

For example, a woman aged 19.9 in a search cell [15,20) with rate .1 moves to
the [20,25) cell with rate .3 after 0.1 years; her calibration contribution
follows the cell, even if her cohort bin straddles age 20. A birth with spacing
.75 makes her ineligible for fertility, removing her fertility-driven
contribution until release; her marriage contribution continues.

## Validation

With option `hammel.check_cells = TRUE`, every search solve compares the
maintained marriage, cohabitation and fertility totals with an explicit scan
of searching women, recomputing each woman's probability from her schedule
(`tests/testthat/test-calibration-intensity.R`). Calibrated runs gave
identical histories before and after the tracker's replacement.

## D082 evidence (previous engine, October 28, 2026)

The checks and timings below were recorded for the original tracker. They
remain the evidence for D082's change from bin-averaged to actual hazards; the
engine they timed has since been replaced.

### Validation and comparison

Compare maintained cell totals against explicit individual scans at boundaries,
in gaps and after pool/state changes, births, spacing releases and restarts.
Check exact boundary classification, queue size, no RNG use, both intensity
channels, narrow bands, age-bin heterogeneity and independent simulations.
Keep the exact male draw and incremental male utility-sum tests. Fixed-seed
comparisons apply to unchanged ungrouped paths and equivalent constant-hazard
cases; broader grouped comparisons are descriptive because the approximation
changes. Use optimized builds of the pre-change commit and working tree with
identical stored baselines, schedules, years and seed sets. Report elapsed times,
workload counts, input provenance and limitations. No empirical rates are fetched
or scientific acceptance thresholds adopted.

### Performance trial results

Compared baseline commit `93efcc1a6c41c2eb2f23ff597edf17089bd36c47` with the
D082 working tree on October 28, 2026. Both were installed into separate
libraries with GCC 16.1.1, C++17 and `-O2` (no `-O0`), using R 4.6.1 on
Linux aarch64, Apple M1 Pro. This is an offline synthetic workload, not a rerun
of the user's empirical report: 50,000 initial persons, maximum initial age 50,
100 simulated years, 75/25 A/B groups in grouped cases, matching mean age gap
0/half-weight distance 5 with A:B effect -2, and dissolution enabled in every
case. Schedules and initial-population seed are fully specified in
[the benchmark script](../validation/milestone4/benchmark-calibration-intensity.R).

The same saved models were used for seeds 42, 71 and 105. Runs were sequential,
with version order old/new, new/old, old/new across those seeds; each process
warmed native code with a one-year ungrouped run before timed simulations.
Garbage collection was outside each timed interval. No builds or tests ran
concurrently with measured simulations. Timings are medians of three runs;
this small sample describes local performance, not a significance claim.

| Workload | Old seconds | New seconds | Change |
| --- | ---: | ---: | ---: |
| Ungrouped rate searches | 0.535 | 0.539 | +0.7% |
| Grouped shared age schedule | 4.454 | 4.555 | +2.3% |
| Grouped previous-union schedules | 4.748 | 4.364 | -8.1% |
| Grouped on-demand searches | 2.097 | 2.170 | +3.5% |

The initial all-old-then-all-new pass suggested a larger improvement, but its
unchanged ungrouped control was also faster. The alternating-version pass above
is the reported comparison to reduce ordering effects. The new approach has
small overhead in simpler cases and avoids repeated category/rate integration
in the previous-status case. Do not extrapolate these percentages to other
population sizes, schedule detail, matching formulas or machines.

Full per-run elapsed times, births, deaths, union episodes, dissolutions,
living A share and directional pairing ratios are retained in
[the benchmark CSV](../validation/milestone4/calibration-intensity-benchmark.csv). Grouped
histories differ because actual hazards replace bin-averaged calibration
hazards; workload counts therefore need not be identical. Three replications
are not a scientific validation of the resulting demographic distributions.
The ungrouped 100-year history at seed 42 was compared in full and is identical
between old/new builds, beyond merely matching event totals.

For the grouped previous-status case at seed 42, Linux process high-water RSS
(`/proc/self/status`, `VmHWM`) was 171.2 MiB old versus 183.3 MiB new: about
12.1 MiB additional peak memory. These include R, all saved input models and
returned histories, not just engine state. Counts differ slightly across the
grouped runs. The tracker adds O(persons) contribution storage and O(available
women) bounded queue storage; it does not add an age/status/duration tensor.

Reproduction after installing old/new sources into separate libraries:

```sh
Rscript validation/milestone4/benchmark-calibration-intensity.R make OLDLIB INPUT.rds
Rscript validation/milestone4/benchmark-calibration-intensity.R compare OLDLIB INPUT.rds NEWLIB OUTPUT.rds
Rscript validation/milestone4/benchmark-calibration-intensity.R memory OLDLIB INPUT.rds OLD_MEMORY.rds
Rscript validation/milestone4/benchmark-calibration-intensity.R memory NEWLIB INPUT.rds NEW_MEMORY.rds
```

The recorded input RDS MD5 was `1ec7e19db1b6baf7b40b7fb410ee57bb`; the new
`src/simulate.cpp` MD5 was `3140b5c733e304c8184b3ae7b2b497c3`.


### Matching association recovery

A subsequent 20-replication recovery check with A:B utility -3 recovered
conditional cross-group log odds ratios -5.97493 (MCSE .05378) for group-only
matching and -5.98931 (MCSE .04935) for age-dependent matching, against target
-6. Age coefficients also recovered their inputs. Both cases included
previous-status search schedules, dissolution and remarriage. Choice-set models
used independently reconstructed available pools and recorded weight offsets.
See [testing evidence](testing.md#matching-parameter-recovery-with-the-new-tracker-october-28-2026)
for inputs, methods, raw results and the distinction from pooled cross-tab odds.
