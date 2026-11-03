# Current-hazard contributions for matching calibration

User-approved implementation and performance trial, October 28, 2026. Replace
age-band-averaged female intensities with actual current hazards and an internal
boundary-update queue. D083 now extends the same mechanism to duration-specific
remarriage without introducing another calibration dimension.

## Scientific behavior

Each living unpartnered woman contributes her current female search hazard to
her group/birth-cohort calibration cell. Rate-mode marriage searches select the
woman's group/previous-union schedule on its selected age/duration clock after
proportional multipliers. The search origin is independent of fertility age.
Fertility-driven searches use each unpartnered eligible woman's cell-derived
hazard (D101: cell intensity divided by eligible members), formed from the
fertility cells at search time, and contribute zero during birth spacing. Keep separate totals
for marriage searches and fertility-driven partner/father searches.

Sum individual contributions into the existing group/cohort cells. Derive male
mirrored targets and expected flows from these same totals. Male-group target
normalization, delta, matching utilities, kinship conventions, and the solver
are unchanged. Quarter-year cohort bins still approximate matching utilities
and cutoff eligibility for calibration; their center-based inclusion rule is
retained. Actual female hazards, rather than bin-averaged hazards, are the
intentional approximation change. Individual event scheduling integrates actual schedules on the selected
age/duration search clock; D082 bookkeeping does not redraw demographic events. Different calibrated weights
can change partner identities and downstream histories; identical seeds need
not reproduce previous grouped outcomes. One/no-group paths and paths with
all exogamy prohibited do not create the tracker or consume additional draws.

## Updates and event timing

Maintain the current two hazard contributions per counted woman and at most
one pending bookkeeping entry for her next relevant rate change. Boundaries
include search age/duration and fertility age band starts and ends, gaps and terminal support, and
birth-spacing release. Rate schedules are piecewise constant, so no update is
needed between changes. This is not a simulation step or a demographic event.
The queue advances through every boundary up to a demographic event/search time
before calibration is used. At an exact boundary the new half-open band applies.
No threshold is drawn and no demographic event is rescheduled by a bookkeeping
update. Use absolute calendar boundaries to avoid subtraction roundoff at ties.

Entering the unpartnered pool inserts contributions and a boundary entry;
forming a union or dying removes both. Dissolution/widowhood inserts the newly
selected previous-union schedule. Birth updates fertility eligibility while
retaining marriage-search intensity. Replace old pending entries when state
changes, keeping queue size bounded by counted women rather than all past
changes. Restarts rebuild contributions and spacing from supplied histories.
The same update interface now accepts a separate search clock origin for D083
duration schedules. It provides a common update point for future state-dependent
rates, without adding parity, region or custom state transitions.

For example, a woman aged 19.9 with search hazards .1 below age 20 and .3
thereafter contributes .1 now and .3 after 0.1 years, even if her cohort bin
straddles age 20. A birth with spacing .75 temporarily sets her fertility-driven
contribution to zero; her marriage hazard continues to follow its schedule.
At spacing release use the fertility rate at the attained age, including any
age boundaries crossed during exclusion.

## Validation and comparison

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

## Performance trial results

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


## Matching association recovery

A subsequent 20-replication recovery check with A:B utility -3 recovered
conditional cross-group log odds ratios -5.97493 (MCSE .05378) for group-only
matching and -5.98931 (MCSE .04935) for age-dependent matching, against target
-6. Age coefficients also recovered their inputs. Both cases included
previous-status search schedules, dissolution and remarriage. Choice-set models
used independently reconstructed available pools and recorded weight offsets.
See [testing evidence](testing.md#matching-parameter-recovery-with-the-new-tracker-october-28-2026)
for inputs, methods, raw results and the distinction from pooled cross-tab odds.
