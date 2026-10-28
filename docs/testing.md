# Testing and milestone validation

Updated October 27, 2026. Milestones 1–3 are accepted as complete (D029,
D050, D077). Software test history and the user's demographic assessment
are distinct evidence; recording acceptance does not independently rerun either.
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

On October 23, 2026, after D072 replaced candidate sampling with exact
full-pool draws and removed `sample_size`, `.n_eligible`, and `.utility_rank`,
`TZ=UTC Rscript -e 'devtools::test(reporter = "summary")'` passed with no
failures or warnings and the live database test skipped. Exact-draw checks
compare draw frequencies with exact full-pool conditional-logit probabilities
(native and custom formulas, group weights and restrictions, men near the age
cutoff, excluded kin, the enumeration fallback, empty pools). Tests that
previously checked sampling semantics now check that an eligible man anywhere
in the pool is found. Draw-by-draw identity is no longer expected between
native and custom formulas or across very steep utilities, which use random
numbers differently; those comparisons are distributional. This verifies
software behavior, not demographic parameter recovery.

| Test file | Main evidence |
| --- | --- |
| `test-groups.R` | Group rate lookup, analytical survival, coverage errors, fallback warnings, maternal inheritance, cross-group unions, restarts, RNG equivalence, and weighted initialization |
| `test-group-demography.R` | Exact group exposure/counts, cross-group parent and spouse attribution, conservation, analytical group TFR/life tables, population sizes, validation, and typed empty tables |
| `test-group-prop-hazard.R` | Multiplier validation, explicit-schedule equivalence, compensation order, omitted/unit/zero factors, warnings, RNG preservation, restarts, numeric codes, empty populations, and overflow |
| `test-group-inheritance.R` | All assignment rules/fallbacks, passing and ancestry probabilities, custom RNG, mixed identity rates, future-group coverage, restarts, and random-parent defaults and explicit maternal inheritance |
| `test-ancestry.R` | Starting identities, supplied proportions and IDs, parental averaging across generations, unknown shares, restarts/deaths, numeric labels, empty tables, and RNG independence |
| `test-group-matching.R` | Symmetric pair offsets and reverse-pair rejection, analytical probabilities, fixed factor contrasts/interactions, native/custom score agreement, father overrides, validation, and unchanged RNG |
| `test-prevent-exogamy.R` | Restrictions on the whole pool, compatible men always found, failure/cancellation, direction/symmetry, father overrides, retained partners, and unchanged disabled RNG |
| `test-detailed-unions.R` | Spouse-ID joins, repeated episodes, deaths, optional group types, empty tables, unchanged inputs, and name collisions |
| `test-weighted-choice.R` | Analytical weights and delta target adjustments, binned versus exact calibration, three groups, zero/unattainable targets, incremental sums versus full rebuild, mid-run identities, strict endogamy, unmet-target detection, narrow and one-year rate bands, recorded log weights, absorbed-formula warnings, unit weights for one/no groups, both modes, fathers, restarts, diagnostics, warning, validation |
| `test-parent-info.R` | Parent-ID joins, ages at birth, deceased/unknown parents, optional groups, empty tables, unchanged inputs, and name collisions |
| `test-model-api.R` | Construction, replacement, schemas, history and metadata retention |
| `test-baseline-population.R` | Reproducibility, bounds, known cases, distributions, empty populations |
| `test-ffab.R` | Female-birth probability validation, defaults, endpoint outcomes, RNG isolation, inheritance, independent runs, and restart settings |
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

On October 26, 2026, after D074 replaced imbalance absorption with named
`delta` target adjustments, `TZ=UTC Rscript -e 'devtools::test(reporter =
"summary")'` passed without failures or test warnings; the live database test
was skipped. Analytical checks cover the 75/25, beta -3 case with delta_B =
0.12274682836 (AB/BA = 3, log cross-product -6), signed and common-shift
adjustments, extreme finite values, three groups and group-stratified hazards.
Integration checks use nonzero deltas in both initiation modes, restarts and
separate father matching. These are software checks, not demographic acceptance.

After D075 removed method selection and made `group_exogamy` symmetric,
`TZ=UTC Rscript -e 'devtools::test(reporter = "summary")'` passed on October
26, 2026 with no failures or test warnings; the live database test was skipped.
Checks cover single-colon symmetric vector/table effects, rejected reverse
pairs and old syntax, native/R utilities, separate father delta, strict
endogamy/unit weights, both initiation modes, reproducibility and restarts.
The accepted milestone 3 report uses the current symmetric, weighted-choice
API and exact full-pool draws.

On October 27, 2026, after D076 changed the default inheritance to
`random_parent` with maternal fallback, `TZ=UTC Rscript -e
'devtools::test(reporter = "summary")'` passed without failures or test warnings;
the live database test was skipped. Checks cover default/explicit/missing-field
equivalence, missing-father fallback, ungrouped RNG equivalence, analytical
random-parent probabilities, and explicit maternal assignment.

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
input data. With a fixed baseline, MCSEs condition on that baseline;
independently generated baselines per replication include baseline variability.
The milestone 3 inheritance assessment uses the latter design.

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
attributes. The validation directory is excluded from package builds. Preserve
completed report sources, rendered artifacts and inputs; some local generated
files are Git-ignored. Use a documented archive when sharing the report.

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
This identifies the inspected artifact, not an independent replication. No
external archive was created as part of that acceptance; preserve the artifact
and its local input data when reproducing the historical assessment.

Limitations remain: initial fertility transients, sparse oldest-age cells and
normal-interval reliability, and fertility compensation after rate changes.
The report defers correction of the restart issue beyond milestone 1. It does
not establish accuracy for all rate schedules or supply a controlled executable
SOCSIM comparison. Further studies and archive/publication work require their
own scope; these are not outstanding conditions on the user's acceptance.

## Current milestone 2 software coverage

| Test files | Main evidence |
| --- | --- |
| `test-union-formation.R`, `test-conditional-logit.R`, `test-age-gap-matching.R` | Constructors, defaults, named coefficients, analytical utilities/probabilities, exact-draw frequencies against full-pool probabilities, native/R formula agreement and age-gap helper |
| `test-on-demand-unions.R`, `test-union-eligibility.R`, `test-union-kinship.R` | Birth-time unions, age/kin screening over the whole pool, exclusivity, widowhood and restarts |
| `test-kinship.R` | Shared recorded-parent traversal, distance witnesses, incomplete ancestry and cycles |
| `test-union-counts.R` | Exact union-year splitting, dissolution/death censoring, duration and period boundaries, repeat episodes, zero exposure, empty grids, group attribution, pooled conservation, endpoint and input validation |
| `test-union-dissolution.R` | Annual duration-hazard survival, one risk per union, input/restart duration, birth/death ties, widowhood invalidation, repeat episodes, availability, group factors/coverage, exposure accounting and seed behavior |
| `test-rate-unions.R`, `test-status-fertility.R` | Search hazards, retries, status-specific fertility, spacing and event invalidation |
| `test-nonmarital-fathers.R`, `test-cancel-births.R` | Separate father matching, parentage without unions, failed-match cancellation and retry semantics |
| `test-union-diagnostics.R` | Diagnostics validation; no candidate-set columns; unchanged histories |

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
it does not recover the actual choice set sampled at that time or apply pair-specific kin
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


## Accepted milestone 3 demographic assessment

The user accepted milestone 3 as complete October 27, 2026 (D077), after
reviewing the source and rendered report and making revisions. The accepted
[report](../validation/milestone3.qmd) includes:

- Group-stratified fertility, mortality and female union-search rates, and
  proportional-hazard alternatives including status-specific fertility factors.
- Weighted group matching with symmetric pair effects and delta, even/uneven
  group-size comparisons, and wife-side conditional-logit recovery with recorded
  log-weight offsets; demographic scenarios vary fertility, mortality and age gaps.
- Built-in inheritance and ancestry comparisons, independent starting populations
  in the inheritance replications, and custom assignment to an emergent AB identity.
- Runtime comparisons for ungrouped and two-, three-, and four-group populations,
  explicit schedules versus multipliers, and native versus custom inheritance.
- The testthat suite executed as part of rendering; live database tests require
  their separate explicit opt-in.

The report concludes that measured behavior meets the user's expectations,
while retaining limitations of small/imbalanced partner pools and the absence
of marriage-squeeze feedback. Other standing limitations, including conditional
fertility transition bias (D045), sparse exposure cells, calibration binning and
kinship omissions, remain. No new acceptance thresholds or universal accuracy
claims are inferred. Baseline regeneration and reducing maximum starting age
were changed together in the inheritance assessment; their effects were not
isolated. Earlier fixed-baseline results remain conditional on their baselines.

This documentation update does not independently rerender the report, download
rates, or rerun simulations/tests, and does not assert every review suggestion
was implemented. The source with the revised conclusion is recorded at commit
`7514916bb8cf0c06faef24dad7fa013304093779`. This is a source reference,
not a final rendered-artifact hash or a report-plus-documentation acceptance
commit. The user handles final rendering; no final HTML hash is assigned here.
Preserve both source and HTML. Reproduction requires the input preparation and
software environment recorded in the report; rendering installs the package,
accesses HMD, runs substantial parallel simulations, and executes the test suite.


## Milestone 4 dissolution implementation checks

On October 27, 2026, `TZ=UTC Rscript -e 'devtools::test(reporter = "summary")'`
passed with the optional duration-dependent dissolution process (D078), without
failures or test warnings; the live database test was skipped. Targeted
`devtools::test(filter = "union-dissolution", reporter = "summary")` checks
constructor/edited-process validation and provenance, constant-hazard survival
against `exp(-rate * duration)`, one risk per couple, duration bands and gaps,
attained-duration restarts, retained dissolution through births, exact
birth/dissolution and death ties, widowhood invalidation in both initiation
modes, repeated formation and spacing, group-specific hazards, future-group
coverage, empty/zero-duration runs, and seed isolation. Exposure accounting
checks confirm that union closure shortens partnered person-years and later
deaths preserve dissolved histories.

These are software and analytical checks, not milestone 4 scientific acceptance.
Age-specific previous-union search schedules are now implemented under D080;
duration-dependent remarriage and cohabitation remain deferred. D045's existing
fertility-transition limitation remains; no new demographic rates or acceptance
thresholds are introduced.

The subsequent birth-status clarification in D078 is checked with
`devtools::test(filter = "union-dissolution|lexis-unions|group-demography",
reporter = "summary")`, which passed without failures or test warnings.
Births exactly at dissolution ends count as partnered while exposure remains
half-open; tests cover period boundaries, overlapping histories, unchanged
death-ending classification and overall-count conservation.

After adding `calculate_dissolution_rates()` on October 27, 2026, the full
`devtools::test(reporter = "summary")` suite passed without failures or test
warnings; the opt-in live database test was skipped. Analytical union-year
checks cover exact calendar/duration splitting, dissolution versus death
endings, censoring, boundary events, repeated episodes, group attribution,
zero-exposure rates, empty grids and pooled accounting identities.

D079 changes dissolution schedules, multipliers and union observations to
symmetric `group_pair` strata. The full suite found only an exact-equality
assertion differing by floating-point rounding between a multiplied rate and
a decimal literal; it was changed to a numerical-tolerance comparison. The
subsequent `devtools::test(filter = "union-counts|union-dissolution", reporter =
"summary")` passed. All other full-suite tests passed, with the opt-in live
database test skipped. Checks include reversed pair equivalence, nonunit
multipliers versus explicit rates, pair coverage (including future identities),
canonical numeric/literal labels, pair-pooled accounting, empty grids and
rejection of ambiguous or overlapping specifications.

## Age-specific remarriage and marriage observations (D080, October 28, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite without test failures or warnings. The opt-in live
HMD/HFD test was skipped; no rates were downloaded. UTC avoids the container's
unavailable system timezone service. The first full run identified one remaining
test using the removed Lexis formation column; it was migrated to the new helper
before the successful full run.

`test-remarriage.R` checks complete category/group schedules, invalid labels and
overlap, recovery from unordered histories, latest-ending replacement, dissolution
and widowhood rescheduling, analytical constant-hazard incidence for all three
categories, reproducibility with shared/equal schedules, explicit group multiplier
equivalence, and calibrated intensities against explicit effective group rates.
`test-marriage-counts.R` checks exact exposure and starts, both sexes, latest
category, pre-window history, half-open time/age boundaries, touching and
zero-duration episodes, grouping and pooling, empty output types and invalid
histories. Existing Lexis/dissolution tests retain birth/death exposure checks
and now use the separate marriage helper for starts.

The shipped milestone 2 example now calls `calculate_marriage_rates()` and
pools category counts/exposure before calculating a combined marriage rate.
Its source parses successfully; the complete empirical example was not run
because it downloads rates. Accepted milestone 2/3 report sources and rendered
artifacts are preserved. Their former Lexis marriage columns are intentionally
removed in D080, so rerender historical assessments at their recorded versions
or migrate those analysis chunks explicitly. The pre-D080 repository head was
`a847d49313c68fec9b3506ba790bfa4f86ee7c5a`.

These checks validate implementation and analytical cases, not milestone 4
scientific acceptance. Duration-dependent remarriage and transition-aware
fertility compensation remain deferred.


The subsequent D080 pooling option passed
`TZ=UTC Rscript -e 'devtools::document(); devtools::test(filter="marriage-counts|remarriage|union-dissolution", reporter="summary")'`
without failures or warnings. Tests cover pooled counts/exposure and rates
against summed previous-category cells, independently with and without group
stratification, empty grouped/ungrouped outputs, and invalid option values.
`by_previous_union` defaults to FALSE; stratified checks explicitly opt in.

## Observation helper names (D081, October 28, 2026)

After renaming the public helpers to `calculate_demog_rates()`,
`calculate_marriage_rates()` and `calculate_dissolution_rates()`,
`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed without failures or warnings; the opt-in live database test was skipped.
Namespace checks confirm the new exports and removal of the old names.
The shipped example and current milestone 4 script parse with the new calls.
Arguments, output columns and demographic calculations are unchanged.

## Actual female hazard contributions (D082, October 28, 2026)

`devtools::test()` passed the full suite using `TZ=UTC`, with no failures or
warnings; the opt-in live database test was skipped. Targeted validation first
passed with `filter = "calibration-intensity|weighted-choice|remarriage|rate-unions|union-dissolution"`.
The native code was built with `devtools::load_all(debug = FALSE, recompile = TRUE)`
for these checks and separately installed with optimized flags for benchmarking.

`test-calibration-intensity.R` compares every maintained cell/channel with an
independent individual scan after age boundaries, narrow bands, gaps, terminal
support, spacing releases, repeated readiness changes, removal and reactivation
with another previous category. It checks RNG preservation, exact zero after
all positive contributors leave, and at most one pending boundary per active
woman. Existing matching, dissolution, fertility, restart and independent-run
checks passed; snapshot calibration tests now use actual individual hazards,
including different risks inside the same quarter-year bin.

The offline optimized old/new comparison used the same saved models, 50,000
starting persons and 100 years, with three seeds and alternating version order.
The ungrouped seed-42 history was identical in full. Median grouped timings
changed +2.3% for shared-rate searches, -8.1% for previous-union schedules and
+3.5% for on-demand searches. Grouped previous-status peak process RSS rose
about 12.1 MiB. See [current-hazard calibration](calibration-intensity.md) for
inputs, runtime environment, raw results, reproduction commands and limits.
These are performance/implementation checks, not milestone 4 acceptance or
validation of duration-dependent remarriage, which remains unimplemented.

## Matching parameter recovery with the new tracker (October 28, 2026)

Ran
`TZ=UTC Rscript validation/validate-calibration-log-odds.R /tmp/hammel-calibration-newlib validation/calibration-log-odds 20 1500`
using the same optimized D082 build as the performance comparison. The script
uses offline synthetic schedules, 20 independently seeded initial populations
of 20,000 persons (max age 50, 75/25 A/B), 60 years, duration-dependent
*dissolution* and age-specific remarriage searches by previous status. Search
rates differ across previous categories and B has a 1.2 search multiplier.
Duration-dependent remarriage is not used. Both matching cases use A:B = -3,
so the conditional cross-group log odds ratio target is -6. The age-dependent
case additionally uses preferred male-minus-female gap 2 and half-weight
distance 5 (linear .11090355, quadratic -.02772589).

For each case/replication, sample 1,500 union events and 50 uniformly sampled
eligible alternate men per event, plus the chosen man. Independently rebuild
the available male pool immediately before formation from full person/union
histories, enforcing living status, age >=15, and unpartnered status. Kin
screening is explicitly disabled in these simulations, avoiding a mismatch
between simulation and reconstructed choice sets. Fit wife-side conditional
logit with the recorded per-event male-group log weights as offsets, separate
AB/BA coefficients, and age/age-squared terms in the age-dependent case.
Uniform sampled alternatives have the same set-selection probability for
each possible chosen member, so no additional alternative-sampling offset is
required. Efron and exact conditional likelihoods coincide with one event per
stratum. Replication SEs below are sd(estimates)/sqrt(20), including simulation,
baseline and choice-subsampling variation; they are not within-fit SEs.

| Case | Target conditional log OR | Mean recovered | MCSE | Pooled cross-tab mean | Pooled MCSE |
| --- | ---: | ---: | ---: | ---: | ---: |
| Group only | -6 | -5.97493 | .05378 | -5.98500 | .01721 |
| Age dependent | -6 | -5.98931 | .04935 | -5.96394 | .01439 |

Directional coefficient means were -3.00624 (AB) and -2.96868 (BA) for group
only, and -3.01641 and -2.97290 for age-dependent matching. Recovered age terms
were .1117834 (MCSE .0017798) and -.0276966 (MCSE .0002403). The input
conditional parameters lie within the corresponding approximate Monte Carlo
uncertainty intervals. These results support parameter recovery for these
workloads; no universal acceptance threshold or milestone completion is inferred.

A pooled cross-tab is a different estimand: changing age distributions, pools
and weights can move its cross-product away from the conditional target.
In particular, do not claim the age-dependent pooled estimate is within two
MCSEs of -6; preservation applies to the conditional utility association, not
an exact constraint on the final pooled table.

Reproduction: [validation script](../validation/validate-calibration-log-odds.R).
Evidence: [replication estimates](../validation/calibration-log-odds-replications.csv),
[summary](../validation/calibration-log-odds-summary.csv), and
[estimates/session metadata](../validation/calibration-log-odds.rds). No package
code changed for this recovery check. The run completed without warnings.
