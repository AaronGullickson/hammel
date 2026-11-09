# Testing and milestone validation

Updated November 6, 2026. Milestones 1–9 are accepted as complete (D029,
D050, D077, D093, D105, D116, D127, D130, D141). Software test history and the user's demographic assessment
are distinct evidence; recording acceptance does not independently rerun either.
Earlier implementation-specific runs and counts remain available in Git history.

## Noun-based specification constructors (November 9, 2026)

D147 renames the matching, age-gap matching, weight, and group-inheritance constructors, retaining `create_baseline_population()`. A parsed-expression comparison against the pre-rename working tree confirmed that executable R definitions differ only by those four name replacements, including updated calls and diagnostic messages. Export checks confirmed exactly those replacements, with 43 exports overall and no aliases under the old names.

A before/after simulation using baseline generation, maternal group inheritance, age-gap matching, and mortality weights produced identical population results and random-number states. All four renamed constructor help pages passed Rd validation and their examples executed successfully. A second `devtools::document()` pass completed without the transitional unresolved-name warnings from the initial namespace rename.

The full offline `devtools::test()` run passed 6,496 expectations with zero failures in 108.6 seconds. The live database test was skipped. One existing environment warning concerned HMDHFDplus/lubridate's `timedatectl` access to the system bus in the sandbox. `git diff --check` passed. Native sources, user-maintained README/vignettes, and historical validation sources/artifacts were untouched; old constructor calls in those historical or user-maintained sources require updating before execution against the renamed API.

## Shared fertility-summary and group-declaration help (November 9, 2026)

D146 combines TFR/PATFR and group construction/declaration into two paired help topics. `devtools::document()` completed successfully, removing the superseded `calculate_patfr.Rd` and `set_groups.Rd` pages. Both shared topics passed Rd validation, alias checks, argument-coverage checks, text rendering, and execution of all examples. The TFR examples exercise both explicit PATFR and automatic parity dispatch; the group examples describe a dimension and declare it on a model.

Parsed executable R expressions and exports are unchanged from the preceding commit. Both matching-model help pages remain unchanged. `git diff --check` passed; no simulation, download, or full-suite rerun was needed for this documentation-only change.

## Shared observed-rate help topic (November 9, 2026)

D145 combines the seven observed-rate help pages into `calculate_rates`. `devtools::document()` completed successfully and removed the seven superseded generated pages. The shared topic passed `tools::checkRd()`, has aliases for all seven public functions, and documents the union of their arguments exactly once. Examples cover all seven functions and executed successfully after `pkgload::load_all()`, using supplied person/union/transition histories with no simulation or data download. Text help rendered successfully through `tools::Rd2txt()`.

Comparing parsed R expressions against the preceding commit confirmed that all executable definitions are unchanged, and `NAMESPACE` is unchanged. `git diff --check` passed. The simulation suite was not rerun for this documentation-only change.

## Shared process help topics (November 9, 2026)

D144 pairs each demographic process constructor with its convenience setter on one help topic, with separate mortality and fertility pages and a standalone generic `set_process()` page. `devtools::document()` completed successfully. All six affected topics passed Rd validation and alias checks; every available example on those topics executed successfully after `pkgload::load_all()`. Union formation and group transitions have no examples in these topics. Comparing parsed R expressions against the preceding commit confirmed that all executable definitions are unchanged, and `NAMESPACE` is unchanged. `git diff --check` passed. The simulation suite was not rerun for this documentation-only change.

## Within-file ordering and section headings (November 9, 2026)

Following the user's instruction under D143, each hand-maintained R file now starts with a purpose comment, followed by exported functions and then internal constants/helpers under `# SECTION NAME ------` headings. Generated wrappers remain untouched. Before/after parsing confirmed that every function and constant is unchanged and remains in its original file; all exported functions precede helpers. Explanatory comments and roxygen blocks were retained, and section headings were updated where necessary.

`devtools::document()` produced no help-page or export changes. The full offline `devtools::test()` run passed 6,496 expectations with zero failures in 108.8 seconds, with the live database test skipped. The same environment warning about `timedatectl` system-bus access occurred while loading HMDHFDplus/lubridate. The final purpose comments and heading adjustments were checked by parsing again and `git diff --check`; they do not change executable R expressions.

## R source organization (November 9, 2026)

D143 reorganizes the R files according to the user-agreed [source organization guide](r-source-organization.md). Before/after parsing confirmed that all 226 function definitions, six shared constants, and package metadata are identical as R expressions. Every original source comment was retained. The resulting 34 R files match the agreed map; `simulate()` is unchanged, and `model_schedules.R`, `hammel_model.R`, generated wrappers, and native sources are unchanged.

`devtools::document()` completed successfully. All 43 exports are unchanged, and every generated help page is identical except for its source-file pointer. Explicit roxygen `@order` tags preserve the original order of functions and parameters in the shared mortality/fertility and process-setter topics.

The full offline `devtools::test()` run passed 6,496 expectations with zero failures in 108 seconds. The opt-in live database test was skipped. One environment warning occurred while loading HMDHFDplus/lubridate: `timedatectl` could not access the system bus in the sandbox, as in the preceding recorded run. No rates were downloaded and no credentials were exposed. `git diff --check` passed. This mechanical reorganization does not change demographic behavior or rerun the user-owned milestone reports.

## Fixed newborn assignment (November 7, 2026)

D142 adds native fixed-category inheritance. `devtools::test(filter =
"group-inheritance")` passed 149 expectations with no failures or warnings.
The full offline `devtools::test()` run passed 6,496 expectations with no
failures; the live database test was skipped as intended. One environment
warning came from HMDHFDplus loading `timedatectl`, whose system bus access
was denied in the sandbox. No live rates were downloaded.

The added cases verify fixed assignment for either sex and an unknown father,
unchanged founders and ancestry, numeric categories, calendar switches,
invalid/undeclared values, reproducibility, and exact event-history/RNG agreement
with an equivalent deterministic custom rule. The implementation avoids the
R callback; no runtime speedup was benchmarked.

## Milestone 9 accepted assessment and closeout

The user accepted milestone 9 as complete on November 6, 2026 (D141),
following the user-authored [assessment](../validation/milestone9.qmd), which
the user maintains and renders. The report's settings and findings remain
authoritative; this closeout neither reruns nor renders it. Code at
acceptance: `fc4b9dbf29dba9782cec9265704f9e87fef4022c` plus closeout
documentation.

Software checks at closeout (November 6, 2026, local; GitHub CI did not run,
because the account's Actions minutes were exhausted):

- The full offline suite (`devtools::test()`) passed all 551 tests with no
  failures or warnings, also with `options(hammel.check_cells = TRUE)`.
- `R CMD check` (`devtools::check()`): 0 errors, 0 warnings, 1 note about
  non-portable compiler flags from the local R configuration (not the package).
- Coverage (`covr::package_coverage(type = "tests")`): 97.8%.
- With one dimension, results were checked bit-identical to the preceding
  stage at each stage (builds with fused multiply-adds disabled), and
  combined-dimension analysis, matching and household moves have exact
  worked-example tests; see [multiple group dimensions](multiple-groups.md)
  and [household moves](group-transitions.md#household-moves-d140).
- Sex symmetry of intermarriage with a second dimension, within combinations
  and within separate markets, and calibration cost, are recorded in
  [multiple group dimensions](multiple-groups.md) with their scripts in
  `validation/milestone9/`.

Limitations retained at acceptance are listed in D141.

## Milestone 8 accepted assessment and closeout

The user accepted milestone 8 as complete on November 5, 2026 (D130),
following the user-authored [assessment](../validation/milestone8.qmd), its
[rendered report](../validation/milestone8.html), and explicit acceptance. The
report's settings and numerical findings remain authoritative; this
documentation closeout neither reruns nor renders it.

Software checks at closeout (November 5, 2026, local, not a GitHub CI rerun):

- The full offline suite (`devtools::test()`) passed all 500 tests with no
  failures, also with `options(hammel.check_cells = TRUE)`, which rescans every
  process's cells, intensity tree and matching-calibration totals at each
  event or search solve.
- Coverage (`covr::package_coverage(type = "tests")`) was 98.11%, with
  `R/weight_model.R` fully covered.
- Weight-model tests (`test-weight-models.R`) check within-cell ratio recovery
  for each process within four standard errors, exact reproduction of
  unweighted runs with zero coefficients or constant custom functions,
  cancellation warnings, NA handling, the custom-function contract for every
  process, weight schedules (including a switch combined with an in-place rate
  update) and recovery of rejection sampling after an extreme-weight member
  leaves.
- Step 1 equivalence, recovery and timing checks are recorded in
  [cell scheduling](binned-group-hazards.md#step-1-validation-and-performance-november-4-2026),
  and weight-model costs in [its cost table](binned-group-hazards.md#cost-november-4-2026).

The D082 calibration-tracker checks below describe the engine before milestone
8; `test-calibration-intensity.R` now runs a simulation-level comparison of
cell-derived calibration totals with explicit scans.

## Milestone 7 accepted assessment and closeout

The user accepted milestone 7 as complete on November 2, 2026 (D127), following
the user-authored [assessment](../validation/milestone7.qmd), its
[rendered report](../validation/milestone7.html), and explicit acceptance.
The report concludes that its expectations were met for calendar-varying
rate/probability and model scenarios. The report's settings and numerical
findings remain authoritative; no new thresholds or validation results are
inferred. The code/report revision at acceptance is
`3b65d2202b006ef6f3e232c98595814b1544f9fe`.

Acceptance retains the current engine, documented model and observation
limitations, and approximate demonstration inputs. Milestone 8 later replaced
the engine's scheduling (D128) and added weight models (D129).
Reports, rendered artifacts and active intermediates are preserved unchanged;
this documentation-only closeout neither reruns nor renders them. The earlier
software-check records below are distinct from scientific acceptance.

## Optional Poisson rate intervals (D124, November 1, 2026)

`TZ=America/Los_Angeles Rscript -e 'devtools::test(stop_on_failure=TRUE)'`
passed the full offline suite: 4,318 assertions, no failures or warnings, one
opt-in live database test skipped. The 182 assertions in `test-rate-confidence.R`
cover all seven public rate calculators, comparison with `stats::poisson.test()`
at 95% and 80%, unchanged NULL/default output and attributes, group/type/parity
strata, positive and zero event counts, zero exposure including endpoint events,
empty grouped tables, invalid levels, and observed-birth intervals alongside
fractional unrealized fertility. Public help was regenerated and documentation
links and `git diff --check` checked. No rates were downloaded or reports edited.
Exactness refers to the Poisson counting model; these tests do not establish
coverage under all demographic dependence or milestone 7 scientific acceptance.

## Milestone 7 calendar-input implementation checks

October 31, 2026 (D118–D120, D122–D123): the full offline `devtools::test(stop_on_failure = TRUE)`
suite passed with 4,136 assertions, no failures or warnings, and one opt-in live
DB test skipped. The run used `TZ=America/Los_Angeles` so optional database-reader
dependencies did not need to query the unavailable sandbox system bus for a time
zone. No rates were downloaded and no credentials were used.

New calendar-rate checks cover malformed bounds, overlapping risks, calendar
coverage and endpoints, fractional years, analytical integrated mortality,
search/formation, dissolution, conversion and directed transition hazards,
retained duration origins, birth spacing and persistent weights, exact
unrealized-intensity accounting, changing parity/union-type fertility schedules,
asynchronous process changes, calibration/cell consistency, exact calendar/event
ties, reproducibility, RNG isolation, future group/pair coverage and empty runs.
The fractional-year example in `inst/examples/milestone7-rates.R` executed
successfully. Generated public help, local documentation links and
`git diff --check` were also checked.

D119 adds offline multi-year HMD/HFD and native HMD shape checks: single-year
schema preservation, chronological vector selection, missing/duplicate years,
download counts, unchanged rate units, all six shapes, infant/child/open age
bands, raw and parsed source periods, shortened endpoints, auxiliary annual
coverage, HFD parity partitions and selected missing-cell counts, provenance,
and direct simulation with prepared calendars. Synthetic local HMD parsing
exercises the installed reader; automated web requests use mocks and the
opt-in live test was not run.

D120 adds 77 assertions for optional HMD open-age pooling: ratios of summed
counts/exposure, sex/period separation, unchanged younger rates, missing replaced
tail Mx, positive pooled exposure, zero-exposure individual cells, zero deaths,
invalid counts and cutoffs, native grouped ages, raw/parsed source periods,
source alignment independent of row order, download counts, safe retained
provenance and analytical survival using a pooled schedule. The default getter
still rejects missing rates, now with year/sex/age context. The pooled web path
was checked with mocked readers; no new live pooling downloads were made.

D122 adds 60 assertions in ten calendar-probability cases: valid/edited/future
probability slices; age/calendar overlap and coverage; group/category coverage;
fractional type switches; non-extreme type shares; unchanged aggregate search
hazards; exact tied searches; retained ongoing episodes and conversion/dissolution
clocks; current age/previous-category selection; zero/empty runs, restarts and
RNG isolation; and shared/asynchronous rate/probability boundaries. A grouped
matching regression changes probabilities before the first realized search and
compares subsequent unions and logged role weights with the equivalent static
new-probability reference, checking immediate calibration refresh. Existing
static cohabitation tests also pass. No rate inputs were downloaded.

D123 adds 112 assertions for whole-model matching/inheritance calendars:
all four argument roles; direct models and normalized list-column tables;
invalid, edited and future entries; overlap and run coverage; exact search and
birth boundaries; dynamic NULL fallback versus fixed explicit models; native
and custom matching formulas; within-group restrictions and weighted calibration
in both directions; hierarchy changes, custom callbacks and numeric/character
future identity registries; unchanged earlier identities and ancestry; restarts;
zero/empty runs; shared rate/probability/model boundaries; reproducibility and
RNG isolation. Controlled candidate pools test event-time selection, and static
reference runs check retained event dates and refreshed role weights. The user's
milestone 7 report was neither edited nor rendered.

This is software verification, not scientific acceptance of milestone 7 or a
new demographic validation report. Calendar boundaries redraw person/union
candidates and require work proportional to the population; redundant boundaries
do not promise fixed-seed agreement with static inputs. Probability-only and model-only
boundaries retain all pending clocks and consume no demographic randomness.
Other scalar process arguments remain static. See
[the rate contract](time-varying-rates.md) and
[the model contract](time-varying-models.md).

## Milestone 6 accepted assessment and active final render

The user accepted milestone 6 on October 30, 2026 (D116) after report and code
review. The [report source](../validation/milestone6.qmd) and
[rendered assessment](../validation/milestone6.html) are the user's scientific
assessment of age/duration transition recovery, sex/union-type strata,
fertility/mortality and formation-group summaries, population shares and
performance. The final render is active at this closeout. Preserve the source,
HTML and all intermediates; this documentation update neither reruns nor
certifies completion of that render or its final numerical results.

Software checks and prior developmental measurements below are separate
records. Acceptance retains finite-sample and sparse-cell limits, complete
history requirements, duration fallback to age without prehistory,
non-checkpoint restarts and the documented life-table aggregation caveat.
Household moves remain deferred (D114), with unresolved regional applications
retained under milestone 9's multiple group columns (D121, D131). Calendar-varying rates
are implemented under D118, with matching and inheritance schedules under D123.
Milestone 7 was subsequently accepted complete on November 2, 2026 (D127).

## Milestone 5 accepted assessment and closeout

The user accepted milestone 5 on October 29, 2026 (D105) after revising and
reviewing the [report source](../validation/milestone5.qmd). The final
[rendered assessment](../validation/milestone5.html) was being produced when
closeout was requested. This documentation update does not rerun simulations,
download HFD/HMD data, execute tests or certify final-render completion. Preserve
report sources, artifacts and active intermediates.

The assessment covers parity-specific conditional fertility, cell-accounted
spacing and revisited baseline/restart/union-type checks, inherited lognormal
weights, childbearing dispersion, mother–daughter association and performance.
Its limitations include sparse or structurally infeasible cells, unrealized
expected births (an accounting diagnostic rather than actual births),
history completeness and selection-driven weight distributions. Accepted
completion does not introduce numerical thresholds or universal rate-recovery
claims. Individual group transitions (milestone 6) are now implemented
under [their contract](group-transitions.md), accepted October 30 (D116).

Earlier dated verification records below describe the versions tested at the
time, including superseded corrections and intermediate proposal status.
Current execution includes D101–D104 and D112; accepted milestone status includes
D116. The October 29 milestone 5 closeout checked local Markdown links and
whitespace; that closeout made no engine or report edits.
Supporting evidence and reproduction-version notes are in [validation/milestone5/](../validation/milestone5/README.md).

## Milestone 4 accepted assessment

The user accepted milestone 4 on October 29, 2026 (D093), after final report
corrections and adding duration-dependent conversion. The
[report source](../validation/milestone4.qmd) and
[rendered assessment](../validation/milestone4.html) are the accepted scientific
assessment. The author's final HTML render was active at documentation closeout;
this edit does not independently rerender it, download HMD rates or rerun tests.
The previously reviewed rendered report contained 322 test cases: 321 passed,
one opt-in live database case skipped, no failures/errors/warnings. Treat the
final report's own results as the record for that render.

Supplementary benchmark and conditional-choice artifacts are retained under
[validation/milestone4/](../validation/milestone4/README.md), with their original
run provenance. Paths in reproduction commands below reflect that move.
The final report supersedes the removed exploratory `milestone4_tests.R`.
[Validation contents](../validation/README.md) lists retained inputs/evidence.
Earlier dated entries below distinguish implementation checks from scientific
acceptance at the time; D093 records the subsequent user acceptance. Completion
retains sparse-cell, initialization, fertility-transition and calibration limits.

## Run the suite

October 30 external-review follow-up: `prepare_union_group_observations()` now
uses native indexed partner histories for both grouped dissolution and
conversion. On synthetic data with 40,000 people, 20,000 unions and 8,000
moves, the old adapter took 8.016 seconds and the new one 0.032 seconds,
producing identical adapted union histories. Full grouped dissolution and
conversion calls took 0.299 and 0.211 seconds, respectively, including
validation. These are local development timings, not performance guarantees.
The full offline `devtools::test()` suite passed without failures or warnings;
the opt-in live database test was skipped. New regressions cover simultaneous
partner moves, typed numeric group pairs, conversion exposure, and contextual
errors for missing groups and overlapping union histories. No event-generation
or RNG behavior changed.

October 30 PR #7 coverage follow-up: added tests for invalid transition units,
history references/life dates, unsupported processes, unknown/nonmoving identity
lookups, legacy populations without transition components, and union exposure
when only an unrelated person moves. Targeted `devtools::test()` checks and a
full `covr::package_coverage()` run passed. On the local ARM/GCC environment,
overall line coverage increased from 97.5539% to 97.8234%; both
`R/group_transition.R` and `R/calculate_transition_rates.R` reached 100%.
GitHub's compiler/platform can produce slightly different totals. No engine
behavior, coverage exclusions or thresholds changed; the native consistency
error guard remains uncovered rather than forcing invalid internal state.

October 30 transition-rate counter: the full offline `devtools::test()` suite
passed after adding `calculate_transition_rates()`. The live database test was
skipped and the existing environment `timedatectl` warning remained.
`test-transition-rate.R` checks hand-calculated age/duration exposure with
repeated moves and death, directional denominators, sex/type pooling, exact
formation/conversion/dissolution boundaries, calendar endpoints, empty typed
outputs, invalid inputs, and recovery of simulated competing annual hazards.
A synthetic 50,000-person, 100,000-move period/duration/sex/type calculation
took 1.63 seconds including population validation, recovering 100,000 events
and 5,000,000 person-years (two groups, so one destination per origin).
This timing is a development measurement, not a performance guarantee.

October 30, 2026 group-transition implementation: `devtools::test()` passed
the full offline suite with no failures. The live database case was skipped
without opt-in; the local HMD reader test emitted an environment warning from
`timedatectl` (status 1). Focused `devtools::test(filter = "group-transition")`
checks cover analytical competing destination probabilities, age and duration
support, duration across restarts, sex/union-type strata, initially empty
destinations, malformed schedules, typed empty transition histories, ancestry
retention, mortality after moving, matching-pool and fertility-cell updates,
pair-dependent dissolution and cohabitation conversion, and historical group
exposure/event attribution. Cell consistency diagnostics are enabled in the
pool/fertility integration case. These are software and analytical checks;
the user subsequently accepted milestone 6 under D116.

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
limitations are recorded in [data provenance](../validation/data/usa_2000_2001/usa_2000_rates.README.md).
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
duration-dependent remarriage was not covered by these checks (now D083 below).
Cohabitation is still deferred and D045's fertility-transition limitation remains;
no new demographic rates or acceptance thresholds are introduced.

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

The shipped milestone 2 example now calls `calculate_union_rates()` and
pools category counts/exposure before calculating a combined marriage rate.
Its source parses successfully; the complete empirical example was not run
because it downloads rates. Accepted milestone 2/3 report sources and rendered
artifacts are preserved. Their former Lexis marriage columns are intentionally
removed in D080, so rerender historical assessments at their recorded versions
or migrate those analysis chunks explicitly. The pre-D080 repository head was
`a847d49313c68fec9b3506ba790bfa4f86ee7c5a`.

These checks validate implementation and analytical cases, not milestone 4
scientific acceptance. Duration-dependent remarriage was outside that check;
D083 below covers its implementation. Transition-aware fertility compensation
remains deferred.


The subsequent D080 pooling option passed
`TZ=UTC Rscript -e 'devtools::document(); devtools::test(filter="marriage-counts|remarriage|union-dissolution", reporter="summary")'`
without failures or warnings. Tests cover pooled counts/exposure and rates
against summed previous-category cells, independently with and without group
stratification, empty grouped/ungrouped outputs, and invalid option values.
Previous-union strata are off by default; stratified checks explicitly request them.

## Observation helper names (D081, October 28, 2026)

After renaming the public helpers to `calculate_demog_rates()`,
`calculate_union_rates()` and `calculate_dissolution_rates()`,
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
`TZ=UTC Rscript validation/milestone4/validate-calibration-log-odds.R /tmp/hammel-calibration-newlib validation/milestone4/calibration-log-odds 20 1500`
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

Reproduction: [validation script](../validation/milestone4/validate-calibration-log-odds.R).
Evidence: [replication estimates](../validation/milestone4/calibration-log-odds-replications.csv),
[summary](../validation/milestone4/calibration-log-odds-summary.csv), and
[estimates/session metadata](../validation/milestone4/calibration-log-odds.rds). No package
code changed for this recovery check. The run completed without warnings.

## Duration-specific remarriage (D083, October 28, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite without failures or warnings; the live database test was
skipped. Subsequent tests of births and failed searches preserving the duration
origin passed with `devtools::test(filter = "duration-remarriage", reporter = "summary")`.
The final duration test file passed all 51 expectations. No rates were fetched.

`test-duration-remarriage.R` covers the single-table age/duration format,
exclusive clock selection, all-category coverage, consistent clock per category
across groups, invalid/partial/overlapping bands, required age clock for `none`,
optional all-NA duration columns, duration gaps and analytical cumulative hazard,
attained-duration restarts, unordered histories, age independence, repeated
ending resets, widowhood origin, and preservation through births/failed searches.
A grouped duration schedule and its mathematically equivalent age schedule
produce identical full histories and weight diagnostics. Tracker tests use
separate birth and search origins to verify duration boundaries independently
of fertility ages and spacing. Missing previous-ending origins fail explicitly
in deterministic calibration checks. Existing age-only grouped simulation
(seed 42, 50,000 initial persons, 100 years, prior D082 benchmark input) was
compared in full against the installed pre-D083 build and was identical.

Additional parameter-recovery check:

```sh
TZ=UTC Rscript validation/milestone4/validate-calibration-log-odds.R /tmp/hammel-duration-lib validation/milestone4/duration-remarriage-log-odds 10 1500 duration
```

This uses ten independent 20,000-person baselines per matching case and the
same 60-year synthetic setup and choice-set estimation method described above.
`none` uses the original age schedule. Dissolved/widowed schedules now use
durations [0,1), [1,3), [3,8), [8,Inf), with rates 0/.15/.10/.05 for dissolved
and 0/.07/.04/.02 for widowed. Thus actual duration boundaries, delays,
dissolution, widowhood, repeated remarriage and fertility all operate together.
B retains its 1.2 search multiplier. Kinship screening is explicitly disabled
for exact eligibility reconstruction. No convergence warnings were emitted.

| Matching case | Target conditional log OR | Mean recovered | MCSE |
| --- | ---: | ---: | ---: |
| Group only | -6 | -5.91620 | .07840 |
| Age dependent | -6 | -6.05778 | .07244 |

Mean age coefficients were .1087908 (MCSE .0023092) and -.0280502
(MCSE .0004512), against inputs .1109035 and -.0277259. Conditional combined
log odds ratios are consistent with their target at this Monte Carlo precision.
These checks support the implementation in these cases; ten replications do
not establish universal scientific accuracy or milestone 4 acceptance. Pooled
cross-tab means are retained in the raw results and are not exact target
constraints under changing pools/age composition.

Evidence: [replication estimates](../validation/milestone4/duration-remarriage-log-odds-replications.csv),
[summary](../validation/milestone4/duration-remarriage-log-odds-summary.csv), and
[estimates/session metadata](../validation/milestone4/duration-remarriage-log-odds.rds).
The public marriage-rate observation helper remains age-based; duration-band
exposure summaries require a separate interface specification.

## Duration marriage observations (D084, October 28, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite without failures or warnings; the opt-in live database
test was skipped. The earlier targeted `filter = "marriage"` checks also passed.
No rates were downloaded and no simulation mechanisms were changed for D084.

`test-marriage-duration-rates.R` verifies exact duration exposure and starts,
automatic previous-category stratification even with an explicit FALSE flag,
different numbers of age/duration bands, unused boundary NA columns, finite
last-duration Inf extension, pre-window history, resetting after each ending,
half-open age/duration/period event boundaries, death censoring, all-age duration
exposure, grouping/type preservation, and conservation of exposure and counts
against compatible age-only totals. Touching/zero-duration episodes, empty
grouped/ungrouped results, invalid breaks and unchanged NULL output are covered.

A 3,000-woman synthetic simulation with duration-specific hazards .1, .3 and .2
on [0,1), [1,2), [2,Inf) recovers those female hazards with the new exposure
helper under the analytical six-standard-error software check. Available men
are sufficient and mortality/fertility are zero, so realized formations reflect
search hazards in that test. This does not establish empirical validity or a
new milestone acceptance threshold.


## Typed cohabitation (D086, October 28, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite without failures or warnings. The opt-in live database
check was skipped; no rates were downloaded. New cohabitation checks cover
probability gaps and strata, independent matching preferences, exclusive male
availability, linked conversion episodes, exact event ties, restart history,
legal previous-union recovery and its toggle, typed fertility/dissolution,
pair multipliers, and fertility multiplier lists. Analytical conversion versus
dissolution competition and calibration contribution sums at probability band
boundaries are included. A saved pre-feature grouped age/previous-status run
(50,000 initial people, 100 years, seed 42) retains identical complete output
when cohabitation is disabled.

`validation/milestone4/validate-cohabitation-matching.R` was run for ten independently
initialized 20,000-person populations over 60 years. Probability varied by age,
group and previous union; conversion, typed fertility and typed dissolution
were active. Conditional-choice fits used 1,500 entry events per type per
replication, sampled available alternatives, age terms and applied group-weight
offsets. Conversion episodes were excluded from new-partner choices. Combined
conditional log odds averaged -2.0402 (MCSE .04280) for cohabitation against -2,
and -6.0298 (MCSE .03625) for marriage against -6. These support parameter
recovery in this synthetic case, not empirical validity or milestone acceptance.

Evidence: [script](../validation/milestone4/validate-cohabitation-matching.R),
[replications](../validation/milestone4/cohabitation-matching-replications.csv),
[summary](../validation/milestone4/cohabitation-matching-summary.csv), and
[estimates/session metadata](../validation/milestone4/cohabitation-matching.rds).


## Union and marriage observation split (D087–D088, October 28, 2026)

`TZ=UTC Rscript -e 'devtools::test(filter="marriage|cohabitation|lexis|group-demography", reporter="summary")'`
passed without failures or warnings after regenerating documentation and native
exports. This covers both renamed union entry observations and new marriage
entry observations. Exact cases verify shared exposure and count/rate
conservation for type-specific union rates, group/previous/duration combinations,
conversion exclusion from union entry, conversion inclusion in marriage entry,
cohabiting exposure, half-open age/period boundaries, death censoring,
dissolution, direct and legacy marriage equivalence, zero exposure and empty
histories. No simulation scheduling behavior changed in this split, and no
empirical downloads were used. This is targeted verification, not a new full-suite
or scientific milestone acceptance claim.


D089 dissolution type observation: `devtools::test(filter =
"dissolution|marriage-entry", reporter = "summary")` passed without failures or
warnings. Tests cover separate type exposure/events, conservation against pooled
totals, crossing with symmetric group pairs, conversion censoring and duration
reset, empty inputs, zero-exposure NA rates, invalid type columns and logical
flags. Documentation was regenerated; no native engine changes were needed.


## Demographic current-type strata (D090, October 28, 2026)

`devtools::test(filter = "demog|lexis|cohabitation|dissolution|marriage",
reporter = "summary")` passed without failures or warnings after documentation
and native exports were regenerated. Exact typed cases verify disjoint exposure,
conservation against pooled births/deaths/person-years with and without groups,
births at conversion and dissolution boundaries, left-limit death type,
discovered single types, empty typed histories, and rejection of ambiguous
overlapping typed episodes. Migrated untyped cases retain merged overlap
coverage, event boundaries and downstream TFR/life-table grouping. Pooled default
observations remain covered by the demographic tests. This was targeted offline
verification; no simulation scheduling or matching behavior changed.


D091 conversion observation: `devtools::test(filter = "cohab|dissolution",
reporter = "summary")` passed without failures or warnings after documentation
and native exports were regenerated. Exact histories verify duration and woman's
age, independent axis-cell boundary assignment, exposure/count conservation,
woman-group labels (including numeric labels), conversion-only events,
competing death/dissolution censoring, ongoing exposure, pre-window histories,
finite last-break extension, zero exposure, empty tables and input validation.
Existing dissolution tests verify unchanged behavior of the shared accumulator.
This was offline targeted verification; no empirical rate downloads were used.


D092 extraction: `devtools::test(filter =
"unpartnered-births|parent-info|demog-union-types", reporter = "summary")`
passed without failures or warnings after documentation and native exports were
regenerated. Exact histories cover marriage/cohabitation exclusion, formation,
dissolution/conversion/death boundaries, repeated parental pairs, overlapping
coverage, unknown parents, untyped unions, preserved custom columns/group types,
row order and input independence, empty schemas and invalid references. No
simulation scheduling or empirical inputs were changed.


## Milestone 4 documentation and evidence closeout (October 29, 2026)

All three scripts relocated to `validation/milestone4/` parse. The eleven CSV/RDS
result files are byte-for-byte identical to their pre-move Git versions.
Relative Markdown file targets in changed documentation and both new validation
guides resolve, and `git diff --check` passes. No simulation, database request,
report render or test suite was started for this documentation-only cleanup.
The author's active render and final report source were left unchanged.


## Conditional HFD parity preparation (D094, October 29, 2026)

`TZ=UTC Rscript -e 'devtools::test(filter =
"hfd-parity|database-processes|status-fertility", reporter = "summary")'`
passed without failures or warnings. Synthetic local HFD text and mocked web
readers verify conditional product mi, revision/year selection, birth-order to
current-parity mapping including [4,Inf), retained annual values, omission of
open age endpoints, missing-versus-zero distinction, duplicate/invalid bands,
invalid request flags, absent products/years, RDS reuse and credential-free
provenance. Existing total-rate/database and status-fertility tests passed.
Parity tables are rejected for simulation until engine support is implemented.
No live HFD data or credentials were used and no engine changes were made.
The user's validation/milestone5.qmd was left unchanged.


## Upper-open fertility age bands (D095, October 29, 2026)

`TZ=UTC Rscript -e 'devtools::test(reporter="summary")'` passed the full suite
without failures or warnings; the opt-in live database test was skipped.
Targeted fertility-adjustment, demographic and database checks also passed.
Analytical tests cover the interior-limit adjustment f/(1-b*f), zero spacing,
zero fertility and infeasible eligibility. Execution tests verify births above
age 55 and complete fixed-seed history equality after adding a zero upper tail.
TFR tests verify truncated finite-range results, warnings for nonzero omitted
female open-band rates, no such warning for zero/unknown rates, strata and
empty results. Offline HFD fixtures verify 12- -> [12,13), 55+ -> [55,Inf),
parity bounds, unchanged annual values, missing cells and provenance. No HFD
rates were downloaded. The user's milestone5 report was not modified.


HFD conditional header correction (October 29, 2026): the user's real USA mi
file identifies the highest-order rate as m5px, not m5x. The importer and both
mock/local-text fixtures now use m1x–m4x and m5px; mapping remains parity 0–4+.
`devtools::test(filter = "hfd-parity|database-processes", reporter = "summary")`
passed without failures or warnings. No live request or credentials were used.
The earlier synthetic fixtures incorrectly mirrored the assumed m5x header;
the user's download exposed that schema error.


D096 zero-fill preparation (October 29, 2026):
`devtools::test(filter = "hfd-parity|database-processes", reporter = "summary")`
passed without failures or warnings. Tests cover missing cells mapped to zero,
unchanged known/explicit-zero rates, all-missing logical columns from text
parsing, endpoint retention, selected-year-only provenance counts and unchanged
total-rate validation. No live data was downloaded. Historical D094 checks of
NA preservation describe the superseded import policy.


## PATFR summaries and automatic dispatch (D097, October 29, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(filter =
"patfr|demography|lexis", reporter = "summary")'` passed without failures or
warnings. An independent hand-worked one-year HFD midpoint example yields 12/55
births per starting woman. Tests also verify first-birth survival under the table
probabilities, constant-hazard exposure conservation, recurrent births in the
open parity state, invariant constant one-/five-year schedules, fractional age
steps, broader parity bands, all existing TFR strata, female-only selection,
process inputs, typed empty outputs, metadata, included missing-cell propagation,
upper-open age warnings, and rejection of incomplete/overlapping age/parity
coverage or invalid progression probabilities. Direct and delegated calls agree.
No empirical HFD download or comparison with its full published PATFR was run;
source endpoint/zero-fill conventions can change that comparison. This is a
summary helper, not parity-dependent engine execution. The user's milestone5
report was left unchanged.


## Parity-dependent fertility execution (D099, October 29, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite without failures or test warnings; the opt-in live database
case was skipped. Targeted integration checks and final parity tests passed.
Coverage includes complete group/type parity partitions, edited-model validation,
finite age overlap, first-birth exemption, unadjusted preceding-parity integrals,
upper-age interior limits, top-category entry versus repeated births, broad
parity bands without expanding to every possible parity, actual maternal history
including deceased children, ignored custom parity columns, unknown fathers,
cancellations, spacing, restarts, typed fertility/conversion, group multiplier
versus explicit-schedule equivalence, and grouped zero-spacing controls.

A saved pre-feature 50,000-person grouped/previous-status run over 100 years at
seed 42 retained identical complete output without parity columns. No empirical
rates were downloaded or SOCSIM reference simulations run. SOCSIM's source was
inspected at the documented pinned commit; no source was copied/adapted.

The independent [spacing diagnostic](../validation/milestone5/check-parity-spacing.R)
uses 10 seeds, 3,000 initial women aged 20 at parity zero, no mortality and ten
years of follow-up. With zero spacing, observed means for parities 0,1,2+ are
.199760 (MCSE .001479), .299198 (.002280), .403456 (.001528), versus inputs
.2,.3,.4. With .75-year spacing, they are .199760 (.001479), .266465 (.001925),
.333460 (.000637). Exposure is reconstructed independently from birth histories.
This initialized-cohort check demonstrates a rate-recovery limitation of the
selected SOCSIM-style approximation, rather than successful scientific parity-
rate validation. Do not claim exact eligibility adjustment or milestone 5
acceptance from the passing software tests. The [parity specification](parity-fertility.md)
records the limitation and remaining scope. [Replications](../validation/milestone5/parity-spacing-replications.csv)
and [summary](../validation/milestone5/parity-spacing-summary.csv) retain the evidence.
The user's validation/milestone5.qmd remains unmodified by this work.

## Separate mortality/fertility observations and parity ASFRs (D100)

`devtools::document(); devtools::test(reporter = "summary")` passed the full
suite with no failures. The live database test remained opt-in/skipped. The
local HMD reader test emitted an environment warning because `timedatectl`
returned status 1; this was not a demographic assertion failure. The final
focused `devtools::test(filter = "fertility-rate-parity", reporter = "summary")`
passed all 29 assertions, including the subsequently added simultaneous-birth
and final-endpoint checks.

Checks cover process-specific output columns, female-only fertility, pre-window
births including deceased children, ignored custom parity columns, pre-birth
parity assignment, exact parity exposure, broad/open parity bands, zero-exposure
NA rates, complete grouped/type/period/age pooling identities, conversion and
dissolution birth-first ties, invalid boundaries, typed empty outputs, and
successive parity progression for simultaneous recorded births. Existing pooled
Lexis tests continue exercising the shared internal kernel. No simulation
behavior or fertility compensation was changed. The user-authored milestone 5
report was not edited by this work.


## Executable SOCSIM parity/spacing reference (October 29, 2026)

Built an unmodified scratch copy of MPIDR/rsocsim commit
`125450cff968bb3ff762fa99d4c926a8136062a2` (GPLv3), installed only under /tmp.
No engine source copied into Hammel, network rates downloaded, or Hammel engine
changes made. Forty independent fresh-process reference jobs completed: ten
seeds per birth-interval/experiment combination. Conditional maternal parity and
monthly exposure were reconstructed independently from raw output. Saved
[replications, summaries and reproduction scripts](../validation/milestone5/socsim/README.md)
include initialized cohorts and 100-year generational runs. All 20 generational
histories passed explicit unity-fertility-multiplier checks; output categories
remained single. The cohort exposure/event accounting is checked in the runner.

The initialized-cohort higher-parity shortfall is also present in SOCSIM. At
age [18,19) in years [80,100), nine-month-spacing means for parities 1,2,3,4+
are .272932,.284820,.186698,.054263 versus inputs .3,.4,.4,.4; all ten
replications contribute. No-spacing controls are .302024,.413259,.428620,.435636.
Both monthly-engine control offsets and composition sensitivity remain; the
reference does not imply exact equivalence with Hammel or recovery at every
parity. Earlier exploratory variable-multiplier runs were excluded, not pooled
with controlled results. See [the collaboration investigation](parity-spacing-investigation.md)
for user-reported empirical-schedule evidence and unresolved calibration choices.
The user-authored milestone5.qmd was not edited by this work.

## Cell-accounted fertility (D101, October 29, 2026)

`tests/testthat/test-fertility-cells.R` covers `spacing_method` validation and
defaults, metadata, exact unrealized-intensity cases (1-year age cells and an
actual-parity cell inside a pooled band), zero spacing,
unstratified recovery from a no-history start, pooled upper-parity recovery,
groups, union types, cohabitation, dissolution, matching, restarts and
cancellation. The grouped cases run with `options(hammel.check_cells = TRUE)`,
which rescans cell counts at every event and calibration fertility totals at
every fertility-driven search. A separate test covers the unreachable-cell
warning (earliest eligible ages, partly reachable cells, zero spacing,
unstratified schedules, delayed starts, zero intermediate parities, groups).
Tests of `calculate_fertility_rates(unrealized = TRUE)` cover exact
attribution, NA for breaks that cut a period or age year, restarts and
overall schedules split by union type, pooling actual parities into observed
bands, rejection of midpoint populations, and restored initialized-cohort
rates. Synthetic parity fixtures elsewhere wrap construction in the test helper
`quiet_unreachable()`, which silences only that advisory warning. Midpoint fixed-seed histories were compared in
full before and after the `BoundaryQueue` extraction and the cell integration
(grouped on-demand, grouped rate mode with cohabitation and dissolution,
grouped parity-specific) and are identical apart from the new
`spacing_method` metadata field. The full suite passed.

Midpoint compensation was later removed (D101 became the only method); the
midpoint fixed-seed comparison above and the method comparisons below run at
commit `89d451b`. Offline engine evidence, 10 seeds per method:
`validation/milestone5/check-cell-spacing.R` (recovery across cohort,
generational, restart, union-type and group scenarios),
`check-cell-resolution.R` (supplied wide, single-year and per-parity bands
give identical results; finer supplied bands give finer cells) and
`benchmark-cell-spacing.R LIB` (timing of an optimized installed build in LIB;
`load_all()` builds use `-O0` and are not representative). Results and the revised cell-resolution
decision are in [cell-accounted fertility](fertility-cell-accounting.md#implementation-and-validation).
No scientific acceptance threshold is set.

## Persistent fertility heterogeneity (D102, October 29, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test(reporter="summary")'`
passed the full suite with no test failures; the live database test stayed
opt-in/skipped. Final targeted verification after the native optional-argument
cleanup passed `hetfert|fertility-cells|simulation|calibration-intensity` without
warnings. The heterogeneity file has 62 assertions, including subsequent
extreme-weight removal and failed-initialization RNG checks.

Checks cover constructor/edited-model validation, independent mean-one draws
at the requested CV, zero-year/empty results, retained supplied weights,
strict-positive female validation, unused male entries, seed restoration and
unseeded advancement, exact alpha endpoints, intermediate normal innovations,
fixed-intensity weighted sampling, spacing/re-entry, restarts, deceased women,
and grouped parity/type/cohabitation/dissolution with explicit weighted
calibration scans. A dominant 1e100 weight can leave eligibility without losing
a remaining 1e-100 weight; no arbitrary clipping is introduced.

An independently built HEAD control in /tmp used 2,000 initial persons, two
groups, ten years of grouped on-demand formation (seed 92) and rate/cohabitation
formation with dissolution (seed 93), from baseline seed 91. Complete population
histories and metadata were identical to current execution with hetfert NULL.
CV-zero unit weights also retained identical output to the uniform path after
removing only the new weight column/specification metadata.

These are software and distribution checks, not an empirical acceptance claim
or a guarantee of a stationary inherited weight distribution under reproductive
selection. No online rates were downloaded; the user's milestone 5 report was
not edited.

## Disabled heterogeneity lifecycle and validation follow-up (D104)

`TZ=UTC Rscript -e 'devtools::test(filter="hetfert|simulation|fertility-cells", reporter="summary")'`
passed without failures or warnings after the D104 change. The enabled ->
disabled -> enabled regression creates daughters during the disabled interval,
checks that disabled results omit the weight column, verifies the original
population is unchanged, and confirms positive fresh weights on re-enabling.
Zero-year and empty disabled runs also omit the column. Validator checks cover
missing sex, logical all-NA placeholders in male-only/empty histories, and
incomplete female weights with explicit redraw guidance. Help and schema now
consistently describe weighted selection and the disabled-output rule. Native
weighted sampling and inheritance algorithms were not changed in this follow-up.

October 30 detailed-union attribution correction: `devtools::test(filter =
"detailed-unions|group-transition")` passed without warnings or failures.
Partner groups now refer to union formation; regression cases cover later
and return moves, exact-time priority, numeric and character labels, repeated
unions, empty histories, unknown endpoints and unchanged inputs. Existing
simulation outputs can be reanalyzed using their retained transition tables.

For baseline comparisons, an attached all-zero transition process is not
fixed-seed equivalent to omitting the process: outgoing schedule draws and
partner rescheduling can change RNG consumption without changing the model's
zero-transition distribution. Sharing a seed across scenarios therefore does
not guarantee identical paths; compare replicated demographic quantities.

## Detailed and pooled previous unions (D126, November 2, 2026)

`TZ=UTC Rscript -e 'devtools::document(); devtools::test()'` passed 6,404
assertions with zero failures or warnings and one skipped opt-in live database
test. No live rates were downloaded. Detailed inputs and outputs now use
none/separated/widowed; pooled inputs and outputs use none/dissolved. Regression
checks compare pooled inputs with explicit detailed copies for both ending
causes, age and duration clocks, grouped calendar schedules, cohabitation
probabilities and fixed-seed execution. Incomplete/mixed alternatives and
inconsistent effective clocks are rejected.

Calculator checks sum detailed events/exposure independently for every pooled
cell and compare rate and exact Poisson limits with those totals, including
marriage conversions, union-type strata, periods, sex, numeric groups, age and
duration bands, empty cells and the requirement for previous-union strata.
Documentation and help were regenerated. These are software checks, not a new
demographic validation or milestone acceptance. User-authored validation reports
and their data files were not edited.

The user-selected follow-up renames the pooling flag to `pool_dissolved` and
requires previous-union stratification rather than implicitly enabling it.
After regenerating help, `devtools::test(filter =
"previous-union-pooling|marriage-entry|marriage-duration|marriage-counts")`
passed 2,339 assertions without failures, warnings or skips. Tests include the
explicit-flag requirement and duration breaks' existing implicit stratification.

## Calendar review follow-up (November 2, 2026)

`TZ=UTC Rscript -e 'devtools::test()'` passed 6,430 assertions with no failures
or warnings and one skipped opt-in live database test. No validation report
was edited or rendered and no external rates were downloaded.

Analytical regression cases use a mother excluded by birth spacing throughout
the run and fertility changes at offsets .5, 1.25 and 2.2. Both integer and
fractional run starts recover the exact expected unrealized totals .15, .35,
.32 and .12 in run-anchored annual/partial-final windows. Only the supplied
initial birth counts as observed; no simulated births occur. Additional cases
accept complete in-run mortality, fertility, search, probability, conversion and
dissolution group/pair schedules despite incomplete unreachable periods. They
retain errors for active gaps, active incomplete age support, zero-run coverage
and malformed supplied rows. Inheritance-only boundaries retain matching
outcomes and demographic clocks in a fixed-seed comparison without births.

Engine output now binds all unrealized epoch pieces in one call rather than
repeatedly copying a growing table. Calendar boundaries avoid probability
reloads unless probabilities change, and skip intensity-tracker rebuilding
when only inheritance changes. These allocation/work reductions are established
by code inspection; no runtime speedup or empirical benchmark is claimed.

R calendar segmentation and time-coverage checking are shared by rate/model
coverage, pooled previous-category expansion and contextual dissolution
validation. Most rate validators already used `validate_calendar_rates()` before
this review; their domain-specific cross-period clock checks remain necessary.
The stale model-schedule status in AGENTS and the detailed dissolution category
in the dissolution specification were corrected. No demographic defaults or
calendar-gap rules were changed.
