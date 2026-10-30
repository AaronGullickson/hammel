# Current project handoff

Updated October 29, 2026. This is a dated status record, not authorization for
unspecified features or a standing pause instruction. The current user request
controls work. Start with [the documentation index](README.md),
[scientific contract](scientific-contract.md) and [decisions](decisions.md).

## Accepted milestones and active render

Milestones 1–5 are complete by user acceptance: mortality/fertility (D029),
union formation (D050), groups (D077), dissolution/remarriage/cohabitation
(D093), and parity-specific fertility/cell accounting/heterogeneity (D105).
Acceptance dates are October 15, 19, 27, 29 and 29, 2026 respectively.

Preserve all five user-authored report sources and rendered artifacts, their
input provenance, and explanatory comments. Milestone 1 must be reproduced at
its recorded historical commit because current execution requires a union
process. The final milestone 5 render was active when the user requested this
closeout. Do not delete or move its intermediates while rendering; do not
rerender, download rates or rerun simulations merely to update documentation.
No final-render result or new test count is claimed by this closeout.

## Next milestone: 6, group transitions

The user selected group transitions next, ahead of calendar-time inputs in
milestone 7 (D103/D105). Use "group transitions" for life-course group changes;
cohabitation-to-marriage conversion retains its existing name. Groups should
cover most categorical characteristics, with multiple types of groups an
optional direction for more complex applications. The former separate
individual-characteristics milestone is removed; milestones 8–10 retain their
numbers. See [the roadmap](roadmap.md).

Before implementation, specify transition inputs, age/duration clocks,
eligibility, timing/ties, historical group observation, ancestry implications,
newborn assignment, and invalidation/update of pending events, fertility cells
and matching contributions. A multiple-group schema/API is not selected.
Fixed newborn group assignment remains part of this discussion, not an
implemented method; custom mixed-parent assignment already exists. Completion
of milestone 5 does not authorize unresolved transition behavior.

## Current engine and fertility

- Independent C++17/Rcpp continuous-time engine with R seed isolation, tibble
  histories, exact full-pool partner draws and known-kin screening.
- Births use D101 cell accounting exclusively: intensity is input rate times
  all living cell members; assign attempts to an eligible mother. Age cells
  are at most one year, split at input edges, with actual-parity cells for
  parity schedules. Zero-rate bands form no cells; sufficiently low-rate
  upper-open bands remain pooled under the documented .001 tolerance.
- No-eligible-member intensity is dropped and reported in unrealized_births.
  This is expected lost intensity, not actual births. Adding it to measured
  rates is an accounting diagnostic. Sparse/structurally infeasible young
  high-parity cells may not recover their schedules.
- HFD conditional mi rates map birth order 5+ (m5px) to current parity 4+.
  Endpoint and missing-cell conventions are explicit (D094–D096). PATFR uses
  HFD midpoint-table rules and calculate_tfr delegates on parity bounds (D097).
- Parity and last-birth dates come from recorded maternal histories, including
  deceased children and unknown fathers, with no overriding parity column or
  imputed prehistory (D098–D099). Successful births alone advance parity.
- D102 het_fert(cv, alpha) specifies positive mean-one lognormal fresh weights
  and maternal log-scale inheritance with sqrt(1-alpha^2) normal innovation.
  Cell intensity is fixed; eligible mother selection uses relative weights.
  Matching calibration uses the same weighted individual hazards.
- Existing female persons$hetfert values are retained in enabled runs; absent
  columns are drawn under the simulation seed. Output always includes numeric
  weights when enabled, including zero-year runs. Newly generated males get NA.
  Alpha one copies maternal weights exactly; alpha zero draws independently.
  Reproductive selection can change cohort weight distributions.
- Under D104, disabled runs omit the weight column even if it was supplied.
  Input populations are unchanged. Subsequent re-enabling draws fresh weights
  for all recorded women, rather than keeping a partially populated column.
- D028/D095's interior-limit adjustment and D099's spacing correction were
  superseded by D101. The former spacing_method argument was removed.
  Midpoint comparison diagnostics require recorded commit 89d451b.

## Other implemented processes

Female rate-mode searches optionally stratify by previous_union none/dissolved/
widowed, with age or duration per previous category (D080/D083). On-demand
formation remains birth-triggered. Dissolution uses one female-carried duration
clock with symmetric group-pair schedules/multipliers. Optional rate-mode
cohabitation supports cohab_prob, separate matching, type-specific fertility
and dissolution, and age/duration union_convert_rate. Conversion links episodes
through previous_cohab_id and resets marriage dissolution duration. Exact ties
prioritize death, birth, dissolution, conversion, search.

Newborn identity and ancestry are separate; random_parent with maternal
fallback is the default. Custom inheritance is supported. Cohabitation endings
recover earlier legal status by default, with a toggle to count previous unions
instead. Unpartnered father matching can use a separate model. Failed matches
cancel births by default; pending spacing and weighted-cell state remain
consistent after cancellations.

## Observation APIs and evidence

- calculate_fertility_rates: female exposure/births, optional parity_breaks,
  by_group/by_union_type, and unrealized-intensity accounting.
- calculate_mortality_rates: deaths/exposure for both sexes, optional group/type.
- calculate_union_rates: entries from unpartnered exposure, excluding conversion;
  optional previous_union and duration clocks, group and entry type.
- calculate_marriage_rates: direct marriage plus conversion over nonmarried
  exposure, optional previous_union and group; age only.
- calculate_dissolution_rates: union-years/dissolutions by optional group_pair/type.
- calculate_cohab_convert_rate: conversion over cohab union-years, age or duration.
- get_unpartnered_births: known-parent births outside partnerships, preserving
  repeated father-selection events.

The accepted [milestone 5 report](../validation/milestone5.qmd) revisits baseline,
restart and union-type fertility, then assesses parity and heterogeneity,
childbearing dispersion, mother–daughter association and runtime. Software
checks are recorded separately in [testing](testing.md). Supporting scripts
and their version restrictions live in [validation/milestone5/](../validation/milestone5/README.md).
The pinned SOCSIM source audits and narrow executable parity/spacing and
heterogeneity checks are reference evidence, not broad engine equivalence.

Retain sparse-cell, infeasible-input, history completeness, reproductive
selection, failed-search and approximate matching-bin limitations. Restarts
are not exact checkpoints. Poor-pool waiting (D044), former-parent continuity
(D042), general two-sex squeeze models, custom process execution and calendar-
time schedules remain future work. The postponed African historical benchmark
has not been selected. Scientific acceptance remains distinct from passing
software checks or universal accuracy.
