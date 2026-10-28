# Current project handoff

Updated October 28, 2026. This is a dated status record, not authorization for
additional features. Start with the [documentation index](README.md); the
current user request controls work.

## Implemented state

- Call-local C++17/Rcpp continuous event scheduling, R RNG isolation, tibble
  interfaces, and recorded population/union histories.
- Mortality and fertility with shared hazard sampling, per-schedule midpoint
  compensation, and actual-birth spacing preserved across restarts.
- Mandatory union formation with on-demand or female rate-based initiation,
  exact full-pool conditional-logit partner selection, minimum age and known-kin
  screening, exclusive ongoing unions, widowhood, and status-event invalidation.
- Overall or partnered/unpartnered fertility schedules. Rate-mode nonmarital
  births select fathers without forming unions; a separate `father_matching`
  model may override general matching. Failed matches cancel births by default.
- Group-specific schedules and proportional hazard multipliers, including
  separate fertility factors by union status, plus group-stratified observations.
- Weighted choice calibrated to mirrored male-group per-capita targets with
  optional delta adjustments. Symmetric pair effects use `A:B`; directed pair
  utilities, `::`, candidate sample sizes and method selection are removed.
  Default age-gap matching is `create_age_gap_matching_model(0, 5)`.
- `diagnostics = TRUE` records applied `.log_weight_<group>` offsets on new
  unions. There are no `.n_eligible` or `.utility_rank` columns.
- Separate ancestry proportions and newborn identity. Built-in and custom
  inheritance support possible future identities and custom mixed categories.
  Default inheritance is `random_parent`, with maternal fallback (D076).
- Shared native kinship screening and `extract_kin()`, parent/spouse extraction,
  and Lexis summaries with group and overall/partnered/unpartnered observations.

See [design](design.md), [scientific contract](scientific-contract.md),
[group inheritance](group-inheritance.md), and
[weighted choice](weighted-choice-matching.md) for current rules.

## Validation and ownership

Milestones 1, 2, and 3 were accepted by the user on October 15, October 19,
and October 27, 2026 respectively (D029, D050, D077). Preserve all three
validation source reports, rendered artifacts, and the user's explanatory
comments. Milestone 1 is historical and must be reproduced at its recorded
commit because current execution requires a union process. See [testing](testing.md).

Milestone 3 covers group-specific and multiplied schedules, weighted matching
and delta, wife-side conditional-logit recovery with log-weight offsets,
ancestry and built-in/custom inheritance, starting-population sensitivity,
and performance. Inheritance runs regenerate baselines per replication and
reduce maximum starting age to 50; both changes were made together. Earlier
fixed-baseline MCSEs condition on their shared populations. The experiments did
not establish a calibration bug or a systematic 3:1 group-size effect.

The user revised the report and added its conclusion, then requested closeout
while rerendering. The revised source is at
`7514916bb8cf0c06faef24dad7fa013304093779`. This documentation update does not
independently rerender reports, download rates, or rerun simulations/tests.
No final HTML hash or combined acceptance commit is recorded. Do not interfere
with the user's rendering or modify the report as part of this closeout.
Acceptance is the user's scientific assessment and does not imply that every
review suggestion was implemented or establish universal accuracy guarantees.

## Active milestone 4 investigation

The user created `milestone4-union-dissolution` and requested an audit of SOCSIM
union dissolution and remarriage. The
[static audit](socsim-dissolution-audit.md) traces the same pinned clean rsocsim
checkout used by earlier audits. It finds duration-based divorce candidates
for either spouse, ordinary marriage-search schedules selected by divorced or
widowed status, optional duration-based remarriage, inherited schedule defaults,
partner rescheduling, cohabitation histories and conditional former-partner
exclusion. No reference simulation was run or source copied.

After this audit, the user selected duration-dependent annual hazards with
one female-carried clock per union and death/birth/dissolution exact tie
priority (D078). D079 replaces woman-only grouping with symmetric group-pair
schedules and multipliers using `group_pair` labels such as A:A, A:B, B:B.
The optional dissolution process is the first implementation: `dur_start`/
`dur_end` bands, one retained clock per union, status invalidation and availability
updates, plus latest-run dissolution counts. `calculate_union_counts()` also
returns period/duration dissolution rates with union-year exposure, optionally
by symmetric group pair; each episode contributes once. The user checked the
dissolution demonstrations and moved on to age-specific remarriage (D080).
Optional `previous_union` on search tables uses `none`, `dissolved`, `widowed`,
with all categories required per group; absent columns share rates. Calibration
sums category-specific intensities within group/age bins. State is recovered
from complete exclusive histories and updated on endings.
`calculate_marriage_counts()` now returns unpartnered exposure, starts and rates
by sex, optionally by person group. `by_previous_union = TRUE` adds previous
categories; the default FALSE pools counts and exposure and omits that column. The user requested moving
`unions`/`rate_marriage` out of Lexis outputs entirely; current-status birth,
death and exposure summaries remain. Duration-dependent remarriage and
cohabitation stay deferred. Milestone 4 is not accepted as complete. See
[union dissolution](union-dissolution.md) and [remarriage](remarriage.md).

The reported grouped rate-mode slowdown was largely debug compilation:
`devtools::load_all()` defaults to `debug = TRUE` (`-O0`). The user's optimized
installed runs took 0.568 seconds without groups and 5.302 seconds with two
groups and dissolution, for 50,000 initial persons over 100 years. Those are
user-supplied timings, not a controlled performance study. For development
timing use `devtools::load_all(debug = FALSE, recompile = TRUE)`.

## Retained limitations and next scope

Transition-aware fertility compensation (D045), poor-pool waiting (D044), and
previous-father continuity (D042) remain deferred. Sparse exposure cells,
initialization transients, approximate calibration ages and omitted calibration
kinship remain documented. Exact partner draws still screen individual age and
known kinship. Small or sex-imbalanced markets may yield unexpected directional
ratios, and fixed delta does not guarantee pairwise ratios over time.

The current female-driven model does not represent marriage-squeeze feedback.
Full two-sex models remain a future alternative. The [roadmap](roadmap.md)
records active dissolution/remarriage work as milestone 4 and later parity,
time-varying inputs, group transitions, summaries, migration and characteristics.
None is authorized by closing milestone 3.

The user deferred fixed newborn assignment (for example, minimum educational
category at birth) to the group-transitions discussion. No fixed inheritance
method was added. Custom inheritance already permits mixed-parent births to
receive a specified category. Multiple simultaneous identities, arbitrary state
callbacks, exact checkpoints, and initialization correction remain future work.
The African historical benchmark remains postponed.
