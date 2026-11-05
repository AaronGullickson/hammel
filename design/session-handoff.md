# Current project handoff

Updated November 5, 2026. This is a dated status record, not authorization for
unspecified features or a standing pause instruction. The current user request
controls work. Start with [the documentation index](README.md),
[scientific contract](scientific-contract.md) and [decisions](decisions.md).

## Accepted milestones and active render

Milestones 1–8 are complete by user acceptance: mortality/fertility (D029,
October 15), union formation (D050, October 19), groups (D077, October 27),
dissolution/remarriage/cohabitation (D093, October 29), parity-specific
fertility/cell accounting/heterogeneity (D105, October 29), and individual
group transitions (D116, October 30), calendar-varying inputs (D127,
November 2) and cell scheduling with weight models (D130, November 5), all in
2026.

Preserve all eight user-authored report sources and rendered artifacts, input
provenance and explanatory comments. Milestone 1 must be reproduced at its
recorded historical commit because current execution requires a union process.
Milestone 6 had an active final render at its acceptance. Preserve
`validation/milestone6.qmd`, `.html`, and all active render intermediates;
do not move, remove, rerender or interrupt them for documentation work. Actual
replication counts and population settings belong to the report. Acceptance
is explicit, but final-render completion and numerical results are not
independently certified by this documentation update.

## Milestone 8 accepted complete

The user accepted milestone 8 on November 5, 2026 (D130), based on
`validation/milestone8.qmd` and its rendered HTML. Every rate process is
scheduled by cells with one clock per process (D128): no individual
waiting-time clocks remain. Fertility alone splits to single-year and
actual-parity cells; cells stratify by group only when rates do; dissolution
and conversion units are unions; unrealized births are counts; calendar rate
changes update cells in place when the layout allows; and matching calibration
reads the cells (the former calibration tracker is gone). Weight models
(`create_weight_model()`, formulas or custom functions, optionally scheduled)
choose who in a cell experiences an event (D129). See
[cell scheduling and weight models](binned-group-hazards.md), including its
deferred issues. The user plans, after the milestone, to add `set.seed()` to
the articles; `group_prop_hazard` has since been replaced by
`stratify_rates()` (D132). The user then reordered the roadmap (D131): multiple group columns are
next as milestone 9, external immigration is 10, and summaries, diagnostics
and helpers move to 11. Milestone 9's schema and mechanisms are still to
specify.

## Milestone 7 accepted complete

The user approved calendar input semantics and authorized the six annual-rate
inputs first (D118). Rates now accept flat `time_start`/`time_end` bands;
[the contract](time-varying-rates.md) records coverage, fractional years,
boundary precedence, resampling, histories and observations. D123 implements
[matching/inheritance model schedules](time-varying-models.md) for all four
arguments, with a shared `time_start`/`time_end`/`model` interface, dynamic NULL
matching fallbacks and a registry covering all future inheritance identities. D122 adds calendar bands
to the flat `cohab_prob` table, with current probabilities selected at searches
and matching-calibration updates that preserve event clocks at probability-only
boundaries. Other arguments remain static. D119 adds multi-year HMD/HFD getters and native HMD shape selection;
D120 adds optional `open_age` pooling from native deaths/exposure products;
see [database preparation](database-rates.md). The user deferred stationary-year
and historical-run setup to fix the getter first. The user has now accepted
milestone 7 complete (D127), based on `validation/milestone7.qmd` and its retained
HTML assessment. The code/report revision at acceptance is
`3b65d2202b006ef6f3e232c98595814b1544f9fe`. This closeout does not rerender or
independently reproduce the report.

Milestone 7 was completed with the previous engine (D127); milestone 8
then carried out D121's [cell scheduling and weight models](binned-group-hazards.md).
Under D131, multiple independent group columns are milestone 9, with no
selected schema, and summaries/helpers are milestone 11. Fixed newborn
assignment remains unselected; custom mixed-parent assignment already exists.
The standalone regions/internal-migration milestone is removed; regional and
household applications may be considered under milestone 9 (D121, D131).
D114's rate-interpretation, household-membership and conditional-rescheduling
questions are preserved in the roadmap. No household mechanism is implemented.

D126 uses detailed previous-union categories `none`/`separated`/`widowed`, or
pooled `none`/`dissolved`, consistently for search and probability inputs.
Union/marriage calculators use `by_previous_union = TRUE` with optional
`pool_dissolved = TRUE`; confidence limits use combined counts/exposure.
Recent review fixes retain run-anchored unrealized-birth accounting and restrict
execution group-coverage checks to active run slices. The last recorded full
offline suite passed 6,430 assertions with no failures or warnings and one
opt-in live test skipped; the subsequent `withr` test dependency declaration
passed R's dependency scan and 133 focused assertions. These are local checks,
not a claim that GitHub CI was rerun. See [testing](testing.md).

Validation data now live in country/period subfolders; see the
[data guide](../validation/data/README.md). Earlier report sources were not
rewritten for that relocation; preserve historical reproducibility at their
recorded revisions.

## Individual transitions and observations

`group_transition(rates)` and `set_group_transition()` implement person-level
competing hazards, with required `group_origin`, `group_dest`, and `rate`, and
one common age or duration axis across the table. Optional `sex` and
`union_type` stratify risks; absent dimensions apply universally. Omitted pairs
and gaps have zero risk. Duration runs from latest group entry, or birth without
recorded entry history. Moves preserve ancestry, parentage, parity, spacing and
individual fertility weights while updating affected cells, matching pools,
calibration, fertility cells and union hazards. See [the contract](group-transitions.md).

Grouped outputs always include a typed `transitions` tibble, even with no
process or moves. Columns are `person_id`, `time`, `group_origin`, `group_dest`.
Histories persist through restarts and support group exposure/event attribution.
`calculate_transition_rates()` returns directed Lexis counts over origin
person-years by age or duration, optionally sex and union type. Outgoing
destinations share origin exposure; do not sum those denominator copies.
`get_detailed_unions()` reconstructs groups at formation, not final identity.

Grouped population snapshots use sorted births/deaths/arrivals/departures,
rather than repeated per-person lookups. Dissolution and conversion group
histories are built natively using indexed partner changes. Developmental
performance measurements and regression evidence are in [testing](testing.md).
An attached all-zero transition process does not promise fixed-seed equality
with omitting the process; compare replicated demographic outcomes.

## Current engine and fertility

Each simulation owns independent C++17/Rcpp state and executes continuously,
with R seed isolation, ordinary tibble results, exact full-pool partner draws
and known-kin screening. Public dates/durations are years, hazards per
person-year and exposure person-years. There is no monthly execution grid.

Every rate process uses cells and one clock per process (D128). Births use
D101 cell accounting: a cell generates attempts at its input rate times all
members, assigned to an eligible mother uniformly or by persistent hetfert
weights (D102) combined with any weight model (D129). Cells distinguish actual
parity, union type and, when schedules are group-specific, group. A draw of a
cell with no eligible member is counted as an unrealized birth, an accounting
diagnostic rather than an actual birth.
Parity and spacing derive from retained maternal histories, without imputed
prehistory. Preserve sparse/infeasible-cell limitations and the reproductive
selection caveat for inherited weights. Disabled heterogeneity removes the
weight column (D104); re-enabling draws fresh weights if absent.

D028/D095 midpoint compensation and D099's former spacing correction are
superseded. Do not reintroduce them. HFD endpoint, parity and missing-cell
rules remain D094–D100; see [parity fertility](parity-fertility.md),
[cell accounting](fertility-cell-accounting.md) and
[heterogeneity](fertility-heterogeneity-notes.md).

## Other implemented processes and retained limits

Rate-mode searches use optional previous-union age/duration strata;
on-demand searches remain birth-triggered. Dissolution puts each union in
union-duration cells, with symmetric pair effects. Optional cohabitation and
same-partner conversion retain linked union episodes and distinct observation
risk sets. Ancestry shares remain separate from current identity and average
at birth. Consult the individual contracts rather than inferring behavior
from legacy SOCSIM or future roadmap notes.

Keep the accepted reports, source comments and provenance. Use
`devtools::test()` for ordinary offline checks; database tests need explicit
opt-in. Documentation closeout requires no new downloads, simulations or
package installation. Prior software checks, coverage and timings are recorded
in [testing](testing.md); they are distinct from the user's scientific
acceptance and do not establish universal accuracy or SOCSIM equivalence.

Restarts are not exact checkpoints. Sparse-cell feasibility, incomplete
histories, failed searches and matching-calibration approximations remain
limitations. Poor-pool waiting (D044), previous-father continuity (D042),
external immigration (tentative milestone 10), a two-sex marriage model
(tentative milestone 12, D117/D121), custom process execution and the postponed
African historical benchmark remain deferred.
Milestone 12 assigns a planning priority only; model and API selection remain
open in the [two-sex research notes](two-sex-union-models.md).
