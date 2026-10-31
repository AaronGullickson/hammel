# Current project handoff

Updated October 31, 2026. This is a dated status record, not authorization for
unspecified features or a standing pause instruction. The current user request
controls work. Start with [the documentation index](README.md),
[scientific contract](scientific-contract.md) and [decisions](decisions.md).

## Accepted milestones and active render

Milestones 1–6 are complete by user acceptance: mortality/fertility (D029,
October 15), union formation (D050, October 19), groups (D077, October 27),
dissolution/remarriage/cohabitation (D093, October 29), parity-specific
fertility/cell accounting/heterogeneity (D105, October 29), and individual
group transitions (D116, October 30), all in 2026.

Preserve all six user-authored report sources and rendered artifacts, input
provenance and explanatory comments. Milestone 1 must be reproduced at its
recorded historical commit because current execution requires a union process.
The user is compiling the final milestone 6 report at this closeout. Preserve
`validation/milestone6.qmd`, `.html`, and all active render intermediates;
do not move, remove, rerender or interrupt them for documentation work. Actual
replication counts and population settings belong to the report. Acceptance
is explicit, but final-render completion and numerical results are not
independently certified by this documentation update.

## Milestone 7: calendar-varying rates

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
and historical-run setup to fix the getter first. This is implementation status,
not scientific acceptance of milestone 7.

D121 keeps milestone 7 on the current engine and assigns the future
[binned-group hazard/formula-weight change](binned-group-hazards.md) to milestone
8. It remains unimplemented, with detailed formula and mechanism contracts
still to specify. Summaries/helpers move to milestone 9. Multiple independent
group columns are tentative milestone 10, with no selected schema. Fixed newborn
assignment remains unselected; custom mixed-parent assignment already exists.
The standalone regions/internal-migration milestone is removed; regional and
household applications may be considered under tentative milestone 10 (D121).
D114's rate-interpretation, household-membership and conditional-rescheduling
questions are preserved in the roadmap. No household mechanism is implemented.

## Individual transitions and observations

`group_transition(rates)` and `set_group_transition()` implement person-level
competing hazards, with required `group_origin`, `group_dest`, and `rate`, and
one common age or duration axis across the table. Optional `sex` and
`union_type` stratify risks; absent dimensions apply universally. Omitted pairs
and gaps have zero risk. Duration runs from latest group entry, or birth without
recorded entry history. Moves preserve ancestry, parentage, parity, spacing and
individual fertility weights while updating affected clocks, matching pools,
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

Births use D101 cell accounting exclusively: a cell generates attempts at its
input rate times all members, assigned to an eligible mother uniformly or by
persistent hetfert weights (D102). Cells distinguish actual parity and group/
union type. No-eligible-member intensity is dropped and reported as
unrealized births, an accounting diagnostic rather than actual births.
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
on-demand searches remain birth-triggered. Dissolution uses one female-carried
union-duration clock, with symmetric pair effects. Optional cohabitation and
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
external immigration (tentative milestone 11), a two-sex marriage model
(tentative milestone 12, D117/D121), custom process execution and the postponed
African historical benchmark remain deferred.
Milestone 12 assigns a planning priority only; model and API selection remain
open in the [two-sex research notes](two-sex-union-models.md).
