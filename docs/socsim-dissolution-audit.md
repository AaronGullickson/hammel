# SOCSIM union dissolution and remarriage: static source audit

Inspected October 27, 2026 for Hammel milestone 4. Reference:
[MPIDR/rsocsim](https://github.com/MPIDR/rsocsim/tree/125450cff968bb3ff762fa99d4c926a8136062a2),
commit `125450cff968bb3ff762fa99d4c926a8136062a2` (commit subject
“submission of 1.9.18”). Local checkout: `../rsocsim`, confirmed at that commit
with a clean working tree. That path is a convenience, not a Hammel dependency.
Repository `LICENSE.md` contains GPL version 3. No source was copied or adapted;
this document describes behavior independently. This is static inspection of
the ordinary non-ENHANCED engine, not an executable reference experiment or
an audit of every SOCSIM release, extension, or polygamous configuration.

## Main findings

| Question | Inspected behavior |
| --- | --- |
| What is the divorce time scale? | Duration since the current union began, in integer months; not either spouse's age |
| Who can carry the divorce event? | Either married/cohabiting partner, using their own group/sex/status schedule |
| Is remarriage a different event? | No: ordinary marriage search, selecting schedules by current single/divorced/widowed status |
| What is the remarriage time scale? | Age by default; optionally elapsed time since the latest union ended |
| Is there a built-in post-divorce waiting period? | None in these handlers; zero-rate duration bands can supply a wait |
| What happens to pending events? | The other partner's event is removed; both partners receive new candidates after dissolution |
| Can the former couple reunite? | A former-partner exclusion is attempted, but missing-parent early exits and a duplicated guard limit it |
| Is marriage order a standard schedule dimension? | Not in the inspected marriage/divorce rate indexing; history is retained separately |

## Source map

All paths below are relative to the pinned rsocsim checkout. Line references
are to that commit, not current upstream main.

| File and locations | Evidence |
| --- | --- |
| `src/src/defs.h:267–305,575–631` | Relationship statuses, end reasons, person/union history links |
| `src/src/load.cpp:691–704,1062–1082` | Schedule keys and duration-specific directive |
| `src/src/load.cpp:1644–1698,1924–1932,2109–2247` | Input conversion, numerical zero, schedule inheritance |
| `src/src/events.cpp:1510–1565,1568–1672` | Rate modification and integer-month candidate sampler |
| `src/src/events.cpp:1832–1995` | Event competition, marriage time scale, duration-based divorce |
| `src/src/events.cpp:1016–1082,1214–1223` | Dispatch and initial one-queue male availability |
| `src/src/xevents.c:178–340` | Death, divorce and cohabitation breakup updates |
| `src/src/xevents.c:342–568` | Search failure, union creation and prior-history links |
| `src/src/xevents.c:907–1000,1005–1099` | Former-partner screen and ordinary scoring |
| `src/src/utils.c:1447–1511,1574–1590` | Previous legal status and union-count traversal |
| `src/src/io.c:359–386,982–1003` | Reading/writing union dates, end reasons and prior IDs |
| `R/read_results.R:9–17` | R output documentation; not the native scheduling contract |
| `socsim_oversimplified.md`, “Divorce rates” and “Duration specific transition rates” | Documentation corroborating duration-based divorce |

Key primary sources: [candidate scheduling](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp#L1946-L1995),
[dissolution handler](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L273-L340),
[rate defaults](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/load.cpp#L2109-L2247),
and [reference documentation](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/socsim_oversimplified.md).
The local pinned files were inspected directly; online retrieval was unavailable.

## Dissolution eligibility, inputs and timing

`date_and_event()` considers divorce whenever the person is `MARRIED` or
`COHABITING`. It selects a rate block by the initiating person's group, sex
and current status, then passes elapsed time since `last_marriage->date_start`
to `datev()`. There is no age/duration switch in this divorce branch: duration
is always the time scale. The generic parser's `duration_specific` flag does
not change that fact. No union-order or age-gap dimension is present in this
ordinary rate lookup; ENHANCED hooks can alter behavior and are outside scope.

A husband aged 50 and wife aged 35 whose union began five years ago therefore
both consult their respective divorce schedules at duration 60 months. They
do not use the age-50 and age-35 entries as ages. This is an illustration of
the source rule, not a simulated example.

Ordinary divorce inputs use the nonbirth conversion `log(1-p)` and integer
month endpoints. With no modifiers, a constant input p implies a candidate
probability p in the first month. Modifiers multiply the log survival term,
not the raw probability. Explicit zeros are replaced by small positive values;
missing married divorce schedules for group 1 use a near-zero block. These
are numerical conventions, not exact zero hazards. See the
[rate audit](socsim-rate-audit.md) for the shared conversion/sampler trace.

### Two spouse candidates are two competing dissolution risks

The divorce branch is not restricted to women, including in one-queue mode
(the one-queue restriction concerns marriage searches). Both partners can
carry a candidate. If either executes divorce, the handler ends their shared
union and removes the other partner's pending event before rescheduling.
This prevents the old second candidate from dissolving that union twice.

Consequently, putting the complete desired union dissolution schedule on
both sexes compounds the risk rather than duplicating the same observation.
For independent constant candidate streams with p_F = p_M = 0.01, ignoring
other events, boundaries and numerical floors, the probability that at least
one divorce candidate falls in the first month is
`1 - (1 - p_F) * (1 - p_M) = 0.0199`. This is an analytical consequence of
the two-clock code, not measured reference output. The reference documentation
likewise says divorce schedules generally need not be supplied for both sexes.

### Competition and within-month order

The person-level candidate sequence is death, birth (women), marriage, divorce,
then group transitions. Replacement requires a strictly earlier month, so an
equal-month divorce candidate loses to the earlier candidates for that person.
The monthly queue processes different people in random order, not by one global
event-type priority. Dates can be the current month, and rescheduling can
produce another event in that month. See [competing events](socsim-competing-events-audit.md).

## What dissolution changes

The ordinary `divorce(p)` handler operates on `p->last_marriage`:

1. Identify the other partner and set the union's end month to the current month.
2. Remove the other partner's pending event; the focal person's event was
   already removed by the dispatcher.
3. For marriage, record end reason `DIVORCE` (2) and set both statuses to
   `DIVORCED`. For cohabitation, record `BREAKCOHAB` (4) and recover each
   person's previous legal marital status from their own history.
4. Resample and enqueue each partner's next event under their updated status.
5. Increment the monthly divorce-event count once. The same counter includes
   cohabitation breakups in this handler.
6. In one-queue mode, immediately return the male partner to the available
   marriage queue, recording the current month as his queue-entry date.

Both old death and fertility/search candidates may be redrawn, because SOCSIM
can select rates by legal status for these processes. Birth spacing still uses
the latest recorded child's birth month; dissolution does not erase it. Parent
links and past union records remain. No custody, child relocation, pregnancy
state, or former-father preference is implemented by this handler.

The divorce handler operates only on the latest union, whereas death traverses
all ongoing records. Do not infer validated polygamy support from the ordinary
monogamous path. There is no independent separation state here: legal statuses
are single (1), divorced (2), widowed (3), married (4), and cohabiting (5).

## Remarriage schedules and availability

There is no separate remarriage event. After divorce or widowhood,
`date_and_event()` uses `E_MARRIAGE` indexed by the person's new legal status.
This permits separate search schedules for never-married, divorced and widowed
people. Searches use the ordinary pool, eligibility checks and scoring; they
are not guaranteed realized marriages when suitable partners are absent.

Marriage search uses attained age unless its group/sex/status rate block is
marked duration-specific. For a duration-specific block:

- A cohabiting person's clock starts at the current cohabitation's start.
  Their marriage event legalizes that union with its existing partner.
- Other statuses use the latest union record's end date. For a divorced or
  widowed person with an appropriate history, this is time since divorce or
  widowhood. The code dereferences `last_marriage` directly, so this mode
  requires coherent history; it does not infer a date for a history-free person.

For example, a woman divorcing at age 40 can immediately enter the age-40
remarriage schedule, or duration zero if that schedule is duration-specific.
A rate block with no effective risk for the first 12 duration months can model
an intended waiting period, subject to SOCSIM's near-zero convention. There is
no separate mandatory cooldown in the inspected handler. In one-queue mode,
men become candidates immediately and women follow their search schedules;
in two-queue mode, both sexes enter the pool through ordinary search/queue
rules. A failed search preserves queue availability and resamples the searcher's
next event as described in the [formation audit](socsim-union-formation-audit.md).

Successful remarriage appends a new record with prior-union links for both
partners, sets a new start month and ongoing end reason, changes their statuses,
and reschedules both. The basic schedule arrays do not distinguish second
from third marriage. `get_marity()` counts records, including cohabitations;
the suggested same-union-count scoring bonus in ordinary `score3()` is commented
out, not active remarriage-specific matching.

### Missing schedules can silently reuse first-marriage schedules

For group 1, missing marriage schedules for divorced people inherit the
single schedule of the same sex; widowed schedules inherit divorced schedules.
The loader deliberately skips defaulting the married marriage schedule, which
would otherwise permit an extra union. Missing cohabiting marriage schedules
inherit the preceding married entry, which can remain NULL. Higher groups
inherit missing schedules from the preceding group at the same sex/status.
Duration-specific flags propagate along inherited entries.

Divorce schedules for group-1 married men and women are initialized separately
to the near-zero block when absent: female schedules are not automatically
copied to men. Cohabitation breakup schedules can inherit the married divorce
schedule. Thus omitting a remarriage schedule need not disable remarriage, and
omitting divorce for one sex need not duplicate the other sex's risk.
These legacy defaults are source findings, not proposed Hammel defaults.

## Widowhood, cohabitation and former partners

Death ends an active marriage with reason `DEATH` (3), makes the surviving
spouse `WIDOWED`, and reschedules their event. A surviving man returns to the
one-queue pool. Cohabitation death uses `COCOHABDIES` (6) and restores previous
legal status; legalization ends the cohabitation with `MARAFTERCOHAB` (5) and
creates a new marriage record. The linked history preserves both episodes.

After cohabitation ends, `get_previous_mstatus()` walks backward to recover
single/divorced/widowed status, or married if another ongoing record is found.
The subsequent duration-specific search still uses the latest record's end,
which is the cohabitation end, not necessarily the earlier divorce/widowhood
that supplied the recovered status. Status and elapsed-time origin therefore
need not refer to the same historical episode.

`check_incest()` also attempts to prohibit former partners. It returns early
if either mother is unknown, bypassing that exclusion. Its history guard tests
`p1->last_marriage` twice rather than checking both people. When reached, it
walks the male partner's union history and rejects a candidate previously
recorded as his wife, without filtering union type or end reason. This is an
incomplete and conditional exclusion, not a reliable universal prohibition on
reconciliation. The direct control flow is verified by inspection; its
population-level effects have not been tested.

Native union output retains start/end months, end reasons and each partner's
prior-union ID. R's `read_omar()` documentation lists divorce and death reasons
but not all cohabitation codes; use the native constants for this audit.

## Implications and open choices for Hammel

After the audit, the user selected duration-dependent annual hazards and one
female-carried clock per union (D078), superseding the earlier female-age-based
direction. Group effects now use symmetric unordered pairs (D079), while the
woman still carries the single union clock. Birth precedes dissolution at an
exact same-person tie. See the
[dissolution specification](union-dissolution.md). SOCSIM's duration baseline
informed this decision; its monthly inputs and two-spouse option were not adopted.
Hammel now supports age- or duration-specific search schedules by previous
union in one table (D080/D083); cohabitation remains deferred.

Remaining design areas beyond the accepted dissolution core include:

- Remarriage strata: first versus subsequent union, divorced versus widowed,
  their rate coverage/defaults, and age versus elapsed-duration schedules.
- Whether cohabitation belongs in milestone 4; its formation, conversion,
  breakup and legal-history implications need a separate explicit contract.
- Eligibility after dissolution, including reconciliation with a former partner
  and any waiting interval; nothing in this audit selects a restriction.
- State invalidation for each former partner, search/father pools, fertility
  status, same-time ties, spacing, initial histories and restarts. On-demand
  mode also needs an explicit relationship between birth-triggered reformation
  and any separate remarriage schedule.
- Observation rules: count one dissolution per union, with sex-specific exposure
  and event attribution if requested, preserve end reasons and histories, and
  distinguish input search rates from realized remarriages.

## Suggested validation, not yet executed

For a reference experiment, use the pinned source and explicit schedules to
avoid inherited defaults. Start with known monogamous unions and disable other
risks as far as the reference's numerical conventions permit. Test:

1. One nonzero spouse schedule versus both: compare constant-block candidate
   survival with the analytical combined risk above across replications.
2. Equal union duration at different spouse ages, then equal spouse ages at
   different union durations: establish the divorce time origin operationally.
3. Divorce versus spouse death: verify divorced/widowed histories, closure once,
   invalidation of the other person's old event, and one event count per union.
4. Explicit divorced/widowed search rates versus omitted rates; age versus
   duration schedules and an initial zero-risk duration interval.
5. Remarriage with an available candidate versus none; one- versus two-queue
   availability and repeat same-month events.
6. Former-partner screening with known/unknown mothers and asymmetric history
   availability; cohabitation breakup, legalization and death, including an
   earlier divorce followed by cohabitation.
7. Read/write/restart preservation of dates and prior-union links, distinguishing
   resampled event schedules from exact stochastic continuation.

These cases also suggest Hammel invariants once its scientific specification
is agreed. No reference simulations, fitted comparisons or scientific acceptance
thresholds are claimed here. Static findings do not establish SOCSIM accuracy.
