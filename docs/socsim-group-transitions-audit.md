# SOCSIM group transitions audit

October 29, 2026. Static audit of MPIDR/rsocsim 1.9.18, commit
`125450cff968bb3ff762fa99d4c926a8136062a2`, inspected in the clean local
`../rsocsim` checkout. Repository metadata declares GPL-3 and the root
`LICENSE.md` contains GPLv3. No source was copied or adapted into Hammel.
This audit covers the ordinary non-ENHANCED engine. No reference simulation
was executed for this audit; source defects and their possible consequences
below have not been measured in population experiments. No Hammel transition
semantics or API are selected by these findings.

Pinned sources (line numbers below refer to this commit):

- [Event scheduling](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp)
- [Event execution and household assembly](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c)
- [Input loading and defaults](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/load.cpp)
- [Multipliers and target calculation](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/utils.c)
- [Initialization and history I/O](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/io.c)
- [Guidebook](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/socsim_oversimplified.md#transition-rates)

## Inputs and eligibility

SOCSIM stores one current integer group per person. A transition changes that
group while retaining the person's identity, parents, and union records. It is
an internal move between modeled groups; the execution path neither creates
nor removes a person from the total population. Groups can represent locations
or other categories. The core transition path does not maintain Hammel-style
ancestry proportions or separate multiple categorical dimensions.

`load.cpp` lines 1200–1216 parse `transit <origin> <sex> <marital_status>
<destination>`, followed by ordinary rate blocks. The destination can establish
the existence of a group. Transition schedules are indexed by current origin,
destination, sex and marital status. Both sexes can transition; there is no
general adult-only, unmarried-only, or minimum-age gate in the transition
candidate loop. Eligibility is determined by whether the relevant schedule
exists and its rates. In `fill_rate_gaps()` the ordinary non-birth defaults
loop only over marriage through death (lines 2177–2245), excluding transitions.
Missing transition schedules remain absent rather than inheriting another
group's or marital status's schedule. The guidebook agrees, although its
example misspells `transit` as `trasit`; the parser accepts `transit`.

`add_rate_block()` (lines 1643–1698) interprets transition input as a monthly
probability p, storing log(1-p). Explicit zero is replaced by 1e-12, so an
explicit zero block differs from an absent schedule. These inputs must not be
read as annual person-year hazards. For example, a constant annual hazard 0.12
has an equivalent monthly probability 1-exp(-0.12/12), approximately 0.00995017,
before other modifiers and monthly event competition.

## Clocks and competing destinations

`events.cpp::date_and_event()` (lines 1998–2036) independently samples a
candidate for each destination with a non-null schedule, using `datev()`.
It competes these with the person's death, birth, marriage and divorce
candidates, retaining only one selected event. The destination is encoded in
the event number. No preliminary draw chooses a destination, and outgoing
probabilities are not normalized to sum to one.

The default clock is attained age. `duration_specific transit <origin> <sex>
<status> <destination>` selects duration for that particular stratum
(`load.cpp` lines 1062–1081). Duration is time since the person's **latest
group transition**, irrespective of which transition initiated the current
episode. Without history, duration is age: entry into the initial group is
implicitly placed at birth. Marriage changes do not reset this group clock.
Age- and duration-based schedules can coexist across destinations/strata;
there is no joint age-by-duration surface on this path.

Worked example: someone born in month 0 moves A to B in month 240. At month
300 their B-to-C age schedule is evaluated at 300 months and a duration
schedule at 60 months. Without the transition record, the duration schedule
would instead be evaluated at 300 months. A later move to any group starts a
new group-duration episode.

Candidate comparisons use strict `<`. Death is initialized first, followed
by birth, marriage, divorce, then destinations in increasing group-number
order. Thus a transition loses an equal-month tie with any earlier event;
between tied destinations the lower group number wins. This is same-person
candidate priority, not a global ordering of all people in a month. See the
[competing-events audit](socsim-competing-events-audit.md) for monthly dispatch.

## State changes, histories and rescheduling

`xevents.c::transit()` (lines 1267–1355) determines the destination from the
selected event and assembles the affected people. For each whose group differs
from that destination, it:

1. Decrements the old group population count and increments the new one.
2. Prepends a record containing month, origin and destination to that person's
   transition history, then changes the current group.
3. Removes the person from the marriage queue if present, removes their old
   event-queue entry if they are a following household member, and calls
   `date_and_event()` to install a fresh next event.

Consequently changed people immediately use destination-group mortality,
fertility, marriage, divorce and outgoing transition schedules when their next
event is sampled. Birth histories and union dates remain intact: group change
does not itself reset birth spacing, parity or union duration. Stored individual
multipliers are retained. The function does not dissolve a union or revise an
unchanged spouse's pending event. Household members already in the destination
receive no new history record and are not rescheduled by this loop.

`io.c::write_otx()` (lines 1054–1085) writes one row per individual change:
person ID, month, origin, destination and sequence. Latest is 0, then -1, -2,
and so on. `read_otx()` (lines 537–569) restores the linked history, relying on
the ordered rows and valid person references. Missing initial history is
permitted (`events.cpp` lines 388–391, 469–476); it therefore matters for
duration initialization. Current group alone cannot reconstruct historical
exposure. Newborn inheritance uses parents' current groups at birth or a fixed
group (`xevents.c::birth()`, lines 73–105), not their birth identities.

The aggregate `crnt_month_events[NONTRANSIT + destination]` counter increments
once per initiating event, even when several household members change group.
Individual `.otx` rows and aggregate event totals therefore count different
things under household migration.

## Optional household moves and source caveats

`hhmigration` defaults to zero (`events.cpp` line 223). When enabled,
`assemble_household()` (lines 1357–1459) starts with the initiator and adds:

- For a female initiator, all ongoing husbands/partners and her qualifying
  maternal children.
- For a male initiator, ongoing wives/partners and their qualifying children,
  plus qualifying children of wives/partners whose unions ended through death
  (including cohabitant death).

This uses recorded union and maternal-child links, not a residence or custody
table. Followers move to the initiator's destination regardless of their own
transition rates. A follower already there is skipped. Ongoing cohabitations
pass the same `reason_end == M_NULL` check as marriages.

Three important static caveats need controlled reference tests:

- **Child-age sign:** `add_minor_children()` (lines 1461–1498) computes
  `(birthdate - current_month) / 12 > 21`. For an already-born child this
  quantity is nonpositive, so the supposed upper-age exclusion does not work.
  Living adult offspring can pass. The other test is `last_marriage == NULL`,
  meaning never partnered in the retained history, rather than currently
  unmarried. A never-partnered 40-year-old can pass; a divorced 18-year-old
  cannot. These are direct predicate implications, not measured run results.
- **Marriage-pool removal:** every changed person is removed from the marriage
  queue, but `transit()` has no explicit re-entry for available males in
  one-queue mode. Other event paths explicitly insert such men. This suggests
  a risk that a moving available man disappears from the candidate pool;
  a controlled run is needed to establish persistence and matching effects.
- **Self-transition:** the candidate loop does not exclude destination equal
  to origin. If such a schedule wins, the initiator's update/rescheduling body
  is skipped because their group already equals the destination, although the
  initiating event counter still increments. With household mode, other members
  in different groups can still change. Avoid interpreting a self schedule as
  an ordinary harmless competing event; queue consequences need testing.

These quirks are findings about the pinned version, not requirements for Hammel.
Under D109, Hammel's core scope is individual transitions; linked household
moves are a stretch goal whose membership and follower rules remain open.
ENHANCED builds have before/after hooks and can alter behavior; those extensions
were not audited.

## Heterogeneity and transition targets

`utils.c::usetransmult()` (lines 212–252) returns the schedule lambda unchanged
unless `hettrans` is enabled, then multiplies by the person's persistent
`tmult`. `initialize_person()` initializes `tmult` to 1 (`io.c` line 203).
The ordinary path provides no fresh random distribution here; enabling the
switch alone does not create heterogeneous multipliers. `modify_rates()` also
applies the indexed `rate_factors` (`events.cpp` lines 1539–1542).

`transittarget <origin> <destination> <number>` is parsed in `load.cpp` lines
722–733. `fill_rate_gaps()` (lines 2292–2318) computes a target/expected
`transit_adjust`. However, a source-wide search finds that array only declared,
initialized, assigned and logged: the ordinary transition hazard path does
**not read it**. Thus the calculated target adjustment is not wired into that
path in this version. Do not assume the directive enforces a flow target.

Hammel decision D108 rejects transition-count target adjustment: supplied
hazards will be used directly, with realized counts determined by exposure,
competition and random variation.

The expectation helper (`utils.c` lines 83–132) itself projects current living
origin members across segment months at their initial marital status, without
competing exits, deaths or household effects. It indexes by age even for
duration schedules, sums negative lambdas from `fill_in_rates()`, and truncates
to an integer. Its count initialization starts at sex index 1 while the later
sum includes `MALE` (0), leaving a further uninitialized-count concern. These
are additional reasons to test this facility independently before relying on
it; no numerical target-recovery claim is made here.

## Questions for Hammel and suggested validation

The reference supplies useful alternatives but does not settle milestone 6.
The user selected annual person-year hazards for Hammel's transition inputs
(D106), continuing its existing rate convention. An annual probability p can
be converted with `-log(1-p)` under a constant-hazard, no-competing-exit
interpretation; observed cause-specific probabilities need additional care.
Individual transition execution is now implemented in Hammel; see the
[transition contract](group-transitions.md).

D107/D110 select directed competing schedules and one common age or duration
band type across the transition rate table. Required columns are
`group_origin`, `group_dest`, and `rate`, alongside the bands. Optional `sex`
and `union_type` columns stratify schedules; omitted dimensions apply
universally. The October 30 implementation request and follow-up decisions
select latest-entry duration (age without history), zero risk for missing pairs
and gaps, and unchanged ancestry. The implemented individual process updates
events, pools, calibration, fertility cells, pair hazards and historical group
observations. Household moves remain a stretch goal; fixed newborn assignment
and multiple independent group dimensions remain separate future work.
Scientific milestone acceptance is pending.

Useful future checks include constant one-way hazards against analytical
survival; two outgoing destinations against competing-risk probabilities;
reversible transitions; age versus duration with known initial histories;
equal-time event priority; zero/absent/self schedules; history reconstruction
and restart; group totals and person-years summing to pooled totals; fertility,
matching and union risks immediately after moves; and household followers with
different origins, ages and union histories if that mechanism is selected.
SOCSIM-specific checks should isolate the adult-child predicate, one-queue
male availability, self-transition queue retention, aggregate versus individual
counts, and target adjustment. Changed algorithms require statistical rather
than fixed-seed equivalence comparisons.
