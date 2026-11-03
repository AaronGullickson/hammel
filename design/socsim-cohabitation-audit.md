# SOCSIM cohabitation: static source audit

Inspected October 28, 2026 at MPIDR/rsocsim commit
`125450cff968bb3ff762fa99d4c926a8136062a2` (subject “submission of 1.9.18”).
The local `../rsocsim` checkout was confirmed clean at that commit, with origin
`https://github.com/MPIDR/rsocsim.git`. Its `LICENSE.md` contains GPL version 3.
No code was copied or adapted. This is an independent behavioral description of
the ordinary, non-ENHANCED engine, not a reference simulation or a claim about
all SOCSIM releases, enhancements, or polygamous arrangements. The local path
is a convenience, not a Hammel dependency. Primary links below pin the commit.

## Main findings

SOCSIM treats cohabitation as a distinct current status within its existing
marriage/union machinery. New partners are found through the same search and
matching process. After a successful match, a separate probability determines
whether the episode begins as cohabitation or legal marriage. Cohabiting couples
can then legalize with each other, separate, or experience a partner death.

| Component | Inspected ordinary behavior |
| --- | --- |
| Current status | `COHABITING = 5`, alongside single/divorced/widowed/married |
| Entry | Ordinary marriage search followed by `cohab_prob`, not a separate cohabitation search clock |
| Type probability | Group/sex/age table for the initiating person, falling back to the other partner only if the initiating table is absent |
| Matching | Same partner pool, eligibility and preferences; type is drawn after matching |
| Conversion | `marriage` schedule for cohabiting status legalizes with the current partner; age-specific by default, optionally cohabitation-duration-specific |
| Breakup | `divorce` schedule for cohabiting status, always indexed by current episode duration |
| Fertility | Separate cohabiting-status schedules possible; ongoing cohabitant supplies father |
| After breakup or cohabitant death | Recover the individual's earlier legal single/divorced/widowed status, rather than necessarily becoming divorced/widowed |
| Records | Shared marriage-history table; conversion closes the cohabitation and appends a legal-marriage episode with the same partners |

## Source map

All paths and lines refer to the pinned checkout.

| Source | Evidence |
| --- | --- |
| `src/src/defs.h:267–305,416,575–631` | Status/end-reason constants, cohabitation probability arrays, person and shared union structures |
| `src/src/load.cpp:598–614,706–714,1621–1641` | Probability-table input, age endpoints and raw probability values |
| `src/src/utils.c:272–284` | Age-table lookup, including boundary convention |
| `src/src/xevents.c:342–430` | Ordinary search versus conversion with existing partner; event counting |
| `src/src/xevents.c:435–568` | Conditional type draw, precedence, new records, status updates, rescheduling |
| `src/src/events.cpp:1832–1995` | Status-specific scheduling, conversion clock and dissolution clock |
| `src/src/xevents.c:178–340` | Partner death, breakup, previous-status recovery and pool updates |
| `src/src/utils.c:1447–1511,1574–1590` | Previous legal status and record-count traversal |
| `src/src/xevents.c:35–66` | Father assignment from an ongoing union, including cohabitation |
| `src/src/events.cpp:1865–1914,1040–1064` | Cohabiting fertility schedules and optional birth-triggered legalization |
| `src/src/load.cpp:2130–2247` | Rate defaults by parity, status and group |
| `src/src/io.c:359–386,982–1003`; `R/read_results.R:9–17` | Record input/output and incomplete R end-reason documentation |
| `socsim_oversimplified.md:1267–1311` | Documented default inheritance; code remains the authority for exceptions |

## 1. Formation is a type draw after the ordinary match

The loader accepts `cohab_prob GROUP SEX` followed by age-table entries. Those
entries store raw probabilities; they are not waiting-time rates and do not use
the nonbirth `log(1-p)` conversion. Endpoints are expressed as years/months and
stored in months. `table_lookup()` advances only while `upper_age < age`, so
an age exactly at an upper endpoint still uses that entry. There is no linear
interpolation in this lookup. `add_rate_table()` warns for a negative value but
does not enforce an upper probability bound; do not interpret that permissive
parser as a valid-input recommendation. Explicit 0/1 values act as 0/1 in the
conditional type draw, unlike the near-zero conventions for event schedules.

For an ordinary new match, `new_marriage(p, spouse)` checks, in order:

1. If the initiating person is already cohabiting, legalize directly; do not
   draw another cohabitation type.
2. Otherwise, if their own group/sex `cohab_probs` table exists, look it up at
   their attained age and draw cohabitation with that probability.
3. Only if that table is absent, use the other partner's table at their age.
4. If neither table exists, form a legal marriage without this type draw.

Thus two supplied partner tables are not combined or averaged. A supplied
initiator probability of zero still prevents fallback to the other table.
For example, if a female initiator has .2 and her male partner's table has .8,
this code uses .2; it uses .8 only when her table is absent. In the one-queue
system women initiate searches, so this typically means the woman's table.
The table has no previous-marital-status dimension. No inheritance of missing
`cohab_probs` by preceding group was found in the reconciliation code; the
partner-table fallback above is the observed fallback.

Because the type draw occurs after spouse selection, the ordinary engine has
no independent cohabitation matching preferences or competing cohabitation
search rate. The input marriage-search schedule controls entry into either
type, conditional on finding an acceptable partner. Type probabilities can
then affect realized type-specific formation counts and future availability.
The shared matching restrictions remain active for new partners.

Primary sources: [formation and conversion](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L342-L568),
[table parser](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/load.cpp#L706-L714),
[lookup](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/utils.c#L272-L284).

## 2. Legalization uses the existing partner and a new episode

`marriage(p)` special-cases `COHABITING`: it retrieves the current partner,
removes that partner's pending event and calls `new_marriage()` directly.
There is no new search or reconsideration of other partners. The old episode
gets the current end month and `MARAFTERCOHAB = 5`. A new legal-marriage record
starts in the same month, with the same wife and husband and prior links to
the cohabitation. Both statuses become `MARRIED`, and both event candidates
are resampled under that status.

The conversion candidate uses `E_MARRIAGE` for cohabiting status, indexed by
person group and sex. By default it is age-specific. With `duration_specific`
for that block, its time origin is the current cohabitation's `date_start`.
This is distinct from divorced/widowed duration-search clocks, which use the
latest episode's ending. One-queue mode rejects a winning male marriage
candidate, so female-only conversion schedules are needed there as for searches;
it is not an automatic two-partner conversion risk in that mode. Two-queue
settings require their own rate and eligibility assessment.

A cohabitation starting in 2003 and legalized in 2005 therefore produces two
linked records. A divorce in 2006 uses marital duration since 2005, not total
partnership duration since 2003: the divorce scheduler reads the new latest
record's start. These dates illustrate source rules, not simulated results.

Primary sources: [same-partner conversion and records](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L365-L568),
[conversion time scale](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp#L1946-L1995).

## 3. Breakup and death restore legal history

Cohabiting partners can carry `E_DIVORCE` candidates from their own group/sex/
cohabiting schedule. Its clock is always current episode duration. As with legal
divorce, supplying both partners with risk creates competing dissolution risks;
it does not provide two observations of a single shared clock.

Breakup closes the shared episode with `BREAKCOHAB = 4`, removes the other
partner's pending event, recovers each partner's previous legal status, and
reschedules both. The one-queue divorce handler returns the man to availability.
Partner death instead closes cohabitation with `COCOHABDIES = 6`, restores the
survivor's earlier legal status, and resamples the survivor's event.

`get_previous_mstatus()` starts behind the current cohabitation and walks backward
through that individual's linked history. Earlier legal divorce gives `DIVORCED`,
earlier spouse death gives `WIDOWED`, and no earlier legal marriage gives `SINGLE`.
It skips prior cohabitation endings while continuing to look for legal history;
another ongoing union returns `MARRIED` (polygamous interpretation is not
validated here). This is derived from history rather than a separately stored
pre-cohabitation legal-status field.

Consequently:

- A never-legally-married person who loses a cohabitant returns to `SINGLE`,
  not `WIDOWED`.
- A previously divorced person whose cohabitation ends returns to `DIVORCED`,
  not a new legal-divorce classification caused by that breakup.
- Subsequent duration-specific searching nevertheless starts at the most
  recent cohabitation ending. For a divorce in 2000, cohabitation in 2003 and
  breakup in 2006, status is divorced but search duration starts in 2006.

That last rule separates the origin of legal status from the origin of the
search clock. Hammel currently classifies previous *unions*; adopting cohabitation
therefore requires an explicit decision about legal-marriage versus union history.

Primary sources: [death and dissolution](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L178-L340),
[previous legal status](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/utils.c#L1447-L1511).

## 4. Fertility, defaults and birth-triggered conversion

Fertility scheduling indexes `birth_rate_set[group][mstatus][parity]`, so
cohabiting women can have their own fertility schedule and spacing behavior
uses the latest actual birth. The ordinary birth handler assigns the current
cohabitant as father whenever the latest union has ongoing end reason `M_NULL`,
just as for a husband. The newborn's group rules then use these parent identities.
Mortality and group-transition inputs can also distinguish current cohabiting
status through their ordinary status index; no special cohabitant-death hazard
is added apart from the person's mortality candidate.

Missing rates use legacy defaults, which should not be adopted implicitly:

- In group 1, cohabiting fertility at parity zero inherits married parity zero.
  Missing higher cohabiting parities inherit the preceding cohabiting parity;
  they do not automatically copy the whole married parity curve.
- Group-1 missing cohabiting mortality and divorce entries inherit married
  entries of the same sex. The cohabiting marriage/conversion entry also
  inherits the married entry, but the loader deliberately skips defaulting
  the married marriage-search entry; that entry can remain NULL. Thus an
  omitted cohabiting conversion schedule need not provide active conversion.
- Higher groups inherit missing rate entries from the preceding group at the
  same status (and parity where applicable), rather than automatically from
  their own group's married schedule. Inherited duration flags propagate.

With `option_marriage_after_childbirth`, the birth-event dispatcher calls
`marriage(p)` before birth whenever status is not `MARRIED`. This includes
cohabiting mothers and therefore directly legalizes with the current partner,
even without a scheduled conversion event. This is an optional extra pathway,
not the baseline cohabitation conversion clock. Hammel's birth-triggered
formation does not imply adopting this rule.

Primary sources: [father assignment](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L35-L66),
[rate defaults](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/load.cpp#L2130-L2247),
[birth-triggered conversion](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp#L1040-L1064).

## 5. Output and static implementation caveats

The shared `struct marriage` has partners, start/end months, end reason and
prior links, but no explicit union-type field. Ended cohabitations are identified
by codes 4 (breakup), 5 (legalized) or 6 (partner death). Ongoing unions of either
type use `M_NULL = 16`; interpreting an ongoing episode requires current person
status in the ordinary monogamous case. R's `read_omar()` documentation describes
only legal divorce/death codes and is incomplete for these cohabitation reasons.
Conversion preserves both episodes and their link; `get_marity()` counts records,
so cohabitation plus legal marriage contributes two, not one distinct couple.

Two source-level discrepancies deserve controlled experiments before treating
SOCSIM behavior as a template:

1. **Male availability after cohabitant death.** The death handler explicitly
   returns a surviving man to the one-queue marriage pool in the legal-marriage
   branch. Its cohabitation branch restores status and reschedules, but contains
   no corresponding pool insertion. Formation previously removed the man from
   the pool. This appears to omit immediate re-entry after cohabitant death;
   segment/restart `new_events_for_all()` can reinsert men with single/divorced/
   widowed status. This is a static control-flow finding, not a measured
   population-level effect or a demonstrated behavior of other SOCSIM versions.
2. **Event counters mix/omit categories.** Successful new matches increment
   `E_MARRIAGE` even when the new record is cohabitation. The conversion branch
   appends a legal-marriage record but bypasses that increment. `E_DIVORCE`
   includes cohabitation breakup. Under distribution matching, age-gap tracking
   in `new_marriage()` also runs for cohabitation starts and conversions. These
   counters are not clean legal-marriage/divorce incidence measures.

One-queue male marriage candidates, competing partner risks, implicit rate
inheritance and polygamous history/status ambiguity also limit generalization.
These observations are not fixes to the reference implementation. The earlier
[dissolution audit](socsim-dissolution-audit.md) supplies the shared monthly
probability, event competition, former-partner-screening and restart caveats.

Primary sources: [death-handler branches](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L216-L263),
[pool initialization](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp#L1214-L1223),
[formation counting](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L365-L411),
[history output](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/io.c#L982-L1003).

## Implications for Hammel and subsequent design decisions

The structural lesson is to separate being partnered, union type, and previous
legal/union history. The original audit raised the following design areas;
subsequent accepted direction is recorded under D085 below:

- One formation process plus a conditional type probability (SOCSIM's approach),
  versus separately specified competing entry hazards. Whether type assignment
  precedes matching determines whether type-specific matching is possible.
- Same-partner conversion as its own transition, including its age/duration
  clock, single carrier, and whether total partnership duration or current
  type duration drives later dissolution.
- Which prior category follows cohabitation breakup/death: previous union type,
  previous legal-marriage history, or an additional cohabitation-specific axis.
- Separate fertility/dissolution schedules by union type, and explicit coverage
  rules rather than inherited empirical or numerical defaults.
- Explicit type on history records, and conversion represented as linked
  episodes or a type change within one continuous partnership.
- Availability, father assignment, spacing, queue invalidation, restarts and
  observation accounting across conversion and all endings. No conversion
  should accidentally release/reselect partners or double-count partnered time.
- Whether births can cause conversion. Current on-demand formation does not
  decide this for already cohabiting couples.

The audit itself adopted no cohabitation mechanism. Subsequent user decisions
(D085) select type before matching, separate optional matching roles, same-partner
conversion without previous-status strata, and default recovery of legal status
with an optional toggle counting cohabitation endings as dissolved/widowed.
D086 subsequently implements the authorized rate-mode subset in Hammel; see
[cohabitation](cohabitation.md). Reference queue/counter/default quirks are not
adopted. This audit itself remains a static reference trace.
The user subsequently accepted milestone 4, including optional cohabitation,
on October 29, 2026 (D093), based on the Hammel validation report. That does not
turn this static audit into an executable SOCSIM comparison.

## Suggested reference experiments, not yet run

1. With a fixed eligible pool, use cohabitation probabilities 0, 1 and an
   intermediate value; verify type assignment conditional on successful matches,
   initiator/partner precedence and absent-table behavior.
2. Begin with cohabitants and explicit conversion/breakup schedules; distinguish
   age from cohabitation duration, both-sex versus one-sex risks, and one- versus
   two-queue mode. Verify same partners, old reason 5 and new start date.
3. Give never-married, divorced and widowed people cohabitations, then cause
   breakup or partner death. Check recovered legal status and the subsequent
   duration-search origin separately.
4. Compare fertility by current type and father assignment; toggle the optional
   birth-triggered conversion pathway while scheduled conversion is absent.
5. After a woman's death in cohabitation, verify immediate male pool membership,
   subsequent matching and behavior after segment/restart reconstruction.
6. Compare episode-based legal marriage/cohabitation counts with monthly event
   counters; check conversion duration resets and read/write history recovery.

None of these reference experiments was executed for this audit. No fitted
comparison, scientific acceptance threshold or universal validity claim follows
from the static trace.
