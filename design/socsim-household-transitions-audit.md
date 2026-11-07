# SOCSIM household transition membership

Audited November 6, 2026 at the user's request. This expands the household
section of the [group transition audit](socsim-group-transitions-audit.md).
It establishes reference behavior, not a household design for Hammel.

## Reference and evidence

MPIDR/rsocsim 1.9.18, commit
`125450cff968bb3ff762fa99d4c926a8136062a2`, clean local checkout at
`../rsocsim`, confirmed November 6, 2026. The package declares GPL-3 and
`LICENSE.md` contains GPLv3. The ordinary non-ENHANCED code is examined;
enhancement hooks may change behavior. No reference implementation was copied
or adapted into Hammel's engine.

Sources, with line numbers at that commit:

- [Household assembly and child selection, xevents.c:1357–1498](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L1357).
- [Transition execution, xevents.c:1267–1355](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L1267).
- [Individual destination clocks, events.cpp:1998–2036](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp#L1998).
- [Household switch, load.cpp:1112–1118](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/load.cpp#L1112);
  `events.cpp:223` initializes it to zero.
- [Person and union structures, defs.h:575–632](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/defs.h#L575);
  union-ending constants are at lines 301–306.
- `xevents.c:108–146` maintains child/sibling links at birth;
  `xevents.c:177–255` closes ongoing unions when a partner dies.

Static tracing covers membership, triggering, state updates and accounting.
In addition, [13 executable selection fixtures](../validation/socsim-households/README.md)
passed on November 6. They compile the two selection functions and relevant
structures verbatim from the pinned source in a temporary directory. They
verify returned moving parties, not the full simulation, event queues, input
loading, or population-level migration outcomes. No full reference simulation
was run for this focused audit.

## What SOCSIM calls a household

With `hhmigration 1`, a person's transition event calls
`assemble_household()` to build a temporary linked list of accompanying
people. With the switch disabled, the list contains the initiator only.

There is no household identifier, household membership history, residence
table, or custody rule in this assembly path. The list is reconstructed from
union and maternal-child links at the event. It is not a partition of the
population into fixed, disjoint households. Different initiators can produce
different and overlapping parties.

| Initiator | Partners added | Children considered |
| --- | --- | --- |
| Female | Every ongoing husband/partner in her union-history chain | Her own maternal children, regardless of their fathers or her union history |
| Male | Every ongoing wife/partner in his union-history chain | Children of each ongoing wife/partner, plus children of partners whose unions ended in death or cohabitant death |

An ongoing union means `reason_end == M_NULL`; it covers cohabitation as well
as marriage. The code traverses all union episodes, rather than choosing only
one latest episode. A conversion's closed cohabitation episode does not qualify
as ongoing or death-ended, so that episode alone adds neither a partner nor
children. Its ongoing marriage episode can qualify normally.

The man's child selection does not test paternity: a current partner's child
by another man can accompany him. Conversely, his own child through a woman
who was never his recorded partner, or whose union with him ended in divorce
or cohabitation breakup, is not selected through that relationship. A female
initiator's selection of her own children does not depend on why her unions
ended. A dead partner is not added as a migrant; her living qualifying children
can be added through the death-ended union.

Assembly is not recursive. Adding a spouse does not call that spouse's own
household assembler. Adding a child does not bring that child's parents,
siblings or descendants automatically. Grandparents, other relatives and
unrelated coresidents have no general membership rule here. The spouse loop
does not independently test a partner's death date; ordinary death execution
is expected to close the union. The child loop does test death dates.

## Actual child predicate and confirmed defect

`add_minor_children()` starts at the mother's latest-born child and follows
older maternal siblings. It considers only children whose `deathdate == 0`.
For those children, its purported age limit is the integer-valued expression
`(birthdate - current_month) / 12 > 21`. For an already-born child the
numerator is nonpositive. Thus this condition never excludes an older child.

The other inclusion condition is `last_marriage == NULL`. This means no
recorded union history, not currently unmarried. A previous marriage or
cohabitation excludes the child even after its ending. The fixture with a
living never-partnered 40-year-old and a divorced 18-year-old selects the former
and excludes the latter. This verifies the function-level effect of the sign
error; it does not measure how often it affects an actual simulated population.

The comments describe minor/unmarried children, but the implemented behavior
is living children of any attained age with no recorded union. Correcting the
sign alone would also leave an integer-year threshold requiring review;
neither the defective predicate nor an inferred correction is a Hammel rule.

## Triggering rates and different origins

Household mode changes who follows an event; it does not replace individual
transition clocks with a separately defined household clock. Each person's
origin/destination, sex and marital-status schedules can produce a transition
candidate. When one wins, all selected followers are sent to that destination
without consulting their own transition schedules or sampling a destination
for them. There is no female-only or household-head-only initiation gate in
this path. Users would need to restrict input schedules to impose one.

There is no same-origin test in household assembly. A man moving A to B can
bring his partner from C to B and her child from D to B. A follower already in
B stays there without a transition-history record or resampling. The fixture
with different member groups confirms that origins do not affect assembly;
the destination updates and skipping rule are established by source tracing.

Consequently, supplied risks are initiation risks under household mode, not
automatically the total individual migration risks measured after accompanying
moves. For example, while two partners always follow each other to the same
destination, each can move through either partner's initiation. This is an
interpretation of the source behavior; no numerical rate-recovery experiment
was performed. SOCSIM's monthly timing and event competition still apply.

## Updates and accounting

For every selected person whose group differs from the destination,
`transit()` updates group population counts, prepends a dated origin/destination
history, changes the group, removes marriage-queue membership if present, and
samples a fresh next event. It also deletes a follower's prior event-queue
entry; the initiator's event has already been dispatched. An unchanged person
is not rescheduled. Followers execute ordinary state updates, not a new
recursive household-transition call.

Birth histories and union dates are preserved. The move does not itself
dissolve unions or erase parity or birth-spacing information. Fresh event
selection uses the new group and the person's own current characteristics.

The aggregate monthly transition counter increments once for the initiating
event. Individual `.otx` histories record every actual person's group change.
One initiating event can therefore generate several individual transitions,
potentially from several origins. There is no shared party identifier or
initiator/follower field in those ordinary transition-history rows.

The earlier audit's marriage-queue re-entry and self-transition concerns
remain source-only findings. This fixture harness does not resolve them.

## Consequences for Hammel

SOCSIM provides a concrete family-following rule, but it should not be read as
a general household model or as a calibrated household-rate system. A Hammel
proposal still needs user decisions about membership, coresidence/custody,
child ages and prior unions, initiators and rate denominators, different
origins, and recording linked moves. No such model is adopted by this audit.

The fixture outcomes are useful regression examples for discussing a future
contract: maternal versus paternal children, former-partner and deceased-partner
children, child-initiated moves, adult children, prior child partnerships,
cohabitation/conversion, and members already at the destination. A future
implementation must separately validate rate accounting and event invalidation.
