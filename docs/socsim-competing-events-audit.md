# SOCSIM competing events: initial source trace

Inspected October 12, 2026 at rsocsim commit `125450cff968bb3ff762fa99d4c926a8136062a2`. This is static source inspection, not executable validation. No implementation decision for Hammel follows automatically from these findings.

## Candidate selection for one person

`src/src/events.cpp`, `date_and_event()` (lines 1791–2052), samples a candidate death month first, then candidates for birth (for females), marriage if eligible, divorce if partnered, and group transitions. Each later candidate replaces the selected event only if its month is strictly earlier (`m1 < m`). Only one next event per person is retained in the event queue.

Consequently equal-month candidates are resolved by evaluation order: death, birth, marriage, divorce, then transitions in group iteration order. This is a priority among that person's candidate events, not a global order for all events in a month. Candidate waiting times have already been truncated to integer months by `datev()`, so differences within a month do not resolve ties.

## Processing a month

`process_month()` (around lines 947–1085) repeatedly chooses a random person from the current month's event queue, removes that entry, and dispatches the person's selected event. It continues until the queue is empty. Events affecting different people therefore have randomized sequential execution, not a universal death-before-birth ordering.

## Rescheduling and interactions

`src/src/xevents.c` contains handlers that update state and sample new candidate events:

- Birth reschedules both the mother and newborn (around lines 155–167).
- Death removes a surviving spouse's queued event and reschedules that spouse after widowhood (around lines 222–257).
- New marriage reschedules both partners (around lines 556–564).

Because candidate scheduling can return the current month and the monthly loop drains the queue, the same person can undergo another event in that month after rescheduling. One retained next event does not imply a limit of one realized event per person per month. Birth spacing separately constrains subsequent births.

For example, equal-month death and birth candidates for a mother select death. However, a birth followed by a newly sampled death in the same month is possible. A spouse's death or marriage can also invalidate and replace another person's pending candidate before it executes.

## Remaining work

## Monthly steps versus continuation boundaries

Within a segment, the loop in `events.cpp` (around lines 590–611) generates events once, then advances month by month while retaining future queue entries. At segment boundaries, it calls `initialize_segment_vars()` (around lines 791–798), which clears event queues (around lines 2285–2295); the next segment calls `new_events_for_all()` again. Thus scheduled dates are retained across ordinary month steps but regenerated across segments, even if a segment is only one month long.

`io.c`, `write_person()` (around lines 944–960), writes the selected next-event type but not its scheduled future date or the event queue. A file-based restart therefore does not restore that scheduling state; startup generates fresh candidates. This is not an exact checkpoint of the stochastic simulation.

Verify these behaviors with controlled executable cases. Audit all handlers, eligibility and matching rules, segment changes, numerical floors, and birth-spacing adjustments before asserting a complete scientific contract. In particular, do not equate integer-month tie priority with continuous-time competing-hazard semantics.
