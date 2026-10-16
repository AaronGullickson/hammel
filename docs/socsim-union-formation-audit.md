# SOCSIM union formation: static source trace

Inspected October 15, 2026 in local `../rsocsim`, commit
`125450cff968bb3ff762fa99d4c926a8136062a2` of MPIDR/rsocsim.
This describes that checkout's non-ENHANCED paths, not all SOCSIM releases.
Repository `LICENSE.md` contains GPL version 3. No source was copied or adapted
into Hammel; this is behavioral inspection, not executable validation.

## Initiation

`src/src/events.cpp:1946–1977`, `date_and_event()`, samples a union event from
rates indexed by group, sex, and current relationship status, using age or
configured status duration. It competes with other event dates; dates are
integer months. In two-queue mode, someone already on the partner queue does
not get another scheduled search. In one-queue mode, searches can be resampled
for someone on the queue; scheduled male union events are treated as errors.

`events.cpp:1214–1220` puts available men on the male queue in one-queue mode.
`xevents.c:342–444`, `marriage()`, searches the opposite-sex queue. A failed
search puts the initiator on their own queue if absent and reschedules their
next event. Thus two-queue search failure leaves someone available to be
selected by another initiator; in one-queue mode women can retry later.

`events.cpp:1050–1066` implements `marriage_after_childbirth` by calling the
union handler BEFORE executing the birth if the woman is not MARRIED. Birth
proceeds after a failed search. The condition includes cohabiting women;
`marriage()` converts their existing cohabitation rather than finding a new
partner. The option does not itself suppress separately configured union rates.

## Candidate sampling and matching

`xevents.c:693–762`, `create_working_mqueue()`, visits the whole opposite-sex
queue or samples distinct queue indices without replacement, THEN screens
eligibility. The default `marriage_sample_size` is 1000 (`events.cpp:232`);
zero requests the whole pool through the conditional at `xevents.c:725`.
Sampling therefore does not guarantee the configured number of eligible
candidates, and a sampled search can fail despite eligible people elsewhere.

After screening, only maximum-scoring candidates are retained. `random_spouse2()`
(`xevents.c:573–611`) selects uniformly among those ties. Negative scores are
accepted. There is no default softmax choice or low-score outside alternative.

`score3()` (`xevents.c:1005–1099`) supports:

- PREFERENCE: piecewise linear penalties around a preferred male-minus-female
  age gap, with a configurable ratio of slopes on either side. Gap is in months.
- DISTRIBUTION: target age-gap-bin frequency minus accumulated observed
  frequency (denominator `1 + marriage_tots`), favoring underrepresented bins.
  Gaps are integer years; counts are updated by the woman's group.

## Eligibility and caveats

`marriage_allowable()` (`xevents.c:807–868`) checks male minimum age (default
15), age-gap bounds (bypassed in DISTRIBUTION mode), endogamy, and kinship.
Availability is largely provided by queue membership, not this function.
Endogamy can be a strict within/across-group rule or stochastic candidate
acceptance (`xevents.c:1103–1203`); it is not a guarantee of realized proportions.

`check_incest()` (`xevents.c:907–1001`) compares each person's six recorded
parents/grandparents against the other's. This detects shared parents or
grandparents, including half siblings and first cousins, and some other close
relationships through shared ancestry. It is not a complete ancestor/descendant
traversal. If EITHER mother is unknown it immediately accepts, skipping both
ancestry screening and the subsequent former-partner check. Even with known
mothers, father/daughter and grandparent/grandchild are not directly tested by
comparing ego identity against the other's ancestors. Do not equate source
comments promising kin exclusions with comprehensive coverage.

The former-partner traversal also has a duplicated `p1->last_marriage` guard.
Treat this as an apparent defect needing targeted tests, not a verified outcome
for every history. ENHANCED builds can replace screening/scoring functions.

## Union creation, fertility, and death

`new_marriage()` (`xevents.c:447–568`) records both partners, prior-union links,
start month and ongoing status, updates statuses, and reschedules both partners.
It distinguishes marriage/cohabitation and can convert between them.

Birth hazards select schedules by group, relationship status, and parity
(`events.cpp:1884–1886`); spacing is also enforced. At birth, an ongoing partner
is father; optional random-father selection is a separate fallback
(`xevents.c:50–72`). Status changes consequently trigger new fertility candidates.

Death (`xevents.c:178–260`) removes queue membership, ends ongoing union
records, updates survivors' status, and reschedules their events. In one-queue
mode widowed men are put back in the partner queue. Cohabitation has distinct
end/status handling. SOCSIM also has non-death dissolution, outside Hammel's
milestone 2 scope.

## Consequences for Hammel

This trace informs open decisions without adopting legacy rules. Hammel's
conditional-logit choice deliberately differs from maximum-score selection.
Candidate sampling before versus after eligibility, incomplete-ancestry policy,
failed searches/retries, and birth-time classification need explicit choices.
SOCSIM's status-specific fertility is evidence of state-dependent rate selection,
not proof that Hammel's midpoint compensation is valid across status changes.

## Further kin extraction (additional static inspection October 15, 2026)

The same pinned checkout separates union prohibition from broader kin extraction.
The native `check_incest()` has fixed two-generation arrays; it has no generic
user-selected cousin-depth traversal. `defs.h:584–588` stores mother/father,
next-elder maternal/paternal sibling, and last-born-child pointers. Native
`utils.c`, including `child_census()` around line 657, contains hard-coded
traversal/tally routines extending to aunts/uncles and great-aunts/uncles;
these do not automatically extend the union screen.

The rsocsim R layer has a separate `R/retrieve_kin.R`, `retrieve_kin()`, operating
on output tables. It builds parent and parent-to-children lookups, composes
parent/children/spouse operations, and returns named lists of relative IDs.
Default categories include parents, grandparents, great-grandparents, siblings,
spouse, children and grandchildren. Optional categories include aunts/uncles,
great-aunts/uncles, first cousins, nieces/nephews, and in-laws, with optional
sex-specific outputs. Siblings are children of either parent minus ego, so
half siblings are included without a separate full/half classification.
First cousins are children of aunts/uncles. This is a fixed category interface,
not a generic arbitrary cousin-degree API. These findings are source inspection,
not executable validation of every category, missing-history case, or overlap.
No code was copied. Keep this R-layer function distinct from the SOCSIM engine.

## Previous nonmarital fathers (additional static inspection October 15, 2026)

In the same pinned non-ENHANCED checkout, `xevents.c:50–69` assigns the current
ongoing husband/cohabitant or invokes `find_random_father(child)` when enabled.
`find_random_father()` at lines 1509–1663 scans event-queue month buckets from a
random starting bucket, screens potential fathers with `marriage_allowable()`,
and selects from the first qualifying bucket. It does not consult the mother's
previous children or their fathers, and it deliberately does not require a
potential father to be currently unmarried. It does not use ordinary spouse
scoring. `marriage()` and `score3()` likewise contain no preference for a
previous child's father in the inspected ordinary paths. The mother/child/sibling
pointers store parentage but do not themselves implement preferential reuse.
These are static source findings, not executable validation or claims about
custom ENHANCED variants. No code was copied. Hammel intentionally differs by
using conditional-logit matching restricted to living unpartnered men.
