# Individual group transitions

Implemented contract, October 30, 2026. Authorized for implementation under D106–D111
and the user's October 30 implementation request. Household moves remain a
stretch goal and are excluded here.

`group_transition(rates)` constructs an optional process;
`set_group_transition(model, rates)` attaches it. Required columns are
`group_origin`, `group_dest`, `rate`, and either `age_start`/`age_end` or
`dur_start`/`dur_end`. One band type applies to the entire table. Rates are
finite nonnegative hazards per person-year. Bands are half-open, start at
nonnegative finite values, have increasing ends (possibly Inf), and cannot
overlap within a directed pair/sex/union-type stratum. Self-transitions and
calendar-time schedules are rejected. Input tables and provenance are retained.

Optional `sex` uses m/f. Optional `union_type` follows existing conventions:
unpartnered/partnered, or unpartnered/cohab/marriage, without mixing partnered
with the two explicit partnered types. Omitted columns apply universally;
missing pairs, categories and gaps have zero risk. Existing union types select
actual cohabitation/marriage when explicit types are supplied. Without typed
unions an ongoing union is marriage. Group labels from transition inputs join
the possible identity registry, so an initially empty destination is allowed;
other group-specific processes must cover it under their existing rules.

Each living person samples one candidate for each applicable outgoing pair,
competing with other person events. Duration is time since the latest group
entry, falling back to birth time without transition history (user confirmed).
Age is attained age. Changing union type resamples applicable transition risks
for both partners without resetting group duration. No count-target adjustment
is applied. With no transition process, no transition RNG draws occur.
Same-person exact ties put transitions after death, birth, dissolution,
conversion and search; tied destinations follow their first appearance in the rate table.

At execution, only the initiator changes group. Record the change, remove old
pool/calibration contributions, refresh fertility membership and applicable
hazards, and reinstall contributions under the new group. Refresh an ongoing
union's pair-dependent dissolution clock and the woman's group-dependent
conversion clock when applicable; retain union start, birth spacing, parity,
parents and individual fertility weights. Ancestry proportions stay unchanged
(user confirmed). Newborn inheritance continues to use current parental
identities and the existing ancestry rules.

Grouped populations include a `transitions` tibble even with no process or no
moves. Columns are `person_id`, `time`, `group_origin`, `group_dest`; labels
retain the persons group type. Histories are ordered by time, preserve repeated
moves, must chain to each person's current group, fall within the person's
life and known endpoint, and are retained on restart. Duration resumes from
the latest retained move; pending clocks are redrawn as for other processes.
A transition is effective at its recorded time. Historical observations split
exposure at moves rather than attribute an entire life to final identity.
Births/deaths and union events preceding a same-time transition use the origin;
population snapshots at that time use the destination. No missing prehistory
is imputed.

Worked examples: A-to-B hazard 0.2 and A-to-C hazard 0.3 give time to either
move with total hazard 0.5, and destination shares 0.4/0.6 without other exits.
A duration schedule zero on [0,2) prevents another move for two years after
entry, while marriage changes do not restart that interval. A woman moving
A-to-B at time 3 contributes three years to A and two to B over [0,5), if she
is alive throughout; her fertility after 3 uses B and her ancestry is retained.

Validation: analytical competing hazards, age/duration support, sex/type
selection, absent schedules and malformed input, empty histories, repeated
moves/restarts, RNG isolation, cell/calibration accounting, changed partner
union hazards and historical group exposure/event totals. Software checks do
not establish milestone scientific acceptance.

Example (illustrative rates, not demographic defaults):

```r
rates <- data.frame(
  group_origin = c("A", "A", "B"), group_dest = c("B", "C", "A"),
  age_start = 0, age_end = Inf, rate = c(0.2, 0.3, 0.1)
)
model <- model |> set_group_transition(rates)
result <- simulate(model, years = 5, seed = 42)
result$transitions
```

The model must already contain grouped population, mortality, fertility and
union-formation specifications. Replacing age columns with `dur_start` and
`dur_end` selects duration for the whole transition table.

## Observed transition rates

`calculate_transition_rates(population, time_breaks, age_breaks = NULL,
dur_breaks = NULL, by_sex = FALSE, by_union_type = FALSE)` returns a tibble
with directed `group_origin`/`group_dest`, calendar bounds, the selected age or
duration bounds, `person_years`, `transitions`, and `rate_transition`.
Supply exactly one age or duration axis. Finite final axis breaks gain Inf.
Optional `sex` and `union_type` columns appear only when requested; default
counts pool those dimensions before dividing. With typed union histories use
unpartnered/cohab/marriage, otherwise unpartnered/partnered. Nonoverlapping
union histories are required for union-type disaggregation.

The group registry contains identities in persons and retained transition
histories, preserving label types. Return every distinct directed pair in
that registry, including zero-event cells. Unvisited identities absent from
both tables cannot be inferred from population output alone. Fewer than two
known groups gives a typed empty result. A known observation endpoint and a
grouped population are required; absent transition histories mean no recorded
moves, without imputing missing prehistory.

Exposure is time alive in the origin group, split at calendar, age/duration,
group-change and requested union-type boundaries. A-to-B and A-to-C share the
same origin exposure; do not sum that duplicated denominator across
destinations. Competing moves and death end an origin episode. Duration
resets on each entry, including return moves, and uses age before the first
recorded transition. Union changes do not reset group duration. Events use
the exiting origin's age/duration and union type after any earlier same-time
union event (transition priority is last). Exposure is half-open; an event on
an internal Lexis boundary enters the cell beginning there and events at the
final calendar boundary are excluded. An event can therefore have zero
exposure in its cell. The rate is transitions/person-years, zero with positive
exposure and no moves, NA at zero exposure, and is not a competing-risk
probability.

Worked case: a living person in A from birth at 0 moves to B at 2 and back to
A at 5, observed through 6. A has 3 person-years and one A-to-B move; B has
3 person-years and one B-to-A move. For duration [0,2) and [2,Inf), A has
3 and 0 years, respectively, with its outgoing move at duration 2 in the
second band; B has 2 and 1 years, with its outgoing move at duration 3 in the
second band. Union-type components and sex components sum to pooled exposure
and counts within each directed pair. Validation uses hand-calculated repeated
moves, union changes and deaths, boundary cases, empty histories and
analytical hazard recovery; no simulation state or randomness is changed.
