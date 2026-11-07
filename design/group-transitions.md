# Individual group transitions

Accepted milestone 6 contract, October 30, 2026 (D116). Implementation follows
D106–D115 and the user's October 30 request. The accepted
[report](../validation/milestone6.qmd) was undergoing its final render at
acceptance; preserve all active intermediates. Household moves, deferred by
D114 and excluded from milestone 6, were added in milestone 9 (D140); see
[household moves](#household-moves-d140) below.

Milestone 9 (D134, stage 3) made transitions per group dimension; this page
describes the current behavior.

`group_transition(rates, group = NULL, weights = NULL)` constructs an optional
process for one declared group dimension; `set_group_transition(model, rates,
group = NULL, weights = NULL)` attaches it. Each dimension has at most one
transition process: setting one again for the same dimension replaces it and
keeps the others. `group = NULL` means the model's only dimension; with several
declared it is an error naming them. A process stored before the declaration
exists takes the only dimension at simulation. Required columns are
`origin`, `destination` (categories of the process's dimension), `rate`, and
either `age_start`/`age_end` or `dur_start`/`dur_end`. The former
`group_origin`/`group_dest` names are errors. One band type applies to the entire table. Rates are
finite nonnegative hazards per person-year. Bands are half-open, start at
nonnegative finite values, have increasing ends (possibly Inf), and cannot
overlap within a directed pair/sex/union-type stratum at a given calendar time.
Self-transitions are rejected. Optional calendar bands are supported under
[the rate contract](time-varying-rates.md) (D118), without resetting group
duration. Input tables and provenance are retained.

Columns named by other declared dimensions select schedules by those
categories (education transitions by `region`); a person's risks change when
another dimension changes. The process's own dimension cannot be such a column
(its category is `origin`), and an undeclared column is an error at
simulation. As for sex and union type, uncovered strata have zero risk.

Optional `sex` uses m/f. Optional `union_type` follows existing conventions:
unpartnered/partnered, or unpartnered/cohab/marriage, without mixing partnered
with the two explicit partnered types. Omitted columns apply universally;
missing pairs, categories and gaps have zero risk. Existing union types select
actual cohabitation/marriage when explicit types are supplied. Without typed
unions an ongoing union is marriage. Categories from a dimension's transition
inputs join that dimension's inferred categories, so an initially empty
destination is allowed; other group-specific processes must cover it under
their existing rules.

Each outgoing schedule slot of a stratum (origin category and any other
dimensions the table uses) is its own process with cells and one clock
(D128); a living person is a member of each applicable slot's cells, competing
with other processes. Duration is time since the latest entry into the
person's current category of the process's dimension, falling back to birth
time without transition history in that dimension (user confirmed); moves in
other dimensions do not reset it.
Age is attained age. Changing union type moves both partners to their applicable
transition cells without resetting group duration. No count-target adjustment
is applied. With no transition process, no transition RNG draws occur.
Exact ties put transitions after death, birth, dissolution, conversion and
search; transitions in different dimensions follow the declaration order of
`set_groups()`, and tied destinations follow their first appearance in the
rate table. Each dimension's process has its own weight model; the weight
variable `destination` (formerly `group_dest`) is the destination of the risk,
and custom weight functions receive `origin` and `destination` in `bin` and a
`destination` member column.

At execution, only the initiator changes category, in one dimension. Record the change, remove old
pool/calibration contributions, refresh fertility membership and applicable
hazards, and reinstall contributions under the new group. Refresh an ongoing
union's pair-dependent dissolution cell and its group-dependent conversion
cell when applicable; retain union start, birth spacing, parity,
parents and individual fertility weights. Ancestry proportions stay unchanged
(user confirmed). Newborn inheritance continues to use current parental
identities and the existing ancestry rules.

Grouped populations include a `transitions` tibble even with no process or no
moves. Columns are `person_id`, `time`, `group` (the dimension), `origin`,
`destination`; `group`, `origin` and `destination` are character in every
dimension, since one table holds several dimensions' labels. A supplied
history may omit `group` when one dimension is declared. Histories preserve
repeated moves; each person's history in each dimension must chain to their
current category, fall within the person's life and known endpoint, and is
retained on restart. Same-time moves in different dimensions are allowed.
Duration resumes from the latest retained move in each dimension; cells are
rebuilt as for other processes.
A transition is effective at its recorded time. Historical observations split
exposure at moves rather than attribute an entire life to final identity.
Births/deaths and union events preceding a same-time transition use the origin;
population snapshots at that time use the destination. No missing prehistory
is imputed.

`get_detailed_unions()` reconstructs `<group>_w` and `<group>_h` for every
dimension at each union's start, before a same-time transition, rather than joining final identities.
Later moves therefore do not relabel the pairings used to assess formation.

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
not independently establish scientific acceptance. Milestone 6 acceptance
comes from the user's explicit approval recorded in D116, retaining the
report's limitations.

Example (illustrative rates, not demographic defaults):

```r
rates <- data.frame(
  origin = c("A", "A", "B"), destination = c("B", "C", "A"),
  age_start = 0, age_end = Inf, rate = c(0.2, 0.3, 0.1)
)
model <- model |> set_group_transition(rates, group = "race")
result <- simulate(model, years = 5, seed = 42)
result$transitions
```

The model must already contain grouped population, mortality, fertility and
union-formation specifications. Replacing age columns with `dur_start` and
`dur_end` selects duration for the whole transition table.

## Household moves (D140)

*User decisions, November 6, 2026; informed by the
[SOCSIM household audit](socsim-household-transitions-audit.md).*
**Status: implemented November 6, 2026.**

`set_group_transition(model, rates, group, weights, household = NULL)`, where
`household` is `NULL` (each person moves individually, the default) or a rule
from `household(age_independence = 18)`. With a rule, a move in that
process's dimension carries the mover's household along. The rule is its own
object so it can serve other features later; each process takes its own. Each transition process has
its own flag, so region moves can be household moves while other dimensions'
transitions stay individual. There is no household identifier or residence
record: the moving party is assembled from current union and parent links at
the move, as in SOCSIM, with the rules below in place of SOCSIM's.

**Dependents.** A person is a dependent while they are alive, younger than
`age_independence` (default 18, in years), not in a union (marriage or
cohabitation), and have a living parent to be attached to. A dependent is
attached to their mother while she is alive, otherwise to their father while
he is alive. A dependent with neither parent alive is independent; there is no
attachment to grandparents, stepparents or other relatives. A person under
the cutoff who leaves a union is a dependent again if a parent is alive.

**Who carries a clock (one per household).** In a household process:

- a woman in a union carries the household's clock, with her own
  characteristics (age, union type, other dimensions) selecting her schedule;
- a man in a union has no clock while the union lasts;
- dependents have no clock;
- everyone else, including single adults, adults past the cutoff living with
  a parent, and independent people under the cutoff, carries their own clock.

So each person's move rate is the input rate at the carrier's age, and
observed individual migration rates can be used as input without double
counting moves by partners and children. Children move at their mother's (or
father's) rate, not at a rate for their own age. This holds when households
live together: a partner or child already in the destination does not move, so
in split households (for example, matching that ignores region, or newborns
taking a random parent's region) followers move less often than the input
rate. Keeping households together, with a region prohibition in matching and
mothers passing on their region, gives every person the input rate (checked:
with a constant 0.05, women, men and children of every age moved at
0.049–0.052; with region ignored in matching and random-parent inheritance,
partnered men aged 18–40 moved at 0.039).

**The moving party.** When a carrier moves from origin to destination in the
dimension, the party is the carrier, their partner (if any), and the
dependents attached to any member of the party, recursively: a teenage
dependent who has a child brings that child along. Each member whose category in that
dimension differs from the destination moves there, including a partner or
child in a different origin; members already at the destination stay and are
not recorded. Each mover's time in the dimension restarts at the move; ages,
union durations, parity, birth spacing and ancestry are unchanged. All members
are updated at once before any other event: rate cells, fertility cells,
matching pools and calibration, and union strata.

**Re-attachment.** When a dependent's attachment changes so that they now
depend on a parent living elsewhere in a household dimension, they join that
parent at once, with their own dependents. This happens when a mother dies
and her child re-attaches to the father, and when a person under the age of
independence leaves a union (separation or widowhood) and is a dependent of a
parent again. The move is recorded with the parent as `initiator_id`.

**Reaching the cutoff.** A dependent has no clock until `age_independence`; then
they have their own clock at the ordinary rate for their age. Becoming a
dependent or independent at other times follows births, deaths and union
changes.

**Records.** The transitions history gains `move_id` (shared by all rows of
one move) and `initiator_id` (the carrier who moved). Individual processes
fill them too (one row per move, initiator = the person). Supplied histories
may omit them; they are then NA.

**Limitation.** Household processes require age bands: with duration bands
there is no shared time at which dependents reach their own clock.
`calculate_transition_rates()` counts every individual move (followers
included) over each person's exposure, so it measures individual rates; moves
counted once per household are the rows with `person_id == initiator_id`.

**Worked examples** (one region dimension; household process; ages in years):

1. A woman W (34) and her husband H (36), their child C (10), and H's son S
   (12) by his first wife, who died. W carries the clock; H, C and S have
   none. W moves N to S: four rows with one `move_id` and `initiator_id` W.
2. W's daughter D (19, unpartnered) lives with them: she is past the cutoff,
   carries her own clock, and stays when W moves.
3. C turns 18: from then C carries a clock at the input rate for 18, still in
   the region where the household left them.
4. H already lives in S (the couple formed across regions): W moves N to S;
   C and S move, H does not, so three rows.
5. W and H separate: H carries his own clock again; C stays attached to W and
   S to H.
6. W dies: C is attached to H (C's father). If H also dies, C and S are
   independent and carry their own clocks.
7. A single woman G (40) has a daughter T (17, unpartnered, so a dependent)
   with a baby B. G moves: T and B both move.
8. A mother in N and the separated father in S; their child (8) lives with the
   mother. The mother dies: the child moves to S at that moment, recorded with
   the father as initiator.

*Found in testing (November 6, 2026):* a long run first showed dependents
apart from the parent they were attached to: a teenage mother's baby left
behind (the party was not built recursively) and children of a dead mother
staying while the father lived elsewhere. The user chose the recursive party
and re-attachment moves; with both, households that start together stay
together through births, unions, separations, deaths and moves (tested).

## Observed transition rates

`calculate_transition_rates(population, time_breaks, age_breaks = NULL,
dur_breaks = NULL, by = NULL, conf.level = NULL, pool_dissolved = FALSE,
group = NULL)` counts one dimension's transitions (`group`, required when
several are declared) and returns a tibble
with directed `origin`/`destination`, calendar bounds, the selected age or
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

Attaching a transition process with all-zero rates is distributionally inert
but does not promise fixed-seed equality with omitting the process. Applicable
zero-rate schedules still consume hazard draws, and union-state changes can
reschedule partners' candidates while the process is enabled. Use independent
replications and statistical comparisons; a shared seed does not preserve
identical demographic paths across these configurations. This does not change
the no-extra-transition-draws guarantee when the process is absent.
