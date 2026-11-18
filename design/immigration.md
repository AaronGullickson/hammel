# External immigration design

Updated November 18, 2026. External immigration is implemented: individual
and household arrivals with imported histories (phases 1 to 3), cloning
streams, the automatic `birthplace` group, count and crude-rate intensities,
template cuts, weight schedules and `calculate_immigration_rates()`. This is
the current specification, following D177–D182 as amended and D184–D198;
the [validation report](../validation/external_migration.qmd) awaits the
user's assessment. Behavior depending on the [open decisions](#open-decisions)
requires further user decisions before implementation. See the
[roadmap](roadmap.md#post-release-immigration-from-outside-the-system) for
scope. No emigration mechanism is selected here.

## Purpose and intended workflow

Add people from outside the simulated population during a continuous-time
run. Users should be able to describe different origin populations and
migration flows without constructing every migrant record manually or
restarting the destination simulation for each arrival.

The intended workflow is:

1. Prepare a source `hammel_population`, potentially from a separate burn-in
   simulation. A template need not be stable or represent the whole origin
   population; users control the supplied source.
2. Use the template as supplied. Truncating source histories (D180) is
   deferred by D184; see [template truncation](#template-truncation-deferred).
3. Define one or more streams, each referencing a template and specifying an
   absolute arrival intensity. Templates may change through calendar schedules.
4. Optionally describe selective migration from the template with a weight
   model (D184).
5. Run the destination model, record actual arrivals and preserve sufficient
   migration history for rate strata, classifications and later analysis.

Population preparation, stream/model construction and simulation execution
remain separate responsibilities. The public interface is described in
[interface](#interface).

## Decision summary

| Topic | Current position | Record |
| --- | --- | --- |
| Interface | One `immigration()` process of named `immigration_stream()` components; `set_immigration()` replaces all streams | D184 |
| Arrival flow | Stochastic intensity as a count per year or a crude rate per destination person-year; no controlled/exact-count option | D177, D194 |
| Origins | Multiple streams with separate intensities | D178 |
| Source | Sample an existing population with replacement; preserve age, sex and groups; assign fresh IDs | D179 |
| Template schedules | Cuts of one template every `step` years on a shared calendar (D196), not a list of templates | D180, D184, D196 |
| Imported histories | Arrivals bring their birth, union and transition histories, with dates shifted by the age offset; referenced non-members go in `external_relatives` | D184 |
| Eligibility | No one arrives without their current partner; cloning streams import no relationships, so anyone alive can be cloned | D184, D188 |
| Phasing | Individual arrivals without histories, then imported histories, then household parties | D184 |
| Phase 1 contract | Unpartnered eligibility, per-stream clocks, arrival columns, observation from arrival, group compatibility, `migrant_assignment()`, ancestry fallback | D185 |
| Phase 2 contract | Imported unions, children and transitions; flattened template external relatives; `censor_time`; stored `pre_migration`; population validation; kin extraction excludes external relatives by default | D186 |
| Selective migration | A weight model, or a calendar schedule of them, selects among a stream's eligible template people; intensity is unchanged | D180, D184, D197 |
| Migration histories | Separate personal/family facts from reusable, named classifications | D181 |
| Rate schedules | Stratify by personal immigration status and stream; generation rules stay outside the required rate-table schema | D182 |
| Birthplace | Automatic `birthplace` group (`"native"` or stream), replacing `stream`; optional heritage ancestry; never the only dimension; no transitions | D190, D198 |
| Generation | User-defined groups with custom inheritance and migrant assignment | D192, D193 |
| Immigration rates | `calculate_immigration_rates()`: crude rates by period, `by` strata and stream; intervals adjusted for household dispersion | D194 |
| Phase 3 contract | Household parties under D140's rule; carriers sampled; intensity in persons converted by expected party size; ongoing unions arrive live and are left-truncated | D187 |
| Households | Arrivals are parties; household parties are phase 3 | D179–D180, D184, D187 |
| Cloning | Household streams are the default; `household = NULL` clones age, sex and groups without history | D188 |

## Arrival process and streams

User decisions, November 12 and 17, 2026 (D177, D184, D194). The intensity
tables of D194 were implemented November 17, 2026; see
[testing](testing.md#intensity-tables-november-17-2026).

**Intensity tables (D194).** A stream's `intensity` is a table, as other
processes take rate tables, with exactly one of two value columns:

- `count`: an absolute intensity in arrivals per year (D177). The stream's
  hazard is the count, whatever the destination's size.
- `rate`: a crude immigration rate in arrivals per person-year lived in the
  destination. The stream's hazard at time t is the rate times N(t), the
  number of living people in the destination, counting immigrants once they
  have arrived and never external relatives. Arrivals, births and deaths all
  change N(t), so arrivals feed back into later arrivals, as births do.

A table with both columns or neither is an error; a stream has one kind.
The argument keeps the name `intensity` to mark that it differs from other
rate schedules: it has no age, sex or group strata (who arrives is decided
by the template and weight model), and other columns are an error. A bare
number is no longer accepted; the error points to `tibble(count = ...)`.

**Calendar columns** follow D118, as for other rate tables. Without
`time_start` and `time_end`, the value applies throughout the run. With
them, intervals are half-open `[time_start, time_end)` in absolute calendar
years, ends may be `Inf`, values are constant within an interval, intervals
must not overlap and must cover the run, and an explicit zero pauses a
stream. Boundary changes take effect before same-time events.

**Arrivals as a point process.** With counts, arrivals form a
piecewise-constant Poisson process: over a period the expected count is the
integral of the count, so 100 per year for 2000–2010 then 50 for 2010–2020
gives an expected 1,500 over 2000–2020, with Poisson spread. With rates,
the expected count is the integral of rate × N(t); for a destination near
10,000 people, a rate of 0.004 gives about 40 arrivals a year. Household
streams convert either kind from persons to parties through E[size] (D187).

**Engine.** A count stream keeps its exponential clock, redrawn at
boundaries. A rate stream has one Exp(1) cumulative-hazard threshold: its
hazard is constant between changes in N(t), so the engine accumulates
rate × N × elapsed time at each change and fires the arrival when the total
reaches the threshold, as rate cells do (D128). Exact ties keep D185's
order, with arrivals after group transitions.

Allow multiple streams (D178), potentially representing countries or world
regions of origin, each with its own intensity and source template.
Streams may mix counts and rates. Controlled/exact-count immigration is
excluded.

**Crude immigration rates.** `calculate_immigration_rates(population,
time_breaks, by = NULL, streams = FALSE, conf_level = NULL)` reports, for
each period, `arrivals`, `person_years` lived by the destination population
(immigrants from arrival, everyone else from birth or the start, never
external relatives) and `rate_immigration` = arrivals / person-years, with
Poisson confidence intervals like the other calculators. As elsewhere, `by`
stratifies both: `by = "region"` gives arrivals into a region over
person-years lived there, with arrivals classified by their category at
arrival. `by` accepts `"sex"` and declared dimensions except `birthplace`,
which is an error pointing to `streams`. `streams = TRUE` splits each
stratum's arrivals by stream over the stratum's person-years, so stream
rates sum to the stratum's rate, matching the crude-rate definition of
`rate` intensities. Household arrivals are compound Poisson, so with
`conf_level` each row also reports `dispersion`, the sum of squared
household sizes over its arrivals (households identified by shared arrival
time and stream), and confidence intervals are exact Poisson intervals for
arrivals and person-years both divided by it; individual arrivals have
dispersion 1 and unchanged intervals. Without `conf_level`, `dispersion` is
not reported.

## Interface

User decision, November 16, 2026 (D184). Immigration is one process per
model, made of named streams, following the `set_groups()` pattern:

```r
model |>
  set_immigration(
    mexico = immigration_stream(tibble(count = 100), mexico_2000),
    philippines = immigration_stream(
      tibble(rate = 0.002),
      philippines_2000,
      weights = weight_model(...)
    )
  )
```

- `immigration_stream(intensity, template, weights, assignment, household)`
  describes one stream. Intensity comes first, as rates do in rate
  processes; it is a table with a `count` or `rate` column (D194). It is a
  component, like `group()` or `matching_model()`, not a process.
- `immigration(...)` is the process, holding named streams; argument names
  are stream names. It attaches with `set_process()` like any process.
- `set_immigration(model, ...)` constructs and attaches the process in one
  call, parallel to `set_mortality()`.
- Each call defines the complete set of streams, replacing any earlier set,
  as `set_groups()` replaces the whole declaration. `set_immigration(model)`
  with no streams removes immigration.
- Every stream must be named. Names must be unique and non-empty, and must
  not collide with the label for people who have not immigrated, because
  names identify streams in arrival records, D182 rate strata and D181
  classifications.

Each stream has its own clock (D185), and people who have not immigrated are
`"native"` in the `birthplace` group (D190).

## Population templates and selection

Draw migrant templates with replacement from an existing `hammel_population`
(D179). Preserve age, sex and group identifiers, and assign fresh destination
person IDs on every draw. Preserve age by translating birth time from the
source population's endpoint to the arrival date. The template is a reusable
snapshot and is not depleted by migration.

Selective migration uses a weight model, as other processes do (D184,
replacing D180's predictive-selection proposal). A stream's eligible
template people form one selection set: the stream intensity determines when
an arrival occurs, and the stream's weight model chooses who arrives, with
probability proportional to weight, as D129 selects a member within a rate
cell. Normalized relative weights preserve the absolute stream intensity:
multiplying every eligible person's weight by the same positive constant
changes neither the flow nor the selection distribution. Without a weight
model, eligible people are drawn uniformly. Weights only reweight people
present in the template; they cannot create characteristics the template
lacks. Weights are not probabilities of migrating or individual hazards.

Weight-model inputs are evaluated on template records at the template
endpoint, so age is the age at arrival. Which `weight_model()` variables are
available for template people, whether weight models may follow calendar
schedules here as elsewhere, and the result when every eligible weight is
zero or no one is eligible still need decisions.

Allow each stream's template to change over calendar time using time bounds
and a list column of population objects (D180). Column naming, validation,
coverage and boundary semantics remain open. Selecting the active template
when an arrival occurs, with changes affecting future arrivals only, is the
proposed event-time rule.

Optional household immigration is under discussion. Reuse of household
transition membership rules, household selection probabilities, within-party
ID/link reconstruction and relatives outside the party need specification.
The user explicitly deferred whether intensity means households per year or
migrants per year. D177's individual-arrival example does not settle the
person-count distribution for household arrivals.

## Integration with the destination simulation

Arrivals must compete chronologically with the existing continuous-time
processes. An immigration event adds new destination members rather than
selecting an existing resident to experience an event. Whether streams have
separate clocks or share a clock with stream selection is not decided.

An arriving person must enter all applicable destination process cells and
matching pools with consistent relationships and demographic state. A
template stores records and run metadata, not pending clocks, so nothing is
transplanted: entering the destination changes cell membership and redraws
affected process clocks, as a birth does. Destination rates and models
determine subsequent events. History-dependent state is initialized from
imported records under D184; see
[arriving parties and imported histories](#arriving-parties-and-imported-histories).
Importing supporting histories is distinct from treating every retained source
relative as an arriving member of the destination population.

Migration records and outputs should remain ordinary, inspectable R objects,
with public tables following the package's tibble convention. Their schema,
ID mapping/provenance and restart representation remain open. Each run must
own its state and must not mutate a reusable source template as it samples.

## Arriving parties and imported histories

User decisions, November 16, 2026 (D184). These describe phases 2 and 3 of
the [implementation phases](#implementation-phases); phase 1 imports no
history.

**Parties.** Every arrival is a party. While arrivals are individuals, a
party is one person. Household parties (phase 3) are expected to follow
D140's membership rule: a carrier, their current partner and their
dependents, recursively. One rule covers both: relationships inside the
party arrive as links between destination members, and people outside the
party whom a member's history references are recorded in
`external_relatives`. Arriving household children are members, not
external relatives.

**Eligibility.** No one arrives without their current partner, whether
married or cohabiting: household parties bring the partner along. No
destination record therefore has a partner whose future is not simulated.
Cloning streams (D188) import no relationships, so anyone alive can be
cloned.

**Imported histories.** Each party member brings their recorded births,
unions (including ended unions) and group transitions, including events
before arrival. An imported event is pre-migration when it occurred at or
before the person's arrival time: an event at the template's endpoint
arrives dated at the arrival time itself. A union is classified by its
start time. Unions
and transitions record this in a stored `pre_migration` column (D186).
Pre-migration events are not destination events and contribute no
destination exposure; calculation and observation functions must exclude
them.

**Date translation.** Let the offset be the arrival time minus the template
endpoint. The offset is added to every imported date of party members and
external relatives: births, deaths, union starts and ends, and transitions.
This extends D179's age preservation: age and every elapsed duration at the
template endpoint are unchanged at arrival. A migrant enters exactly as at
the template endpoint; nothing happens to them between that endpoint and
arrival, however much destination time has passed.

**Derived state.** Parity, last-birth time, previous-union status, time
since the last union ended, legal previous-union status, current-union
duration and group duration are derived from the imported records, as for a
starting population. Group duration continues from the latest source
transition in each dimension, or from birth without one; arrival does not
restart it.

**External relatives.** `external_relatives` is a tibble with the same
columns as `persons`, plus a stored `censor_time` (D186). It
holds people referenced by a party member's imported history who are not
destination members: previous partners, alive or dead, and children not in
the party, whether deceased or left behind. Under the eligibility rule no
current partner appears there.

- IDs are fresh destination IDs from the same ID space as `persons`, so
  union and parent references resolve without knowing which table holds
  the person.
- Group values are the person's categories at the template endpoint, or at
  death if earlier, without transition history, parallel to a `persons`
  row.
- `censor_time` is the death time if they died by the template endpoint,
  and the party's arrival time otherwise. Nothing after it is known.
- External relatives are never simulated. They join no rate cells or
  matching pools and contribute no exposure, events or population counts.
- Their own relationships are not imported. A child's other parent who was
  not a previous partner of the migrant has an unknown (`NA`) parent ID.

Under D140's membership rule every household party member is linked to
another member by a union or parenthood, so kin screening within a party
needs no external parents. Multigenerational parties, such as a carrier
with their own parent or adult sibling, would need the linking parents as
external relatives; they are not specified.

**Limitations.** The template is a distribution, not a population: parties
drawn separately are unrelated in the destination, including repeated draws
of the same source person. A previous partner's group is their status at
the template endpoint or death, so partner characteristics of pre-migration
unions, such as those from `augment_unions()`, may be wrong; pre-migration
partnering should be analyzed on the template population.

**Templates without histories.** A template without union or parent records,
such as one from `generate_uniform_population()`, needs no special handling:
migrants arrive unpartnered, at parity zero and without external relatives.

## Implementation phases

Proposed order (D184):

1. **Individual arrivals without histories.** Streams and intensity,
   sampling, fresh IDs, age translation, group validation against the
   destination, entry into cells and matching pools, arrival records, and
   arrival-aware accounting and observation. Every migrant arrives
   unrelated, never partnered and at parity zero. This is scaffolding and
   not a releasable immigration feature. Source people with a current
   partner are not eligible (D185).
2. **Imported histories.** Histories, date translation, derived state,
   carried-over group duration and `external_relatives`.
3. **Household parties.** Party sampling under D140's membership rule;
   arriving children become links instead of external relatives (D187).

Weight-model selection is part of phase 1 (D185). Template schedules,
D182 rate strata and D181 classifications follow once arrivals and
observation work.

## Phase 1 specification

User decisions, November 16, 2026 (D185). Implemented November 16, 2026;
see [testing](testing.md#external-immigration-phase-1-november-16-2026).

**Scope.** Phase 1 adds `immigration()`, `immigration_stream()`,
`set_immigration()` and `migrant_assignment()`; constant or calendar-varying
intensity; optional weight models; group assignment for dimensions the
template lacks; arrival records; and arrival-aware observation. It excludes
imported histories and external relatives (phase 2), household parties
(phase 3), calendar-scheduled templates, weight models and assignments,
D182 rate strata and D181 classifications.

**Eligibility.** A template person was eligible when alive at the
template's `end_time` and without a current partner. Under D188, phase 1's
individual arrivals became cloning streams, where anyone alive is
eligible.

**Arrival times.** Each stream has its own exponential clock at its current
intensity (D128's one clock per process). At an intensity boundary the
stream's clock is redrawn at the new intensity; by memorylessness this
preserves the piecewise-constant Poisson process. The intensity column is
`intensity`. Same-time events are ordered death, birth, dissolution,
conversion, search, group transition, then arrival. A model without streams
makes no immigration random draws, so existing seeded runs are unchanged.

**Selection.** Conditional on an arrival, one eligible template person is
selected with probability proportional to the stream's weight model, or
uniformly without one. In phase 1 weight models may use age at the template
endpoint (equal to age at arrival), sex and the groups the template
supplies; assigned groups are not available, because assignment follows
selection. Template and weights are fixed for the run, so selection
probabilities are computed once per stream. A stream with no eligible
person, or with zero total eligible weight, is an error before the run
starts.

**The arriving person.** The destination record has a fresh ID, the
selected person's sex, and birth time equal to the arrival time minus the
selected person's age at `end_time`. Parent IDs are unknown (`NA`). They
arrive with no union history, parity zero, no last birth and no transition
history; group duration therefore falls back to birth time under the
existing rule. They then enter every applicable rate cell, fertility cell,
matching pool and calibration count, and the affected process clocks are
redrawn.

**Groups.**

- A dimension declared by both template and destination must declare the
  same categories; otherwise it is an error. The migrant takes the
  template value.
- Template dimensions the destination does not declare are ignored.
- A destination dimension the template lacks needs a
  `migrant_assignment()` for that dimension in the stream's `assignment`
  argument, a list named by dimension; without one it is an error. An
  assignment for a dimension the template declares is an error: assignment
  does not override template values.

`migrant_assignment()` methods:

| Method | Assigns |
| --- | --- |
| `"fixed"` | One category to every migrant on the stream |
| `"probability"` | A category drawn from fixed probabilities over the dimension's categories |
| `"custom"` | Under D193, the category a user function returns for each arriving member, called at arrival |

Argument names and validation for each method are specified with the R
interface draft.

**Ancestry.** For each destination dimension that tracks ancestry, a
migrant takes the selected person's ancestry shares when the template
tracks ancestry for that dimension. Otherwise, because the template lacks
the dimension or does not track its ancestry, the migrant has 100% ancestry
in their category for that dimension, whether supplied by the template or
assigned.

**Arrival records.** `persons` gains `arrival_time` and `stream`, `NA` for
people who did not immigrate. The columns appear only when the model has
immigration streams or the input population already has them, so restarts
retain them and simulations without immigration are unchanged.

**Observation.** A person is observed from the later of birth and arrival.
Exposure, events and population counts before arrival are excluded, so a
migrant's own birth is not a destination birth, and their age-specific
exposure begins at their age at arrival. Every calculation and observation
function applies this rule; a population without arrival columns is
treated as having no arrivals.

## Phase 2 specification

User decisions, November 16, 2026 (D186). Implemented November 16, 2026;
see [testing](testing.md#external-immigration-phase-2-november-16-2026). Phase 2 added imported histories to individual arrivals under D184;
under D188 they apply to every member of a household party, and cloning
streams import none. Selection, groups and ancestry of the arriving person
are unchanged from [phase 1](#phase-1-specification).

**What an arrival imports.** For the selected template person P:

- **Unions.** Every template union in which P is a partner. By eligibility
  all have ended. Each becomes a destination union with a fresh
  `union_id`, dates shifted by the offset, its type and end reason, and
  `previous_cohab_id` remapped to the new ID of the linked cohabitation.
  P's partner in each union is imported as an external relative.
- **Children.** Every template person whose mother or father is P. Each
  becomes an external relative whose parent ID for P is P's new ID. The
  other parent's ID is the new ID of that parent when they are a previous
  partner of P imported above, and unknown (`NA`) otherwise.
- **Transitions.** P's template transitions in dimensions the template
  shares with the destination, with times shifted, fresh `move_id`s and
  `initiator_id` set to P's new ID. Transitions in template-only
  dimensions are dropped.
- **P's own parents** are not imported; P's parent IDs are unknown.

Template people are looked up in both the template's `persons` and its
`external_relatives`; where they came from in the template does not matter
in the destination. The template's own `arrival_time`, `stream` and
`pre_migration` values are ignored: arrival records describe arrivals into
the destination only.

**External relatives.** Imported partners and children form the
destination's `external_relatives`, a top-level element of
`hammel_population`, with fresh person IDs from the same ID space as
`persons`, shifted dates, sex, the parent IDs above, and group values as of
the template endpoint or their death. In groups the template lacks, their
values are `NA`: assignment applies only to arriving migrants, and external
relatives never lived in the destination. `censor_time` is:

- the death time, when they died by the template endpoint;
- their template `censor_time`, shifted, when they were already an external
  relative in the template;
- otherwise P's arrival time.

**`pre_migration`.** When a population has arrival records, `unions` and
`transitions` gain a logical `pre_migration` column. It is TRUE for a union
that started before either partner's arrival and for a transition before
the person's arrival; FALSE otherwise, including every destination event.
A birth is pre-migration when the child was born before the mother's
arrival; no birth flag is stored.

**Derived state at arrival.** Parity, last-birth time, previous-union
status and its end time, legal previous-union status and group duration in
shared groups are derived from the imported records, including children
who are external relatives. In assigned groups, which the template lacks,
group duration runs from the arrival time: arrival is entry into the
assigned category. Each such entry is recorded as a `transitions` row at
the arrival time with `origin = NA`, the assigned category as
`destination` and `pre_migration = FALSE`. Observation and restarts read
these rows as entries, never as moves, so transition counts and rates
exclude them. These are computed before the run for
every eligible template person and set when the migrant arrives. On
restart, the same quantities are derived from `persons`, `unions`,
`transitions` and `external_relatives` together.

**Validation of populations.** `check_population()`, and so
`validate_population()`, templates and restarts, requires:

1. `external_relatives` has the required `persons` columns plus
   `censor_time`.
2. Person IDs are unique across `persons` and `external_relatives`.
3. Union partner, mother and father IDs in either table resolve in one of
   the two tables.
4. Each external relative is referenced by a union or parent link of a
   `persons` row.
5. `censor_time` equals `death_time` for those who died; otherwise it lies
   between `birth_time` and `end_time`.
6. No external relative is in a union without an end time.
7. Group values use declared categories; external relatives may have `NA`
   group values, and have no ancestry or transition rows.
8. `pre_migration` is present on `unions` and `transitions` exactly when
   `persons` has arrival records, and has no missing values.
9. A transition may have `origin = NA` only as an entry: the person's
   first row in that group. Arrival entries are at the person's
   `arrival_time`; entries imported from a template (its own arrival
   entries) are pre-migration.

**Observation.** External relatives are never counted: no exposure, births,
deaths, population counts or rate cells. Observation still reads their
records to build members' histories, so a migrant mother's parity includes
children who stayed behind and her union status reflects pre-migration
unions. Exposure and events before a person's arrival are excluded, as in
phase 1; union-level calculations exclude pre-migration unions. Functions
that look up partners or parents, such as `augment_unions()`, resolve IDs
in both tables. `extract_kin()` and `count_kin()` exclude external
relatives by default; `include_external = TRUE` includes them.

## Phase 3 specification

User decisions, November 16, 2026 (D187). Implemented November 16, 2026;
see [testing](testing.md#external-immigration-phase-3-november-16-2026).
Phase 3 adds household parties under D140's membership rule.

**Interface.** `immigration_stream()` gains `household`, which defaults to
`household()` (D188), the same object `set_group_transition()` takes.
`NULL` makes a [cloning stream](#cloning-streams). Household and cloning
streams may be mixed in one model.

**Households in the template.** Households are read from the template's
living `persons` at its endpoint with D140's rule. A dependent is alive,
unpartnered and younger than `age_independence`, attached to their mother
while she is alive, otherwise to their father while he is alive. A
household is a carrier, the carrier's current partner and the dependents
of either, recursively. The carrier is the woman in a current union, or an
unpartnered person who is not a dependent, including an independent
orphan. Every living template person belongs to exactly one household.
Future, more complex household rules, such as multigenerational households,
must keep membership exclusive. Template external relatives are never
household members. A dependent whose mother is a template external
relative attaches to their father if he is a living member of the
template; otherwise they are independent.

**Selection.** A stream with a household rule samples carriers, with
replacement. The stream's weight model applies to carriers'
characteristics, so a household's selection probability is its carrier's
normalized weight. Under D189, carriers also carry their household's
`partnered` (0/1), `household_size`, `n_dependents` and `partner_age` (`NA`
when unpartnered) for weights; `custom` assignment sees every member
(D193). A formula term with a
missing value contributes nothing to the log weight, as in the engine's
weight models, so `partner_age` acts only on couples. Cloning streams can't
use these variables.

**Intensity.** `intensity` stays in persons per year. Let E[size] be the
expected party size under the stream's selection probabilities: the sum of
each household's size times its probability. The stream's clock runs at
`intensity / E[size]` parties per year, so expected arrivals are
`intensity` persons per year. A calendar-varying intensity is converted
band by band with the same constant. Persons arrive in clumps, so the
count over T years is compound Poisson, with mean `intensity * T` and a
variance larger than the mean: `(intensity / E[size]) * T * E[size^2]`.
For a cloning stream, E[size] is 1.

**Arrival.** Every member of a party arrives at the same time with the
same stream and the same date offset. A `"fixed"` or `"probability"`
assignment gets one draw per party, so members share it. A `"custom"`
assignment (D193) is called at arrival with one row per member, holding
their own `person_id`, `age`, `sex`, shared groups and `arrival_time`, plus
`role` (`"carrier"`, `"partner"` or `"dependent"`), a `household` number
and the household's D189 variables. Rows come in household order, carrier
first. It returns one category per member; random choices are made inside
the function with R's random numbers, which follow the run's seed, so a
household stays together only when the function gives its members the same
category. Shared dimensions and ancestry are copied for each member as in
phases 1 and 2; a template column is the fastest way to give members
person-level categories.

**Imported histories.** Each member imports their histories as in phase
2, with these changes:

- A union between two members becomes a destination union between them.
  This includes the current union, which arrives live.
- A child who is a member has their member parents' new IDs. A child's
  parent who is not a member is an external relative when they are a
  previous partner of a member, and unknown (`NA`) otherwise.
- Anyone a member's history references who is not a member becomes an
  external relative. Each is created once per party, even when several
  members reference them, such as a dead mother of a dependent who was
  also the partner's previous wife.
- A carrier's and partner's own parents are not imported, as in phase 2.
- `censor_time` of a living external relative is the party's arrival time.

**Ongoing unions.** A current union arrives live, with duration from its
original start, shifted by the offset. It is pre-migration
(`pre_migration = TRUE`), so its start is not a destination union
formation. From the arrival time it is simulated like any destination
union: it enters dissolution cells, and conversion cells when cohabiting,
in the duration band its original start implies; the woman's fertility
uses its union type; and neither partner searches. Previous-union state
comes from each partner's ended unions, as in phase 2.

**Observation of unions.** Union calculations left-truncate unions at the
arrival time instead of dropping pre-migration unions. A union contributes
union-years, dissolutions and conversions only after both partners have
arrived. Formation counts exclude unions that started at or before arrival.
A
pre-migration union that ended before arrival therefore still contributes
nothing, as in phase 2.

**Party records.** Members of a party share `arrival_time` and `stream`;
no party identifier is stored.

**Validation of populations.** The nine phase 2 rules apply. Rule 6 still
holds: a live union always has two members.

**Worked example.** A template ends at year 100 with four living
households: three unpartnered adults and one couple. The woman, born 70,
married a man born 66 at 90. Their children were born at 92 (died 93), 95
and 99. The man has a son born 85 from a previous marriage that ended with
his first wife's death at 89. With `age_independence = 18`, the son (15)
attaches to his father, since his mother is dead, and the children born at
95 and 99 attach to their mother. Their household has five members. With
equal weights, E[size] is (1 + 1 + 1 + 5) / 4 = 2, so an intensity of 100
persons per year runs the clock at 50 parties per year.

When this household arrives at 2035.4 (offset 1935.4), all five arrive at
2035.4 with fresh IDs. The couple's union starts at 2025.4, is live and
pre-migration, and enters dissolution cells at duration 10. The first wife
and the daughter who died are external relatives with censor times at
their deaths, 2024.4 and 2028.4. The first marriage is a pre-migration
union between the man and his first wife, ended at 2024.4. The son's
parents are the man and the first wife. The woman has parity three. The
couple's union contributes union-years from 2035.4 at duration 10; it adds
no destination union formation.

## Cloning streams

User decisions, November 16, 2026 (D188). Household parties are the default
and the recommended way to import migrants. A stream with `household = NULL`
clones template people instead, for uses such as labor migration from a
restricted template of young adults.

**Eligibility.** Every person alive at the template's `end_time` can be
drawn, whether partnered or not.

**What a clone keeps.** Age, sex, shared group categories and ancestry, as in
phase 1. Nothing else: a clone's history begins at arrival. They arrive
never partnered (`previous_union = "none"`), at parity zero with no last
birth, with unknown parents, and with no imported unions, transitions or
external relatives. Group duration starts at arrival in every dimension,
recorded as an entry row (`origin = NA`) per dimension.

**Selection and intensity.** Weights see every eligible person, and
`custom` assignment the arriving person, without household variables,
`role` or `household`. Each arrival is one person, so the clock runs at
`intensity`.

**Warning.** Constructing a cloning stream from a template with any unions
or transitions warns (kind `no_effect`) that these histories will not come
with migrants, and that a `household()` rule imports them.

## Birthplace group

User decisions, November 16 and 17, 2026 (D190). Implemented November 17,
2026; see [testing](testing.md#birthplace-group-november-17-2026). The
dimension is named `birthplace` because `origin` and `destination` are
transition columns.

**Always on.** Every model with an immigration process has a group
dimension `birthplace`, created automatically and marked `automatic`. Its
categories are `"native"`, any categories the population already has, and
the stream names in declaration order. It replaces the former `stream`
column of `persons`: an immigrant's `birthplace` is their stream's name,
and everyone else's is `"native"`. `arrival_time` remains the arrival
record. Stream names can't be `"native"` or `"birthplace_ancestry"`, and a
model with immigration can't declare its own `birthplace` group. Populations
that already have the automatic dimension keep it in later runs, with or
without immigration, including when `set_groups()` redeclares the user's
dimensions without it.

**Category: personal status.** Newborns are always `"native"` (fixed
inheritance, not configurable). An arrival's category is their stream's
name; the template's own `birthplace`, if any, is ignored, and
`assignment` can't set it. The category therefore records personal
immigration status and stream, which D182 rate strata, weight models,
matching and `by = "birthplace"` observation use like any group. Each
arrival gets an entry row in `birthplace`, so group duration in it starts
at arrival.

**No transitions.** Birthplace is fixed at birth or arrival, so the
dimension can't have a transition process (D198). Naming it in
`set_group_transition(group = )` is an error at `run_simulation()`, and a
transition without `group` never falls back to it. A model without
immigration may declare its own `birthplace` group and move people in it.

**Only-dimension defaults.** The automatic dimension never counts as a
model's only dimension: transitions without `group`, bare matching settings,
default matching pools, transition rates without `group`, `transitions`
input without a `group` column and observation by a single dimension all
ignore it. It is used only when named, and never by transitions.

**Ancestry: heritage.** `immigration(..., birthplace_ancestry = TRUE)` (and
`set_immigration()`) tracks ancestry in `birthplace`; the default is
`FALSE`, as for `group()`. Shares follow the existing ancestry rule: an
immigrant has all their ancestry in their stream, and a newborn has the
average of the parents' shares, with an unknown father's half in the
unknown share. Weight models can use `ancestry_birthplace_<category>`. The
template's own birthplace ancestry is ignored.

**Starting populations and restarts.** A population without `birthplace`
starts with every person `"native"` and, with ancestry tracking, all their
ancestry native. `"native"` means not an immigrant during the simulation,
as expected of a synthetic population; users who know more supply the
`birthplace` column (and ancestry). A population with the dimension keeps
its categories and shares; streams added on restart add categories, and
categories of streams no longer present are kept for the people who
arrived through them. With immigration and `birthplace_ancestry = FALSE`,
existing birthplace ancestry is dropped; turning it on for a population
without shares starts everyone at 100% of their category, as for other
groups. A restart without immigration keeps the dimension's ancestry
setting.

**External relatives** have `NA` birthplace, as in groups templates lack.

**Worked example.** Streams `asia` and `latin_america`; the starting
population is native. An Asian immigrant woman has `birthplace = "asia"`
and ancestry 100% asia. Her daughter with a native man is `"native"`, 50%
asia and 50% native. The daughter's son by a native man is `"native"`, 25%
asia and 75% native. His child with a Latin American immigrant woman is
`"native"`, 12.5% asia, 37.5% native and 50% latin_america. A child of the
immigrant woman with an unknown father is `"native"`, 50% asia and 50%
unknown. A rate table with `birthplace` rows `native`, `asia` and
`latin_america` gives the immigrants their own schedules from arrival and
their descendants the native schedule; a weight model on
`ancestry_birthplace_asia` lets heritage act within those cells.

## Migration facts, family history and classifications

Keep migration history separate from derived generation/origin classifications
(D181). One generation number cannot describe every mixed-stream pedigree.

| Layer | Meaning |
| --- | --- |
| Personal migration facts | Whether a person actually immigrated into the destination, when, and through which stream |
| Family migration history | Connections to migration histories through parent–child relationships |
| Derived classifications | Named variables under a specified rule, usable in demographic mechanisms and analysis |

For example, a locally born child may have a mother who immigrated through
stream A and a locally born father whose parents immigrated through stream B.
The child has an immigrant parent connected to A and immigrant grandparents
connected to B. Preserve enough information to distinguish those connections;
do not silently select one origin or a universal generation rule.

An immigrant arriving as a child has their own arrival history. A source
relative retained only to support an arriving person's reproductive history
has not thereby immigrated. Source genealogy alone does not establish a
relative's migration history with respect to the destination.

Migration facts, group membership and the package's ancestry proportions are
distinct. Classification may use several streams without requiring a group
category for each stream. Mappings among these concepts must be explicit.

## Concise access in custom demographic models

Users should be able to define a migration classification once and reuse its
named outputs in partner matching, group inheritance, weight models and
analysis (D181). Custom functions should receive derived values in their
documented inputs rather than having to traverse pedigrees. This is a new
interface requirement, not existing callback functionality.

Illustrative outputs include an immigrant indicator, personal arrival stream,
a generation summary under a chosen rule, and origin-specific generation
classifications. These are examples, not selected field names, default
generation rules or a mandatory fixed set of variables. A simple model may
need one summary, while a mixed-origin model may need several outputs.

Specify classification inputs, output names and types, unknown-history
behavior, and initialization/update timing, especially at birth and
immigration. Distinguish unknown history from a known absence of connection
to a stream. Custom mechanisms must retain existing state-consistency and
event-invalidation requirements; unrestricted state-changing callbacks are
not authorized.

## Rate stratification

Support rate schedules stratified by personal immigration status and arrival
stream (D182), without requiring a generational classification in the rate
table. Illustrative categories distinguish people who have not themselves
immigrated, immigrants through stream A, and immigrants through stream B.
These are conceptual categories, not selected column names or encodings.

A destination-born child of stream-A immigrants belongs to the category of
people who have not themselves immigrated. Their family migration history
still connects them to stream A for D181 classifications and custom models.
Arrival-stream status records personal immigration, not inherited origin.
Under D190, starting populations default to `"native"`, meaning not an
immigrant during the simulation; users supply `origin` when they know
more.

Status/stream-specific schedules set aggregate rate targets for those strata.
Derived generational variables in normalized weight models instead influence
who experiences events within rate cells, without changing cell intensities.
They do not independently impose generation-specific aggregate rates.

Whether to use an ordinary declared group or a built-in migration stratifier
remains open, as do process coverage, schedule validation and category names.
General generation-specific rate stratification is not selected. This design
does not change the current implementation's requirement that categorical
rate strata use declared group dimensions.

## Template cuts

User decisions, November 17, 2026 (D196). Implemented November 17, 2026;
see [testing](testing.md#template-cuts-november-17-2026).

**Interface.** `immigration_stream(..., step = NULL)`. With `step = NULL`
(default), migrants are drawn from the template's living population at its
`end_time`, as before. With a positive number of years, the run is divided
into intervals `[start + k * step, start + (k + 1) * step)`, the last ending
at the run's end, and an arrival in interval k is drawn from the template as
it stood at the cut `c_k = start + k * step`, the start of its interval. The
template and the destination share one calendar.

**Coverage.** Every cut must lie within the template's simulated history:
from the start of its first recorded run (`unrealized_since`) to its
`end_time`. Before its first run a template has only its starting
population. A cut outside that range, or a template without a recorded run,
is an error naming the cuts and the template's simulated period.

**The template at a cut.** An internal helper truncates the template to its
state at time c:

- People born after c are removed, as are the template's own immigrants who
  had not arrived by c, with unions and external relatives only they
  reference. Deaths after c are undone.
- Unions starting after c are removed; unions ending after c are live at c
  (no end time or reason). A cohabitation converted after c is a live
  cohabitation, and the marriage it became is removed.
- Transitions after c are removed, and each person's groups are their
  categories at c from their transition history.
- Template external relatives born after c are removed; those alive at c
  are censored at c (`censor_time` no later than c, deaths after c undone).
- Ancestry is kept for remaining people, and the truncated `end_time` is c.

Parity, previous-union state, households, eligibility, weights, custom
assignment, imported histories and external relatives then follow from the
truncated template exactly as for `step = NULL`, with c in place of the
template's `end_time`: a migrant arrives as they were at c, with every date
shifted by the arrival time minus c (less than `step`). Template IDs are
the same at every cut, and custom weight functions' `bin` gains `cut`.

**Weight schedules (D197).** A stream's `weights` may be a calendar
schedule in the D123 format. Slices then start at every cut and every
schedule boundary in the run, each with its latest cut's template and the
weight model active at its start; a custom function's `bin` holds the
interval's `time_start` and `time_end`.

**Engine.** Each cut's parties, selection probabilities, assignment
probabilities and imported state are prepared before the run. Cut times
inside the run are calendar boundaries: at each, the stream switches to the
cut's parties and redraws its arrival clock (memoryless). Household streams
divide intensity by each cut's E[size], so the person intensity holds as
household composition changes. Intensity tables are unaffected.

**Cost.** Setup work and memory grow with the number of cuts; at report
scale (50 yearly cuts of a template of about 25,000) setup is estimated at
tens of seconds, to be measured.

**Worked example.** A template simulated from 0 to 200 serves a destination
run from 100 to 150 with `step = 10`. Arrivals in [100, 110) are drawn from
the template as of 100, those in [110, 120) as of 110, and so on to [140,
150) from 140. A woman aged 30 in the template at 120, married at 115 and
mother of a child born at 118, arriving at 124.5 is 30 at arrival, born at
94.5; her marriage starts at 119.5 and arrives live, and her child, aged 2
at the cut, arrives aged 2 if a dependent. Her template births after 120
are not part of her history. A template simulated only from 50 to 140
still serves this run, whose last cut is 140, but not a run from 100 to
160, which needs a cut at 150; nor does a template whose first run starts
after 100.

## Template truncation (deferred)

Truncating a template to its state at a cut time (D196) is internal and
distinct from the following. D180 recorded a requested helper to truncate
templates to living people and
their descendants. D184 defers it: the user is not certain it is worth its
complexity, and migrants no longer depend on it because external relatives
carry the records parity and previous-union status need. Templates are used
as supplied. If revisited, truncation would need retained records for
deceased previous partners of living people and for deceased ancestors
within the kinship-screening depth, and summaries preserving migration
classifications. See the
[burn-in helper discussion](roadmap.md#living-population-extraction-after-burn-in).

Observation is arrival-aware without truncation (D185, D186): imported
events at or before arrival are never destination events or exposure.

## Open decisions

The following choices are open; listing an option does not select it.
Source sampling (D189, D195), imported state and observation (D186),
households (D187), template and weight schedules (D196, D197), intensity
(D194), rate strata (D190) and the interface are settled.

| Area | Decisions still needed |
| --- | --- |
| Calendar inputs | Calendar schedules of `migrant_assignment()` |
| Households | Multigenerational household rules beyond D140's |
| Migration history | Source provenance, history depth and restart persistence |
| Duration of residence | Use in formulas and matching models; custom weight functions already see `arrival_time` (D195) |
| Classifications | Recency shares by generation depth, awaiting quantitative characteristics; generation itself is user-defined (D192, D193) |
| Template preparation | A truncation helper, deferred by D184 |
| Emigration | No mechanism selected |

No implementation milestone or scientific acceptance threshold is established
here.

## Worked examples for review

These illustrate accepted meanings and proposed validation cases, not executed
results or additional demographic defaults.

**Age and replacement.** A source snapshot ends in 2000 and contains a person
born in 1975. Sampling that person for an arrival in 2050 gives age 25 and
destination birth time 2025. Sampling the same source person again in 2051
gives another age-25 person with birth time 2026 and a different destination
ID. Neither draw removes or ages the source record. Imported histories shift
by the same 50-year offset (D184).

**Multiple flows.** A stream with intensity 100 and another with intensity
40 have expected individual arrival counts of 100 and 40 in a one-year
interval. Those are expected flows, not quotas or expected final population
increments: arrivals may subsequently die or have children. These flows do
not depend on the destination's current population size.

**Imported history and external relatives.** A template ends at year 100.
A woman in it was born at 70 and married a man born at 66 at 88. Their
children were born at 90 (died 92) and 95 (alive), and her husband died at
98. She arrives at 2035.4, an offset of 1935.4, in a household stream
with `household(age_independence = 0)`, under which her living child is
independent, so she arrives alone. In the
destination she is born at 2005.4 (age 30). Her union runs from 2023.4 to
2033.4, ended by death, and is pre-migration. She arrives widowed, 2 years
after that union ended, at parity two, with last birth at 2030.4. Her
husband and children are external relatives born at 2001.4, 2025.4 and
2030.4. The husband and first child have death and censor times 2033.4
and 2027.4; the second child has no death time and censor time 2035.4. The
union's `pre_migration` is TRUE.
None of these records adds destination births, deaths or exposure, and her
destination exposure begins at 2035.4.

**Phase 1 arrival.** The destination declares `ethnicity` (tracking
ancestry) and `region`. A stream's template ends at 2000 and declares
`ethnicity`, with ancestry, and `education`; its `assignment` gives `region`
a `"probability"` rule of 0.6 west and 0.4 south. An unpartnered woman born
in 1975 is selected for an arrival at 2050.3. Her destination record has a
fresh ID, birth time 2025.3 (age 25), unknown parents, `arrival_time`
2050.3 and the stream's name. She keeps her template ethnicity and its
ancestry shares; education is dropped; her region is west with probability
0.6. If the destination tracked region ancestry, she would have 100% of her
assigned region. A life table over 2050–2060 counts her exposure from
2050.3, beginning at age 25; her birth is not a destination birth.

**Personal status versus family origin.** An immigrant through stream A and
a locally born person with immigrant parents from stream B have a child in
the destination. The child has not personally immigrated and does not inherit
stream A or B as a personal arrival stream. Their family history supports
connections to both streams. A named classification can summarize these
connections for custom inheritance or matching under a user-specified rule;
no universal numeric generation is assigned here.

## Validation plan

Validation must follow the completed specification. For individual arrivals
with a valid, nonempty source, the initial analytical and accounting cases are:

| Case | Expected check |
| --- | --- |
| Constant intensity lambda for T years | Across replications, arrival counts follow Poisson(lambda * T), and waiting times follow Exp(lambda) |
| Zero intensity or zero duration | No arrivals |
| Replacement draws | Repeated source selections get distinct destination IDs; the original template remains unchanged |
| Destination accounting | Final living population equals initial living population plus actual arrivals plus destination births minus destination deaths |
| Imported histories | Historical source births/deaths are excluded from destination event counts and pre-entry destination exposure |
| Imported state | Arrivals' parity, last birth, previous-union status and end, and group duration equal their template-endpoint values, with dates shifted by the offset |
| External relatives | References resolve; external relatives never enter cells, pools, exposure, event counts or population totals |
| Phase 1 selection | Across replications, selection frequencies match normalized weights; scaling all weights leaves them unchanged |
| Phase 1 arrivals | Age at arrival equals template age at `end_time`; parents unknown; never partnered; parity zero |
| Phase 1 groups | Shared dimensions copy template values; ignored dimensions are absent; assigned categories follow fixed, probability or custom rules across replications |
| Phase 1 ancestry | Template shares are copied where tracked; otherwise 100% of the migrant's category |
| Phase 1 observation | No exposure, events or population count before arrival; a migrant's own birth is not a destination birth |
| No streams | Seeded output identical to a run without the immigration code, with no arrival columns |
| Piecewise intensity | Expected arrivals equal the integral of intensity; 100 per year for 10 years then 50 for 10 years gives 1,500 |
| Rate intensity | With no births or deaths and individual arrivals, the destination grows as a pure-birth (Yule) process: E[N(T)] = N(0) e^(rate T); `calculate_immigration_rates()` recovers the input rate within Poisson error |
| Intensity tables | One of `count` or `rate`; bare numbers, both, neither and extra columns are errors; D118 coverage rules apply |
| Phase 2 import | Unions, partners and children arrive with shifted dates and fresh IDs; a child's other parent is known only when a previous partner |
| Phase 2 derived state | Imported parity selects parity-specific fertility; imported last births keep spacing; imported widowhood selects remarriage schedules; group duration continues from the template or starts at arrival |
| Phase 2 observation | External relatives add no counts; pre-migration unions add no union-years; entry rows are not counted as transitions |
| Phase 2 restarts | IDs stay unique; parity from external children persists |
| Phase 3 households | Every living template person is in exactly one household; a party is its carrier, partner and dependents, recursively |
| Phase 3 selection | Across replications, household frequencies match normalized carrier weights |
| Phase 3 intensity | Expected persons equal `intensity * T`; parties equal `intensity * T / E[size]`; person-count variance matches the compound Poisson value |
| Phase 3 arrival | Members share arrival time, stream and assigned categories; internal unions and parent links connect destination members; each external relative is created once per party |
| Phase 3 ongoing unions | A live union enters dissolution and conversion cells at its imported duration; partners do not search; union calculations count it from arrival and exclude its formation |
| Phase 3 restarts | Live imported unions and party links persist |
| Cloning streams | Anyone alive can be drawn; clones arrive never partnered at parity zero with no imported records, an entry row per dimension, and unknown parents; templates with unions or transitions warn |

Additional tests should cover multiple streams, age preservation, template
schedule boundaries, future group categories, mixed-stream pedigrees,
classification inputs in custom models, unknown histories, and classification
and external-relative preservation through restart. Source and destination
should be reusable across independent seeded runs under the existing RNG
contract.
Statistical test methods must be justified when implementation begins; no
acceptance thresholds or successful validation results are claimed here.
