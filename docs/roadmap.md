# Development roadmap

Updated October 29, 2026. Milestones 1–4 are complete by user acceptance
(D029, D050, D077, D093). The group validation report is the accepted milestone 3
assessment. The user proposed the following sequence on
October 20, 2026. These are planning priorities, not authorization to implement
unspecified demographic rules or a selected public API.

| Milestone | Proposed scope | Status |
| --- | --- | --- |
| 4 | Union dissolution and age/duration remarriage; possibly cohabitation | In progress (D078–D083); cohabitation undecided |
| 5 | Parity-specific fertility and fertility heterogeneity | Proposed |
| 6 | Time-varying rates, matching, and inheritance within a simulation run | Proposed |
| 7 | Group conversions | Provisional; distinguish identity from characteristics and region |
| 8 | Polished summaries, diagnostic checks, and post-simulation helpers | Proposed |
| 9 | Immigration from outside the system, tracking generation | Tentative |
| 10 | Internal migration between regions, with unions restricted to the same region | Tentative |
| 11 | Individual characteristics with custom inheritance and life-course change, affecting rates and union formation | Tentative |

## Milestone 3: groups (complete)

Accepted October 27, 2026 (D077); see the
[validation report](../validation/milestone3.qmd) and [testing](testing.md).

Group schedules and weighted baseline sampling originated under D052;
see [groups](groups.md). Random-parent inheritance with maternal fallback
is now the default (D076). Group
observations (D053) and proportional input-rate multipliers, including separate
partnered/unpartnered fertility factors (D054), are also implemented. Ancestry
proportions are tracked separately with parental averaging (D055); see
[ancestry](ancestry.md). Group matching formula terms and symmetric pair effects
are implemented (D056/D075); calibrated male-group weights, delta adjustments,
and exact full-pool selection follow D067–D075. Directed pair utilities and
method selection are removed. See [group matching](group-matching.md) and
[weighted-choice matching](weighted-choice-matching.md).
Built-in/custom newborn identity assignment is implemented under D058; see
[group inheritance](group-inheritance.md). Multiple simultaneous identities
remain future work. Custom rules can already assign discrete mixed identities
such as AB with independently supplied AB demographic schedules.

The user selected the following directions for group support:

- Include intergroup unions in formula-based matching, alongside age difference.
- Allow proportional rate differences by group relative to a baseline schedule,
  avoiding a complete duplicate rate table for each group. An illustrative
  fertility multiplier is 1.1 times the baseline, not a demographic default.
- Support flexible child group assignment when parents have different groups.
  Parental, ancestry, hierarchy/passing, and custom rules are now implemented
  under D058, with random-parent assignment as the default (D076).
- Consider tracking each person's proportional ancestry by group separately
  from assigned group identity.

### Accepted organization (October 20, 2026; D051)

Group identity and optional ancestry data belong in the population. Group effects
belong in the relevant demographic processes, including group-specific rate
schedules and proportional multipliers, and group predictors in matching. There
will be no separate general group process for milestone 3. Fertility will hold
an inheritance model object that specifies newborn group assignment and, where
supported, ancestry inheritance.

Inheritance rules now use character-string modes in an object attached to
fertility. Built-ins cover parental, ancestry, and hierarchy assignment with
passing; custom functions support more complex rules and declared identities
(D058). Ancestry remains separate from identity. Multiple-identity representation
and specialized model constructors remain unselected.

Mortality, fertility, and rate-based searches now accept group multipliers;
fertility factors apply before compensation and may differ by union status
(D054). Matching now supports female/male group factors with explicit fixed
levels, numeric same-group indicators, and additive pair effects (D056).
Ancestry initialization and unknown-parent treatment are
now specified under D055: starting identities supply pure proportions unless
shares are provided, and missing parental contributions enter an unknown share.
These are expected genealogical contributions, independent of discrete group
identity. They do not select a majority-ancestry identity classification.

## Milestone 4 complete: dissolution, remarriage and cohabitation

Active investigation began October 27, 2026 at the user's request on branch
`milestone4-union-dissolution`. The first
[SOCSIM dissolution/remarriage audit](socsim-dissolution-audit.md) is a pinned
static source trace, not executable validation or a new implementation decision.
The user accepted milestone 4 on October 29, 2026 (D093), based on the
[validation report](../validation/milestone4.qmd) and rendered assessment.

The user selected duration-dependent dissolution on October 27, 2026 (D078),
superseding the earlier female-age-based direction. Use annual hazards and one
clock per union carried by the woman, with shared schedules or symmetric
group-pair schedules and multipliers (D079). Exact same-person ties prioritize death, birth, then
dissolution. See the [dissolution specification](union-dissolution.md).

The first implementation uses a separate optional `union_dissolution` process
with `dur_start`/`dur_end` input bands. Leaving it absent retains death-only
endings. See [union dissolution](union-dissolution.md) for the interface and
worked examples. The accepted scope retains the limitations described below.

Consider allowing dissolution risk to vary with couple characteristics in
addition to union duration. The user's example is higher divorce risk for couples
with large age differences. This is an extension to investigate, not an accepted
effect or coefficient: the functional form, supported characteristics, and how
effects relate to baseline duration-specific rates remain to be specified.

Female search rates now optionally distinguish `previous_union`
as `none`, `dissolved`, or `widowed` (D080), with all three strata required
within each scheduled group. Separate marriage observations use the same
categories. Age or duration since the latest ending can be selected per
previous category in the same table (D083); `none` remains age-specific. See
[remarriage](remarriage.md). Optional cohabitation is now implemented in rate mode under D086, with typed
matching/rates and linked same-partner conversion. Marriage-order dependence
and richer couple covariates remain deferred. Acceptance is the user’s
assessment of the report, distinct from software checks alone.

## Milestone 5 design: parity-specific fertility and heterogeneity

The user selected optional half-open `parity_start`/`parity_end` bands with Inf
for the open category (D094). The HFD helper prepares conditional rates in
that format via `by_parity = TRUE`; D099 implements parity-specific execution
with history-derived counts and SOCSIM-style spacing compensation. See [parity fertility](parity-fertility.md).

Add fertility schedules that depend on parity as well as the existing schedule
axes. The user selected initialization from recorded maternal birth histories
(D098), with zero for women lacking recorded births and no overriding parity
column. D099 selects preceding-parity spacing compensation, birth-driven schedule
changes and calibration updates. The approximation requires demographic
assessment. D100 adds parity observations via `calculate_fertility_rates()`;
heterogeneity remains later work.

Consider person-level fertility differences to represent persistent high- and
low-fertility individuals and increase outcome variance beyond the current
model. The user reports SOCSIM has a mechanism for this; its precise behavior
has not been audited for this task. Distribution, persistence, association with
other attributes, inheritance and calibration to aggregate rates are unselected.
No heterogeneity parameter or demographic default is adopted here. See the
[distribution and inheritance discussion](fertility-heterogeneity-notes.md) and
[pinned SOCSIM heterogeneity audit](socsim-fertility-heterogeneity-audit.md).

## Proposed milestone 6: time-varying inputs within a run

Allow demographic schedules to change during a single simulation, rather than
requiring population restarts. Specify calendar-time inputs, change boundaries,
which processes may change, event ordering at boundaries, and invalidation of
pending events. Preserve birth spacing from actual births and address the
interaction with fertility compensation after rate changes (D045). No rate-change
API, interpolation rule, or demographic default is selected.

### Input representation proposal (October 21, 2026)

The user requested retaining the following approach for milestone 6, allowing
inputs to change without restarting the simulation:

- Keep rate schedules as flat tibbles, adding `time_start` and `time_end`
  alongside age, sex, group, union status, and other supported rate strata.
  This preserves ordinary filtering, joins, plotting, and inspection.
- For structured inputs that are not already rate tibbles, use a schedule
  tibble with `time_start`, `time_end`, and a list column containing a complete
  model object for each interval. Matching models and group inheritance models
  are the current examples. Changing the whole matching model permits changes
  to formulas, coefficients, and group exogamy effects without making each
  internal parameter separately time-varying.
- Keep the time schedule outside each individual model object. A shared
  scheduling helper could validate intervals and select objects, while rates
  retain their flat representation. Constructor names and public arguments
  remain to be designed; no generic framework beyond these needs is selected.

For example, a matching schedule could have the following shape (illustrative
data structure, not an implemented process interface):

```r
tibble::tibble(
  time_start = c(2000, 2050),
  time_end = c(2050, Inf),
  matching = list(matching_model_before, matching_model_after)
)
```

An inheritance schedule would similarly hold complete inheritance models in
a list column. Shared half-open calendar intervals `[time_start, time_end)`
are a proposed convention. Selecting matching at search/father-selection
execution and inheritance at birth, without rewriting existing unions or
identities, is a proposal requiring confirmation in the scientific specification.
Before implementation, resolve coverage, gaps, overlaps, boundary ties,
validation of all scheduled models, and interactions with pending events.
Rate sampling must account for future rate boundaries rather than assume the
rate at scheduling time persists indefinitely. This note records a future
design direction, not implementation authorization or accepted timing rules.

## Provisional milestone 7: group conversions

Consider changes in group identity during life. The user is uncertain about the
scope: many applications should instead use separate individual characteristics
(milestone 11), while location should use a dedicated region mechanism
(milestone 10). Do not treat groups as a universal state variable. Conversion
rules, timing, ancestry implications, and effects on pending events remain open.
On October 27, 2026, the user deferred fixed newborn group assignment (for
example, all children starting at a minimum education category) to this
discussion. No `fixed` inheritance method or API was adopted. A custom
inheritance function already handles assigning mixed-parent births to one
specified category.

## Proposed milestone 8: summaries, diagnostics, and helpers

The user also requested future model diagnostic tools (October 20, 2026),
including checks that help identify unintentionally omitted group-pair settings.
No diagnostic interface or implementation is selected.

Add custom print/summary displays for `hammel_model` and simulation results
(the user's proposed `hammel_simulation` interface). Aim for a clear, polished
console presentation inspired by usethis, making configured models and results
easy to inspect. Specific contents and formatting remain to be designed.

Simulation currently returns a `hammel_population` with a `simulation` metadata
component. The summary proposal does not itself select a new result class or
authorize changing that return contract. Keep underlying tables and settings
available as ordinary inspectable R objects.

### Post-simulation helpers

Expand the existing extraction/summary functions into a coherent user toolkit
for useful quantities from simulation histories. Existing functions cover Lexis
counts, fertility, life tables, population sizes and recorded kin. Identify
additional quantities and exposure/observation definitions from actual workflows.
Status-specific exposure, fertility and mortality summaries use
`calculate_fertility_rates(by_union_type = TRUE)`; formation counts and rates by
previous union use `calculate_union_rates()` (D080).
See [demographic analysis](demographic-analysis.md) for their observation rules.

### Living-population extraction after burn-in

On October 29, 2026, the user requested a future helper to truncate a simulated
population to living individuals, so a burn-in run can provide a new starting
population without retaining the full initialization history. The desired result
preserves correct reproductive histories and parent–child relationships, while
allowing the deeper ancestral kinship history to be discarded. This is a
milestone 8 proposal; no function name or implementation is selected.

Before implementation, specify how to preserve births to deceased children,
parity and last-birth dates after removing deceased person records. Simply
filtering living persons would lose some recorded births and conflict with
history-derived parity initialization (D098). Also specify handling of removed
parent IDs, retained parent–child links, ongoing unions, previous-union histories,
ancestry, observation endpoints and restart validation. The helper must state
which kinship information has been discarded; it cannot imply complete pedigree
or historical exposure after truncation. The representation of minimal history
needed for correct restarting remains to be designed.

## Tentative milestone 9: immigration from outside the system

Add immigrants within a running simulation, rather than requiring a restart.
The user envisages specifying numbers and an age/sex distribution and tracking
individual immigrant generational status. Timing (batch dates versus arrivals),
meaning of supplied numbers, attribute/kin/union initialization, ancestry and
generation definitions, and interactions with group assignment remain open.
No emigration process is implied by the immigration request.

## Tentative milestone 10: internal migration and regions

Represent location with a dedicated `region` column, separately from group
identity. Allow internal migration using rates between regions, and restrict
union formation to partners within the same region. The region schema, migration
clock and rate strata, individual versus household movement, and treatment of
existing unions when partners move remain to be specified. No migration default
or rule about relocating existing partners is selected.

## Tentative milestone 11: individual characteristics

Consider characteristics that can affect rates and matching without functioning
as definitive group identities. Education is the user's example. Users would
provide custom functions for inheritance and life-course change. Characteristics
would be exposed to modify demographic rates and union formation. Initialization
and callback contracts remain to be specified. Define permitted changes and event invalidation, not arbitrary state
mutation. No distinction between categorical/continuous attributes or transition
model has been fixed by these notes.

## Two-sex union formation and union API

Discussed October 24, 2026, without an assigned milestone; see
[two-sex union models](two-sex-union-models.md) and D073. Possible
extensions: a participation function making a woman's chance of forming a
union depend on market conditions (overlapping D044); a Schoen harmonic-mean
model with a market-level union clock, possibly modified along Dagsvik's
lines; possibly Choo–Siow. These would let composition affect union timing
and capture marriage squeezes, at the cost of not necessarily reproducing
input rates; users may want either class. Implemented rate-mode cohabitation (milestone 4) and
nonmarital father selection must be designed consistently with any two-sex
model. The union-formation API stays unchanged until a single working
implementation incorporates all union-formation components. No model,
parameterization or API has been selected.

## Preparing input schedules

Provide more helpers for building required input rates, including retrieving or
using model life tables and splitting overall fertility into partnered and
unpartnered schedules. No source, model-life-table family, fitting method or
split rule has been selected. Overall fertility alone does not identify two
status-specific schedules: decomposition needs additional information or
assumptions, such as age-specific status proportions and relative fertility.
Keep empirical rates, inferred rates and user assumptions inspectable with
provenance. Avoid automatic downloads as part of simulation execution.

## Population feedback

Consider allowing fertility or other demographic rates to respond to population
size or age composition. Discussed October 17, 2026 as future work, without an
assigned milestone. No feedback rule, demographic response, or public API has
been selected.

The event queue can coexist with population summaries: a living-person count
could be updated at births and deaths, while age composition could be calculated
from living persons' birth dates when needed. The main scheduling issue is that
a population-wide rate change makes pending events sampled under the previous
rates obsolete. Immediate feedback after every birth or death could therefore
require frequent rescheduling across the affected population.

A proposed approach is periodic feedback, with the user defining the period
length in years. Schedule feedback boundaries explicitly, so they are processed
even when no demographic event occurs during a period. At each boundary, assess
the current population, calculate any rate adjustment, and invalidate and
reschedule affected individuals' pending events from that boundary before
executing later events. Hold the adjustment fixed until the next boundary.
For fertility changes, reconsider each affected woman's next competing event,
preserving her current state and spacing restrictions from actual births.
This would retain continuous demographic event times; the feedback period is
a response interval, not a simulation execution step.

Before implementation, specify the population signal and reference population,
response function, any lag, first feedback time, boundary tie ordering, and
which processes and individuals are affected. The response interval is a
demographic assumption as well as a computational choice. Interaction with
fertility compensation after rate changes remains an open issue (D045).
Any user-defined feedback mechanism must obey state consistency and event
invalidation rules; arbitrary callbacks must not mutate simulation state.

For a precedent, see Mike Murphy's
[Bringing behaviour back into micro-simulation: Feedback mechanisms in demographic models](https://www.demogr.mpg.de/Papers/workshops/010221_paper10.pdf)
(2001 workshop paper). Murphy extended SOCSIM with population-pressure feedback
through nuptiality, evaluating population size at 25-year intervals. This is a
reference for investigation, not an adopted response rule or period length.

## Previously deferred items

- Previous-father continuity at subsequent births or union formation (D042).
- Declining poor nonempty candidate pools and waiting (D044).
- Transition-aware fertility compensation for status/rate changes (D045).
- Multiple concurrent unions remain unassigned. Parity-specific fertility is
  proposed for milestone 5; non-death dissolution is complete in milestone 4.

R formula-fallback performance comparisons are now included in the user-authored
milestone 2 validation report, accepted by the user on October 19, 2026 (D050).
Initialization and exact checkpoints remain separate future topics. Calendar-varying
schedules now fall under proposed milestone 6 above.
