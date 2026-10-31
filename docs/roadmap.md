# Development roadmap

Updated October 31, 2026. Milestones 1–6 are complete by user acceptance
(D029, D050, D077, D093, D105, D116). **Milestone 7 calendar-varying rates
are implemented under D118; scientific assessment and model schedules remain
pending.** Finish milestone 7 with the current engine. The user selected
binned-group hazards and formula weights for milestone 8 (D121).
Future entries record planning priorities, not authorization to implement unspecified demographic behavior
or an unselected public API.

| Milestone | Scope | Status |
| --- | --- | --- |
| 4 | Union dissolution, age/duration remarriage, and cohabitation | Complete (D093) |
| 5 | Parity-specific fertility and fertility heterogeneity | Complete (D105); limitations retained |
| 6 | Group transitions, including categorical characteristics represented by groups | Complete October 30, 2026 (D116); final render active at acceptance |
| 7 | Time-varying rates, matching, and inheritance within a simulation run | Rates implemented (D118); scientific assessment pending; matching/inheritance schedules deferred |
| 8 | Binned-group hazards for all rate processes; custom formula weights on all rates | User-selected direction (D121); specification/API and implementation pending |
| 9 | Polished summaries, diagnostic checks, and post-simulation helpers | Proposed; moved from former milestone 8 (D121) |
| 10 | More than one group column, potentially including regions and internal migration | Tentative (D121); schema and mechanism contracts unselected |
| 11 | Immigration from outside the system, tracking generation | Tentative; moved from former milestone 9 (D121) |
| 12 | Two-sex marriage model | Tentative (D117/D121); model and API unselected |

Milestone 6's scope is individual group transitions. Household moves are
deferred beyond milestone 6 (D114, superseding D109). D121 removes the separate
regions milestone and places multiple group columns at tentative milestone 10;
regions/internal migration and household moves remain possible applications
requiring their own specifications. D106–D108 select annual hazards, directed competing
age- or duration-specific schedules with optional strata,
and no transition-count target adjustment. The implemented contract is in [group transitions](group-transitions.md).
D110 requires `group_origin`, `group_dest`, `rate`, and one common age or
duration band type across the transition table, with optional `sex` and
`union_type` columns; absent stratification columns apply universally.
Grouped simulation output must always include a `transitions` history tibble,
even with no configured transition process or no observed moves (D111).

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

## Milestone 5 complete: parity-specific fertility and heterogeneity

The user selected optional half-open `parity_start`/`parity_end` bands with Inf
for the open category (D094). The HFD helper prepares conditional rates via
`by_parity = TRUE`. Parity is initialized from recorded maternal birth histories
(D098), with zero for women lacking recorded births and no overriding parity
column. D100 adds observed parity-specific rates through
`calculate_fertility_rates(parity_breaks = ...)`.

D101 cell-accounted fertility replaces D099's preceding-parity spacing
approximation: cells generate births at input rate times all their members,
then select an eligible mother. This also handles composition changes during
start-up, restarts and union-type transitions. Cells without eligible women
report unrealized intensity. See [parity fertility](parity-fertility.md) and
[cell accounting](fertility-cell-accounting.md) for the contract and limitations.

D102 implements persistent positive lognormal selection weights through
`set_fertility(..., hetfert = het_fert(cv, alpha))`. Dispersion and maternal
inheritance have separate inputs; supplied person weights persist, otherwise
simulation draws them. Weights redistribute births within cells without
changing their target intensities. Reproductive selection can change cohort
weight distributions; no generational normalization is imposed. See the
[heterogeneity specification](fertility-heterogeneity-notes.md) and
[SOCSIM audit](socsim-fertility-heterogeneity-audit.md). The user accepted
[the milestone 5 report](../validation/milestone5.qmd) on
October 29, 2026 (D105). Unrealized-intensity/sparse-cell limitations, incomplete
birth histories and reproductive selection remain explicit; completion does
not imply universal rate recovery or stationary cohort weight distributions.

## Milestone 6: group transitions (complete)

Group membership can change during life. Use **group transitions** for
this feature, reserving conversion terminology for cohabitation becoming
marriage. Most categorical characteristics can use groups and group-specific
rates, matching and inheritance rather than a separate general individual-
characteristics milestone. Education is one example: people could start in a
minimum category and transition to higher categories over life.

Multiple group columns that represent different dimensions simultaneously are
now tentative milestone 10 (D121). No schema or public API is selected.
How multiple group types interact in rate stratification, matching, newborn assignment, transitions and
ancestry remains to be specified. The roadmap does not introduce an arbitrary
attribute-mutation callback or a continuous-characteristic model.

Individual transition schedules, age/duration dependence, eligibility, timing,
event/tracker updates and retained ancestry are implemented under the
[group-transition contract](group-transitions.md). Household moves remain deferred
beyond this milestone (D114), with their unresolved questions retained under
[tentative milestone 10](#tentative-milestone-10-multiple-group-columns) (D121).
The user accepted milestone 6 on October 30, 2026 (D116), based on the [report](../validation/milestone6.qmd); its final render
was active at acceptance and must be preserved.
The former standalone regions/internal-migration milestone is removed; those
applications may be handled through multiple group columns in tentative
milestone 10. Household and geographic rules remain unselected.

Fixed newborn group assignment (for example all children starting in the minimum
education category) was deferred to this discussion on October 27, 2026. No
`fixed` inheritance method or API has yet been adopted. A custom inheritance
function already handles assigning mixed-parent births to a specified category.

## Milestone 7: time-varying inputs within a run

Calendar-varying rates are implemented under the user's D118 authorization;
scientific assessment remains pending. Flat rate tables use optional half-open
`time_start`/`time_end` bands, including fractional years. The engine updates
pending risks, fertility cells and calibration while preserving histories,
spacing and duration origins. See [the rate contract](time-varying-rates.md).

The user selected annual rates first. Matching and inheritance model schedules
are deferred to the second step, along with their interface implementation.
D122 also implements flat calendar-banded `cohab_prob`, using the same interval
and coverage rules as rates. Probability-only boundaries update type-selection
and calibration without changing pending event dates. Other non-rate arguments,
including multipliers, remain static.
The representation and execution-time model selection below were approved in
D118, but model-schedule constructors and public arguments remain to be designed.
Flat rate representation and timing are specified by the implemented D118 contract.

### Input representation (October 21 proposal; October 30 D118 selection)

The user requested retaining the following approach for time-varying inputs,
now scheduled for milestone 7, allowing inputs to change without restarting
the simulation:

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
are accepted under D118. Selecting matching at search/father-selection
execution and inheritance at birth, without rewriting existing unions or
identities, is also approved; model-schedule implementation is deferred.
Before that step, specify the public API and validation of all scheduled models
and apply the approved coverage and boundary rules.
Rate sampling now accounts for future rate boundaries under D118. The
whole-model scheduling interface remains to be designed; the user explicitly
deferred its implementation to the second step.

## Milestone 8: binned-group hazards and formula weights

The user selected this engine change on October 31, 2026 (D121), after milestone
7 is finished with its current engine. Extend cell-based event scheduling to all
six existing rate processes and allow custom formula weights to distribute
individual or union risk within cells while preserving input cell intensities.
Use single-year age or duration subdivisions, retain fractional input boundaries
and all other existing rate strata, and preserve fertility's actual-parity rule.

See [the discussion and design direction](binned-group-hazards.md) for equations,
worked examples, the distinction between rate strata and relative weights,
composition-dependent normalization, possible redundancy warnings, future
user-defined characteristics, process-specific caveats and proposed validation.
Formula APIs, safe evaluation, eligibility/failure rules and scheduling details
still require specification before implementation. This milestone does not
redefine the remaining milestone 7 work or imply its scientific acceptance.

## Proposed milestone 9: summaries, diagnostics, and helpers

Moved from former milestone 8 by D121; its scope and unresolved API choices
remain as described below.

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
milestone 9 proposal; no function name or implementation is selected.

Before implementation, specify how to preserve births to deceased children,
parity and last-birth dates after removing deceased person records. Simply
filtering living persons would lose some recorded births and conflict with
history-derived parity initialization (D098). Also specify handling of removed
parent IDs, retained parent–child links, ongoing unions, previous-union histories,
ancestry, observation endpoints and restart validation. The helper must state
which kinship information has been discarded; it cannot imply complete pedigree
or historical exposure after truncation. The representation of minimal history
needed for correct restarting remains to be designed.

## Tentative milestone 10: multiple group columns

The user assigned the possibility of more than one group column to tentative
milestone 10 on October 31, 2026 (D121). Different columns could represent
identity, education or location simultaneously. Representation, supported rate
strata, matching predictors/restrictions, dimension-specific histories,
transitions, newborn assignment and ancestry interactions remain to be specified.
No naming convention, default dimension or independent-dimension assumption is
selected.

The former separate milestone for internal migration and regions is removed.
Location may be represented as one group dimension, with internal migration as
transitions in that dimension and geographic matching restrictions where
explicitly specified. This replaces the former roadmap's dedicated-region-column
direction; it does not implement or settle any geographic or household behavior.
The following D114 questions remain relevant to a future regional application.

### Household moves: deferred design questions (D114)

D114 deferred household moves from milestone 6 to the former regions milestone.
D121 removes that standalone milestone while retaining multiple group columns
as tentative milestone 10. These concerns are retained for a possible regional
application; no household mechanism is authorized by this placement. The
following questions still require a concrete scientific specification.

- **Meaning of input rates.** An individual's total move hazard includes
  initiation by that person and moves induced by others. For two partners who
  always move together to the same destination, independent initiation hazards
  a and b give each partner total hazard a+b while that arrangement holds.
  Observed individual migration rates generally include accompanying moves;
  treating them as initiation rates for every member would overpredict moves.
  Setting one partner's initiation risk and children's risk to zero can avoid
  double counting in a restricted household arrangement, but imposes an
  initiator rule and becomes awkward with age-specific rates, mixed origins,
  household changes and different destination choices. A household-level clock
  is a candidate alternative; no rate unit, carrier, calibration or API is
  selected. Individual data do not automatically identify household hazards.
- **Conditional rescheduling.** A move triggered by another member can
  legitimately resample a person's pending hazards from the current state,
  conditional on survival/event-free history up to that time. The algorithm
  must integrate forward from attained age and retained durations, not restart
  age, union duration or birth spacing. In integrated-hazard terms a fresh
  exponential threshold has the correct conditional law when the modeled
  state captures the relevant dependence; age/duration-dependent waiting times
  themselves are not generally memoryless. Preserve persistent heterogeneity.
  Unmodeled shared frailty or dependence between moving and another event
  would need an explicit model; resampling does not supply it. See
  [Anderson's integrated-intensity simulation method](https://arxiv.org/abs/0708.0370).
- **Household membership and eligibility.** Partner and parent links do not
  establish coresidence or custody. Specify cohabiting partners, children of
  earlier unions, adult children, multiple residences, and members joining or
  leaving. Formation, dissolution, death and eligibility-age boundaries change
  who can initiate a joint move and who can follow, and may invalidate clocks.
- **Origins, destinations and partial moves.** Partners can start in different
  regions; an A-to-B initiation may cause a C-to-B follower move, or no change
  for a follower already in B. Specify whether following is mandatory, whether
  members may stay, and what happens to existing unions when partners separate
  geographically. Residence duration resets only for actual changes of region.
- **Atomic state updates.** Apply one linked move as a single state change
  before scheduling further events. Update all affected fertility cells,
  matching pools/calibration, region restrictions and union hazards; preserve
  parity, birth histories and ancestry. Invalidate displaced candidates so an
  old follower clock cannot execute a second move or trigger an accidental
  cascade. Fertility uses cell accounting, not an individual birth-clock redraw.
- **Histories, estimands and validation.** Consider a shared move identifier,
  initiator and follower roles alongside individual origin/destination rows.
  Separate initiating household events from individual moves. Individual Lexis
  rates count all moves, so they need not recover input initiation hazards;
  household hazards require a household exposure denominator. Joint moves also
  correlate individual outcomes. Validate conditional rescheduling against
  analytical unchanged-hazard cases, duration preservation and reset rules,
  household accounting, linked-event invalidation and the chosen input-rate
  interpretation before assessing demographic performance.

## Tentative milestone 11: immigration from outside the system

Moved from former milestone 9 by D121. Add immigrants within a running
simulation, rather than requiring a restart.
The user envisages specifying numbers and an age/sex distribution and tracking
individual immigrant generational status. Timing (batch dates versus arrivals),
meaning of supplied numbers, attribute/kin/union initialization, ancestry and
generation definitions, and interactions with group assignment remain open.
No emigration process is implied by the immigration request.

## Tentative milestone 12: two-sex marriage model

Originally assigned tentative milestone 11 on October 30, 2026 (D117), then
moved to tentative milestone 12 on October 31 (D121), following the October 24
discussion; see
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
cell-intensity updates, eligibility and matching calibration must be specified
for any feedback mechanism (D101/D082).
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
- Multiple concurrent unions remain unassigned. Parity-specific fertility is
  complete under milestone 5; non-death dissolution is complete in milestone 4.

R formula-fallback performance comparisons are now included in the user-authored
milestone 2 validation report, accepted by the user on October 19, 2026 (D050).
Initialization and exact checkpoints remain separate future topics. Calendar-varying
schedules now fall under planned milestone 7 above.
