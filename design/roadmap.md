# Development roadmap

Updated November 9, 2026. Milestones 1–9 are complete by user acceptance
(D029, D050, D077, D093, D105, D116, D127, D130, D141). **Milestone 9,
multiple group dimensions, was accepted complete on November 6, 2026.** **Milestone 8 was accepted
complete on November 5, 2026:** cell scheduling for every rate process (D128)
and weight models (D129), retaining documented limitations. On November 5 the
user reordered later milestones (D131). On November 6, D136 superseded that
sequence: finish milestone 9, then prepare the first version bump and public
release. External immigration and two-sex marriage models are post-release
projects. Rich console output, diagnostic bullets, condition classes, and
polished summaries are deferred to a future output-presentation project,
likely after the first release milestone (D151). Other diagnostic checks and
helpers remain proposed; their release scope has not yet been selected.
Future entries record planning priorities, not authorization to implement unspecified demographic behavior
or an unselected public API.

| Milestone | Scope | Status |
| --- | --- | --- |
| 4 | Union dissolution, age/duration remarriage, and cohabitation | Complete (D093) |
| 5 | Parity-specific fertility and fertility heterogeneity | Complete (D105); limitations retained |
| 6 | Group transitions, including categorical characteristics represented by groups | Complete October 30, 2026 (D116); final render active at acceptance |
| 7 | Time-varying rates, matching, and inheritance within a simulation run | Complete November 2, 2026 (D127); current engine and documented limitations retained |
| 8 | Cell scheduling for all rate processes; formula and custom weight models | Complete November 5, 2026 (D130); documented limitations retained |
| 9 | More than one group column, potentially including regions and internal migration | Complete November 6, 2026 (D141); documented limitations retained |
| Post-release | Immigration from outside the system, tracking generation (formerly 10) | Deferred until after first public release (D136) |
| Release preparation | First version bump and public release after milestone 9 | User-selected next focus (D136); version number and release criteria unselected |
| Deferred, likely post-release | Rich console output, diagnostic bullets, condition classes, and polished summaries | Deferred output-presentation project (D151); scope, timing, and dependencies unselected |
| Unassigned | Model diagnostic checks and post-simulation helpers (formerly 11) | Proposed; scope and timing to select during release planning |
| Post-release | Two-sex marriage model (formerly 12) | Deferred until after first public release (D136); model and API unselected |

Milestone 6's scope is individual group transitions. Household moves are
deferred beyond milestone 6 (D114, superseding D109). D121 removes the separate
regions milestone and places multiple group columns at milestone 9 (D131);
regions/internal migration are applications of multiple group columns, and
household moves are implemented for group transitions in milestone 9 (D140). D106–D108 select annual hazards, directed competing
age- or duration-specific schedules with optional strata,
and no transition-count target adjustment. The implemented contract is in [group transitions](group-transitions.md).
D110 requires `group_origin`, `group_dest`, `rate`, and one common age or
duration band type across the transition table, with optional `sex` and
`union_type` columns; absent stratification columns apply universally.
Grouped simulation output must always include a `transitions` history tibble,
even with no configured transition process or no observed moves (D111).

## GitHub development workflow (D137)

Use GitHub milestones to group related issues and PRs, and Projects to track
work status and priority. Split a milestone into bounded issues and reviewable
PRs as appropriate; no one-PR-per-milestone requirement applies. Link issues
and PRs to the relevant repository specifications and validation evidence.
Component merges do not replace the user's milestone acceptance. Existing
GitHub milestones 10 and 12 represent post-release work under D136; milestone
11 (cleanup, documentation, helpers and polish) provides a starting point for
release planning without selecting its full scope as release requirements.

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
as `none`, `separated`, or `widowed` (D080), with all three strata required
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
count unrealized births (D128). See [parity fertility](parity-fertility.md) and
[cell accounting](fertility-cell-accounting.md) for the contract and limitations.

D102 implements persistent positive lognormal selection weights, now internal
only (D176): users cannot enable them, pending a general quantitative-
characteristics feature. Dispersion and maternal
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
now milestone 9 (D121, D131). No schema or public API is selected.
How multiple group types interact in rate stratification, matching, newborn assignment, transitions and
ancestry remains to be specified. The roadmap does not introduce an arbitrary
attribute-mutation callback or a continuous-characteristic model.

Individual transition schedules, age/duration dependence, eligibility, timing,
event/tracker updates and retained ancestry are implemented under the
[group-transition contract](group-transitions.md). Household moves remain deferred
beyond this milestone (D114), with their unresolved questions retained under
[milestone 9](#milestone-9-multiple-group-columns) (D121, D131).
The user accepted milestone 6 on October 30, 2026 (D116), based on the [report](../validation/milestone6.qmd); its final render
was active at acceptance and must be preserved.
The former standalone regions/internal-migration milestone is removed; those
applications may be handled through multiple group columns in milestone 9. Household and geographic rules remain unselected.

Fixed newborn group assignment was selected November 7, 2026 (D142):
`group_inheritance("fixed", fixed_value = ...)` assigns every newborn to a
specified category, such as the minimum education category. See
[group inheritance](group-inheritance.md).

## Milestone 7: time-varying inputs within a run

Accepted complete November 2, 2026 (D127), based on the user-authored
[validation report](../validation/milestone7.qmd) and explicit acceptance.
Calendar-varying rates are implemented under the user's D118 authorization. Flat rate tables use optional half-open
`time_start`/`time_end` bands, including fractional years. The engine updates
pending risks, fertility cells and calibration while preserving histories,
spacing and duration origins. See [the rate contract](time-varying-rates.md).

The user selected annual rates first. D122 also implements flat calendar-
banded `cohab_prob`, using the same interval and coverage rules as rates.
D123 implements complete model schedules for `matching`, `father_matching`,
`cohab_matching` and fertility's `inheritance`: direct objects remain supported,
or supply a tibble with `time_start`, `time_end` and list-column `model`.

Select matching at partner/father choice and inheritance at birth, without
rewriting existing unions or identities. Model/probability-only boundaries retain
pending event dates. All supplied models and their possible identities validate
upfront; future inheritance outputs share one simulation-wide group registry.
See [calendar-varying models](time-varying-models.md) for the implemented API,
coverage rules, callbacks, fallback roles and validation. The schedule remains
outside each individual model, as selected in D118. Other non-rate arguments,
including multipliers, spacing, `ffab` and persistent fertility heterogeneity,
remain static.

## Milestone 8: cell scheduling and weight models (complete)

The user selected this engine change on October 31, 2026 (D121) and accepted it
complete on November 5, 2026 (D130). Every rate process is scheduled by cells
with one clock per process (D128), and weight models choose who in a cell
experiences an event (D129). See [the specification](binned-group-hazards.md)
and [the validation report](../validation/milestone8.qmd). Deferred items
(eligibility rules beyond spacing, unnormalized weights, user-defined
characteristics, weights within bands) are listed in
[its deferred issues](binned-group-hazards.md#deferred-issues).

## Milestone 9: multiple group columns

The user assigned the possibility of more than one group column to tentative
milestone 10 on October 31, 2026 (D121) and moved it to milestone 9, the next
milestone, on November 5, 2026 (D131). Different columns could represent
identity, education or location simultaneously. The design was decided on
November 5, 2026 (D134) and is recorded in
[multiple group dimensions](multiple-groups.md): declared dimensions with
their own names, rate strata, inheritance, ancestry, transitions, analysis
splits, weight-model variables and per-dimension matching settings. Steps 1
and 2 are implemented; matching may combine effects from several dimensions,
calibrated over their combinations. Cross-dimension matching interactions
are not offered.

The former separate milestone for internal migration and regions is removed.
Location may be represented as one group dimension, with internal migration as
transitions in that dimension and geographic matching restrictions where
explicitly specified. This replaces the former roadmap's dedicated-region-column
direction; it does not implement or settle any geographic or household behavior.
The following D114 questions remain relevant to a future regional application.

### Household moves: deferred design questions (D114)

*Resolved for group transitions by D140 (November 6, 2026): see
[household moves](group-transitions.md#household-moves-d140). The questions
below are retained as the record of what that specification answers; any
household rules beyond D140 (grandparents, stepparent custody, coresidence)
remain open.*

D114 deferred household moves from milestone 6 to the former regions milestone.
D121 removes that standalone milestone while retaining multiple group columns
as milestone 9. These concerns are retained for a possible regional
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

## Post-release: immigration from outside the system

Formerly milestone 10 under D131; deferred until after the first public release
by D136. Add immigrants within a running
simulation, rather than requiring a restart.
Under D177, the user selected absolute arrival intensity in migrants per year,
with stochastic arrivals and no controlled/exact-count option. A constant
intensity gives exponential interarrival times and Poisson period counts;
resident population changes do not themselves change that intensity.
Under D178, support multiple streams, potentially representing countries or
world regions of origin, with separate arrival intensities. Stream-specific
migrant characteristics and how origin is retained or mapped to groups still
need specification; no clock organization is selected yet.
Under D179, draw migrant templates with replacement from an existing
`hammel_population`, preserving age, sex and group identifiers and assigning
new person IDs. The source is not depleted. Optional household sampling is
under discussion; selection, intensity units and relationship/history
handling are not yet selected.
D180 requests calendar schedules of population templates using time bounds
and a list column. Selective migration through a predictive selection model
is under discussion as an alternative to tuning the source population's
distributions; no weighting or prediction contract is selected yet.
Calendar-varying intensity inputs and their boundary handling still need a
contract. The user also envisages an age/sex distribution and tracking
individual immigrant generational status. Migrant composition,
attribute/kin/union initialization, ancestry and generation definitions,
interactions with group assignment, entry observations and the public API
remain open. This is a modeling decision, not implementation authorization.
Under D181, retain personal/family migration history separately from reusable,
named generation/origin classifications. Custom matching, inheritance and
weight models should be able to consume concise derived variables without
traversing pedigrees. Definitions, schemas, update timing and history retention
through truncation remain to be specified. Stream origin, group membership
and ancestry are distinct. See [external immigration](immigration.md) for the
consolidated accepted direction, examples and unresolved contracts.
Under D182, support rate stratification by personal immigration status and
arrival stream. Generational classifications remain available to custom
models; normalized weights do not replace generation-specific aggregate rate
targets. Group-based versus built-in stratifier representation and rate-table
contracts remain open.
No emigration process is implied by the immigration request.

## Proposed summaries, diagnostics, and helpers

Moved from former milestone 8 by D121, then from milestone 9 to 11 by D131;
D136 makes release preparation the next focus after milestone 9. The scope
and timing of these helpers remain to be selected; they are not automatically
requirements for the first public release.

On November 9 the user deferred rich diagnostics and summaries to a future
project for improving output presentation, likely after the first release
milestone (D151). Consider tidyverse-style diagnostic bullets, condition
classes, and polished model/result print and summary displays together when
that project resumes. The interface, contents, migration approach, and use of
`cli` or `rlang` remain unselected; no dependency or bulk diagnostic rewrite
is authorized. These improvements are not first-release requirements. This
deferral does not select a schedule or scope for other model checks or
post-simulation analysis helpers.

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

### Unweighted or default intervals in model schedules

Recorded November 4, 2026, during milestone 8 (D129), for milestone 11 (D131). Calendar model schedules (the
`time_start`/`time_end`/`model` tibbles) currently require a complete model in
every interval. For weight schedules, a period without weighting is written as
a formula with zero coefficients. The user is open to allowing `NULL` entries
to mean "no model" for an interval (unweighted, or the default model), but only
if the change applies in parallel to all three schedule types: matching,
inheritance and weight models. Not selected; the meaning of `NULL` for
matching and inheritance (for example, falling back to the general or default
model) needs specification first.

### Post-simulation helpers

Expand the existing extraction/summary functions into a coherent user toolkit
for useful quantities from simulation histories. Existing functions cover Lexis
counts, fertility, life tables, population sizes and recorded kin. Identify
additional quantities and exposure/observation definitions from actual workflows.
Status-specific exposure, fertility and mortality summaries use
`calculate_fertility_rates(by = "union_type")`; formation counts and rates by
previous union use `calculate_union_rates()` (D080).
See [demographic analysis](demographic-analysis.md) for their observation rules.

### Living-population extraction after burn-in

On October 29, 2026, the user requested a future helper to truncate a simulated
population to living individuals, so a burn-in run can provide a new starting
population without retaining the full initialization history. The desired result
preserves correct reproductive histories and parent–child relationships, while
allowing the deeper ancestral kinship history to be discarded. This is a
milestone 11 proposal; no function name or implementation is selected.

Revisited November 12, 2026 (D180) for both burn-in starters and migrant
templates. The user selected retention of living people and their descendants,
including deceased descendants. Retain the descendant closure so recorded
births to living women, including births of children who died, remain available
for parity and last-birth reconstruction. Deceased ancestors are not retained
solely because of their relationship to living people. A related proposal would
model selection from migrant templates instead of modifying their demographic
distributions through tuning helpers.

Retention of deceased descendants addresses the loss of recorded births that
would otherwise conflict with history-derived parity initialization (D098);
it does not impute unknown births. Before implementation, specify handling of removed
parent IDs, retained parent–child links, ongoing unions, previous-union histories,
ancestry, observation endpoints and restart validation. The helper must state
which kinship information has been discarded; it cannot imply complete pedigree
or historical exposure after truncation. The representation of minimal history
needed for correct restarting remains to be designed.

## Post-release: two-sex marriage model

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

## First public release (D136)

After milestone 9 is completed and accepted, focus on the first version bump
and public release. The version number, distribution channel and concrete
release criteria remain to be selected. Release planning should assess the
existing API, documentation/examples, package metadata and licensing, checks
and retained validation evidence, and clearly describe supported behavior and
limitations. This is a planning direction, not authorization to publish or a
claim of comprehensive SOCSIM equivalence. Immigration and two-sex marriage
models are not prerequisites for this release.

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
schedules are implemented under milestone 7, accepted complete (D127).
