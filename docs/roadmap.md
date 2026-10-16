# Directions beyond milestone 2

Recorded October 15, 2026 from the user's planning notes. Milestone 2 remains
in validation. The user selected **groups as the theme for milestone 3**;
other items below are prospective work without a milestone assignment. These
notes preserve intentions and open questions, not authorization to implement
unspecified demographic rules or a selected public API.

## Analysis helpers

Expand the existing extraction/summary functions into a coherent user toolkit
for useful quantities from simulation histories. Existing functions cover Lexis
counts, fertility, life tables, population sizes and recorded kin. Identify
additional quantities and exposure/observation definitions from actual workflows.
Status-specific exposure and fertility summaries are relevant to milestone 2
validation; do not assume they are implemented by this roadmap.

## Preparing input schedules

Provide more helpers for building required input rates, including retrieving or
using model life tables and splitting overall fertility into partnered and
unpartnered schedules. No source, model-life-table family, fitting method or
split rule has been selected. Overall fertility alone does not identify two
status-specific schedules: decomposition needs additional information or
assumptions, such as age-specific status proportions and relative fertility.
Keep empirical rates, inferred rates and user assumptions inspectable with
provenance. Avoid automatic downloads as part of simulation execution.

## Milestone 3: groups

The user selected the following directions for group support:

- Include intergroup unions in formula-based matching, alongside age difference.
- Allow proportional rate differences by group relative to a baseline schedule,
  avoiding a complete duplicate rate table for each group. An illustrative
  fertility multiplier is 1.1 times the baseline, not a demographic default.
- Support flexible child group assignment when parents have different groups,
  with more user control than the user has found in SOCSIM. The rules have not
  been selected; do not presume maternal, paternal or probabilistic assignment.
- Consider tracking each person's proportional ancestry by group separately
  from assigned group identity.

Before implementation specify which processes can use multipliers, how they
combine with status schedules and birth-spacing compensation, and the formula
predictors/contrasts for each parent's group. For ancestry proportions clarify
the reference ancestors, missing-parent treatment, initialization and whether
proportions represent expected genealogical contributions or another quantity.
A discrete group identity need not equal a majority-ancestry classification.
These are open specification questions, not accepted restrictions on the design.

## Fertility heterogeneity

Consider person-level fertility differences to represent persistent high- and
low-fertility individuals and increase outcome variance beyond the current
model. The user reports SOCSIM has a mechanism for this; its precise behavior
has not been audited for this task. Distribution, persistence, association with
other attributes, inheritance and calibration to aggregate rates are unselected.
No heterogeneity parameter or demographic default is adopted here.

## Characteristics

Consider characteristics that can affect rates and matching without functioning
as definitive group identities. Education is the user's example. Users would
need to specify intergenerational inheritance/initialization and life-course
models. Define permitted changes and event invalidation, not arbitrary state
mutation. No distinction between categorical/continuous attributes or transition
model has been fixed by these notes.

## Immigration

Add immigrants within a running simulation, rather than requiring a restart.
The user envisages specifying numbers and an age/sex distribution and tracking
individual immigrant generational status. Timing (batch dates versus arrivals),
meaning of supplied numbers, attribute/kin/union initialization, ancestry and
generation definitions, and interactions with group assignment remain open.
No emigration process is implied by the immigration request.

## Previously deferred items

- Previous-father continuity at subsequent births or union formation (D042).
- Declining poor nonempty candidate pools and waiting (D044).
- Transition-aware fertility compensation for status/rate changes (D045).
- Parity-specific fertility, non-death union dissolution, and multiple concurrent
  unions; no implementation milestone has yet been selected.

R formula-fallback performance benchmarking belongs to milestone 2's final
validation report. It is not deferred beyond that milestone. Broader existing
future notes about calendar-varying schedules, initialization and exact
checkpoints remain separate from the user's new priorities above.
