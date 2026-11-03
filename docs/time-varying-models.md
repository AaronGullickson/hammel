# Calendar-varying matching and inheritance models

Implemented milestone 7 contract, October 31, 2026 (D123). The user authorized all four
model arguments—`matching`, `cohab_matching`, `father_matching`, and fertility's
`inheritance`—to accept either a direct model or a schedule tibble with
`time_start`, `time_end`, and a list column named `model`. This extends D118's
approved event-time model selection. It does not introduce the milestone 8
binned-group hazard rewrite. Milestone 7 was subsequently accepted complete
on November 2, 2026 (D127); see the [validation report](../validation/milestone7.qmd).

## Input contract

A direct model keeps its existing static interpretation. Matching entries are
complete objects from `create_matching_model()` or its age-gap helper;
inheritance entries are complete `create_group_inheritance()` objects. Existing
static inheritance method-name aliases and top-level NULL defaults remain
supported. A schedule contains at least one row; base data frames with list
columns are accepted and normalized to tibbles. Extra columns and source
attributes are retained, without assigning them execution semantics.

Both interval columns are numeric absolute years, with finite starts (negative
and fractional years allowed) and greater ends (possibly Inf). Intervals are
half-open `[time_start, time_end)`, must not overlap, and are not interpolated.
Every row must contain a valid model of the correct kind; NULL entries, nested
schedules and method-name strings inside the list column are errors. All rows
are validated, including those outside the requested run. The complete supplied
schedule is retained; private lookup may sort intervals without reordering it.

Execution requires coverage throughout `[start, end)`. A zero-year run requires
one active model at start. Gaps are errors when the run reaches them, even if
no search or birth happens there; no model is carried forward implicitly.
Each argument may have different interval boundaries.

Top-level `father_matching = NULL` and `cohab_matching = NULL` continue to share
whatever general `matching` model is active at the event. An explicitly supplied
static model stays static even when general matching varies. An explicit
schedule supplies every interval itself, without NULL rows for fallback.
Cohabitation matching still requires cohabitation enabled; father matching is
used only for unpartnered births in rate mode. Supplying an unused role does
not bypass validation of its schedule.

## Identity registry and validation

Use one simulation-wide identity registry covering starting groups, ancestry,
transition identities and declared possible outputs from every inheritance
model in the schedule, including future intervals. Every model uses the same
native group indices. Rate coverage, matching factor levels, delta labels and
other existing group checks apply against this full registry before execution.
A future custom inheritance identity must therefore be declared and supported
by the configured grouped processes from the outset. Registry preparation does
not assign that identity to existing people or change ancestry.

## Event-time semantics

- General matching is selected when a direct-marriage search or an on-demand
  birth-triggered partner search executes.
- Cohabitation matching is selected after intended type is drawn at the search.
  Its top-level NULL fallback uses the current general matching model.
- Father matching is selected when choosing a father for an unpartnered birth
  in rate mode. Partnered births continue to use the current partner.
- Inheritance is selected for each successful birth, after sex and parental
  ancestry are known, before scheduling the newborn's risks. Existing identity,
  parent, missing-father and ancestry rules remain in force.

At an exact boundary, the newly active model applies. Model-only changes
consume no demographic random draws and preserve pending event dates, duration
origins, spacing, parity, histories and persistent weights. They do not form or
end unions, reselect existing partners/fathers, reassign existing identities or
change ancestry. Shared rate/probability/model boundaries apply all active inputs
before events; rate boundaries retain D118 rescheduling.

Matching changes include formulas, coefficients, group preferences/restrictions
and delta targets in the complete object. Refresh the active partner-draw and
role-specific weighted-calibration model before the next selection. Switching
between within-group-only and weighted matching must work in either direction.
Screening parameters such as minimum male age and kinship bound remain static.

Inheritance changes include methods, hierarchy order/passing and custom callbacks
in the complete object. The current `hetfert` mechanism is separate and remains
static; identity inheritance schedules do not alter fertility-weight inheritance.
Custom inheritance and matching retain their existing safe callback and RNG
contracts. No inheritance callback is invoked merely because a boundary is
crossed. Deterministic matching formulas may be evaluated to prepare a model's
cached draw/calibration structures, as under the existing matching contract.

## Worked examples

```r
matching_schedule <- tibble::tibble(
  time_start = c(1811, 1900.5),
  time_end = c(1900.5, Inf),
  model = list(matching_early, matching_late)
)
inheritance_schedule <- tibble::tibble(
  time_start = c(1811, 1920),
  time_end = c(1920, Inf),
  model = list(create_group_inheritance("mother"),
               create_group_inheritance("father"))
)
```

These model names/years illustrate the API, not demographic defaults.

- A search scheduled in 1899 for 1901 uses the matching model active in 1901;
  changing preferences in 1900.5 does not move the search date.
- With mother-A/father-B parents, maternal inheritance before 1920 assigns A;
  paternal inheritance at and after 1920 assigns B. Older children retain A and
  both children's ancestry is calculated by the existing parental-average rule.
- If father matching is NULL while general matching changes, subsequent
  unpartnered father draws use the new general preferences. If father matching
  is an explicit static model, those draws retain its preferences.
- Changing hard exogamy restrictions does not invalidate existing cross-group
  unions. A failed future search still follows the existing failure rules.

## Implementation and validation requirements

Model construction stays separate from simulation. Keep full inspectable models
in the public list column, including formula environments and custom functions;
prepare private native structures with the common registry. Cache reusable draw
and calibration structures within the simulation call; no hidden global state.
Static direct-model paths should preserve existing behavior and avoid scheduled-
model overhead. New algorithms must not claim universal fixed-seed equivalence;
use exact comparisons only when random-number consumption is preserved.

Validate malformed intervals/list entries, all scheduled model kinds and future
models, edited objects, coverage/overlap, static defaults, dynamic NULL fallback,
all four event roles, fractional/exact boundaries, current model restrictions,
role calibration/delta/formula changes, expanded future identities, custom callbacks,
unchanged ancestry and older histories, simultaneous/asynchronous input changes,
restarts, zero/empty runs, RNG isolation and independent calls. Use controlled
candidate pools and analytical choice probabilities where possible. Run the
ordinary offline `devtools::test()` suite; no database requests or report renders
are required. Scientific acceptance remains the user's assessment.
