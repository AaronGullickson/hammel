# Milestone 9: multiple group dimensions

Started November 5, 2026. **Status: design decided (D134); stages 1 and 2
implemented: declared dimensions of any name, several at once, in every rate
table, with per-dimension inheritance and ancestry. Transitions, matching and
weight-model group effects, inheritance schedules and `by` group splits still
use a single dimension (stages 3–5).**
This records the user's decisions as they are made. Milestone 9 was assigned
multiple group columns by D131. See the [roadmap](roadmap.md#milestone-9-multiple-group-columns)
for the motivation and the household-move questions retained from D114.

## Purpose

A population can carry several categorical group dimensions at once, such as
race, region and education, each with its own categories, inheritance,
transitions and effects on rates and matching. Today a single `group` column
serves all purposes, and combining dimensions means encoding them in one
label.

## Decisions so far (user, November 5, 2026)

**Groups are always declared.** There is no implicit single `group` column:
a population's group dimensions are those declared with

```r
set_groups(model,
  race   = group(categories = c("A", "B", "AB"), ancestry = TRUE,
                 inheritance = create_group_inheritance("custom", fun = mixed_rule)),
  region = group(inheritance = create_group_inheritance("mother")))
```

A one-dimension model uses the same mechanism, for example
`set_groups(model, group = group())`. Pre-release, no legacy single-group path
is kept.

- **Names.** Each dimension is named by its persons column (`race`, `region`),
  not by a prefix. Any name is allowed except names Hammel already uses as
  variables or columns (for example `sex`, `union_type`, `previous_union`,
  `parity`, `partnered`, `age_difference`, `person_id`, `origin`,
  `destination`), names starting with `partner_`, `same_`, `female_`, `male_`
  or `ancestry_`, and names ending in `_w` or `_h`.
- **Derived names**, for a dimension `race`:

  | Where | Names |
  | --- | --- |
  | Person weight variables | `race`, `partner_race`, `ancestry_race_A`, `ancestry_race_unknown` |
  | Union weight variables | `race_w`, `race_h`, `same_race`, `ancestry_race_A_w` |
  | Matching formulas | `female_race`, `male_race`, `same_race` |
  | `by` in person rates | `"race"` |
  | `by` in dissolution and conversion | `"race_w"`, `"race_h"` |
  | Transitions | `origin`, `destination` (the transition's own dimension), replacing `group_origin`, `group_dest` and the weight variable `group_dest`; other dimensions by name; `calculate_transition_rates(group = "education")` selects the dimension |

- **Ancestry storage.** The population's `ancestry` component is a named list
  of tibbles, one per tracked dimension, each with `person_id`, one column per
  category and `unknown`.
- **Where the declaration lives.** `set_groups()` stores it with the model's
  population, before or after `set_population()`; it is validated whenever both
  are present and at simulation. Simulation output carries it, so analysis and
  restarts need nothing further. `hammel_population(..., groups = ...)` accepts
  the same declaration for analyzing a population directly, without a model.
- **Undeclared columns** are ordinary custom columns: kept, `NA` for newborns,
  never changed and not used by any process, weight model or `by`. Simulation
  errors if a rate table has a column matching an undeclared person column,
  asking the user to declare it.
- **`group()`** describes one dimension:
  - `categories`: optional complete list in order. The default infers the
    population's values in first-occurrence order plus the dimension's
    transition destinations. When given, every population value and transition
    destination must be in it; extra categories may start empty. Anything that
    creates a member (custom inheritance, transitions) must use a listed
    category. Group-specific rate tables cover every category. This replaces
    the `groups` argument of `create_group_inheritance()`.
  - The first category is the reference level for group terms in matching and
    weight formulas (coefficient names such as `raceB`). This replaces
    `group_levels` in `create_matching_model()` and `create_weight_model()`;
    coefficient names are checked when the model meets the declaration.
  - `ancestry`: whether to track ancestry shares for this dimension, default
    `FALSE`. The `"ancestry"` inheritance method and ancestry variables in
    weight models require `ancestry = TRUE` for their dimension, and are
    errors otherwise. Only tracked dimensions appear in the population's
    ancestry output.
  - `inheritance`: the newborn rule for this dimension, moved here from
    `set_fertility()`, including calendar schedules (D123). Default: random
    parent with maternal fallback.
- **Validation against the declaration.** Everything that names a dimension or
  a category outside `set_groups()` is checked against it, at simulation and
  earlier when the declaration is already present: matching settings
  (`group_exogamy`, `prevent_exogamy`, group formula terms, per dimension),
  rate-table group columns and coverage, transition tables, weight-model
  variables and `by`. Unknown dimensions are errors naming them.
- **Transitions** stay a separate rate process, one call per dimension:
  `set_group_transition(model, rates, group = "region", weights = NULL)`.
  Calling again for the same dimension replaces it; each dimension has its own
  weights. Tables use `origin`/`destination` (categories of that dimension), age or duration bands, optional `sex` and `union_type`, and may
  be stratified by other declared dimensions (education transitions by
  `region`); a person's risks change when another dimension changes. Their own
  dimension cannot be a stratum. Duration is time since entering the current
  category of that dimension only. The population's `transitions` history
  gains a column naming the dimension. Same-instant transitions in different
  dimensions follow the declaration order. Transitions were kept out of
  `group()` because they share the rate-process machinery (calendar bands,
  `stratify_rates()`, weights, cells) and are changed more often than
  declarations.
- **Matching** settings are per dimension. Formula predictors take the
  dimension's name: `female_race`, `male_race`, `same_race`, with the first
  category as reference. `group_exogamy`, `prevent_exogamy` and `delta` become
  named lists keyed by dimension, for example
  `group_exogamy = list(race = c("A:B" = -1), region = c("N:S" = -2))`.
  Effects from different dimensions add; a candidate excluded in any dimension
  is ineligible; per-dimension `delta` values add for each combination of
  categories. Calibration works jointly over the cross-classification of the
  dimensions the matching model uses (dimensions it does not mention stay
  out), with a benchmark of its cost. Cross-dimension interactions (race
  exogamy that differs by region) are not offered; a custom function may
  allow them later.
- **Staging (user direction).** Step 1 gets everything else about multiple
  dimensions working, with group effects in matching limited to one
  dimension, calibrated as now. Step 2 adds joint calibration across
  dimensions.
- **Baseline populations.** `create_baseline_population()` renames
  `group_weights` to `groups`, which takes either a named list of
  per-dimension weights (drawn independently), such as
  `list(race = c(A = .6, B = .4), region = c(N = .5, S = .5))`, or a data frame
  with one column per dimension and `weight`, giving the joint distribution.
  The population carries a default declaration (each dimension with `group()`
  defaults), so it can be analyzed directly; a later `set_groups()` replaces
  that declaration. Ancestry turned on there initializes shares from each
  person's category.
- **`stratify_rates()`** (D132) must accept dimension names, so `by` may name
  any column that is not reserved; the stratum check moves to where the table
  meets the declaration.

## Order of work

Step 1, each stage leaving the test suite passing:

1. Declaration: `set_groups()`, `group()`, storage on the population,
   validation, baseline `groups`, inheritance moved out of `set_fertility()`.
2. Rates and cells: group columns in all rate tables, `stratify_rates()`
   accepting dimension names, coverage per dimension, engine cells crossing
   dimensions.
3. Inheritance per dimension, ancestry per flagged dimension, transitions per
   dimension (histories naming the dimension; `origin`/`destination`).
4. Observation: `by` with dimension names, union-level `race_w`/`race_h`,
   transition rates by dimension.
5. Weight-model variables; matching with group effects in one dimension.
6. Tests, articles and documents.

Step 2: joint matching calibration across dimensions.

### Stage 2 as implemented (November 5, 2026)

User decisions: dissolution and conversion rates use directional
`<dim>_w`/`<dim>_h` columns (the wife's and husband's categories); the
symmetric `group_pair` is gone, and symmetric rates list both directions.
Search rates and cohabitation probabilities use the woman's categories under
the dimension's own name.

- Each person carries a category per dimension in the engine. Each rate
  process stratifies by the dimensions its table has columns for, combined in
  declaration order (wife before husband); cells are only as fine as the rates
  (D128). Every combination of the used dimensions' categories must be
  covered. Rows for categories absent from an inferred registry are unused;
  with declared categories they are errors.
- Setters cannot know the declaration, so any non-schema column counts as a
  possible stratum in their overlap checks; at simulation only declared
  dimensions stratify, and overlapping bands within a real stratum are errors.
  `stratify_rates()` accepts any name except band, rate and calendar columns.
- Newborns get each dimension's category from its own rule, in declaration
  order (one dimension reproduces stage 1 draws exactly). Ancestry is a named
  list of tibbles, one per tracked dimension, with columns named by category
  (so tracked dimensions cannot have categories named `unknown` or
  `person_id`); custom inheritance rules receive that dimension's shares.
- Unrealized births are reported with a column per dimension.
- With one dimension, results are bit-identical to stage 1 when both builds
  disable fused multiply-adds; ordinary optimized builds can differ in the last
  bit because compilers fuse differently around changed code.
- With several dimensions, transitions, matching group effects, weight-model
  group variables, inheritance schedules and `by` group splits stop with a
  message; with one dimension of any name they work, and `by` uses its name.

### Stage 1 as implemented (November 5, 2026)

`group()`, `set_groups()`, `hammel_population(groups =)` and baseline `groups`
(named list of weights or a joint data frame) work as specified, with the
declaration stored on the population (`population$groups`) and on the model
(`model$groups`, which wins). Inheritance and categories moved into `group()`
(`set_fertility(inheritance =)` and `create_group_inheritance(groups =)` are
gone); ancestry is tracked only with `ancestry = TRUE`. An undeclared `group`
column is hidden from the engine during simulation and restored with `NA` for
newborns; rate tables using an undeclared persons column are errors, which
replaces the former alphabetical fallback for ungrouped populations. Until
stage 2, a declaration may contain only one dimension, named `group`, so the
engine and rate tables are unchanged. The broader documentation pass is
stage 6; until then the other design documents still describe an implicit
`group` column.
