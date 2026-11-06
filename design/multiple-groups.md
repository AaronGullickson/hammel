# Milestone 9: multiple group dimensions

Started November 5, 2026. **Status: design decided (D134); stages 1–5
implemented: declared dimensions of any name, several at once, in every rate
table, with per-dimension inheritance (including calendar schedules),
ancestry, transitions, analysis splits and weight-model variables, and
matching group effects in any one dimension. Remaining: stage 6 (tests,
articles and documents) and step 2 (joint matching calibration across
dimensions).**
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

### Stage 5 as implemented (November 6, 2026)

- **Matching.** Predictors are `female_<dim>`, `male_<dim>` and `same_<dim>`;
  `group_exogamy`, `prevent_exogamy` and `delta` are lists keyed by dimension,
  or a bare value for the only declared dimension (a convenience, like
  `group = NULL` for transitions). `group_levels` is removed: levels are the
  dimension's categories, the first being the reference. A model's group
  effects may use one dimension, and every matching model of a run must use
  the same one; using two is an error (joint effects are step 2). That
  dimension becomes the engine's matching dimension (pools, pair effects,
  restrictions, `delta`, calibration). Without group effects and with several
  dimensions there is none. Coefficient names of group predictors are checked
  at simulation, when the categories are known. At simulation each public
  model is translated to the internal one-dimension form (`female_group`,
  explicit levels), so the engine path is unchanged.
- **Weight models.** Variables are named by dimension: `race`,
  `partner_race`, `ancestry_race_A`, `ancestry_race_unknown` (person and
  transition processes), `destination` (transitions, in their own
  dimension), and `race_w`, `race_h`, `same_race`, `ancestry_race_A_w`,
  `ancestry_race_A_h` (union processes), for any declared dimension and in
  any combination. `group_levels` is removed; levels are the categories.
  Names that could be group variables are resolved and their coefficients
  checked at simulation (also for zero-year runs); other names still fail
  when the model is set. In the engine a group variable's factor carries its
  dimension's index.
- **Custom weight functions.** Members have a column per dimension
  (`race`, `partner_race`; `race_w`, `race_h`, `same_race`); without declared
  groups there are no group columns. `bin` describes the cell's strata by
  dimension (`race`, or `race_w`/`race_h` for union processes; NA when the
  rates do not stratify by it; `origin`/`destination` and other dimensions for
  transitions). This fixes the stale stage 2 `group_pair` label for
  dissolution.
- **Bug fixed.** A calendar change in fertility rates sized the fertility
  strata by the matching dimension's categories instead of the fertility
  table's own strata (equal with one dimension).
- With one dimension named `group`, results (persons, unions, ancestry,
  transitions, custom-function calls) are bit-identical to stage 4 for a model
  using matching pair effects, restrictions, `delta`, categorical predictors,
  calibration, and formula and custom weights on every process, when both
  builds disable fused multiply-adds. Undeclared categories follow first
  appearance, so a reference level that came from `group_levels` order may
  now need declared `categories`.

### Stage 4 as implemented (November 2026)

- `by` in every calculator accepts any declared dimensions by name, several at
  once, and union rates accept any combination of `<dim>_w`/`<dim>_h`.
  Results have one column per dimension; the grid is the full cross of each
  dimension's observed categories (first-appearance order, as before).
- Implementation: the listed dimensions form one combined internal category
  with a merged history of their moves (same-instant moves in declaration
  order); results are split back, and union rates sum over partner-dimension
  combinations not requested. With one dimension the code path is unchanged
  and results are identical to stage 3.
- The observation engine now reads a history with several moves at one
  instant (state after a time follows all of them), and takes an optional map
  from combined category to one dimension: `calculate_transition_rates(...,
  group = "race", by = "region")` counts only race moves, with race duration
  unaffected by region moves.
- `by` errors list the declared dimension names rather than the internal
  `group`.
- Summaries (`calculate_tfr()`, `calculate_life_table()`,
  `calculate_patfr()`) now keep every non-band, non-measure column as a
  stratum. Previously only `group` was kept, so a rate table split by a
  dimension with any other name (possible since stage 2) was pooled or failed.

### Stage 3 as implemented (November 2026)

- `set_group_transition(model, rates, group = NULL, weights = NULL)` keeps one
  process per dimension in `model$processes$group_transition`, a list keyed by
  dimension. `group = NULL` takes the only declared dimension; with several it
  is an error naming them. A process set before the declaration exists is
  stored under `.default` and resolved at simulation. Setting a dimension
  again replaces its process.
- Tables use `origin`/`destination`; `group_origin`/`group_dest` are errors
  naming the new columns. Other declared dimensions may stratify a table;
  uncovered strata have zero risk, as for sex and union type. The process's own
  dimension cannot be a column.
- The engine runs one transition process per dimension in declaration order,
  each with strata over its origin and the other dimensions its table uses,
  its own weight model, and entry times in that dimension only. Same-instant
  transitions follow declaration order.
- Histories have `person_id`, `time`, `group`, `origin`, `destination`, all
  labels character (one table holds every dimension's labels; previously labels
  kept the persons column's type). A supplied history may omit `group` when one
  dimension is declared. Validation chains each person's history per dimension.
- The transition weight variable `group_dest` is now `destination`; custom
  weight functions see `origin`/`destination` in `bin` and a `destination`
  member column. Group variables in weight models (including `destination`)
  still require one dimension (stage 5).
- Inheritance calendar schedules work for every dimension.
- `calculate_transition_rates(..., group = NULL)` counts one dimension's
  transitions, returning `origin`/`destination`. `get_detailed_unions()`
  replays each dimension's transitions to union start. Analysis code works on a
  one-dimension view of the population.
- With one dimension, results are bit-identical to stage 2 when both builds
  disable fused multiply-adds.
- Known issue for stage 5: custom weight functions' `bin` for mortality,
  fertility, search, conversion and dissolution still describes strata by the
  single dimension (`group`, and a stale `group_pair` for dissolution), which is
  wrong for directional dissolution strata and for several dimensions.

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
