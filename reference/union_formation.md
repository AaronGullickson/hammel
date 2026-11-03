# Specify a union-formation process

Construct or attach a union-formation specification without executing
it. In on-demand mode, at a birth use the woman's current partner, or
choose a partner from all living, unpartnered men using
conditional-logit probabilities from the matching formula to form a
union. Men below `min_male_age` are excluded. If no eligible man exists,
cancel the birth by default; `cancel_birth_if_no_match = FALSE`
preserves it with unknown father. Known kin within `kinship_generations`
on both sides of a common ancestor are excluded. Each person may have
only one ongoing union. Partnered men are not available; surviving
husbands return to the pool at widowhood. Death ends affected unions. In
rate mode unpartnered women face competing search events and use their
unpartnered fertility schedule and may give birth without initiating a
search. Failed searches draw another waiting time; formation and
widowhood switch fertility schedules while preserving spacing. An
overall fertility schedule applies to both statuses, with a simulation
warning. On-demand mode requires an overall fertility schedule.
Unpartnered rate-mode births use the same eligibility and exact draw to
assign a father, using `father_matching` when supplied or otherwise
`matching`, without forming a union or removing him from availability. A
search with no eligible man cancels the birth by default. No
previous-father preference is applied.

## Usage

``` r
union_formation(
  initiation,
  matching = create_age_gap_matching_model(0, 5),
  rates = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  father_matching = NULL,
  cancel_birth_if_no_match = TRUE,
  diagnostics = FALSE,
  group_prop_hazard = NULL,
  cohab_prob = NULL,
  cohab_matching = NULL,
  union_convert_rate = NULL,
  cohabitation_counts_as_previous_union = FALSE
)

set_union_formation(
  model,
  initiation,
  matching = create_age_gap_matching_model(0, 5),
  rates = NULL,
  min_male_age = 15,
  kinship_generations = 2,
  father_matching = NULL,
  cancel_birth_if_no_match = TRUE,
  diagnostics = FALSE,
  group_prop_hazard = NULL,
  cohab_prob = NULL,
  cohab_matching = NULL,
  union_convert_rate = NULL,
  cohabitation_counts_as_previous_union = FALSE
)
```

## Arguments

- initiation:

  Either `"on_demand"` for birth-triggered searches or `"rate"` for
  scheduled searches. No default is selected.

- matching:

  A specification created by
  [`create_matching_model()`](https://aarongullickson.github.io/hammel/reference/create_matching_model.md),
  or a tibble with numeric `time_start`, `time_end`, and list-column
  `model` holding complete matching objects. Model intervals are
  half-open absolute years, cannot overlap, and must cover the requested
  run. The active model is used at partner selection, without changing
  pending search dates or existing unions. Base data frames with list
  columns are also accepted. Defaults to
  `create_age_gap_matching_model(0, 5)`: preferred age gap zero years,
  with selection weight halving five years on either side.

- rates:

  For `"rate"` initiation, an age-specific female search-hazard table
  with `age_start`, `age_end`, and `rate`, in years and per person-year.
  These specify search attempts, not guaranteed realized union
  formation. An optional `group` column selects the woman's group.
  Without it, all groups share the schedule. Group does not restrict
  partner matching. With no population groups, simulation warns and
  selects the first alphabetical group. See
  [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
  for group-label and rate-coverage rules. Nonoverlapping bands may have
  gaps or an open end; gaps and ages outside support have zero hazard.
  Optional `previous_union` supplies either `"none"`, `"separated"`,
  `"widowed"`, or pooled `"none"`, `"dissolved"`. Dissolved pools
  separation and widowhood. Each group and active calendar slice must
  supply one complete set; mixing the sets is not allowed. Without it,
  the schedule applies to all previous union types. Categories describe
  the most recent recorded episode; complete histories are required.
  Optionally include both `dur_start` and `dur_end` for elapsed years
  since the most recent union ended. Each row supplies either age or
  duration bands with both unused boundaries NA. `none` must use age;
  each other category uses one clock consistently across all its bands
  and groups. Age and duration cannot apply together. All-age inputs may
  omit duration columns. Must be NULL for `"on_demand"`.

- min_male_age:

  Minimum male age in continuous years at a new partner search, default
  `15`. Finite nonnegative scalar; zero disables this age exclusion.
  Existing unions are retained. No minimum is added for actual female
  event eligibility. Weighted choice counts unpartnered men and women at
  or above this same cutoff when calibrating its weights.

- kinship_generations:

  Generation bound on both sides of a known common ancestor, default
  `2`. Excludes parents, grandparents, siblings, aunts/uncles,
  nieces/nephews and first cousins, including half relationships, when
  known. Zero imposes no cross-person kin exclusion. Missing ancestry is
  not proof of unrelatedness; only demonstrated relationships are
  screened. Existing unions are retained. See
  [`extract_kin()`](https://aarongullickson.github.io/hammel/reference/extract_kin.md)
  for the shared distance definitions.

- father_matching:

  Optional matching specification for fathers of unpartnered births in
  rate mode, created by
  [`create_matching_model()`](https://aarongullickson.github.io/hammel/reference/create_matching_model.md),
  or a `time_start`/`time_end`/`model` schedule as for `matching`.
  Defaults to `NULL`, using the general `matching` model active at that
  birth. Union searches always use `matching`. Not used by on-demand
  initiation. Process-level age and kinship screening settings are
  shared; NULL inherits the general model's settings.

- cancel_birth_if_no_match:

  Logical, default `TRUE`. Cancel a birth when its on-demand partner
  search or rate-mode nonmarital father selection fails. Set `FALSE` to
  allow births with unknown fathers. Historical missing father IDs are
  retained. A cancellation starts no new spacing interval and draws
  fresh waiting times from the failed attempt. Partner scarcity can
  reduce realized fertility below inputs.

- diagnostics:

  Logical, default `FALSE`. When the general matching model uses
  weighted choice in a grouped population, record `.log_weight_<group>`
  columns on newly formed unions: the log group weights applied at the
  search that formed each union (0 for the reference group and when
  calibration cannot affect choices). Use the column for each
  alternative's group as a conditional-logit offset in wife-side choice
  sets.

- group_prop_hazard:

  Optional named numeric vector of finite nonnegative search-rate
  multipliers by the woman's group. Missing groups use 1.0; zero
  suppresses searches. Applies to shared or group-specific input rates.
  Without population groups, simulation warns and ignores multipliers.
  Requires `initiation = "rate"`; on-demand formation has no search
  hazard. Matching utilities and input rate tables are unchanged.

- cohab_prob:

  Optional age-specific probability table with `age_start`, `age_end`,
  `prob`, optionally `group` and `previous_union`. Enables typed unions
  only in rate mode. Age gaps mean zero. Optional `time_start` and
  `time_end` supply half-open calendar intervals; each active slice must
  satisfy the same age/group/previous-category rules and cover the run.
  Type uses the active probability at the search instant, before
  matching. Existing unions keep their types. NULL retains untyped
  partnerships and their RNG path.

- cohab_matching:

  Optional cohabitation matching model or whole-model schedule; NULL
  shares the general `matching` model active at that search. Marriage,
  cohabitation and unpartnered father searches calibrate weights using
  their own intensities and share the available male pool.

- union_convert_rate:

  Optional age- or duration-specific conversion table, optionally by
  woman's group. Duration uses the current cohabitation start. Uses
  `rate` and either age or duration boundaries, not both clocks, with no
  previous-union axis. NULL gives no conversions and warns at
  simulation. Conversion keeps partners, closes cohab with reason
  `conversion`, and starts a linked marriage with a new
  dissolution-duration origin.

- cohabitation_counts_as_previous_union:

  Logical, default FALSE. FALSE restores previous legal-marriage
  category after a cohabitation ends; TRUE classifies cohab
  breakup/death as separated/widowed. Search duration still starts at
  the latest partnership ending. Only applies to typed unions.

- model:

  A `hammel_model` object.

## Value

`union_formation()` returns a `hammel_union_formation` process;
`set_union_formation()` returns an updated model. Rate tables are
tibbles.

## Details

Partners are drawn exactly from the full conditional-logit distribution
over all eligible men, without a candidate sample. Living, unpartnered
men are grouped by group and quarter-year birth bin; a bin is chosen in
proportion to its size times the largest weight any member could have, a
member is accepted with the ratio of his weight to that maximum, and
only he is screened for age and kinship. For the built-in
linear/quadratic age model this is exact. Custom R formulas are bounded
on a 1/32-year grid; if a drawn man's utility exceeds the bound,
simulation warns that draws were approximate. Existing unions are
retained even if a new-search restriction would exclude the partner.

Matching always uses calibrated male-group weights: a woman chooses
among screened candidates with each man's selection weight multiplied by
a weight for his group. Weights are recalculated at every search, in
continuous time, so that expected male union flows by group match a
mirrored target: each unpartnered man has the per-capita intensity of an
unpartnered woman of his group, adjusted by the matching model's
`delta`. Women's intensity is the search hazard for rate-mode searches
and the unpartnered fertility hazard (zero during birth spacing) for
on-demand searches and nonmarital father selection, which has its own
weights and uses the father model's preferences and delta when one is
supplied. Weights depend only on the man's group, so they do not change
the association between partners' characteristics in the static
group-only illustration. With one group, all group weights are one. When
no exact weights exist, approximate weights are used with a warning.

## Calendar rates

Search `rates` and `union_convert_rate` optionally accept `time_start`
and `time_end` for half-open absolute calendar-year bands, including
fractional years. Both columns are required together, with finite starts
and greater ends (possibly Inf). Every active slice must satisfy its
existing rate rules and cover the requested run. Clock choices stay
consistent across periods; changes precede same-time events and retain
duration origins. `cohab_prob` accepts the same calendar columns, with
probabilities selected at the search instant. Probability-only changes
refresh matching calibration without resampling pending event clocks.
All three matching arguments also accept whole-model schedule tibbles
with `time_start`, `time_end`, and list-column `model`. Model changes
update role-specific preferences and calibration before same-time
choices without resampling event clocks. All supplied models, including
future intervals, are validated; gaps in a run are errors. NULL
father/cohabitation matching dynamically follows general matching. Other
non-rate arguments remain static.
