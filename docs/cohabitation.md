# Cohabitation and marriage within partnerships

Milestone 4 stretch-goal specification, October 28, 2026 (D085/D086).
Implementation is authorized; the user selected the rules below. This is not
milestone scientific acceptance. Cohabitation is supported only in rate mode.

## Inputs and entry

`union_formation()` and its setter accept `cohab_prob`, `cohab_matching`,
`union_convert_rate`, and `cohabitation_counts_as_previous_union = FALSE`.
`cohab_prob = NULL` retains untyped partnerships, without type draws or new
clocks. When enabled, probability tables use `age_start`, `age_end`, `prob`,
optionally `group` and `previous_union`; probabilities are finite in [0,1].
Bands cannot overlap within strata. Gaps and uncovered ages have probability
zero. Group coverage and all-three previous categories follow search schedules.

The overall search hazard still schedules attempts. At the attempt time, draw
intended type before matching using current age/group/previous category.
Matching uses the general model for marriage and `cohab_matching` for cohabitation
(NULL shares the general preferences). Both draw from the same exclusive pool
of unpartnered men. Failure leaves the woman unpartnered, resamples her ordinary
search risk, and draws type anew at the next attempt. Father selection remains
separate through `father_matching`; partnered births use the current partner.

Calibration maintains marriage lambda*(1-p), cohabitation lambda*p, and
unpartnered-father fertility intensity separately, with role-specific weights
and shared pool utility sums when models are identical. Probability age bands
join the contribution boundary queue, without additional state cross-products.

## Fertility and dissolution

Fertility input stratification is renamed from `union_status` to `union_type`.
Accepted complete sets per group are {unpartnered,partnered} or
{unpartnered,cohab,marriage}; no mixed or partial set. Shared partnered rates
apply to both types. Omission means overall fertility as before. The old input
column errors with migration guidance. Demographic observation strata use
`union_type` with `by_union_type = TRUE`. Compensation remains separate per effective fertility schedule and
actual birth spacing persists across all transitions (D045 limitation retained).
Type-specific fertility requires cohabitation mode. On-demand stays overall-only.

Dissolution can omit `union_type` for a common duration hazard or supply both
cohab and marriage schedules within each shared/group-pair stratum. Pair
multipliers apply to either type. There is one female-carried clock per episode,
retained through births. No process means no non-death dissolution.

## Conversion and history

`union_convert_rate` is a separate rate table on union formation: age or years
since the current cohab episode began, optionally by woman's group. It uses one
clock consistently across groups, not both axes or previous-union strata.
Gaps have zero hazard. NULL disables conversion and emits a simulation warning
when cohabitation is enabled. Explicit zero bands disable risk without that
missing-schedule warning. Inputs use rate per cohabiting woman-year.

One female conversion clock is sampled per cohab episode and retained through
births. Same-person tie priority is death, birth, dissolution, conversion, search.
Conversion closes the cohab episode with `end_reason = "conversion"`, then
starts a marriage at exactly the same time with the same partners and a new ID.
`previous_cohab_id` links the new marriage to its cohab episode. Both
partners remain partnered and unavailable throughout; no matching or type draw
occurs. Marriage dissolution duration resets to zero. Fertility changes to the
marriage schedule while preserving the latest actual birth/spacing. No immediate
conversion is automatically triggered by a birth.

Typed union histories have `union_type` cohab/marriage. Legacy untyped records
are treated as marriage when enabling cohabitation. Unchanged untyped mode
retains its original schema and RNG path. A typed cohabitation history requires
cohabitation enabled to continue without silently losing its meaning.
Converted episodes must link the same partners at the same boundary, with one
marriage successor. A conversion counts neither as a new search-driven union
entry nor dissolution. Partnered exposure stays continuous; episode-duration
summaries reset at conversion. Marriage-entry summaries exclude conversion.

## Previous history, restarts and observations

Track legal-marriage category independently from current partnership type.
Under the default toggle, a cohab ending restores none/dissolved/widowed from
legal marriages; with TRUE it assigns dissolved/widowed from that cohab ending.
Every actual partnership ending resets the subsequent search-duration origin,
including a cohab ending when an older legal category is recovered. Conversion
is not an unpartnered ending and does not reset the previous-category cache.
Both partners re-enter availability after genuine endings, including a surviving
man after his cohabitant dies. Keep all histories for restart reconstruction.

Observation helpers infer the classification policy from latest simulation
metadata (default FALSE without metadata). A changed policy reclassifies supplied
history under that policy; exact historical policy changes are not reconstructed.
Conversion entries are excluded from unpartnered marriage-entry counts. Generic
partnered Lexis exposure includes both types and conversion boundaries. Shared
or typed dissolution observations count only dissolution reasons.

## Worked cases and validation

- With p=0/p=1 all new entries are marriage/cohab; absent probability consumes
  no additional draws. Input probability gaps select marriage.
- A birth at a conversion time uses cohab fertility/partner first, then converts;
  a dissolution at the same time wins and cancels conversion.
- Cohab A+B converts at year 3: same pair, cohab closes at 3, marriage begins at
  3, dissolution duration zero, no available-pool insertion.
- A never-married cohabitor who loses a partner returns to none by default, and
  widowed when the toggle is TRUE; both use the most recent ending origin if
  the selected previous category has a duration schedule.
- Restart in an active cohabitation resamples conversion/dissolution at attained
  age/duration and reconstructs legal history, rather than fabricating origins.

Check input coverage, exact clocks and ties, analytical conversion/dissolution
competition, spacing and fathers, type-specific fertility and matching recovery,
contribution sums at probability boundaries, partner exclusivity, linkage,
restart/reproducibility and old untyped controls. Use offline synthetic inputs;
no empirical rates or scientific acceptance thresholds are defaulted.

## Fertility multiplier example

Use list names present in the supplied `union_type` schedule:

```r
set_fertility(model, fertility_rates,
  group_prop_hazard = list(
    unpartnered = c(B = 1.1),
    cohab = c(B = 1.2),
    marriage = c(B = 0.9)
  )
)
```

For shared partnered rates, use `list(unpartnered = c(B = 1.1),
partnered = c(B = 1.2))` instead. Omitted types and omitted groups use 1.
A single named numeric vector applies to every supplied type. A list naming
`partnered` against a three-type table, or `cohab` against a shared partnered
table, errors rather than silently ignoring the factor.

## Union entry observations (D087)

`calculate_union_rates()` replaces `calculate_marriage_rates()` and returns
`rate_union`. By default entries pool direct marriages and cohabitations.
`by_union_type = TRUE` requires typed histories and returns both types, including
zero-entry cells. Each uses the same unpartnered person-years within its age,
period, sex, group and previous-union stratum. Thus type-specific rates sum to
the pooled rate wherever exposure is positive; exposure is repeated, not divided
between future outcomes. Conversion remains excluded. Untyped histories remain
supported for pooled observation but cannot establish a type split.


## Marriage entry observations (D088)

`calculate_marriage_rates(population, time_breaks, age_breaks, by_group = FALSE, by_previous_union = FALSE)`
counts every marriage start, including conversion, over person-years outside
marriage. Cohabitation contributes exposure continuously until conversion,
dissolution, partner death or observation censoring. Direct marriages and
conversion-created marriage episodes each count once per spouse. Age and period
cells are half-open, with optional individual-group strata; empty cells remain
and zero exposure gives an NA rate. Outputs include `person_years`, `marriages`
and `rate_marriage`. Untyped unions are interpreted as marriages, consistent
with legacy history promotion. Previous-union stratification is optional. Cohabiting exposure and conversion
retain the category held before that cohabitation began; endings follow the
recorded cohabitation classification policy. Duration stratification remains
on the partnership-entry helper.

Worked example: a woman is unpartnered from year 0 to 1, cohabits from 1 to 3,
then converts to marriage through year 5. Partnership entry has one event over
one unpartnered year; marriage entry has one event over three nonmarried years.
The conversion boundary ends cohabitation exposure and counts a marriage at 3,
without an additional partnership entry.


## Dissolution observations by type

`calculate_dissolution_rates(..., by_union_type = TRUE)` returns separate cohab
and marriage rows, optionally crossed with `group_pair`. Each episode contributes
its own type's union-years and dissolution events. Conversion ends cohabitation
exposure without counting as dissolution, and the linked marriage starts at
duration zero. Both types and empty cells are retained; zero exposure gives NA.
The default FALSE pools counts and exposure before calculating rates. Splitting
requires a valid `union_type` column, including for bare union tables.


## Demographic rates by current union type

`calculate_demog_rates(..., by_union_type = TRUE)` replaces `include_unions`.
It returns a disjoint partition labeled `union_type`: unpartnered plus the types
discovered in union histories, without overall rows. Untyped histories retain
partnered/unpartnered; typed empty histories have only unpartnered. FALSE keeps
pooled demographic rates without a type column. Exposure/deaths/births conserve
pooled totals. At conversion a tied birth remains cohabiting (birth first),
while exposure moves into marriage at that instant. Dissolution-tied births
retain the ending type; death uses the left-limit type. Ambiguous overlapping
typed histories are rejected; untyped overlap merging remains supported.
