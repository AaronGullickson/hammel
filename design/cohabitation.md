# Cohabitation and marriage within partnerships

Implemented milestone 4 specification (D085/D086), accepted by the user
October 29, 2026 (D093); see the [validation report](../validation/milestone4.qmd). Cohabitation is supported only in rate mode.

## Inputs and entry

`union_formation()` and its setter accept `cohab_prob`, `cohab_matching`,
`union_convert_rate`, and `cohabitation_counts_as_previous_union = FALSE`.
`cohab_prob = NULL` retains untyped partnerships, without type draws or new
clocks. When enabled, probability tables use `age_start`, `age_end`, `prob`,
optionally `group` and `previous_union`; probabilities are finite in [0,1].
Bands cannot overlap within strata. Gaps and uncovered ages have probability
zero. Group coverage follows search schedules. When `previous_union` is used,
supply either `none`/`dissolved` (pooling breakup and widowhood), or
`none`/`separated`/`widowed`, independently per group and active period (D126).

### Calendar-varying type probabilities (D122)

The user authorized `cohab_prob` to use the same flat calendar-band format as
rate inputs. Optional numeric `time_start` and `time_end` are supplied together,
in absolute years, including fractional years. Intervals are half-open, with
finite starts and greater ends (possibly Inf), and no interpolation. Values
remain probabilities in [0,1], not annual hazards. Without calendar columns,
the existing static behavior applies.

Every active calendar slice must satisfy the existing age/group/previous-union
rules. Overlapping probabilities for the same stratum and age are errors.
Calendar coverage is required across the requested run, including a valid
starting slice for zero-year runs. Gaps within the age schedule still select
marriage. Missing calendar periods are not interpreted as probability zero;
use an explicit zero-probability interval when intended.

At a search attempt, the woman's current age, group, previous category and
calendar time select the probability. A change is effective before an attempt
at exactly the same instant. It does not relabel existing episodes or trigger
conversion, matching or a type draw on its own. Typed cohabitation mode stays
enabled throughout a configured table's coverage, even during probability-zero
periods. Other cohabitation and history rules are unchanged.

At a probability-only calendar boundary, refresh marriage/cohabitation
calibration contributions to lambda*(1-p) and lambda*p. Preserve pending person,
union and fertility clocks: changing intended type probabilities does not change
the total search hazard. This bookkeeping consumes no random draws. At a shared
rate/probability boundary, apply both active slices before the existing rate
rescheduling. Future actual events may consume different random draws when type
outcomes differ; no general fixed-seed agreement is implied.

Worked examples (illustrative probabilities, not demographic defaults):

- Probability zero in `[2000,2005.5)` and one in `[2005.5,Inf)` produces direct
  marriages for successful searches before 2005.5 and cohabitations afterward,
  including at 2005.5. Existing marriages stay marriages.
- A woman has a pending search at 2006 when probability changes in 2005.5.
  That search stays at 2006 and uses the new probability. There is no additional
  search at the change date.
- Search hazard 0.2 and probability 0.25 give intended marriage/cohabitation
  intensities 0.15 and 0.05 per unpartnered woman-year. Changing probability to
  0.75 reverses those components without changing their sum of 0.2. Realized
  formations still depend on matching success.

Validate construction and edited specifications, malformed and overlapping
intervals, calendar/group/category coverage, fractional boundaries, exact
search ties, zero-year and empty runs, age gaps, current group/previous category,
preserved existing episodes and clocks, type probabilities and search recovery,
calibration refresh with distinct group probabilities, asynchronous rate changes,
seed isolation, and unchanged static controls. Ordinary checks use synthetic
inputs and `devtools::test()` without live downloads.

The overall search hazard still schedules attempts. At the attempt time, draw
intended type before matching using current age/group/previous category.
Matching uses the general model for marriage and `cohab_matching` for cohabitation
(NULL shares the general preferences active at that instant). Both matching
arguments also accept complete-model calendar schedules under
[D123](time-varying-models.md). Both draw from the same exclusive pool
of unpartnered men. Failure leaves the woman unpartnered in her search cell,
and type is drawn anew at the next attempt. Father selection remains
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
`union_type` with `by_union_type = TRUE`. Fertility cells are separate per effective schedule and
actual birth spacing persists across all transitions; cell accounting covers them (D101).
Type-specific fertility requires cohabitation mode. On-demand stays overall-only.

Dissolution can omit `union_type` for a common duration hazard or supply both
cohab and marriage schedules within each shared/group-pair stratum. Pair
multipliers apply to either type. Each episode is one unit in the dissolution
cells (D128), unaffected by births. No process means no non-death dissolution.

## Conversion and history

`union_convert_rate` is a separate rate table on union formation: age or years
since the current cohab episode began, optionally by woman's group. It uses one
clock consistently across groups, not both axes or previous-union strata.
Gaps have zero hazard. NULL disables conversion and emits a simulation warning
when cohabitation is enabled. Explicit zero bands disable risk without that
missing-schedule warning. Inputs use rate per cohabiting woman-year.

Each cohab episode is one unit in the conversion cells, by the woman's age or
episode duration (D128), unaffected by births. Same-person tie priority is death, birth, dissolution, conversion, search.
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
Under the default toggle, a cohab ending restores none/separated/widowed from
legal marriages; with TRUE it assigns separated/widowed from that cohab ending.
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
- Restart in an active cohabitation rebuilds conversion/dissolution cells at
  attained age/duration and reconstructs legal history, rather than fabricating origins.

Check input coverage, exact clocks and ties, analytical conversion/dissolution
competition, spacing and fathers, type-specific fertility and matching recovery,
contribution sums at probability boundaries, partner exclusivity, linkage,
restart/reproducibility and old untyped controls. Use offline synthetic inputs;
no empirical rates or scientific acceptance thresholds are defaulted.

## Fertility stratification example

Different group effects by union type are built by stratifying each type's
rows and binding them (D132):

```r
fertility_rates <- rbind(
  stratify_rates(fert[fert$union_type == "unpartnered", ], list(group = c(A = 1, B = 1.1))),
  stratify_rates(fert[fert$union_type == "cohab", ], list(group = c(A = 1, B = 1.2))),
  stratify_rates(fert[fert$union_type == "marriage", ], list(group = c(A = 1, B = 0.9))))
set_fertility(model, fertility_rates)
```

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

`calculate_fertility_rates(..., by_union_type = TRUE)` replaces `include_unions`.
It returns a disjoint partition labeled `union_type`: unpartnered plus the types
discovered in union histories, without overall rows. Untyped histories retain
partnered/unpartnered; typed empty histories have only unpartnered. FALSE keeps
pooled demographic rates without a type column. Exposure/deaths/births conserve
pooled totals. At conversion a tied birth remains cohabiting (birth first),
while exposure moves into marriage at that instant. Dissolution-tied births
retain the ending type; death uses the left-limit type. Ambiguous overlapping
typed histories are rejected; untyped overlap merging remains supported.


## Cohabitation conversion observations

`calculate_cohab_convert_rate(population, time_breaks, dur_breaks = NULL,
age_breaks = NULL, by_group = FALSE)` requires typed population histories and
exactly one duration/age axis. Duration is elapsed years since each cohabitation
began; age is the woman's attained age, matching the conversion input schedule.
Each cohabitation contributes union-years once, ending at conversion,
dissolution, either partner's death or observation censoring. Only end reason
conversion contributes an event. Marriage episodes contribute neither exposure
nor events. Rate is conversions per cohabitation union-year; zero exposure gives
NA. Period and axis cells are half-open; the final finite axis break gains Inf.
Optional groups refer to the woman, matching the input conversion schedule.
Output has `union_years`, `conversions`, `rate_cohab_convert` and period/axis
boundaries, plus optional `group`. Retain empty grid cells. Example: a cohab
lasting two years then converting contributes two union-years and one conversion;
a cohab lasting three years then dissolving contributes three union-years and
zero conversions, for a pooled conversion rate of 1/5 per union-year.
