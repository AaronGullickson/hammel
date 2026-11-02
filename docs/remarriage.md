# Age- or duration-specific searches by previous union

Milestone 4 specification, October 28, 2026 (D080/D083). Rate-based union
formation uses one rate data frame with optional `previous_union` strata and
optional duration bands for previously partnered women. These are annual
female search hazards,
not probabilities or guaranteed formation rates. Men remain available for
female-initiated matching; no male search clock is added.

Milestone 4 was accepted October 29, 2026 (D093); see the
[validation report](../validation/milestone4.qmd).

## Schedules and state

The single rate data frame requires `age_start`, `age_end`, and `rate`.
Optional `previous_union` takes exactly `none`, `dissolved`, or `widowed`.
When supplied, `none` and `dissolved` are required within each scheduled group.
An omitted `widowed` category inherits the complete `dissolved` schedule
(D125); an explicit `widowed` category overrides that fallback. Explicit zero
rates are allowed. Without the column, the same age-specific schedule
applies to all unpartnered women.

The fallback is selected independently for each group and active calendar
slice. It copies the dissolved clock, bands and rates, including duration since
the most recent union ended. An explicit widowed schedule is used as a whole:
its age/duration gaps retain zero hazard rather than borrowing dissolved rates.
Effective categories must still use a consistent clock across groups and time.
Public input tables retain the user's two or three categories; execution expands
a private copy. People and observations retain their actual `widowed` status.
This convenience applies to union-search rates, not `cohab_prob` strata.

For example, `none: [15,Inf) at .1` and `dissolved: [20,Inf) at .2`
give widowed women the .2 search hazard from age 20. Adding an explicit
`widowed: [30,Inf) at .05` instead gives zero widowed search hazard below 30
and .05 afterward. A two-category duration schedule for dissolved women also
uses elapsed time since widowhood for widowed women. Validation compares the
fallback with explicitly duplicated schedules under fixed seeds, including
calendar changes, groups, duration clocks and restart histories.

Optionally include both `dur_start` and `dur_end`. Each row supplies exactly
one complete pair of boundaries, with both unused boundaries NA. `none` must
use age. `dissolved` and `widowed` each independently use age or duration, but a
category must use the same clock in all its rows, including across groups.
Age and duration cannot act simultaneously, and bands within a category cannot
switch clocks. Supplying just one duration column, partial bands, both clocks,
neither clock, or actual duration bands without `previous_union` is an error.
All-age schedules may omit duration columns or supply them as all NA.

Starts are finite nonnegative years, ends exceed starts and may be Inf, and
rates are finite nonnegative annual search hazards. Bands must not overlap
within each group/category on the selected clock; gaps have zero risk. Group
multipliers retain their existing form and apply to every category.

Duration is elapsed time since the **most recent union ended**, either by
dissolution or partner death, not duration since the first such event, since
simulation start, or since union formation. A new ending resets the origin.
Restart histories recover both the category and its ending time; no missing
origin is invented. While partnered there is no remarriage search risk.

For example, a dissolution at 2005 and duration schedule `[0,1): 0`,
`[1,5): .2`, `[5,Inf): .1` gives zero risk until 2006, .2 until 2010 and .1
thereafter, irrespective of the woman's age. Restarting in 2008 begins at
attained duration 3; it does not restart the one-year delay. A subsequent union
ending in 2012 resets duration to zero. Individual search waiting times
integrate across the selected bands using the same annual-hazard sampler.

```r
search_rates <- tibble::tribble(
  ~previous_union, ~age_start, ~age_end, ~dur_start, ~dur_end, ~rate,
  "none",          18,         80,       NA,          NA,       .10,
  "dissolved",     NA,         NA,       0,           1,        0,
  "dissolved",     NA,         NA,       1,           Inf,      .15,
  "widowed",       NA,         NA,       0,           2,        0,
  "widowed",       NA,         NA,       2,           Inf,      .08
)
model <- set_union_formation(model, "rate", rates = search_rates)
```

These example rates are illustrative, not empirical estimates or defaults.

`none` means no recorded prior union. After an episode ends, dissolution gives
`dissolved` to each living partner; partner death gives `widowed` to the survivor.
The most recent episode determines the category, not whether a person has ever
experienced either event. A new union removes both partners from the available
pool. Its subsequent ending replaces their previous category. Births and failed
searches leave the category unchanged. Newborns start at `none`.

Initial state and restarts use the supplied union history, independent of table
row order. No public person status column is required. Histories must be complete
for the intended interpretation: no recorded episode does not establish lifetime
never-partnered status when history was omitted. Overlapping episodes cannot
have an unambiguous exclusive previous-union history and are rejected when this
stratification is requested. Touching episodes are allowed. At a shared endpoint,
the ending precedes the next formation; zero-duration episodes follow union ID
order when their starts coincide. Events after a person's death are invalid for
these summaries. Existing pooled observation methods retain their contract.

For example, a woman forms at 20, dissolves at 25, forms again at 28, and is
widowed at 35. Her unpartnered exposure is `none` before 20, `dissolved` from
25 to 28, and `widowed` after 35. The starts at 20 and 28 count respectively
in `none` and `dissolved`. Partnered time is outside all three denominators.
A later dissolution would replace `widowed` with `dissolved`.

## Timing and weighted choice

The woman's previous category selects an age- or duration-specific search schedule.
Union endings update state before the surviving woman is rescheduled. Death,
birth, dissolution, and search retain their existing same-person priority.
Any delay must be expressed as a zero-risk input band; no default waiting
period or former-spouse exclusion is added.

Calibration now maintains each unpartnered woman's actual current hazard
contributions, summed by group/birth-cohort bin (D082). A common boundary-update
queue refreshes marriage-search totals on the selected age/duration clock and
fertility-driven totals at age bands and spacing releases, with state/pool
changes replacing contributions. Previous
categories select individual schedules but are not an extra calibration count
dimension. Matching utilities and mirrored targets retain their definition.
Individual event scheduling and random-number rules are unchanged; actual
hazards replace bin-averaged calibration hazards. See
[current-hazard contributions](calibration-intensity.md). The search clock
origin is separate from the birth-cohort cell used for matching.

## Validation

Check category/group coverage, exclusive clock selection and overlap; duration
reset/restart and age-independent waiting; exact zero-hazard eligibility;
repeated transitions and restart history; equal category rates against shared
rates; analytical constant hazards with adequate candidates; calibration
intensity against an explicit sum over individuals; and observation accounting
against pooled unpartnered exposure and formations. Test half-open time/age
boundaries, histories predating the window, zero-duration and touching episodes,
sex-specific widowhood, grouping, empty histories, and invalid histories.

## Marriage observations

`calculate_union_rates(population, time_breaks, age_breaks, by_group = FALSE,
by_previous_union = FALSE)`
returns `person_years`, `unions`, and `rate_union`, split by period, age and sex. Set
`by_previous_union = TRUE` to add the three previous-union categories; otherwise
counts and exposure are pooled before calculating rates and that column is omitted. All exposure is unpartnered; no additional `union_status`
column is needed. Optional `group` follows each individual's identity, including
men. When stratifying, all three categories and empty cells are retained. Rates with zero exposure
are NA; positive exposure and no starts gives zero. Starts are classified by
category immediately before formation. Both spouses contribute one event and
their own exposure; men's realized rates are outcomes, not search inputs.

Union histories before the observation window determine categories within it.
Cells are half-open in time and age. Touching and zero-duration episodes may
contribute starts without positive exposure in the applicable category. To pool
categories or replications, sum counts and exposure before dividing. A mean of
category rates does not generally equal the pooled rate.

`calculate_fertility_rates(by_union_type = TRUE)` and the separate
`calculate_mortality_rates()` partition their respective counts and exposure by
unpartnered and discovered union types, without an overall row. Formation
counts and rates use `calculate_union_rates()` or `calculate_marriage_rates()`.
This intentionally changes its output schema; use the separate marriage helper
for formation counts and rates.

```r
search_rates <- tibble::tibble(
  previous_union = c("none", "dissolved", "widowed"),
  age_start = 18, age_end = 80, rate = c(.1, .15, .08))
model <- set_union_formation(model, "rate", rates = search_rates)
marriages <- calculate_union_rates(result, seq(2025, 2100, 5), 0:100,
                                      by_group = TRUE, by_previous_union = TRUE)
```

The rates above are illustrative inputs, not empirical estimates or defaults.


Optional `dur_breaks` on `calculate_union_rates()` now measures dissolved
and widowed exposure/events by duration since the latest ending (D084).
Supplying it automatically enables previous-category stratification, including
when `by_previous_union = FALSE`. `none` still uses required `age_breaks`.
The single returned tibble has `age_start`/`age_end` and `dur_start`/`dur_end`,
with the unused pair NA. This selects a different clock per category, not an
age-by-duration surface. Group/sex/period strata remain available.

Both clocks use half-open bands; a finite final duration boundary gains Inf.
Duration resets after each ending and does not reset at observation-window
boundaries. Durations apply regardless of attained age; age_breaks restrict
only `none` rows in this mode. Starts count using category and duration just
before formation. Partnered time adds no exposure, deaths censor it, and ongoing
unpartnered spells censor at the observation endpoint. Zero-duration/touching
spells may add events with zero exposure. Complete exclusive histories are
required. Without dur_breaks the existing age-based pooling behavior and schema
are unchanged.

```r
marriages <- calculate_union_rates(result,
  time_breaks = c(2025, 2100), age_breaks = seq(15, 70, 5),
  dur_breaks = c(0, 1, 3, 5, 10, Inf), by_group = TRUE)
```


D083 software checks passed, including analytical duration hazards, restart and
reset behavior, unchanged age-only histories, and grouped calibration equivalence.
A ten-replication duration-schedule recovery check estimated combined conditional
log odds ratios -5.9162 (MCSE .0784) and -6.0578 (MCSE .0724) for group-only and
age-dependent matching against target -6. See [testing](testing.md#duration-specific-remarriage-d083-october-28-2026)
for inputs, methods, artifacts and limitations. Those checks are implementation
evidence; subsequent user acceptance of milestone 4 is recorded in D093.
