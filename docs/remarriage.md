# Age-specific searches by previous union

Milestone 4 specification, October 28, 2026. The user selected an optional
`previous_union` column on rate-based union-formation schedules. Duration-based
remarriage schedules remain deferred. These are annual female search hazards,
not probabilities or guaranteed formation rates. Men remain available for
female-initiated matching; no male search clock is added.

## Schedules and state

`age_start`, `age_end`, and `rate` retain their existing meaning. Optional
`previous_union` takes exactly `none`, `dissolved`, or `widowed`. When supplied,
all three categories are required within each scheduled group, with explicit
zero rates allowed. Bands must not overlap within group/category; gaps have
zero risk. Without the column, the same schedule applies to all unpartnered
women. Existing woman-group multipliers apply to each category.

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

The current woman's age and previous category select the search schedule.
Union endings update state before the surviving woman is rescheduled. Death,
birth, dissolution, and search retain their existing same-person priority.
There is no new waiting period, duration clock, or former-spouse exclusion.

Calibration maintains unpartnered female counts by group, birth-cohort bin,
and previous category. Rate-mode intensity for a group/age bin is the sum of
category count times that category's age-band-averaged hazard. Because matching
utilities do not depend on previous category, these intensities can be summed
before the existing expected-flow solver. Female denominators and male mirrored
targets retain their existing definition; men are not assigned category-specific
targets. Nonmarital father selection and on-demand formation retain their
fertility-driven intensities. Shared schedules keep the existing calculation
path and random-number policy.

## Validation

Check category/group coverage and overlap; exact zero-hazard eligibility;
repeated transitions and restart history; equal category rates against shared
rates; analytical constant hazards with adequate candidates; calibration
intensity against an explicit sum over category bins; and observation accounting
against pooled unpartnered exposure and formations. Test half-open time/age
boundaries, histories predating the window, zero-duration and touching episodes,
sex-specific widowhood, grouping, empty histories, and invalid histories.

## Marriage observations

`calculate_marriage_rates(population, time_breaks, age_breaks, by_group = FALSE,
by_previous_union = FALSE)`
returns `person_years`, `unions`, and `rate_marriage`, split by period, age and sex. Set
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

Per the user's request, `calculate_demog_rates(include_unions = TRUE)` no longer
returns `unions` or `rate_marriage`. Its overall/partnered/unpartnered exposure,
births, deaths, fertility and mortality rates retain their existing semantics.
This intentionally changes its output schema; use the separate marriage helper
for formation counts and rates.

```r
search_rates <- tibble::tibble(
  previous_union = c("none", "dissolved", "widowed"),
  age_start = 18, age_end = 80, rate = c(.1, .15, .08))
model <- set_union_formation(model, "rate", rates = search_rates)
marriages <- calculate_marriage_rates(result, seq(2025, 2100, 5), 0:100,
                                      by_group = TRUE, by_previous_union = TRUE)
```

The rates above are illustrative inputs, not empirical estimates or defaults.
