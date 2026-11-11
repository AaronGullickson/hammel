# Duration-dependent union dissolution

Core rules accepted by the user October 27, 2026 (D078). This specification
covers the separate optional dissolution process. Age-specific remarriage
schedules, including duration since the latest ending, are described in
[remarriage](remarriage.md) (D080/D083). Optional cohabitation is implemented under [D086](cohabitation.md).

Milestone 4 was accepted October 29, 2026 (D093); see the
[validation report](../validation/milestone4.qmd).

## Risk, units and eligibility

Each ongoing union is one unit at risk in the dissolution process's cells
(D128; until milestone 8 the woman carried a per-union clock). Duration is current time minus the union's recorded `start_time`, in
years. Input `rate` is a continuous hazard per partnered woman-year, equivalently
per union-year under Hammel's exclusive unions. It is not an annual probability.
There is no second male risk and no multiplication by two. Neither partner's
age selects this baseline schedule. Shared schedules apply to all unions;
optional schedules and multipliers select an unordered group pair (D079).
Thus A-wife/B-husband and B-wife/A-husband share the same dissolution risk;
carrier sex does not determine the group stratum.

Only living, currently partnered women carry risk. The man has no extra
candidate. Existing ongoing unions start at their attained duration, not zero.
Each newly formed episode starts at duration zero. Recorded ended unions are
unchanged. The process is optional; absent means death-only endings and no
additional random draws. There is no demographic default dissolution rate.

## Interface and schedules

`union_dissolution(rates, weights = NULL)` constructs a process;
`set_union_dissolution(model, rates, weights = NULL)` attaches it.
Rates are a nonempty data frame with `dur_start`, `dur_end`, and `rate`,
returned/stored as a tibble. Duration bands are half-open `[dur_start, dur_end)`.
Starts and rates are finite and nonnegative; ends exceed starts and may be
`Inf`. Bands cannot overlap within an unordered group pair. Gaps and durations outside supplied
support have zero risk. Input order, custom columns and source provenance are
preserved. Sex, age and union-status columns are rejected. Optional
`time_start`/`time_end` bands follow [the rate contract](time-varying-rates.md)
(D118); calendar changes retain the union-duration origin.

Optional `<group>_w` and `<group>_h` columns, for any declared group dimension
(D134), select schedules by the wife's and husband's categories. Strata are
directional: a symmetric rate lists both `A`/`B` and `B`/`A` rows. Either
column may be omitted to pool over that partner. A stratified table must cover
every combination of the used columns' categories, including future
identities. Missing combinations error, rather than silently imply zero risk.
Within each stratum, gaps still mean zero risk. A bare dimension column or the
former `group_pair` is an error at simulation, as is a column for an
undeclared dimension.

Proportional effects are written with `stratify_rates()` (D132), for example
`stratify_rates(rates, list(group_w = c(A = 1, B = 2)))`. The input table, its
ordering and provenance remain unchanged.

Illustrative input, not a demographic default:

```r
rates <- tibble::tibble(
  dur_start = c(0, 5, 10),
  dur_end = c(5, 10, Inf),
  rate = c(0.04, 0.02, 0.01)
)
model <- set_union_dissolution(model, rates)
```

## Timing and state changes

A union sits in the cell of its group pair, union type and duration band, and
moves at duration band edges. The dissolution process's clock fires at the
total intensity of its cells and picks a union within the chosen cell (D128).
Births and cancelled birth attempts do not affect it; a birth must not reset
union duration. Restarts rebuild cells at attained duration, like other
processes, and are not exact checkpoints.

Exact ties between processes put death before birth, then dissolution, then
conversion and search. In a birth/dissolution tie, birth uses the
still-partnered status; the dissolution executes afterward if the union is
still ongoing. The run covers
`[start, start + years)` and excludes endpoint events.

When dissolution executes:

- End exactly the active episode at that time with `end_reason = "dissolution"`.
- Clear both current-partner links and return the living man to the available
  father/partner pool. Restore the woman's unpartnered calibration counts,
  respecting her actual-birth spacing interval.
- Invalidate and reschedule the woman's pending person event using her attained
  age, current unpartnered fertility schedule and existing search rules. Retain
  her latest actual birth time and all parent and union histories.
- Retain the man's death event: his current hazards do not depend on union
  status. His availability changes immediately, without a second divorce risk.
- Cancel the old episode's dissolution when either partner dies. The surviving
  woman must not execute a stale dissolution after widowhood, in either mode.

Unpartnered rate-mode women use the currently supplied search rates; on-demand
women may form another union at a later birth attempt. Age-specific searches
select the detailed `separated` schedule, or the pooled `dissolved` schedule
when separation and widowhood share rates (D126). No waiting period or
former-partner exclusion is introduced.
Current kinship/matching restrictions still apply. New union episodes have new
IDs and their own duration clocks, even when the partners previously formed a
union. Mortality, identity, and ancestry are not altered by dissolution itself.

## Observations and worked examples

Unions expose the end time and reason using the existing schema. With this
process configured, latest-run metadata includes `dissolutions`, counting one
per executed union ending, not one per former spouse. Existing Lexis partnered
exposure uses the shortened half-open union interval. A birth exactly at an
episode ending with `dissolution` counts as partnered, reflecting birth-first
execution; this adds no exposure beyond the end. Other union-end birth
classification remains unchanged. The union-duration summary below measures
dissolutions against union-year exposure.

- Constant rate 0.1 over all durations: without competing risks, survival for
  another two years is `exp(-0.1 * 2)`, not `1 - 0.1 * 2`.
- Rate zero for `[0,5)`, then 0.2: a new union cannot dissolve before duration
  five; an input union already at duration seven faces rate 0.2 immediately.
- Bands `[0,1)` at 0.1 and `[2,Inf)` at 0.2: a union faces 0.1 in its first
  year, no risk in its second, and 0.2 from duration two; the integrated hazard
  from duration zero reaches 0.3 at duration three.
- With factor `"A:B" = 2`, both wife-A/husband-B and wife-B/husband-A
  couples use twice the baseline; A:A and B:B use the baseline unless specified.
- A birth and pending dissolution at time 1: record the partnered birth with
  that husband, then end the episode at 1. Count this birth as partnered in
  summaries;
  partnered exposure still ends at 1. This exception is specific to dissolution
  end reasons.
- Dissolution at a run endpoint is excluded. Widowhood earlier in the run closes
  the episode with a death reason and prevents the old dissolution executing.

## Union-year counts and dissolution rates

Requested by the user October 27, 2026. `calculate_dissolution_rates(population,
time_breaks, dur_breaks, by = NULL)` returns a tibble with `time_start`,
`time_end`, `dur_start`, `dur_end`, `union_years`, `dissolutions`, and
`rate_dissolution`. There is no duplicated male/female row: one episode is one
unit of observation, and repeated episodes contribute separately.

Calendar and duration bands are half-open. Exposure is the exact overlap of
`[start_time, end_time)` with each calendar period and duration band translated
from the episode's own start. Missing ends are censored at the observation
endpoint. Both dissolution and death end exposure; only reason `dissolution`
counts toward the event numerator. A dissolution at an internal time/duration
boundary belongs to the cell beginning there; one at the final time boundary
is excluded. Zero-duration episodes can contribute an event but no exposure.
Zero exposure gives `NA` for the rate; positive exposure and zero dissolutions
give zero. Duration boundaries use years and append `Inf` when the last supplied
boundary is finite, as `calculate_fertility_rates() and calculate_mortality_rates()` does for age.

Input must be a population object with a known endpoint; observation must not
exceed it. Bare union tables are not accepted (D159).
Do not stack independent replications with repeated union IDs before extraction.
Use one history per call, then combine the summaries. All grid cells are retained,
including typed zero-exposure cells for empty ungrouped histories.

`by` accepts `"group_w"`, `"group_h"` and `"union_type"` (D133). The group
variables are the wife's and husband's current groups, following transition
histories, with a full period/duration grid for every combination of
population groups (including combinations without observed episodes). Rows are
directional; sum both directions for the unordered pairs that dissolution
inputs use.

Worked accounting: U1 starts at 0 and dissolves at 3; U2 starts at 1 and ends
at 4 through death; U3 starts at 2 and remains ongoing at endpoint 5. Over
`[0,5)`, they contribute 3 + 3 + 3 = 9 union-years and one dissolution. With
periods `[0,2)`/`[2,5)` and durations `[0,2)`/`[2,Inf)`, exposures are 3, 0,
3, and 3 for (early/short, early/long, late/short, late/long). Only late/long
contains the dissolution, giving rate 1/3 in that cell. The pooled rate is 1/9.
Death and censoring do not add dissolution events. Splitting periods/durations
must preserve total exposure and event counts over the same covered domain.

## Validation and limits

Check constructor/attachment/execution validation, tibble/provenance retention,
rate-band inversion, analytical duration survival, single risk per couple,
attained-duration restarts, formation and repeat episodes, death invalidation,
birth/dissolution ties, spacing and status fertility, availability/calibration
updates, group coverage/multipliers, seed isolation and independent runs, and
absent-process history/RNG preservation. Use synthetic local inputs only.

D101 cell accounting, accepted in milestone 5 (D105), replaces the former
D045 fertility-transition approximation: dissolution moves the woman between
current fertility cells while retaining her spacing eligibility. Sparse cells
without eligible women still count unrealized births. Cohabitation is
implemented under D086; richer couple covariates remain deferred. Calendar-varying
dissolution is implemented under [D118](time-varying-rates.md).
Duration-based remarriage observations are available via `dur_breaks` under D084. These rules intentionally differ from SOCSIM's
monthly probability inputs and potentially two-spouse candidate clocks, while
sharing its duration-based baseline idea. No scientific validation threshold
or demographic rate is supplied by this implementation.


D086 adds optional union_type cohab/marriage strata within each shared or
group-stratified schedule. Supply both when stratifying; without the column,
each episode uses the same dissolution schedule. Conversion is not dissolution:
it ends the cohab episode and starts marriage with a fresh duration clock.

For calendar-varying tables, both types are required in each nonempty active
calendar interval. Missing-type errors name the absent type, the half-open
calendar interval, and the group pair only when pair-specific rates are used.
For example, cohabitation rates starting in 1990 and marriage rates starting
in 1992 report a missing marriage schedule during `[1990,1992)`. Validation
still examines supplied intervals outside the eventual simulation window.
