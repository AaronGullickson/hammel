# Duration-dependent union dissolution

Core rules accepted by the user October 27, 2026 (D078). This specification
covers the separate optional dissolution process. Age-specific remarriage
schedules are now described in [remarriage](remarriage.md) (D080);
duration-dependent remarriage and cohabitation remain deferred.

## Risk, units and eligibility

Each ongoing union has one dissolution clock, carried by its living female
partner. Duration is current time minus the union's recorded `start_time`, in
years. Input `rate` is a continuous hazard per partnered woman-year, equivalently
per union-year under Hammel's exclusive unions. It is not an annual probability.
There is no second male clock and no multiplication by two. Neither partner's
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

`union_dissolution(rates, group_prop_hazard = NULL)` constructs a process;
`set_union_dissolution(model, rates, group_prop_hazard = NULL)` attaches it.
Rates are a nonempty data frame with `dur_start`, `dur_end`, and `rate`,
returned/stored as a tibble. Duration bands are half-open `[dur_start, dur_end)`.
Starts and rates are finite and nonnegative; ends exceed starts and may be
`Inf`. Bands cannot overlap within an unordered group pair. Gaps and durations outside supplied
support have zero risk. Input order, custom columns and source provenance are
preserved. Sex, age, union-status and calendar-time schedule columns are rejected
rather than interpreted as additional risk strata.

An optional `group_pair` column gives symmetric pair labels, e.g. `A:A`,
`A:B`, `B:B`. Pair order is normalized using character labels in radix order,
so `B:A` selects the same schedule as `A:B`. A pair-stratified table must cover
all unordered pairs from the possible identity registry, including future
identities. Missing pairs error, rather than silently imply zero risk. Within
each pair, gaps still mean zero risk. Shared rates omit `group_pair`. Bare
`group`, `group_w`, and `group_h` schedule columns are rejected. Pair names
require two nonempty colon-free labels separated by a single colon.

`group_prop_hazard` is a named numeric vector with the same unordered pair
names, e.g. `c("A:B" = 1.2)`. Factors are dimensionless, finite and nonnegative;
omitted pairs use 1, zero suppresses dissolution, and specifying both `A:B`
and `B:A` is an error. Factors do not introduce identities. With population
groups, preparation selects each couple's pair and scales an execution-only
copy. Without population groups, warn and ignore multipliers; pair-stratified
rates warn and use the first canonical pair alphabetically. The input table,
its ordering and provenance remain unchanged.

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

The clock uses the shared cumulative-hazard inversion (`Band`/`waiting_time()`),
with duration as its time coordinate. Sample one pending dissolution date per
union episode. Births and cancelled birth attempts retain that date; a birth
must not reset union duration. Restarts resample at attained duration, like
other processes, and are not exact checkpoints.

The pending date competes with the woman's death and birth. For exact
same-person ties, death precedes birth, then dissolution, then search. In a
birth/dissolution tie, birth uses the still-partnered status; the pending
dissolution executes afterward if the union is still ongoing. Across different
people, existing insertion-order tie semantics remain. The run covers
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
may now select the `dissolved` schedule under D080. No waiting period or
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
- Bands `[0,1)` at 0.1 and `[2,Inf)` at 0.2: threshold 0.3 from duration zero
  crosses the first year, a zero-risk gap, and one additional year, reaching
  duration three. From duration 2.5, threshold 0.1 gives a wait of 0.5 years.
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
time_breaks, dur_breaks, by_group = FALSE)` returns a tibble with `time_start`,
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
boundary is finite, as `calculate_demog_rates()` does for age.

Input may be a population object with a known endpoint or a union data frame
with unique `union_id`, finite `start_time`, numeric `end_time` (NA if ongoing),
and a valid `end_reason`. Population observation must not exceed its endpoint;
for bare tables the caller ensures complete follow-up to the requested end.
Do not stack independent replications with repeated union IDs before extraction.
Use one history per call, then combine the summaries. All grid cells are retained,
including typed zero-exposure cells for empty ungrouped histories.

With `by_group = TRUE`, population inputs use both partners' recorded groups
and return `group_pair`, with a full period/duration grid for every unordered
pair of population groups (including pairs without observed episodes). Labels
are canonical character strings in radix order, not numeric group codes.
Bare union tables require `group_pair`; reverse labels are normalized and
combined. Empty grouped inputs return a typed zero-row tibble. This replaces
the earlier woman-only grouping; no single-partner option is added (D079).

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

D045's fertility-transition compensation limitation remains: dissolution adds
another status transition but no history-aware compensation. Cohabitation,
duration-dependent remarriage, couple covariates and calendar-varying dissolution
require further specifications. These rules intentionally differ from SOCSIM's
monthly probability inputs and potentially two-spouse candidate clocks, while
sharing its duration-based baseline idea. No scientific validation threshold
or demographic rate is supplied by this implementation.
