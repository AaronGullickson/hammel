# Calendar-varying rates and cohabitation probabilities

Milestone 7 rate specification, October 30, 2026. The user approved these
semantics and authorized rate implementation. D122 (October 31) adds flat
calendar-banded cohabitation probabilities. Matching, inheritance and other
arguments remain static in this step. Milestone 7 scientific acceptance remains
separate from implementation and software verification.

## Inputs and coverage

Each of the six existing rate tables in [the census](rate-inputs.md) optionally
takes numeric `time_start` and `time_end` columns. Both are required together.
They use the same absolute calendar years as population histories, with finite
starts (negative years allowed), greater ends (possibly `Inf`), and half-open
intervals `[time_start, time_end)`. Rates are constant within supplied calendar
and age/duration bands; there is no interpolation. Static tables omit both
columns and continue to apply throughout the run.

Existing age/duration clocks, supported strata, rate units and eligibility
rules remain in force at each calendar instant. An age/duration band may have
its own calendar boundaries; all simultaneously active rows must form a valid
existing rate table. Overlapping calendar intervals are legal for different
age/duration bands or strata, but overlapping rates for the same risk are
invalid. Calendar gaps in a configured table during the requested run are
errors. Use explicit zero rates to disable a process for an interval. Existing
rules for zero risk in omitted age/duration bands and transition pairs/categories
still apply. Required groups, pairs, sexes, union types and parity partitions
must be covered in each active table under their existing rules. Mortality
continues to require complete age coverage.

Construction validates every supplied calendar slice, including future slices.
Execution additionally checks calendar coverage over the requested run. A
zero-year run requires a valid active table at its starting instant. Different
processes may have different calendar boundaries. Existing clock-selection
restrictions apply across the entire table, including across calendar periods.
Input tables and source attributes are retained; group multipliers remain static.
The static-only parity reachability warning is not applied to calendar-varying
fertility, because earlier rate periods can produce different parity histories.

`cohab_prob` accepts the same optional calendar columns under D122, while
remaining a probability table with `prob` in [0,1]. All active age/group/previous-
union strata must be valid; calendar coverage is required and age gaps mean
probability zero. At a search, select type using the current probability.
Probability changes refresh marriage/cohabitation calibration without changing
the total search hazard or existing union types. Probability-only boundaries
retain pending event dates and consume no random draws. A shared rate boundary
keeps the existing rate rescheduling, with both changes effective first.
See [the cohabitation contract](cohabitation.md#calendar-varying-type-probabilities-d122).

Matching, inheritance, `group_prop_hazard`, birth spacing, `ffab`, `hetfert`,
and previous-union classification settings remain static.

The database getters can prepare these calendars directly from year vectors;
see [HMD/HFD preparation](database-rates.md) (D119). Single-year requests remain
static. HMD also supports its native grouped age/time products.

## Event timing and state

Calendar changes take precedence over demographic events at the same instant.
The run remains `[start, end)`: changes and events at its final endpoint do not
execute. Pending candidates under an old rate are invalidated at changes;
new waiting times use attained age/duration and the active rate table. Sampling
must never execute a future candidate using a superseded rate. The first implementation redraws all pending person and union candidates at
each calendar rate boundary; conditional exponential resampling after survival
to that boundary preserves the hazard interpretation. This requires work
proportional to the population at each boundary. Fixed-seed equality with static
or split-run alternatives is not promised, including redundant rate boundaries.

The simulation retains its population, histories, ancestry, parity, fertility
weights, latest actual birth times and all duration origins at boundaries.
No population restart occurs. New fertility cells use the active input rate
times their current membership, preserving spacing eligibility. Matching
calibration uses current search and fertility intensities. Unrealized fertility
intensity accrued before each change is retained, with observation intervals
split at fertility schedule changes. Observations continue to use recorded
events and exposure; rate changes do not create demographic events.

## Worked examples and validation

For a runnable fractional-year workflow, see
[the milestone 7 rate example](../inst/examples/milestone7-rates.R).
`simulate(model, years = 4.3)` starting at population `end_time = 2015.5`
ends at 2019.8. Calendar bounds may also be fractional, such as
`[2015.5, 2019.8)`, with no rounding to an annual or monthly grid.

- Mortality 0.1 in `[2000, 2002)`, then 0.3 in `[2002, Inf)`: without other
  exits, survival from 2000 to 2004 is `exp(-0.1 * 2 - 0.3 * 2)`.
- Search hazard zero before 2005 and positive thereafter: no searches occur
  before 2005. A widow whose partner died in 2003 has attained search duration
  two at that boundary, not zero.
- A cohabitation begun in 2001 faces a changed duration-based conversion rate
  in 2004 at duration three. Marriage conversion still starts a new episode.
- A woman's last birth at 2004.8 and spacing of 0.75 keep her ineligible until
  2005.55 even if fertility changes in 2005. A cell with one ineligible member
  and rates 0.2 for half a year then 0.4 for half a year accrues 0.3 unrealized
  births over that year, if membership and ineligibility persist.
- An A-to-B transition rate zero before 2005 becomes positive afterward.
  Duration remains time since the latest recorded entry into A. Existing
  ancestry and historical group observations remain unchanged.

Validate malformed intervals, calendar/age/duration overlap and coverage,
future slices, edited processes, asynchronous boundaries, group/type/parity
strata, analytical integrated mortality and competing hazards, conversion and
dissolution duration origins, cell accounting and spacing, calibration updates,
seed isolation, independent runs, restart histories, empty populations and
zero-year runs. Ordinary verification uses synthetic inputs and the offline
`devtools::test()` suite; no database downloads or new scientific acceptance
thresholds are required.
