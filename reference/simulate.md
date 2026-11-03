# Run a demographic simulation

Advance a configured model and return a population that can initialize a
later run. The C++ engine samples candidate waiting times across age
bands (fertility cell clocks for births), executes at continuous times,
and prioritizes death over birth on exact ties. Events execute
chronologically. In on-demand mode all living females face the fertility
schedule. Rate mode uses either overall rates for all women (with a
warning), or separate `union_type` schedules. A union-formation process
is required. With on-demand union formation, use the woman's current
partner, or draw one from all eligible men using conditional-logit
probabilities from the matching formula to form a union at the birth.
Men must meet the configured minimum age (default 15 years). Known kin
are excluded within the configured generation bound on each side of a
common ancestor (default two). Hard group exclusions also apply. If no
eligible man exists, cancel the birth by default. Explicitly set
`cancel_birth_if_no_match = FALSE` to preserve it with unknown father.
Rate initiation adds age- or duration-specific search events for
unpartnered women only. Duration schedules for separated/widowed
categories measure time since the most recent union ended; `none`
remains age-specific. Failed attempts draw another waiting time.
Widowhood cancels pending births and switches to unpartnered fertility;
a new union switches to partnered fertility, preserving birth spacing.
Unpartnered births do not trigger searches in rate mode. Instead, select
a father using the same available male pool, screening and
`father_matching` model (falling back to `matching`), without forming a
union or reserving the man. A failed match cancels the birth unless
explicitly allowed by the union process. Cancellation starts no new
birth-spacing interval. On-demand mode rejects status schedules.
Optional
[`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md)
adds one duration-dependent annual union hazard carried by the woman. It
preserves the pending dissolution through births, ends the episode with
`end_reason = "dissolution"`, and restores unpartnered availability and
fertility/search eligibility. Exact same-person ties select death,
birth, dissolution, then search. Dissolution `group_pair` schedules and
multipliers use symmetric pairings of both partners' groups. Rate-based
formation optionally selects age/duration schedules using
`previous_union`: none, separated, or widowed (or pooled
none/dissolved). Duration origins and categories are recovered from
recorded histories. Matching uses male-group weights recalibrated at
each search; simulation warns once if any search needed approximate
weights. A separate father model has its own preferences and delta
target adjustments. Event initiation is unchanged. Optional cohabitation
uses rate mode only: type is selected from `cohab_prob` before
marriage/cohabitation-specific matching. Conversion retains partners and
creates linked episodes, with same-person ties ordered death, birth,
dissolution, conversion, search. Fertility input strata use
`union_type`.

## Usage

``` r
simulate(model, years, seed = NULL)
```

## Arguments

- model:

  A `hammel_model` containing a population, mortality, fertility, and
  union-formation processes, optionally union dissolution and group
  transitions. The population must have a finite `end_time`.

- years:

  Nonnegative simulation duration in years; fractional years are
  supported, as are fractional population endpoints and calendar-rate
  bands.

- seed:

  Optional integer seed. With a seed, the caller's RNG state is restored
  afterward. With `NULL`, the simulation advances the current RNG.

## Value

A `hammel_population` with updated histories and endpoint, plus a
`simulation` component recording run settings and birth/death counts,
plus `dissolutions` when a dissolution process is configured.

## Details

Optional
[`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)
adds directed competing age- or duration-based moves, with optional sex
and union-type strata. A move updates current group and affected
event/pool/cell state, retaining ancestry. Grouped output always
includes a transitions tibble; restarts preserve the full transition
history.

Optional fertility `hetfert = het_fert(cv, alpha)` weights eligible
mother selection within fixed-intensity cells. Supplied positive female
weights are retained; absent weights are drawn under the simulation
seed, and daughters inherit on the log scale. Output persons includes
`hetfert` whenever enabled, and omits it whenever disabled, even if the
input supplied weights.

Minimum birth spacing is configured in
[`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
or
[`set_fertility()`](https://aarongullickson.github.io/hammel/reference/set_process.md)
and uses recorded maternal births, including on restart. Births are
generated by cell-accounted fertility (D101; see
[`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md)):
each one-year age/parity/union-type/group cell produces births at its
input rate times all its members, assigned to an eligible member
uniformly or proportionally to `hetfert` weights when enabled. Births
that cannot occur because a cell has no eligible member are reported in
`simulation$unrealized_births`. Fertility outside supplied bands is
zero. Mortality must cover all ages for both sexes, ending at `Inf`.

Events are recorded at their sampled times in `[start, start + years)`.
Rates and all times remain in person-year and year units. Newborn group
identity uses fertility's inheritance model (random parent by default,
with maternal fallback), if groups are present; other newborn custom
columns are missing. Grouped rate tables use each person's group (the
woman's group for union searches). Ungrouped tables apply to everyone.
Missing population groups in a grouped rate table are errors. Without
population groups, each grouped table uses its first alphabetical group
with a warning; this does not create group identities. Matching can use
group predictors and pair effects. Existing custom columns and histories
are preserved. Restarting samples new pending events and does not
reproduce an uninterrupted run exactly. When ancestry is tracked, each
newborn receives half of each parent's shares, with a missing parent
contributing its half to unknown ancestry. Initial shares are retained,
including on restart; grouped populations without shares initialize from
recorded identities. No additional random draws are made.

Calendar-banded rates change within the same run at
`time_start`/`time_end` boundaries, before same-time demographic events.
Histories, birth spacing and duration origins persist. Pending
candidates are resampled at boundaries; adding redundant rate boundaries
need not preserve a fixed-seed trajectory. Calendar-banded `cohab_prob`
is evaluated at the search instant; probability- only boundaries refresh
role calibration while retaining pending event dates. Matching and
inheritance arguments may be direct models or calendar schedules with
`time_start`, `time_end` and list-column `model`. The current matching
model is selected at partner/father choice and inheritance at birth.
Model- only boundaries preserve pending events and earlier histories.

## Examples

``` r
population <- create_baseline_population(100, max_age = 50)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
    age_end = Inf, rate = 0.01)) |>
  set_fertility(data.frame(age_start = 15, age_end = 50, rate = 0.1)) |>
  set_union_formation("on_demand",
    create_matching_model(~ age_difference, c(age_difference = 0)))
result <- simulate(model, years = 2, seed = 42)
result$simulation
#> $start_time
#> [1] 0
#> 
#> $end_time
#> [1] 2
#> 
#> $birth_interval
#> [1] 0.75
#> 
#> $ffab
#> [1] 0.4886
#> 
#> $births
#> [1] 2
#> 
#> $deaths
#> [1] 0
#> 
#> $unrealized_births
#> # A tibble: 3 × 6
#>   period_start period_end age_start age_end unrealized_births
#>          <dbl>      <dbl>     <dbl>   <dbl>             <dbl>
#> 1            1          2        22      23           0.00758
#> 2            1          2        23      24           0.0605 
#> 3            1          2        25      26           0.00831
#> # ℹ 1 more variable: years_without_eligible <dbl>
#> 
```
