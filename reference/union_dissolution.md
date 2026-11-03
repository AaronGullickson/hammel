# Specify duration-dependent union dissolution

Construct or attach an optional dissolution process. One continuous
hazard per union is carried by the female partner, using years since the
union began, independently of either partner's age. Rates are per
partnered woman-year (equivalently union-year), not annual
probabilities. Omitting this process retains death-only union endings.

## Usage

``` r
union_dissolution(rates, group_prop_hazard = NULL)

set_union_dissolution(model, rates, group_prop_hazard = NULL)
```

## Arguments

- rates:

  A nonempty data frame with numeric `dur_start`, `dur_end`, and `rate`.
  Bands are half-open `[dur_start, dur_end)` in years. Starts and rates
  must be finite and nonnegative; ends must exceed starts and may be
  `Inf`. Bands cannot overlap within an unordered group pair. Gaps and
  durations outside support have zero hazard. Sex, age and union-status
  columns are not supported. Input order, custom columns and provenance
  are retained. An optional `group_pair` column selects an unordered
  pair such as `"A:A"`, `"A:B"`, or `"B:B"`; `"A:B"` and `"B:A"` select
  the same schedule. Optional `union_type` supplies both `cohab` and
  `marriage` schedules within each pair stratum; without it the same
  duration risk applies to either type. Pair names require nonempty
  colon-free group labels. Pair-stratified tables must cover every
  unordered pair of possible identities, including future ones. Without
  population groups, warn and use the first canonical pair
  alphabetically. Omit `group_pair` for shared rates. `group`, `group_w`
  and `group_h` schedule columns are not supported.

- group_prop_hazard:

  Optional named numeric vector of finite nonnegative rate multipliers
  by unordered group pair, e.g. `c("A:B" = 1.2)` applies in both
  directions. Omitted pairs use 1; zero suppresses risk. Reversed
  duplicate keys are errors. With no population groups, simulation warns
  and ignores multipliers.

- model:

  A `hammel_model` object.

## Value

`union_dissolution()` returns an inspectable `hammel_union_dissolution`
process with a tibble of rates; `set_union_dissolution()` returns an
updated model.

## Details

Births preserve the pending dissolution date. Exact same-person ties
select death, then birth, then dissolution; a tied birth occurs while
partnered. Dissolution ends the episode, restores both partners'
availability, and reschedules the woman under her unpartnered fertility
and existing search rules, retaining birth spacing. Optional
age-specific search schedules by previous union are specified through
[`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md).
Restarts resample at attained union duration rather than reset duration.

## Calendar rates

Optional `time_start`/`time_end` select half-open absolute calendar-year
bands, including fractional years. Active slices retain the existing
strata rules and must cover the requested run. Changes precede same-time
events without resetting the union-duration origin. Other arguments
remain static.

## Examples

``` r
# Illustrative rates, not demographic defaults.
rates <- data.frame(dur_start = c(0, 5), dur_end = c(5, Inf),
                    rate = c(0.04, 0.02))
model <- hammel_model() |> set_union_dissolution(rates)
```
