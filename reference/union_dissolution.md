# Specify duration-dependent union dissolution

Constructs a union dissolution process specification for a simulation
model. The convenience function `set_union_dissolution` is equivalent to
`set_process(model, union_dissolution(...))` and allows for direct
piping from a model.

## Usage

``` r
union_dissolution(rates, weights = NULL)

set_union_dissolution(model, rates, weights = NULL)
```

## Arguments

- rates:

  A nonempty data frame with at least numeric `dur_start`, `dur_end`,
  and `rate`. See details for format and description of stratification
  variables that can be applied.

- weights:

  NULL (default) for unweighted selection, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  using union variables, which chooses which union in a cell dissolves
  without changing the cell's input rate.

- model:

  A `hammel_model` object.

## Value

`union_dissolution()` returns an inspectable `hammel_union_dissolution`
process with a tibble of rates; `set_union_dissolution()` returns an
updated model.

## Details

The primary and only required input to the function is a union
dissolution rate schedule for the `rates` argument. This data frame must
contain `dur_start` and `dur_end` columns to define the half-open age
intervals (inclusive at the start and exclusive at the end) and a `rate`
column for the duration-specific dissolution rates.

Rates are always specified as union dissolutions over union-years lived.
Gaps in the duration band are treated as cases with zero rates.

Columns named by declared group dimensions (see
[`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md)),
such as `race` or `region`, supply separate schedules for each
combination of their categories, but they must be suffixed by `_h` or
`_w` to indicate the group of the husband or wife, respectively (e.g.
`race_h` and `race_w`).

Optional `time_start`/`time_end` columns select half-open calendar
intervals in absolute years, including fractional years. Both columns
are required together; starts are finite and ends exceed starts
(possibly Inf). Active slices must satisfy the strata and coverage rules
and cover the requested simulation run.

If cohabitations are in use in the simulation, the rates can be
stratified by `union_type` to supply both `cohab` and `marriage`
schedules; without it the same duration risk applies to either type. The
[`stratify_rates()`](https://aarongullickson.github.io/hammel/reference/stratify_rates.md)
function provides a quick means to add stratification to an existing
rate schedule.

## Examples

``` r
# Illustrative rates
diss_rates <- data.frame(
  dur_start = c(0, 5), dur_end = c(5, Inf),
  rate = c(0.04, 0.02)
)
model <- hammel_model() |> set_union_dissolution(diss_rates)
```
