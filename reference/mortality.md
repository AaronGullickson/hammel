# Specify an age-specific mortality process

Constructs a mortality process specification for a simulation model. The
convenience function `set_mortality` is equivalent to
`set_process(model, mortality(...))` and allows for direct piping from a
model.

## Usage

``` r
mortality(rates, weights = NULL)

set_mortality(model, rates, weights = NULL)
```

## Arguments

- rates:

  A nonempty data frame defining the mortality rate schedules. Must
  contain at least numeric `age_start`, `age_end`, and `rate`, as well
  as character `sex` columns. See details for format and description of
  stratification variables that can be applied.

- weights:

  NULL (default) for unweighted selection within rate cells, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  specification. Weights choose which member of a cell experiences an
  event without changing the cell's input rate.

- model:

  A `hammel_model` object, used by `set_mortality()`.

## Value

`mortality()` returns an inspectable process inheriting from
`hammel_mortality` and `hammel_process`, with `name`, `rates`, and
`weights`. `set_mortality()` returns an updated `hammel_model`.

## Details

The primary and only required input to the function is a fertility rate
schedule for the `rates` argument. This data frame must contain
`age_start` and `age_end` columns to define the half-open age intervals
(inclusive at the start and exclusive at the end), a `sex` column ("f"
for female and "m" for male) and a `rate` column for the age-specific
fertility rates. The example below shows the basic format:

    #> # A tibble: 44 x 4
    #>    sex   age_start age_end     rate
    #>    <chr>     <int>   <dbl>    <dbl>
    #>  1 f             0       1 0.00628
    #>  2 f             1       5 0.000288
    #>  3 f             5      10 0.000142
    #>  4 f            10      15 0.000159
    #>  5 f            15      20 0.000396
    #>  6 f            20      25 0.000471
    #>  7 f            25      30 0.000536
    #>  8 f            30      35 0.000745
    #>  9 f            35      40 0.00114
    #> 10 f            40      45 0.00174
    #> # i 34 more rows

Rates are always specified as deaths over person-years lived. An
open-ended age group specified with `age_end` of `Inf` must be included.
No gap in the age band is allowed.

Columns named by declared group dimensions (see
[`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md)),
such as `race` or `region`, supply separate schedules for each
combination of their categories; dimensions without a column share the
schedule. Every combination of categories of the dimensions must be
included in the rates if the dimension is used.

Optional `time_start`/`time_end` columns select half-open calendar
intervals in absolute years, including fractional years. Both columns
are required together; starts are finite and ends exceed starts
(possibly Inf). Active slices must satisfy the strata and coverage rules
and cover the requested simulation run.

The
[`stratify_rates()`](https://aarongullickson.github.io/hammel/reference/stratify_rates.md)
function provides a quick means to add stratification to an existing
rate schedule.

## See also

[`set_process()`](https://aarongullickson.github.io/hammel/reference/set_process.md),
[`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
[`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md),
[`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md),
[`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)

## Examples

``` r
mort_process <- mortality(usa_2000_mortality)
model <- hammel_model() |> set_process(mort_process)
# Equivalent convenience setter:
model <- hammel_model() |> set_mortality(usa_2000_mortality)
```
