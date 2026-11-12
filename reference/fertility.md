# Specify an age-specific fertility process

Constructs a fertility process specification for a simulation model. The
convenience function `set_fertility` is equivalent to
`set_process(model, fertility(...))` and allows for direct piping from a
model.

## Usage

``` r
fertility(rates, birth_interval = 0.75, ffab = 0.4886, weights = NULL)

set_fertility(
  model,
  rates,
  birth_interval = 0.75,
  ffab = 0.4886,
  weights = NULL
)
```

## Arguments

- rates:

  A nonempty data frame defining the fertility rate schedules. Must
  contain at least numeric `age_start`, `age_end`, and `rate`. See
  details for description of stratification variables that can be
  applied.

- birth_interval:

  Minimum maternal birth-to-birth interval in years, defaulting to
  `0.75`. Must be finite and nonnegative.

- ffab:

  Probability that a newborn is female, default `0.4886`. Must be a
  finite numeric scalar between 0 and 1.

- weights:

  NULL (default) for unweighted selection within rate cells, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  specification. Weights choose which member of a cell experiences an
  event without changing the cell's input rate.

- model:

  A `hammel_model` object, used by `set_fertility()`.

## Value

`fertility()` returns an inspectable process inheriting from
`hammel_fertility` and `hammel_process`, with `name`, `rates`,
`birth_interval`, `ffab`, and `weights`. `set_fertility()` returns an
updated `hammel_model`.

## Details

The primary and only required input to the function is a fertility rate
schedule for the `rates` argument. This data frame must contain
`age_start` and `age_end` columns to define the half-open age intervals
(inclusive at the start and exclusive at the end) and a `rate` column
for the age-specific fertility rates. The example below shows the basic
format:

    #> # A tibble: 6 x 3
    #>   age_start age_end   rate
    #>       <int>   <int>  <dbl>
    #> 1        15      20 0.0477
    #> 2        20      25 0.110
    #> 3        25      30 0.114
    #> 4        30      35 0.0912
    #> 5        35      40 0.0397
    #> 6        40      45 0.008

Rates are always specified as births over person-years lived for women.
An open-ended age group specified with `age_end` of `Inf` can be
included but is not required. Any age groups without a specified rate
will be assumed to have a zero rate.

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

Fertility optionally accepts `union_type` containing both `"partnered"`
and `"unpartnered"`, or all three of `"unpartnered"`, `"cohab"`,
`"marriage"`. Shared partnered rates cover either partnership type.
Rates use woman-years in the indicated type; the old fertility input
`union_status` is rejected. Without that column, one overall schedule
applies to all women. On-demand union formation requires overall rates;
rate-based formation accepts either form, but using overall fertility
rates with specified union formation rates will lead to independence
between marriage and fertility, which is likely undesirable in most
applications. A warning is printed if this is the case.

Optional `parity_start`/`parity_end` use half-open current maternal
parity bands, covering zero through Inf without gaps or overlap in every
group/type stratum.

The
[`stratify_rates()`](https://aarongullickson.github.io/hammel/reference/stratify_rates.md)
function provides a quick means to add stratification to an existing
rate schedule.

## See also

[`set_process()`](https://aarongullickson.github.io/hammel/reference/set_process.md),
[`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
[`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md),
[`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md)

## Examples

``` r
fert_process <- fertility(usa_2000_fertility)
model <- hammel_model() |> set_process(fert_process)
# Equivalent convenience setter:
model <- hammel_model() |> set_fertility(usa_2000_fertility)
```
