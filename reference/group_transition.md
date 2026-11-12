# Specify individual group transitions

Constructs a group-specific transition process specification for a
simulation model. The convenience function `set_group_tansition` is
equivalent to `set_process(model, group_transition(...))` and allows for
direct piping from a model.

## Usage

``` r
group_transition(rates, group = NULL, weights = NULL, household = NULL)

set_group_transition(
  model,
  rates,
  group = NULL,
  weights = NULL,
  household = NULL
)
```

## Arguments

- rates:

  A nonempty data frame defining the mortality rate schedules. Must
  contain at least `origin`, `destination`, `rate`, and either
  `age_start`/`age_end` or `dur_start`/`dur_end` throughout. See details
  for format and description of stratification variables that can be
  applied.

- group:

  The name of the declared group dimension (see
  [`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md))
  whose categories are affected by this transition rate schedule. `NULL`
  (default) means the model's only dimension, and is an error when
  several are declared.

- weights:

  NULL (default) for unweighted selection, or a
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  using person variables and `destination`, which chooses who in a cell
  moves without changing the cell's input rate.

- household:

  `NULL` (default): each person moves individually. A
  [`household()`](https://aarongullickson.github.io/hammel/reference/household.md):
  a move carries the mover's household with them on a transition. If the
  household is set, then the rates apply to the household collectively
  or individuals if they are not in households. Can only be used with
  age bands.

- model:

  A `hammel_model` object, used by `set_group_transition()`.

## Value

`group_transition()` returns a process with a rates tibble;
`set_group_transition()` returns the updated model.

## Details

The primary and only required input to the function is a transition rate
schedule for the `rates` argument. This data frame must contain either
`age_start` and `age_end` columns or `dur_start` and `dur_end` columns
to define the half-open age/duration intervals (inclusive at the start
and exclusive at the end), an `origin` and `destination` column to
define the origin and destination categories for the transition and a
`rate` column for the transition rates. Transition rate are defined as
the number of transitions from the origin to destination category per
person-year spent in the origin category.

For a set of K categories, there will be K(K-1) possible transition
rates at every age/duration band. Any cases left unspecified will be
assumed to have a zero rate. The rates are not assumed to be symmetric
(e.g.A \> B is not the same as B \< A), but users can provide schedules
that make these rates

Each group dimension has its own transition process: and users must call
`set_group_transition()` separately for each group dimension for which
they desire a process. Columns named by declared group dimensions
besides the focal group supply separate schedules for each combination
of their categories; dimensions without a column share the same
schedule.

Optional `time_start`/`time_end` columns select half-open calendar
intervals in absolute years, including fractional years. Both columns
are required together; starts are finite and ends exceed starts
(possibly Inf). Active slices must satisfy the strata and coverage rules
and cover the requested simulation run.

Transition rates can also be stratified by `union_type` containing
either `"partnered"` and `"unpartnered"`, or all three of
`"unpartnered"`, `"cohab"`, `"marriage"`, or by `sex` ("m" and "f").

## Examples

``` r
# create grouped population
pop_start <- generate_uniform_population(
  1000,
  groups = list(ethnicity = c(A = 1, B = 1))
)

# create transition matrix with more A > B than B > A
ethnic_trans <- tibble::tribble(
 ~age_start, ~age_end, ~origin, ~destination, ~rate,
 0,          Inf,      "A",     "B",          0.05,
 0,          Inf,      "B",     "A",          0.02
)

# create simulation model and run it
model <- hammel_model() |>
  set_population(pop_start) |>
  set_fertility(usa_2000_fertility) |>
  set_mortality(usa_2000_mortality) |>
  set_union_formation("on_demand") |>
  set_groups(ethnicity = group(ancestry = TRUE)) |>
  set_group_transition(ethnic_trans)
pop_end <- run_simulation(model, 100)

# check population sizes - group A will grow relative to group B
calculate_population_size(pop_end, seq(0, 100, 20), by = "ethnicity")
#> # A tibble: 12 × 3
#>     time ethnicity population
#>    <dbl> <chr>          <int>
#>  1     0 B                497
#>  2    20 B                709
#>  3    40 B                767
#>  4    60 B                805
#>  5    80 B                825
#>  6   100 B                797
#>  7     0 A                503
#>  8    20 A                389
#>  9    40 A                351
#> 10    60 A                356
#> 11    80 A                331
#> 12   100 A                343
```
