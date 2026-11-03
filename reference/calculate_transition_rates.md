# Calculate age- or duration-specific group-transition rates

Count individual moves and person-years in their origin group in
half-open Lexis cells. Each outgoing destination shares its origin's
exposure: do not sum denominators across destinations. Duration is years
since the latest group entry, or age before the first recorded
transition; return moves reset duration and union changes do not. Events
use the exiting group's duration. Death and competing moves end
exposure. Rates are cause-specific hazards, not transition
probabilities. Internal-boundary events enter the cell starting there;
the final calendar boundary is excluded.

## Usage

``` r
calculate_transition_rates(
  population,
  time_breaks,
  age_breaks = NULL,
  dur_breaks = NULL,
  by_sex = FALSE,
  by_union_type = FALSE,
  conf.level = NULL
)
```

## Arguments

- population:

  A grouped `hammel_population` with a known endpoint and complete
  person, transition and (when requested) union histories.

- time_breaks:

  Strictly increasing finite calendar boundaries in years.

- age_breaks:

  Optional strictly increasing nonnegative age boundaries.

- dur_breaks:

  Optional strictly increasing nonnegative group-duration boundaries.
  Supply exactly one of age_breaks and dur_breaks. A finite final axis
  boundary gains Inf.

- by_sex:

  Logical, default FALSE. Disaggregate by `sex` (m/f).

- by_union_type:

  Logical, default FALSE. Disaggregate exposure and events by current
  union type. Typed histories use unpartnered/cohab/marriage; untyped
  histories use unpartnered/partnered. Requires exclusive histories.

- conf.level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1. Adds `conf_low` and `conf_high` from
  two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

## Value

A tibble with `time_start`, `time_end`, age or duration bounds,
`group_origin`, `group_dest`, optional `sex` and `union_type`,
`person_years`, `transitions`, and `rate_transition` (NA at zero
exposure). All distinct directed pairs among current and historical
identities are retained, including zero-event cells. Fewer than two
known groups returns a typed empty tibble. Period varies fastest within
axis bands and pair blocks. With non-NULL `conf.level`, also includes
`conf_low` and `conf_high`.

## Examples

``` r
persons <- data.frame(person_id = 1, sex = "f", birth_time = 0,
  death_time = NA_real_, mother_id = NA_real_, father_id = NA_real_, group = "B")
population <- hammel_population(persons, end_time = 5,
  transitions = data.frame(person_id = 1, time = 2,
    group_origin = "A", group_dest = "B"))
calculate_transition_rates(population, c(0, 5), age_breaks = c(0, Inf))
#> # A tibble: 2 × 9
#>   time_start time_end age_start age_end group_origin group_dest person_years
#>        <dbl>    <dbl>     <dbl>   <dbl> <chr>        <chr>             <dbl>
#> 1          0        5         0     Inf B            A                     3
#> 2          0        5         0     Inf A            B                     2
#> # ℹ 2 more variables: transitions <dbl>, rate_transition <dbl>
```
