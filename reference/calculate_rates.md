# Calculate observed demographic event rates

Calculate age-specific mortality, fertility, union formation,
dissolution, cohabitation conversion, and group transition rates from
simulated data by time intervals. Rates can be stratified across a
variety of dimensions.

## Usage

``` r
calculate_mortality_rates(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
)

calculate_fertility_rates(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  parity_breaks = NULL,
  unrealized = FALSE,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
)

calculate_union_rates(
  population,
  time_breaks,
  age_breaks,
  dur_breaks = NULL,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
)

calculate_marriage_rates(
  population,
  time_breaks,
  age_breaks,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE
)

calculate_dissolution_rates(
  population,
  time_breaks,
  dur_breaks,
  by = NULL,
  conf_level = NULL
)

calculate_cohab_convert_rate(
  population,
  time_breaks,
  dur_breaks = NULL,
  age_breaks = NULL,
  by = NULL,
  conf_level = NULL
)

calculate_transition_rates(
  population,
  time_breaks,
  age_breaks = NULL,
  dur_breaks = NULL,
  by = NULL,
  conf_level = NULL,
  use_detailed_previous_union = TRUE,
  group = NULL
)
```

## Arguments

- population:

  A `hammel_population` object returned from
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

- time_breaks:

  Strictly increasing finite calendar boundaries in years or fractions
  of years that define the time intervals of the Lexis surface.

- age_breaks:

  Strictly increasing nonnegative age boundaries in years or fractions
  of years that define the age intervals of the Lexis surface. If the
  last boundary is finite, an open-ended age interval extended to Inf is
  added. Required for mortality, fertility, and marriage. Cohabitation
  conversion and group transitions require exactly one of `age_breaks`
  and `dur_breaks`. Union formation can take both (see details below).

- by:

  Character vector of variables to stratify rates across. Defaults to
  `NULL` for no stratification. Accepted stratification variables differ
  by calculator; see the sections below.

- conf_level:

  NULL (default) to omit confidence intervals, or a finite numeric
  scalar strictly between 0 and 1 indicating the level of confidence to
  use for rate confidence intervals. Adds `conf_low` and `conf_high`
  from two-sided exact
  [`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)
  intervals for observed event counts over exposure. Bounds are NA at
  zero exposure. For fertility, these describe the observed rate, not
  the unrealized-adjusted rate.

- use_detailed_previous_union:

  Logical value that indicates whether stratification by
  `"previous_union"` will use `"separated"` and `"widowed"` or a
  combined `"dissolved"` category.

- parity_breaks:

  For fertility only: NULL to pool parities, or strictly increasing
  nonnegative integer boundaries starting at zero that define categories
  for stratifying by parity.

- unrealized:

  For fertility only: logical, default `FALSE`. When `TRUE`, adds a
  column for `rate_fertility_adjusted` that gives the fertility rate
  when unrealized births (counted births that could not be realized
  because their fertility cell had members but no eligible women).

- dur_breaks:

  Strictly increasing nonnegative duration boundaries in years or
  fractions of years that define the duration intervals of the Lexis
  surface. If the last boundary is finite, an open-ended duration
  interval extended to Inf is added. Required for union dissolution;
  optional for union formation; an alternative to age for cohabitation
  conversion and group transitions.

- group:

  For group transitions only. The declared group dimension whose
  transitions are counted. `NULL` (default) means the only declared
  dimension and will only work when one group is defined.

## Value

A tibble containing calendar bands, age or duration bands where
applicable, requested strata, event counts, exposure, and the observed
rate. Rates are NA at zero exposure; zero events with positive exposure
give zero. Confidence limits are included only when requested.
Person-based exposure uses `person_years` (woman-years for fertility);
union-based exposure uses `union_years`. Each calculator's section
specifies its count/rate columns, and additional output.

## Details

Rates are calculated by counting person-years and events within a given
Lexis triangle from a Lexis surface. Some rates allow for either
`age_breaks` or `dur_breaks` but one must always be provided.
`time_breaks` are also always required but users can specify a single
time interval for the entire simulation.

All breaks defined by the function use half-open intervals that are
inclusive at the lower bound and exclusive at the upper bound.

For more details on each rate, see below.

## Choosing a calculator

|  |  |  |
|----|----|----|
| Calculator | Events | Exposure |
| Mortality | Deaths | Person-years alive |
| Fertility | Maternal births | Woman-years alive |
| Union formation | Partnership starts | Unpartnered person-years |
| Marriage | Marriage starts and conversions | Nonmarried person-years |
| Dissolution | Dissolutions | Union-years |
| Cohabitation conversion | Conversions | Cohabitation union-years |
| Group transitions | Directed group moves | Person-years in origin category |

## Stratifying with `by`

The values `by` accepts in each calculator. Group means the name of a
declared group dimension, such as `"race"` or `"region"`.

Person-level rates count events and exposure for individuals:

|                                |         |        |              |                  |
|--------------------------------|---------|--------|--------------|------------------|
| Function                       | Group   | Sex    | `union_type` | `previous_union` |
| `calculate_mortality_rates()`  | Yes     | Always | Yes          | Yes              |
| `calculate_fertility_rates()`  | Yes     | Women  | Yes          | Yes              |
| `calculate_union_rates()`      | Yes     | Always | Yes          | Yes (a)          |
| `calculate_marriage_rates()`   | Yes     | Always | No           | Yes              |
| `calculate_transition_rates()` | Yes (b) | Yes    | Yes          | Yes              |

Union-level rates count events and exposure for unions:

|                                  |         |              |                  |
|----------------------------------|---------|--------------|------------------|
| Function                         | Group   | `union_type` | `previous_union` |
| `calculate_dissolution_rates()`  | Yes (c) | Yes          | No               |
| `calculate_cohab_convert_rate()` | Yes (c) | No           | No               |

Always: results are always split by sex, so `sex` is not a `by` value.
Women: fertility is reported for women only.

\(a\) Supplying `dur_breaks` adds `previous_union` automatically.

\(b\) Other declared dimensions only, not the one whose transitions are
counted (the `group` argument).

\(c\) The partners' groups: `"<dim>_w"` for the wife's current category,
`"<dim>_h"` for the husband's, or both for the couple's combination,
such as `c("race_w", "race_h")`. A bare dimension name is ambiguous for
a union.

\(d\)

## Mortality

Use `calculate_mortality_rates()`: Reports deaths per person-year lived.

## Fertility

Use `calculate_fertility_rates()`: Reports maternal births per
woman-year, including time excluded by minimum birth spacing. Only
female rows are returned. Current parity is reconstructed from all
recorded maternal births, including pre-window births and deceased
children. Setting `unrealized = TRUE` will calculate a separate
`rate_fertility_adjusted` column that calculates fertility rates as if
unrealized births were realized.

## Union formation

Use `calculate_union_rates()`: Reports union starts (both marriage and
cohabitation) per unpartnered person-year lived. Marriages that occur
because of the conversion of a cohabitation are not counted as events in
this function. To calculate cohabitation and direct marriage separately,
use `by = "union_type`.

The `age_breaks` argument is always required for this function, but
`dur_breaks` can also be used when stratified by `previous_union` to get
duration-based rates for union formation since the previous marriage for
individuals who have been previously married.

## Marriage

Use `calculate_marriage_rates()`: Reports direct marriages and
conversions to marriage from cohabitation over nonmarried person-years.
Unpartnered and cohabiting time both contribute exposure. Married time
contributes none.

## Dissolution

Use `calculate_dissolution_rates()`: Reports union dissolutions over
union-years lived. To calculate dissolution separately by marriage and
cohabitation, use `by = "union_type"`.

Dissolution can only be calculated using `dur_breaks` where duration is
measured since the formation of the current union.

## Cohabitation conversion

Use `calculate_cohab_convert_rate()`: Reports conversions from
cohabitation to marriage over cohabitation union-years.

## Group transitions

use `calculate_transition_rates()`: Reports state transitions from
origin to destination category per person-years lived in origin
category.

When more than one group is in use in a simulation, the `group` argument
must be used to specify which group transition rates are desired.
Transition rates can be stratified by other groups besides the focal
group.

Transition rates can be defined by `age_breaks` or `dur_breaks`. If
`dur_breaks` is used, The duration for the first transition is
equivalent to age.

## See also

[`calculate_population_size()`](https://aarongullickson.github.io/hammel/reference/calculate_population_size.md),
[`calculate_tfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md),
[`calculate_patfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md),
[`calculate_life_table()`](https://aarongullickson.github.io/hammel/reference/calculate_life_table.md)

## Examples

``` r

# run a simulation
pop_start <- generate_uniform_population(1000, end_time = 0)
model <- hammel_model() |>
  set_population(pop_start) |>
  set_fertility(usa_2000_fertility_union) |>
  set_mortality(usa_2000_mortality) |>
  set_union_formation("rate", rate = usa_2000_marriage)
pop_end <- run_simulation(model, 100)

calculate_mortality_rates(
  pop_end,
  time_breaks = seq(0, 100, 20),
  age_breaks = c(0, 1, seq(5, 100, 5))
)
#> # A tibble: 220 × 8
#>    time_start time_end sex   age_start age_end person_years deaths
#>         <dbl>    <dbl> <chr>     <dbl>   <dbl>        <dbl>  <dbl>
#>  1          0       20 m             0       1         129.      2
#>  2         20       40 m             0       1         108.      1
#>  3         40       60 m             0       1         120.      1
#>  4         60       80 m             0       1         121.      2
#>  5         80      100 m             0       1         111.      1
#>  6          0       20 m             1       5         541.      0
#>  7         20       40 m             1       5         425.      0
#>  8         40       60 m             1       5         484.      0
#>  9         60       80 m             1       5         469.      0
#> 10         80      100 m             1       5         459.      0
#> # ℹ 210 more rows
#> # ℹ 1 more variable: rate_mortality <dbl>

calculate_fertility_rates(
  pop_end,
  time_breaks = seq(0, 100, 20),
  age_breaks = seq(15, 40, 5),
  by = "union_type"
)
#> # A tibble: 60 × 9
#>    time_start time_end sex   age_start age_end union_type  person_years births
#>         <dbl>    <dbl> <chr>     <dbl>   <dbl> <chr>              <dbl>  <dbl>
#>  1          0       20 f            15      20 unpartnered         605.     15
#>  2         20       40 f            15      20 unpartnered         517.     17
#>  3         40       60 f            15      20 unpartnered         485.     18
#>  4         60       80 f            15      20 unpartnered         523.     21
#>  5         80      100 f            15      20 unpartnered         512.     16
#>  6          0       20 f            20      25 unpartnered         513.     38
#>  7         20       40 f            20      25 unpartnered         393.     24
#>  8         40       60 f            20      25 unpartnered         355.     26
#>  9         60       80 f            20      25 unpartnered         397.     28
#> 10         80      100 f            20      25 unpartnered         452.     29
#> # ℹ 50 more rows
#> # ℹ 1 more variable: rate_fertility <dbl>
```
