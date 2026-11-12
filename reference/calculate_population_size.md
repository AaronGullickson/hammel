# Calculate the living population at time points

Calculate population counts from a returned simulation result or a
starting population, at chosen time points. Counts can be disaggregated
by stratifying variables.

## Usage

``` r
calculate_population_size(
  population,
  times,
  by = NULL,
  use_detailed_previous_union = TRUE
)
```

## Arguments

- population:

  A `hammel_population` object either returned from
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)
  or built as a baseline starting population.

- times:

  A numeric vector of finite time points in years at which to calculate
  population size.

- by:

  Character vector of variables to disaggregate population counts by.
  Defaults to `NULL` for no disaggregation. Populations counts can be
  disaggregated by `"sex"`, `"union_type"`, `"previous_union"` and any
  declared group dimensions.

- use_detailed_previous_union:

  Logical value that indicates whether disaggregation by
  `"previous_union"` will use `"separated"` and `"widowed"` or a
  combined `"dissolved"` category.

## Value

A tibble with `time`, the `by` columns and `population`.

## Details

The `by` argument can take the standard variables of `sex`,
`union_type`, or `previous_union` as well as any group variables defined
in the `population` input.

The user can input any numeric value for `times` (including decimal
values), but no value may come after the `end_time` of the simulation.

## Examples

``` r
# calculate population size on a generatd population
pop_start <- generate_uniform_population(1000, end_time = 0)
calculate_population_size(pop_start, times = c(-10, -5, 0))
#> # A tibble: 3 × 2
#>    time population
#>   <dbl>      <int>
#> 1   -10        861
#> 2    -5        934
#> 3     0       1000

# run a simulation
model <- hammel_model() |>
  set_population(pop_start) |>
  set_fertility(usa_2000_fertility) |>
  set_mortality(usa_2000_mortality) |>
  set_union_formation("on_demand")
pop_end <- run_simulation(model, 100)

# calculate at 20 year intervals
calculate_population_size(pop_end, times = seq(0, 100, 20))
#> # A tibble: 6 × 2
#>    time population
#>   <dbl>      <int>
#> 1     0       1000
#> 2    20       1100
#> 3    40       1146
#> 4    60       1165
#> 5    80       1245
#> 6   100       1313

# caclulate at 20 year intervals by sex
calculate_population_size(pop_end, times = seq(0, 100, 20), by = "sex")
#> # A tibble: 12 × 3
#>     time sex   population
#>    <dbl> <chr>      <int>
#>  1     0 m            491
#>  2    20 m            511
#>  3    40 m            524
#>  4    60 m            524
#>  5    80 m            580
#>  6   100 m            612
#>  7     0 f            509
#>  8    20 f            589
#>  9    40 f            622
#> 10    60 f            641
#> 11    80 f            665
#> 12   100 f            701
```
