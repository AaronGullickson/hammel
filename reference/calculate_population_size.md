# Calculate the living population at time points

Count persons with `birth_time <= time` and missing `death_time` or
`death_time > time`. Counts include births at the requested time and
exclude deaths at that time, representing the population after recorded
events there. Results preserve the order of the requested times,
including repeated times.

## Usage

``` r
calculate_population_size(population, times, by_group = FALSE)
```

## Arguments

- population:

  A `hammel_population` or persons data frame representing one
  independent history. For bare tables, ensure follow-up is complete
  through every requested time.

- times:

  A numeric vector of finite time points in years. A population object
  cannot be queried beyond its known `end_time`.

- by_group:

  Logical, default `FALSE`. With `TRUE`, require a persons `group`
  column and count separately by group. Groups retain their input type
  and first-occurrence order within each requested time. No pooled group
  row is added. Empty histories return a typed zero-row table.

## Value

A tibble with `time` and `population` columns, plus `group` when
`by_group = TRUE`.

## Examples

``` r
population <- create_baseline_population(100, end_time = 0)
calculate_population_size(population, times = c(-10, -5, 0))
#> # A tibble: 3 × 2
#>    time population
#>   <dbl>      <int>
#> 1   -10         90
#> 2    -5         95
#> 3     0        100
```
