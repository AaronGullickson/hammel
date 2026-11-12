# Extract union data from simulation results with additional detail

Extract the unions data frame from a `hammel_population` object,
augmented with additional information.

## Usage

``` r
augment_unions(population)
```

## Arguments

- population:

  A `hammel_population` object containing person and union histories,
  such as the result of
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

## Value

A separate tibble with all union columns plus `birth_time_w`,
`death_time_w`, `birth_time_h`, and `death_time_h`, plus `<group>_w` and
`<group>_h` for each declared group dimension. Suffixes `_w` and `_h`
identify wives and husbands. Empty union histories retain the enriched
column types.

## Details

The unions record will match wives and husbands to the population's
person table by ID and return their recorded birth times, death times,
and group membership at formation time alongside every union record.

## Examples

``` r
population <- generate_uniform_population(1000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility) |>
  set_union_formation(
    "on_demand",
    matching = age_gap_matching_model(3, 5)
  )
result <- run_simulation(model, years = 50, seed = 42)
result |>
  augment_unions()
#> # A tibble: 347 × 10
#>    union_id wife_id husband_id start_time end_time end_reason    birth_time_w
#>       <dbl>   <dbl>      <dbl>      <dbl>    <dbl> <chr>                <dbl>
#>  1        1     727        118     0.0472     43.9 death_wife           -29.1
#>  2        2     463        275     0.0812     NA   NA                   -26.0
#>  3        3     415        992     0.0871     23.0 death_husband        -36.7
#>  4        4     589        835     0.135      NA   NA                   -14.9
#>  5        5     846        776     0.490      42.6 death_husband        -23.0
#>  6        6     630        201     0.530      37.8 death_husband        -22.3
#>  7        7     880        321     0.555      31.0 death_husband        -20.3
#>  8        8      74        307     0.601      46.7 death_husband        -31.2
#>  9        9     672        177     0.685      48.6 death_husband        -23.8
#> 10       10     857        149     0.723      NA   NA                   -22.2
#> # ℹ 337 more rows
#> # ℹ 3 more variables: death_time_w <dbl>, birth_time_h <dbl>,
#> #   death_time_h <dbl>
```
