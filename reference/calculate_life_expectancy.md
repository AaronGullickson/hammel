# Calculate life expectancy at given ages

A convenience wrapper around
[`calculate_life_table()`](https://aarongullickson.github.io/hammel/reference/calculate_life_table.md)
that returns only remaining life expectancy (`ex`) at the requested
ages, for each stratum of the input rates.

## Usage

``` r
calculate_life_expectancy(rates, age = 0)
```

## Arguments

- rates:

  A table returned by
  [`calculate_mortality_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
  a mortality process object from
  [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
  or a data frame with at least `age_start`, `age_end`, and
  `rate_mortality` or `rate` columns. Must contain an open-ended age
  group specified by an `Inf` value for `age_end`. Input may be further
  stratified by other variables.

- age:

  Ages at which to report life expectancy, in years. Each must be the
  start of an age band in `rates`. Defaults to `0`, life expectancy at
  birth.

## Value

A tibble with the stratifying columns of the input (when supplied),
`age`, and `ex`, with one row per stratum and requested age.

## See also

[`calculate_life_table()`](https://aarongullickson.github.io/hammel/reference/calculate_life_table.md)

## Examples

``` r
calculate_life_expectancy(usa_2000_mortality)
#> # A tibble: 2 × 3
#>   sex     age    ex
#>   <chr> <int> <dbl>
#> 1 f         0  79.4
#> 2 m         0  74.0
calculate_life_expectancy(usa_2000_mortality, age = c(0, 30, 65))
#> # A tibble: 6 × 3
#>   sex     age    ex
#>   <chr> <int> <dbl>
#> 1 f         0  79.4
#> 2 f        30  50.5
#> 3 f        65  19.0
#> 4 m         0  74.0
#> 5 m        30  45.8
#> 6 m        65  16.0
```
