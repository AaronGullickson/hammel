# Extract persons data from simulation results with additional detail

Extract the persons data frame from a `hammel_population` object,
augmented with additional information.

## Usage

``` r
augment_persons(population)
```

## Arguments

- population:

  A `hammel_population`, such as the result of
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

## Value

A separate tibble containing all person columns plus `age_mother` and
`age_father`, measured at that person's birth, and for each declared
group dimension `<group>_mother` and `<group>_father` (such as
`race_mother`). Empty populations retain the enriched column types.

## Details

Mothers and fathers will be matched with children by person ID and their
age at the person's birth will be appended along with group membership.

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
  augment_persons() |>
  tail()
#> # A tibble: 6 × 8
#>   person_id sex   birth_time death_time mother_id father_id age_mother
#>       <dbl> <chr>      <dbl>      <dbl>     <dbl>     <dbl>      <dbl>
#> 1      1788 m           49.2         NA      1300      1352       30.5
#> 2      1789 m           49.5         NA      1335      1223       28.7
#> 3      1790 m           49.8         NA      1291      1277       32.2
#> 4      1791 f           49.9         NA      1510      1485       18.7
#> 5      1792 f           49.9         NA      1395      1368       25.7
#> 6      1793 m           50.0         NA      1360      1348       27.3
#> # ℹ 1 more variable: age_father <dbl>
```
