# Extract births with known parents to unpartnered mothers

Return births where the mother was outside both cohabitation and
marriage, with parent ages and optional groups supplied by
[`augment_persons()`](https://aarongullickson.github.io/hammel/reference/augment_persons.md).

## Usage

``` r
extract_unpartnered_births(population)
```

## Arguments

- population:

  A `hammel_population` with person and union histories, such as the
  result of
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

## Value

A separate tibble containing selected person rows and `age_mother`,
`age_father`, plus `group_mother` and `group_father` when groups are
present. Input order, column types and custom person columns are
preserved. Empty results retain the enriched schema. The population is
not modified.

## Examples

``` r
population <- generate_uniform_population(1000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility_union) |>
  set_union_formation(
    "rate",
    rates = usa_2000_marriage,
    matching = age_gap_matching_model(3, 5)
  )
result <- run_simulation(model, years = 50, seed = 42)
extract_unpartnered_births(result)
#> # A tibble: 277 × 8
#>    person_id sex   birth_time death_time mother_id father_id age_mother
#>        <dbl> <chr>      <dbl>      <dbl>     <dbl>     <dbl>      <dbl>
#>  1      1001 m         0.0324         NA       557       456       26.8
#>  2      1002 m         0.197          NA        20       315       20.2
#>  3      1003 m         0.357          NA       225       691       24.3
#>  4      1005 m         0.418          NA       745       830       29.8
#>  5      1006 f         0.458          NA       884       162       39.6
#>  6      1007 m         0.468          NA       877       182       32.5
#>  7      1008 f         0.658          NA       414       286       28.4
#>  8      1009 m         0.847          NA       397       686       27.2
#>  9      1010 m         0.944          NA       563       833       24.6
#> 10      1011 f         0.955          NA        59        90       23.8
#> # ℹ 267 more rows
#> # ℹ 1 more variable: age_father <dbl>
```
