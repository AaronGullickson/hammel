# Add recorded parent information to person histories

Match mothers and fathers by person ID and append their ages at each
child's birth and their optional recorded groups. All persons, including
founders and deceased persons, retain their original row order and
columns. The input population is not modified and no known endpoint is
required.

## Usage

``` r
get_parent_info(population)
```

## Arguments

- population:

  A validated `hammel_population`, such as the result of
  [`simulate()`](https://aarongullickson.github.io/hammel/reference/simulate.md).

## Value

A separate tibble containing all person columns plus `age_mother` and
`age_father`, measured at that person's birth. When groups are present,
also `group_mother` and `group_father`. Empty populations retain the
enriched column types.

## Details

Ages are in continuous years: the child's birth time minus the parent's
birth time. Unknown parents yield missing ages and groups. Group fields
are omitted when persons has no `group` column; otherwise its type is
preserved. These are recorded identities, not reconstructions of
historical group transitions. Existing columns that conflict with the
added fields cause an error.

## Examples

``` r
population <- create_baseline_population(10)
get_parent_info(population)
#> # A tibble: 10 × 8
#>    person_id sex   birth_time death_time mother_id father_id age_mother
#>        <int> <chr>      <dbl>      <dbl>     <int>     <int>      <dbl>
#>  1         1 f          -41.1         NA        NA        NA         NA
#>  2         2 m          -61.6         NA        NA        NA         NA
#>  3         3 m          -15.1         NA        NA        NA         NA
#>  4         4 m          -24.4         NA        NA        NA         NA
#>  5         5 f          -65.2         NA        NA        NA         NA
#>  6         6 m          -61.2         NA        NA        NA         NA
#>  7         7 f          -69.6         NA        NA        NA         NA
#>  8         8 f          -13.7         NA        NA        NA         NA
#>  9         9 f          -17.0         NA        NA        NA         NA
#> 10        10 m          -57.5         NA        NA        NA         NA
#> # ℹ 1 more variable: age_father <dbl>
```
