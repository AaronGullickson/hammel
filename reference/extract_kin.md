# Extract known kin within generation-distance bounds

Extract kin of a set of ego individuals from simulated data. Kin
relations are governed by generational distance.

## Usage

``` r
extract_kin(population, ego_ids, generations = 2, relative_ids = NULL)
```

## Arguments

- population:

  A `hammel_population` with a finite `end_time`.

- ego_ids:

  Unique person IDs whose known kin should be extracted.

- generations:

  Nonnegative whole number of parent generations on EACH side, default
  `2`. Zero returns no kin between distinct people.

- relative_ids:

  Optional unique relative IDs to consider. `NULL` means all recorded
  persons. Unknown IDs are errors; empty vectors return no rows.

## Value

A tibble with one row per ego and relative: `ego_id`, `relative_id`,
`relationship`, `ego_distance`, `relative_distance`,
`n_common_ancestors`, `ego_lineage`, `relative_lineage`,
`ego_ancestry_complete`, and `relative_ancestry_complete`. Lineage is
`"self"`, `"maternal"`, `"paternal"`, or `"both"` for shortest routes to
the closest common ancestors.

## Details

Find recorded common ancestors within `generations` parent steps of both
an ego and a relative, including each person at distance zero. This
function uses the same native kinship module as union eligibility.
Unknown parents terminate their branches. At depth two, known
parents/children, grandparents/grandchildren, siblings, aunts/uncles,
nieces/nephews, and first cousins are included, including half
relationships. Great-grandparents and first cousins once removed require
depth three.

Returns one row per relative, describing their closest common ancestors:
those with the smallest sum of the two distances within the requested
bounds. `n_common_ancestors` counts them, ordinarily 2 for full kin and
1 for half kin; a count of 1 is only recorded evidence when a parent is
unknown. `relationship` names the distance pair, such as `"sibling"`,
`"aunt_uncle"` or `"first_cousin"`, with full and half kin together. The
ego itself is excluded. Includes deceased and living people in recorded
history; use
[`count_kin()`](https://aarongullickson.github.io/hammel/reference/count_kin.md)
for mean counts of living kin at given times.

## See also

[`count_kin()`](https://aarongullickson.github.io/hammel/reference/count_kin.md)

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
extract_kin(result, max(result$persons$person_id))
#> # A tibble: 16 × 10
#>    ego_id relative_id relationship ego_distance relative_distance
#>     <dbl>       <dbl> <chr>               <int>             <int>
#>  1   1688         490 grandparent             2                 0
#>  2   1688         502 grandparent             2                 0
#>  3   1688         604 grandparent             2                 0
#>  4   1688         710 grandparent             2                 0
#>  5   1688        1142 aunt_uncle              2                 1
#>  6   1688        1183 aunt_uncle              2                 1
#>  7   1688        1302 parent                  1                 0
#>  8   1688        1325 aunt_uncle              2                 1
#>  9   1688        1335 parent                  1                 0
#> 10   1688        1396 aunt_uncle              2                 1
#> 11   1688        1479 first_cousin            2                 2
#> 12   1688        1564 first_cousin            2                 2
#> 13   1688        1573 first_cousin            2                 2
#> 14   1688        1609 first_cousin            2                 2
#> 15   1688        1660 first_cousin            2                 2
#> 16   1688        1682 first_cousin            2                 2
#> # ℹ 5 more variables: n_common_ancestors <int>, ego_lineage <chr>,
#> #   relative_lineage <chr>, ego_ancestry_complete <lgl>,
#> #   relative_ancestry_complete <lgl>
```
