# Calculate the average number of living kin over time

Count, for each person living at each requested time, the living
relatives recorded within `generations` parent steps on each side, by
relationship, and average them over living persons.

## Usage

``` r
count_kin(population, times, generations = 2, by = NULL, relationships = NULL)
```

## Arguments

- population:

  A `hammel_population` with a known `end_time` no earlier than the
  requested times.

- times:

  A numeric vector of finite time points in years.

- generations:

  Nonnegative whole number of parent generations on each side, default
  `2`, which covers parents, children, grandparents, grandchildren,
  siblings, aunts and uncles, nieces and nephews, and first cousins.

- by:

  `NULL` (default) or a character vector of the ego's characteristics at
  each time to count separately: `"sex"` and declared group dimensions,
  such as `"race"`.

- relationships:

  `NULL` (default) for all relationships, or a nonempty character vector
  of unique relationship names within `generations`, such as
  `"first_cousin"`. Retains the usual row order and ego denominators.

## Value

A tibble with `time`, the `by` columns, `relationship`, `persons`
(living egos), `kin` (their living kin of that type), `mean_kin`
(`kin / persons`) and `share_with_kin` (share of egos with at least
one). `mean_kin` and `share_with_kin` are `NA` with no egos. Every time,
stratum and selected relationship within `generations` has a row. Rows
vary by time fastest, then the `by` variables in order, then
relationship.

## Details

Counts use recorded parent links only. A starting population without
recorded parents has no recorded kin, so counts are too low until two or
three generations have been simulated. Allow a burn-in before
interpreting them. Relatives share kin, so egos are not independent.

The following named relationships are available:

- "parent"

- "child"

- "grandparent"

- "grandchild"

- "sibling"

- "aunt_uncle"

- "niece_nephew"

- "first_cousin"

- "great_grandparent" (requires `generations = 3`)

- "great_grandchild" (requires `generations = 3`)

- "great_aunt_uncle" (requires `generations = 3`)

- "great_niece_nephew" (requires `generations = 3`)

- "first_cousin_once_removed" (requires `generations = 3`)

- "second_cousin" (requires `generations = 3`)

## See also

[`extract_kin()`](https://aarongullickson.github.io/hammel/reference/extract_kin.md)

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
result <- run_simulation(model, years = 100, seed = 42)
count_kin(result, 100, relationships = "first_cousin")
#> # A tibble: 1 × 6
#>    time relationship persons   kin mean_kin share_with_kin
#>   <dbl> <chr>          <int> <int>    <dbl>          <dbl>
#> 1   100 first_cousin    1263  7282     5.77          0.839
```
