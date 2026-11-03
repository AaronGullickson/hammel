# Add recorded partner information to union histories

Match wives and husbands to the population's person table by ID and
return their recorded birth times, death times, and optional groups at
formation alongside every union record. Union row order and all existing
columns are retained. The population and its union table are not
modified. No known endpoint is needed.

## Usage

``` r
get_detailed_unions(population)
```

## Arguments

- population:

  A validated `hammel_population` containing person and union histories,
  such as the result of
  [`simulate()`](https://aarongullickson.github.io/hammel/reference/simulate.md).

## Value

A separate tibble with all union columns plus `birth_time_w`,
`death_time_w`, `birth_time_h`, and `death_time_h`. When group is
present, also `group_w` and `group_h`. Suffixes `_w` and `_h` identify
wives and husbands. Empty union histories retain the enriched column
types.

## Details

Dates remain in years. Missing death times denote people alive at the
recorded endpoint. Group columns are omitted when persons has no `group`
column; otherwise their types are preserved. Groups are reconstructed at
each union's start from retained transitions; later moves do not relabel
earlier unions. A transition at the exact formation time follows
formation, so its origin group applies. Without transition history, use
recorded person groups. Existing columns whose names conflict with added
partner fields cause an error.

## Examples

``` r
population <- create_baseline_population(10)
get_detailed_unions(population)
#> # A tibble: 0 × 10
#> # ℹ 10 variables: union_id <int>, wife_id <int>, husband_id <int>,
#> #   start_time <dbl>, end_time <dbl>, end_reason <chr>, birth_time_w <dbl>,
#> #   death_time_w <dbl>, birth_time_h <dbl>, death_time_h <dbl>
```
