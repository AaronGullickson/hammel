# Extract births with known parents to unpartnered mothers

Return births where the mother was outside both cohabitation and
marriage, with parent ages and optional groups supplied by
[`get_parent_info()`](https://aarongullickson.github.io/hammel/reference/get_parent_info.md).
Known mother and father IDs are required. Each birth is retained
separately, including repeated parental pairs; missing-parent records
and founders are omitted. These rows can be used to examine
nonpartnership father matching.

## Usage

``` r
get_unpartnered_births(population)
```

## Arguments

- population:

  A validated `hammel_population` with person and union histories, such
  as the result of
  [`simulate()`](https://aarongullickson.github.io/hammel/reference/simulate.md).

## Value

A separate tibble containing selected person rows and `age_mother`,
`age_father`, plus `group_mother` and `group_father` when groups are
present. Input order, column types and custom person columns are
preserved. Empty results retain the enriched schema. The population is
not modified.

## Details

Histories must be complete for the intended interpretation. Partnership
coverage is half-open: births at formation are partnered. Births exactly
at dissolution or conversion retain the ending partnership (birth-first
ties). Births at death endings are unpartnered unless another episode
covers that time. Overlapping episodes do not duplicate births. Untyped
partnerships are treated as partnered too. No time window is imposed and
no endpoint is needed.

## Examples

``` r
population <- create_baseline_population(10)
get_unpartnered_births(population)
#> # A tibble: 0 × 8
#> # ℹ 8 variables: person_id <int>, sex <chr>, birth_time <dbl>,
#> #   death_time <dbl>, mother_id <int>, father_id <int>, age_mother <dbl>,
#> #   age_father <dbl>
```
