# Add proportional strata to a rate table

Stratify a rate schedule by additional columns, applying proportional
changes to rates by stratum.

## Usage

``` r
stratify_rates(rates, by)
```

## Arguments

- rates:

  A data frame with a numeric `rate` column, such as a mortality,
  fertility, search, conversion, dissolution or transition rate table.

- by:

  A named list. Each name is a column to add; each element is a named
  numeric vector of finite, nonnegative factors, one per category. The
  final rate will be determined by taking the original rate times the
  factor provided.

## Value

A tibble with the new columns (character labels) and modified rates.

## Details

Any of the input rate processes such as fertility, mortality, marriage,
transition, dissolution, and conversion can be stratified by this
function. Stratified rates can also be further stratified by this
function.

Columns in `by` are either built-in strata (`sex`, `union_type`,
`previous_union`, `origin`, `destination`) or group dimensions declared
with
[`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md):
a dimension's name (such as `race`) for person-level rates, or
`<name>_w` and `<name>_h` (the wife's and husband's categories) for
dissolution and conversion rates. Band, rate and calendar columns cannot
be added. Process setters reject built-in strata their process does not
support.

If `by` includes multiple columns, every combination of categories gets
the product of its factors.

## Examples

``` r
# stratify fertility rates by group and partnership status
usa_2000_fertility |> calculate_tfr()
#> # A tibble: 1 × 3
#>   age_start age_end   tfr
#>       <int>   <int> <dbl>
#> 1        15      45  2.05
fertility_rates <- stratify_rates(
  usa_2000_fertility,
  by = list(
    group = c(A = 1.1, B = 1),
    # this is not really a great idea for union type as within union
    # type fertility rates have a totally different age profile, but
    # ok for testing purposes
    union_type = c("unpartnered" = 0.3, "partnered" = 1.3)
  )
)
fertility_rates |> calculate_tfr()
#> # A tibble: 4 × 5
#>   group union_type  age_start age_end   tfr
#>   <chr> <chr>           <int>   <int> <dbl>
#> 1 A     unpartnered        15      45 0.676
#> 2 A     partnered          15      45 2.93 
#> 3 B     unpartnered        15      45 0.615
#> 4 B     partnered          15      45 2.66 

# although this works, a proportional hazards assumption by union type
# is unrealistic because the within union type age profile of fertility is
# very different from the overall age profile. So for this case, its better
# to use already stratified rates by union type and then further stratify
# by group
fertility_rates <- stratify_rates(
  usa_2000_fertility_union,
  by = list(group = c(A = 1.1, B = 1))
)
fertility_rates |> calculate_tfr()
#> # A tibble: 4 × 5
#>   union_type  group age_start age_end   tfr
#>   <chr>       <chr>     <int>   <int> <dbl>
#> 1 partnered   A            15      45  4.61
#> 2 unpartnered A            15      45  1.29
#> 3 partnered   B            15      45  4.20
#> 4 unpartnered B            15      45  1.17
```
