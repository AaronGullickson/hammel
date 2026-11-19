# Age-specific fertility in Niger, 2012

Single-year age-specific fertility estimates for Niger in 2012.

## Usage

``` r
niger_fertility
```

## Format

A tibble with 35 rows and three columns:

- age_start, age_end:

  Half-open single-year age bands from `[15, 16)` through `[49, 50)`.

- rate:

  Births per woman-year.

## Source

United Nations, World Population Prospects 2024, fertility rates by
single age. Source metadata and the rebuild script are under
`data-raw/niger_2012/` in the full package repository.

## Examples

``` r
fertility(niger_fertility)
#> $name
#> [1] "fertility"
#> 
#> $rates
#> # A tibble: 35 × 3
#>    age_start age_end   rate
#>        <int>   <int>  <dbl>
#>  1        15      16 0.0803
#>  2        16      17 0.140 
#>  3        17      18 0.200 
#>  4        18      19 0.253 
#>  5        19      20 0.288 
#>  6        20      21 0.303 
#>  7        21      22 0.311 
#>  8        22      23 0.317 
#>  9        23      24 0.321 
#> 10        24      25 0.325 
#> # ℹ 25 more rows
#> 
#> $birth_interval
#> [1] 0.75
#> 
#> $ffab
#> [1] 0.4886
#> 
#> attr(,"class")
#> [1] "hammel_fertility" "hammel_process"  
```
