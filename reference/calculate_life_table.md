# Construct period life tables from age-specific mortality rates

Constructs a period life table from input tables that provide
age-specific mortality rates. Input data can come from user-provided
input rate schedules, mortality process objects, or estimated rates from
a simulation via
[`calculate_mortality_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
and may be stratified by additional variables.

## Usage

``` r
calculate_life_table(rates, radix = 1)
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

- radix:

  Initial number of life-table survivors. Defaults to 1.

## Value

A tibble containing standard life table parameters as well as `sex` and
additional stratifying columns (when supplied).

## Details

When the supplied `rates` table contains columns other than `age_start`,
`age_end`, and either `rate_mortality` or `rate`, the function will
treat these additional columns as stratifying dimensions and will
stratify the output appropriately. It will also stratify output by time
when `time_start` and `time_end` are included. The function also ignores
additional calculated columns from the output of
[`calculate_mortality_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md),
so the results of that function can be directly piped into this
function. Stratification produces separate life tables by strata, not
multistate life tables.

Constructed life tables do not make any adjustment to the \\a_x\\ or
\\q_x\\ values for the 0 to 1 age group. commonly used in demographic
methods. This is to ensure that results are consistent with input rates
used by the simulation. Users who want to model changes in mortality
within the first year of life can provide fractional year age bands for
the first year of life.

Missing values on rates will result in missing values for lx, Lx, Tx,
and ex for any later age bands. A zero mortality rate in the open age
band will result in infinite life expectancy.

## See also

[`calculate_life_expectancy()`](https://aarongullickson.github.io/hammel/reference/calculate_life_expectancy.md)
for life expectancy at selected ages.

## Examples

``` r
data(usa_2000_mortality)

calculate_life_table(usa_2000_mortality)
#> # A tibble: 44 × 13
#>    sex   age_start age_end       mx     n    ax       qx    lx lx_end       dx
#>    <chr>     <int>   <dbl>    <dbl> <dbl> <dbl>    <dbl> <dbl>  <dbl>    <dbl>
#>  1 f             0       1 0.00628      1 0.499 0.00626  1      0.994 0.00626 
#>  2 f             1       5 0.000288     4 2.00  0.00115  0.994  0.993 0.00115 
#>  3 f             5      10 0.000142     5 2.50  0.000712 0.993  0.992 0.000706
#>  4 f            10      15 0.000159     5 2.50  0.000794 0.992  0.991 0.000787
#>  5 f            15      20 0.000396     5 2.50  0.00198  0.991  0.989 0.00196 
#>  6 f            20      25 0.000471     5 2.50  0.00235  0.989  0.987 0.00233 
#>  7 f            25      30 0.000536     5 2.50  0.00268  0.987  0.984 0.00264 
#>  8 f            30      35 0.000745     5 2.50  0.00372  0.984  0.981 0.00366 
#>  9 f            35      40 0.00114      5 2.50  0.00570  0.981  0.975 0.00559 
#> 10 f            40      45 0.00174      5 2.50  0.00866  0.975  0.966 0.00845 
#> # ℹ 34 more rows
#> # ℹ 3 more variables: Lx <dbl>, Tx <dbl>, ex <dbl>

# stratify by groups to get stratified life tables
stratify_rates(
  usa_2000_mortality,
  by = list(ethnicity = c(A = 1.2, B = 0.9))
) |>
  calculate_life_table()
#> # A tibble: 88 × 14
#>    sex   ethnicity age_start age_end       mx     n    ax       qx    lx lx_end
#>    <chr> <chr>         <int>   <dbl>    <dbl> <dbl> <dbl>    <dbl> <dbl>  <dbl>
#>  1 f     A                 0       1 0.00754      1 0.499 0.00751  1      0.992
#>  2 f     A                 1       5 0.000346     4 2.00  0.00138  0.992  0.991
#>  3 f     A                 5      10 0.000171     5 2.50  0.000854 0.991  0.990
#>  4 f     A                10      15 0.000191     5 2.50  0.000952 0.990  0.989
#>  5 f     A                15      20 0.000475     5 2.50  0.00237  0.989  0.987
#>  6 f     A                20      25 0.000566     5 2.50  0.00282  0.987  0.984
#>  7 f     A                25      30 0.000643     5 2.50  0.00321  0.984  0.981
#>  8 f     A                30      35 0.000894     5 2.50  0.00446  0.981  0.977
#>  9 f     A                35      40 0.00137      5 2.50  0.00684  0.977  0.970
#> 10 f     A                40      45 0.00209      5 2.50  0.0104   0.970  0.960
#> # ℹ 78 more rows
#> # ℹ 4 more variables: dx <dbl>, Lx <dbl>, Tx <dbl>, ex <dbl>
```
