# Construct period life tables from mortality rates

Treat mortality rates as piecewise constant hazards. Closed-interval
survival is exponential, and person-years are integrated under that same
assumption. This does not apply conventional HMD infant adjustments.
Missing rates remain unknown. A zero terminal rate gives infinite
remaining lifetime when terminal survivorship is positive. Require a
final explicit open age interval. Tables are constructed separately for
each period, sex, group, and union status present. Status-specific
tables hold that status fixed; they do not model transitions between
statuses.

## Usage

``` r
calculate_life_table(rates, radix = 1)
```

## Arguments

- rates:

  A table returned by
  [`calculate_mortality_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_mortality_rates.md),
  a mortality process object, or a data frame with `age_start`,
  `age_end`, and `rate_mortality`, `mx`, or `rate` per person-year.

- radix:

  Initial life-table survivors, a positive finite number.

## Value

A tibble with period, sex, group, and union-status identifiers (when
supplied), `age_start`, `age_end`, `mx`, `n`, `ax`, `qx`, `lx`,
`lx_end`, `dx`, `Lx`, `Tx`, and `ex`. Summary counts and exposure
columns are omitted. Life expectancy is `Tx / lx` in years.

## Examples

``` r
calculate_life_table(data.frame(age_start = c(0, 1, 5),
  age_end = c(1, 5, Inf), mx = 0.1))
#> # A tibble: 3 × 12
#>   age_start age_end    mx     n     ax     qx    lx lx_end     dx    Lx    Tx
#>       <dbl>   <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1         0       1   0.1     1  0.492 0.0952 1      0.905 0.0952 0.952 10   
#> 2         1       5   0.1     4  1.87  0.330  0.905  0.607 0.298  2.98   9.05
#> 3         5     Inf   0.1   Inf 10     1      0.607  0     0.607  6.07   6.07
#> # ℹ 1 more variable: ex <dbl>
```
