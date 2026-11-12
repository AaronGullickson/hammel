# Set a demographic process specification

Store a process specification under its name, replacing any existing
process with the same name.

## Usage

``` r
set_process(model, process)
```

## Arguments

- model:

  A
  [hammel_model](https://aarongullickson.github.io/hammel/reference/hammel_model.md)
  object.

- process:

  A `hammel_process` specification.

## Value

An updated
[hammel_model](https://aarongullickson.github.io/hammel/reference/hammel_model.md).
The input model is not modified.

## Details

The available processes are:

- [`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md)

- [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)

- [`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md)

- [`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md)

- [`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)

Use this function to attach an already-created process. The convenience
setters
[`set_mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
[`set_fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
[`set_union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md),
[`set_union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md),
and
[`set_group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)
construct and attach the corresponding process in one call.

## See also

[`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md),
[`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md),
[`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md),
[`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md),
[`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)

## Examples

``` r
mortality_process <- mortality(usa_2000_mortality)
model <- hammel_model() |> set_process(mortality_process)
model$processes
#> $mortality
#> $name
#> [1] "mortality"
#> 
#> $rates
#> # A tibble: 44 × 4
#>    sex   age_start age_end     rate
#>    <chr>     <int>   <dbl>    <dbl>
#>  1 f             0       1 0.00628 
#>  2 f             1       5 0.000288
#>  3 f             5      10 0.000142
#>  4 f            10      15 0.000159
#>  5 f            15      20 0.000396
#>  6 f            20      25 0.000471
#>  7 f            25      30 0.000536
#>  8 f            30      35 0.000745
#>  9 f            35      40 0.00114 
#> 10 f            40      45 0.00174 
#> # ℹ 34 more rows
#> 
#> attr(,"class")
#> [1] "hammel_mortality" "hammel_process"  
#> 
```
