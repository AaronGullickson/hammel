# Set a demographic process specification

Store a process specification under its name, replacing an existing
process with the same name. Other processes are retained. Built-in rate
schedules are validated again when attached. No simulation is executed.

## Usage

``` r
set_process(model, process)

set_mortality(model, rates, group_prop_hazard = NULL)

set_fertility(
  model,
  rates,
  birth_interval = 0.75,
  group_prop_hazard = NULL,
  inheritance = create_group_inheritance(),
  ffab = 0.4886,
  hetfert = NULL
)
```

## Arguments

- model:

  A `hammel_model` object.

- process:

  A `hammel_process` specification.

- rates:

  An age-specific rate table passed to
  [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
  or
  [`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md).

- group_prop_hazard:

  Optional named numeric vector of finite nonnegative rate multipliers
  by group. Missing groups use 1.0; values are not normalized. Scales
  shared or group-specific fertility input rates before cell accounting.
  Without population groups, simulation warns and ignores it. Fertility
  also accepts a named list of vectors for `partnered` and/or
  `unpartnered`, or cohab/marriage when rates include those `union_type`
  strata. List names must occur in the schedules; omitted entries use
  1.0.

- birth_interval:

  Minimum maternal birth-to-birth interval in years, defaulting to
  `0.75`. Must be finite and nonnegative.

- inheritance:

  A group identity model from
  [`create_group_inheritance()`](https://aarongullickson.github.io/hammel/reference/create_group_inheritance.md),
  or a tibble with numeric `time_start`, `time_end` and list-column
  `model` containing complete inheritance models. Intervals are
  half-open absolute years, must not overlap and must cover the
  requested run. Select the active model at each birth; older identities
  and ancestry remain unchanged. Base data frames with list columns are
  accepted. All declared possible identities across scheduled models
  join the simulation-wide registry. Static inputs may also be a method
  name requiring no additional parameters: `"mother"`, `"father"`,
  `"random_parent"`, `"same_sex"`, `"opposite_sex"`, or `"ancestry"`.
  Applies to fertility and defaults to selecting either parent with
  equal probability, with maternal fallback. Assignment is independent
  of recorded ancestry proportions.

- ffab:

  Probability that a newborn is female, default `0.4886`. Must be a
  finite numeric scalar in `[0, 1]`. Applies to fertility only.

- hetfert:

  NULL (default) for uniform mother selection, or a
  [`het_fert()`](https://aarongullickson.github.io/hammel/reference/het_fert.md)
  specification for persistent lognormal selection weights.

## Value

An updated `hammel_model`. The input model is not modified.

## Details

`set_mortality()` and `set_fertility()` construct the corresponding
process and delegate to `set_process()`.

## Examples

``` r
mortality_rates <- data.frame(sex = c("m", "f"), age_start = 0,
  age_end = Inf, rate = c(0.01, 0.008))
fertility_rates <- data.frame(age_start = 15, age_end = 50, rate = 0.1)
model <- hammel_model() |>
  set_mortality(mortality_rates) |>
  set_fertility(fertility_rates)
model$processes
#> $mortality
#> $name
#> [1] "mortality"
#> 
#> $rates
#> # A tibble: 2 × 4
#>   sex   age_start age_end  rate
#>   <chr>     <dbl>   <dbl> <dbl>
#> 1 m             0     Inf 0.01 
#> 2 f             0     Inf 0.008
#> 
#> $rate_unit
#> [1] "person_year"
#> 
#> attr(,"class")
#> [1] "hammel_mortality" "hammel_process"  
#> 
#> $fertility
#> $name
#> [1] "fertility"
#> 
#> $rates
#> # A tibble: 1 × 3
#>   age_start age_end  rate
#>       <dbl>   <dbl> <dbl>
#> 1        15      50   0.1
#> 
#> $rate_unit
#> [1] "person_year"
#> 
#> $inheritance
#> $method
#> [1] "random_parent"
#> 
#> $fun
#> NULL
#> 
#> $groups
#> NULL
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
#> 
#> $birth_interval
#> [1] 0.75
#> 
#> $ffab
#> [1] 0.4886
#> 
#> attr(,"class")
#> [1] "hammel_fertility" "hammel_process"  
#> 
```
