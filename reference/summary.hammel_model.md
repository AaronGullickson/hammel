# Summarize a Hammel model

Summary method for class `"hammel_model"`

## Usage

``` r
# S3 method for class 'hammel_model'
summary(object, ...)

# S3 method for class 'summary.hammel_model'
print(x, ...)
```

## Arguments

- object:

  A `hammel_model` object returned by
  [`hammel_model()`](https://aarongullickson.github.io/hammel/reference/hammel_model.md).

- ...:

  Additional arguments passed to
  [`summary.hammel_population()`](https://aarongullickson.github.io/hammel/reference/summary.hammel_population.md).

- x:

  A `summary.hammel_model` object returned by `summary.hammel_model()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_model` containing:

- population:

  The `summary.hammel_population` of the starting population, or `NULL`
  when no population is set.

- groups:

  A named list with one element per group dimension, or `NULL` when the
  model has none. Each element contains `group` (the
  `summary.hammel_group` of the declared dimension, reporting the
  population's categories when a population is set, or `NULL` when the
  dimension is not declared) and `transition` (the
  `summary.hammel_group_transition` of the dimension's transition
  process, or `NULL` when it has none).

- processes:

  A named list of process summaries, in the order mortality, fertility,
  union formation, union dissolution, for the processes present.

- missing:

  A character vector naming the required elements that are not set: any
  of `"population"`, `"mortality"`, `"fertility"`, and
  `"union_formation"`, plus `"transition:<dimension>"` for each group
  transition process whose dimension is not declared
  (`"transition:.default"` for a process without `group` when the model
  does not declare exactly one dimension). Empty when the model has the
  elements needed to run a simulation.

[`print()`](https://rdrr.io/r/base/print.html) displays each component
summary in turn, followed by whether the model has the elements needed
to run a simulation and, if not, which are missing. It returns its input
invisibly.

## Details

The summary combines the summaries of the starting population, of each
group dimension and its transition process, and of each other process in
the model, and records which elements required by
[`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)
are missing.

## See also

[`hammel_model()`](https://aarongullickson.github.io/hammel/reference/hammel_model.md),
[`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)

## Examples

``` r
model <- hammel_model() |>
  set_population(generate_uniform_population(1000)) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility)
summary(model)
#> 
#> ── Hammel model ────────────────────────────────────────────────────────────────
#> 
#> ── Population ──
#> 
#> • First date: -69.948
#> • End time: 0
#> • Persons alive at end time: 1000
#> • Total births: 1000
#> • Total deaths: 0
#> • Total unions: 0
#> • Total dissolutions: 0
#> 
#> ── Living population by age and sex 
#>                      Male     Age     Female                   
#>         ################# |  65-69  | #########################
#>         ################# |  60-64  | ##################       
#>         ################# |  55-59  | #########################
#>         ################# |  50-54  | ##################       
#>          ################ |  45-49  | ####################     
#>     ##################### |  40-44  | ######################   
#>   ####################### |  35-39  | #####################    
#>         ################# |  30-34  | ################         
#>     ##################### |  25-29  | ####################     
#>     ##################### |  20-24  | #####################    
#>          ################ |  15-19  | ##################       
#>     ##################### |  10-14  | #####################    
#>         ################# |   5-9   | ##################       
#>       ################### |   0-4   | ##################       
#> 
#> ── Mortality process ──
#> 
#> • Stratified by: sex (f, m)
#> • Age bands (years): [0, 1), [1, 5), [5, 10), ..., [95, 100), [100, Inf) (22
#>   bands; showing 5)
#> • Weight model: none
#> 
#> ── Fertility process ──
#> 
#> • Minimum birth interval: 0.75 years
#> • Fraction female at birth: 0.4886
#> • Stratified by: none
#> • Age bands (years): [15, 20), [20, 25), [25, 30), [30, 35), [35, 40), [40, 45)
#> • Weight model: none
#> 
#> ✖  The model needs the following before it can run a simulation:
#> • A union formation process, from `set_union_formation()`

model <- set_union_formation(model, "on_demand")
summary(model)
#> 
#> ── Hammel model ────────────────────────────────────────────────────────────────
#> 
#> ── Population ──
#> 
#> • First date: -69.948
#> • End time: 0
#> • Persons alive at end time: 1000
#> • Total births: 1000
#> • Total deaths: 0
#> • Total unions: 0
#> • Total dissolutions: 0
#> 
#> ── Living population by age and sex 
#>                      Male     Age     Female                   
#>         ################# |  65-69  | #########################
#>         ################# |  60-64  | ##################       
#>         ################# |  55-59  | #########################
#>         ################# |  50-54  | ##################       
#>          ################ |  45-49  | ####################     
#>     ##################### |  40-44  | ######################   
#>   ####################### |  35-39  | #####################    
#>         ################# |  30-34  | ################         
#>     ##################### |  25-29  | ####################     
#>     ##################### |  20-24  | #####################    
#>          ################ |  15-19  | ##################       
#>     ##################### |  10-14  | #####################    
#>         ################# |   5-9   | ##################       
#>       ################### |   0-4   | ##################       
#> 
#> ── Mortality process ──
#> 
#> • Stratified by: sex (f, m)
#> • Age bands (years): [0, 1), [1, 5), [5, 10), ..., [95, 100), [100, Inf) (22
#>   bands; showing 5)
#> • Weight model: none
#> 
#> ── Fertility process ──
#> 
#> • Minimum birth interval: 0.75 years
#> • Fraction female at birth: 0.4886
#> • Stratified by: none
#> • Age bands (years): [15, 20), [20, 25), [25, 30), [30, 35), [35, 40), [40, 45)
#> • Weight model: none
#> 
#> ── Union formation process ──
#> 
#> • Search method: on_demand
#> • Matching model: ~age_difference + I(age_difference^2)
#>   • Coefficients: age_difference = 0, I(age_difference^2) = -0.02773
#>   • Preferred age gap: 0 years
#>   • Half-weight distance: 5 years
#> • Births without a matched father: cancelled
#> • Male eligibility restrictions
#>   • Minimum male age: 15
#>   • Generational proximity: 2
#> • Cohabitation: disabled
#> 
#> ✔  The model has the elements needed to run a simulation.
```
