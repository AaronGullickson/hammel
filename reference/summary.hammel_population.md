# Summarize a Hammel population

Summary method for class `"hammel_population"`

## Usage

``` r
# S3 method for class 'hammel_population'
summary(object, interval_pyramid = 5, ...)

# S3 method for class 'summary.hammel_population'
print(x, ...)
```

## Arguments

- object:

  an object of class `"hammel_population"`, usually the result of a call
  to
  [run_simulation](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

- interval_pyramid:

  length of age intervals used for drawing population pyramid

- ...:

  further arguments passed on from other methods.

- x:

  A `summary.hammel_population` object returned by
  `summary.hammel_population()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_population` containing:

- first_time:

  The earliest recorded birth time, in years.

- end_time:

  The population's reference time, in years.

- persons_alive:

  The number of persons alive at the reference time.

- n_births:

  The number of recorded persons with a birth time, including persons in
  the starting population.

- n_births_u:

  The number of retained unrealized birth records, or `NA` when no
  unrealized birth table is present.

- n_deaths:

  The number of recorded persons with a death time.

- n_unions:

  The number of recorded union episodes, including ended unions.

- n_dissolutions:

  The number of union episodes with end reason `"dissolution"`.

- age_counts:

  A tibble of living male and female counts by five-year age group, with
  columns `age`, `male`, and `female`.

- simulation:

  A list describing the latest simulation run, with `start` and `end`
  times in years and living population counts `pop_size_start` and
  `pop_size_end`, or `NULL` when no simulation metadata is present.

- group_counts:

  A named list of tibbles, one per declared group dimension, or `NULL`
  when no groups are declared. Each tibble contains `category`, `count`,
  and `percent` for the living population at the reference time,
  including declared categories with zero members. Percentages are on a
  0–100 scale and are `NaN` when the dimension's total living count is
  zero.

[`print()`](https://rdrr.io/r/base/print.html) displays the summary,
group counts when present, and a text population pyramid, returning its
input invisibly. Pyramid bars share a common scale across both sexes.

## Examples

``` r

# get summary of starting population
population <- generate_uniform_population(1000)
summary(population)
#> 
#> ── Population ──
#> 
#> • First date: -69.885
#> • End time: 0
#> • Persons alive at end time: 1000
#> • Total births: 1000
#> • Total deaths: 0
#> • Total unions: 0
#> • Total dissolutions: 0
#> 
#> ── Living population by age and sex 
#>                      Male     Age     Female                   
#>       ################### |  65-69  | ################         
#>       ################### |  60-64  | ###########              
#>          ################ |  55-59  | #############            
#>        ################## |  50-54  | ###################      
#>              ############ |  45-49  | ################         
#>         ################# |  40-44  | ##############           
#>        ################## |  35-39  | ############             
#>         ################# |  30-34  | #########################
#>        ################## |  25-29  | ####################     
#>           ############### |  20-24  | ###############          
#>      #################### |  15-19  | ###############          
#>       ################### |  10-14  | ############             
#>        ################## |   5-9   | #################        
#>         ################# |   0-4   | ###################      

# get summary of simulatied population
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility) |>
  set_union_formation("on_demand")
result <- run_simulation(model, years = 50, seed = 42)
summary(result)
#> 
#> ── Population ──
#> 
#> • First date: -69.885
#> • End time: 50
#> • Persons alive at end time: 1094
#> • Total births: 1689
#> • Total unrealized births: 0
#> • Total deaths: 595
#> • Total unions: 352
#> • Total dissolutions: 0
#> • Latest simulation run
#>   • Time period: 0-50
#>   • Starting population size: 1000
#>   • Ending population size: 1094
#> 
#> ── Living population by age and sex 
#>                      Male     Age     Female                   
#>                           | 100-104 | #                        
#>                           |  95-99  | ##                       
#>                        ## |  90-94  | ##                       
#>                      #### |  85-89  | ######                   
#>                    ###### |  80-84  | ############             
#>                ########## |  75-79  | ############             
#>                ########## |  70-74  | ##########               
#>        ################## |  65-69  | ##############           
#>        ################## |  60-64  | #############            
#>          ################ |  55-59  | ################         
#>          ################ |  50-54  | ####################     
#>         ################# |  45-49  | #################        
#>            ############## |  40-44  | ##################       
#>           ############### |  35-39  | #############            
#>        ################## |  30-34  | ###################      
#>        ################## |  25-29  | ##############           
#>          ################ |  20-24  | ##############           
#>       ################### |  15-19  | ################         
#> ######################### |  10-14  | ######################   
#>            ############## |   5-9   | ##############           
#>      #################### |   0-4   | ################         
```
