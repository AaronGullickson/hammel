# Generate a population with uniformly distributed ages

Create a starter population with independent uniformly distributed ages,
randomly assigned sex, unknown parents, and no unions.

## Usage

``` r
generate_uniform_population(
  size,
  max_age = 70,
  end_time = 0,
  min_age = 0,
  female_probability = 0.5,
  groups = NULL
)
```

## Arguments

- size:

  Number of living persons, a nonnegative whole number.

- max_age:

  Maximum initial age in years, finite and nonnegative.

- end_time:

  Reference time in years, finite. Birth times are this value minus the
  sampled ages. Defaults to a relative time origin of zero.

- min_age:

  Minimum initial age in years, finite and nonnegative, and no greater
  than `max_age`. Defaults to zero.

- female_probability:

  Probability of assigning sex "f" (female). The other sex value is "m"
  (male).

- groups:

  Optional group dimensions to sample, independently of age and sex.
  Either a named list with one named numeric vector of weights per
  dimension (categories as names; dimensions are drawn independently),
  or a data frame with one column per dimension and a numeric `weight`
  column, giving the joint distribution of category combinations.
  Weights are finite and nonnegative with at least one positive value,
  and are normalized. `NULL` (the default) adds no groups.

## Value

A validated `hammel_population` object.

## Details

Ages are sample from the uniform distribution between `min_age` and
`max_age` and only the living population is generated. This population
will have no father or mother ids. Sex and group categories follow their
supplied probabilities.

The generated population will not have a realistic age distribution.
Typically, users will run a "burn in" simulation to produce a more
realistic population distribution and kinship history for the desired
simulation.

The `end_time` for the population will define the starting point for the
simulation. By default, this value is set to zero so that time is
measured since the beginning of the simulation. This means that birth
times for the existing population will be negative. Users who are
drawing from input data with time increments from real calendar time,
should reset this starting value as appropriate.

If the `groups` argument is used, the population declares each dimension
with
[`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
defaults, except that its `categories` are those given, in order:

- the weights' names

- each dimension's values in order of first appearance in the table's
  rows.

The first is the reference category, and zero-weight categories are
declared, allowing for groups that may emerge later.
[`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md)
can replace that declaration; categories it omits are kept.

## Examples

``` r
set.seed(42)

# generate 50000 people with a maximum age of 50
population <- generate_uniform_population(size = 50000, max_age = 50)

# generate a population divided by race and region, with unequal division
# of race by region
region_race_dist <- tibble::tribble(
  ~race, ~region, ~weight,
 "W",   "NE",    0.145,
 "W",   "MW",    0.190,
 "W",   "S",     0.257,
 "W",   "W",     0.152,
 "B",   "NE",    0.022,
 "B",   "MW",    0.023,
 "B",   "S",     0.067,
 "B",   "W",     0.019,
 "O",   "NE",    0.021,
 "O",   "MW",    0.014,
 "O",   "S",     0.030,
 "O",   "W",     0.060,
)
population <- generate_uniform_population(50000, groups = region_race_dist)
prop.table(table(population$persons$race, population$persons$region))
#>    
#>          MW      NE       S       W
#>   B 0.02316 0.02238 0.06700 0.01818
#>   O 0.01408 0.02154 0.03044 0.05972
#>   W 0.18980 0.14254 0.25886 0.15230
```
