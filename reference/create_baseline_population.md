# Generate a simple baseline population

Create a starter population with independent uniformly distributed ages,
randomly assigned sex, unknown parents, and no unions. Everyone is alive
at the reference time. This is not a stable or stationary population.

## Usage

``` r
create_baseline_population(
  size,
  max_age = 70,
  end_time = 0,
  female_probability = 0.5,
  group_weights = NULL
)
```

## Arguments

- size:

  Number of persons, a nonnegative whole number.

- max_age:

  Maximum initial age in years, finite and nonnegative.

- end_time:

  Reference time in years, finite. Birth times are this value minus the
  sampled ages. Defaults to a relative time origin of zero.

- female_probability:

  Probability of assigning sex "f" (female). The other sex value is "m"
  (male). This does not enforce exact proportions.

- group_weights:

  Optional named numeric vector of nonnegative, finite weights with at
  least one positive value. Unique, nonempty names give group labels.
  Weights are normalized and each person's group is sampled
  independently of age and sex. Realized group counts are random. `NULL`
  (the default) adds no group column and consumes no extra random draws.

## Value

A validated `hammel_population` object.

## Details

Ages are sampled in continuous years without rounding to simulation
steps. Generation uses and advances R's current random-number state; use
[`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility.

## Examples

``` r
set.seed(42)
population <- create_baseline_population(size = 50000, max_age = 70)
model <- hammel_model() |>
  set_population(population)
```
