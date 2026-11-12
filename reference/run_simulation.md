# Run a demographic simulation

Run a demographic simulation with the provided `hammel_model` for a
given time length.

## Usage

``` r
run_simulation(model, years, seed = NULL, progress = FALSE)
```

## Arguments

- model:

  A `hammel_model` containing a population, mortality, fertility, and
  union-formation processes, optionally union dissolution and group
  transition processes. The population must have a finite `end_time`.

- years:

  Simulation duration in years; fractional years are supported.

- seed:

  Optional integer seed for random number generation. If not specified,
  simulation will use current RNG.

- progress:

  Show a console progress bar based on elapsed simulation years.
  Defaults to `FALSE`.

## Value

A `hammel_population` with updated histories and endpoint, plus a
`simulation` component recording selected settings.

## Details

A bare minimum model for simulation will contain a population and
fertility, mortality, and union formation processes, but may contain
additional processes.

The output of the function is a `hammel_population` which can be set as
the new starting population of the model via
[`set_population()`](https://aarongullickson.github.io/hammel/reference/set_population.md)
to chain simulations together. However, users should note that frequent
restarts are slower than a single simulation with changing rates and
models.

If no seed is provided, the simulation will use the current RNG. Setting
the seed will ensure reproducibility of results and simulation will not
affect existing RNG state.

## Examples

``` r
population <- generate_uniform_population(1000)
model <- hammel_model() |>
  set_population(population) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility) |>
  set_union_formation(
    "on_demand",
    matching = age_gap_matching_model(3, 5)
  )
result <- run_simulation(model, years = 50, seed = 42)
```
