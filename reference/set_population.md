# Set the starting population of a Hammel model

Set the starting population for a `hammel_model` with a
`hammel_population` object.

## Usage

``` r
set_population(model, population)
```

## Arguments

- model:

  A `hammel_model` object.

- population:

  A `hammel_population` object. Its `end_time` supplies the model's
  starting time.

## Value

An updated `hammel_model`.

## Details

While a `hammel_population` object can be built from scratch with
[`hammel_population()`](https://aarongullickson.github.io/hammel/reference/hammel_population.md),
this object will more commonly be created by
[`generate_uniform_population()`](https://aarongullickson.github.io/hammel/reference/generate_uniform_population.md)
or returned from
[`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

## Examples

``` r
persons <- data.frame(
  person_id = 1:2, sex = c("m", "f"), birth_time = c(-25, -23),
  death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
)
population <- hammel_population(persons, end_time = 0)
model <- hammel_model() |>
  set_population(population)
```
