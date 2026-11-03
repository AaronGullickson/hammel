# Set the starting population of a Hammel model

Supply a validated population containing person and union histories and
their reference time. Tables are normalized to tibbles; column contents
and additional components are retained. Existing process configuration
is retained. The input model is not modified. Grouped populations
without ancestry proportions initialize each starting person at 100%
ancestry in their recorded group. Supplied ancestry is retained and
normalized to a tibble in person order; see
[`hammel_population()`](https://aarongullickson.github.io/hammel/reference/hammel_population.md).

## Usage

``` r
set_population(model, population)
```

## Arguments

- model:

  A `hammel_model` object.

- population:

  A `hammel_population` object. Construct it from tables with
  [`hammel_population()`](https://aarongullickson.github.io/hammel/reference/hammel_population.md).
  Its `end_time` supplies the model's starting time.

## Value

An updated `hammel_model`.

## Examples

``` r
persons <- data.frame(person_id = 1:2, sex = c("m", "f"),
  birth_time = c(-25, -23), death_time = NA_real_,
  mother_id = NA_integer_, father_id = NA_integer_)
population <- hammel_population(persons, end_time = 0)
model <- hammel_model() |>
  set_population(population)
```
