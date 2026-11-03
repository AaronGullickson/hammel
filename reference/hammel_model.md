# Create an empty Hammel model specification

Set the starting population and demographic processes using pipeable
setters. Construction and configuration do not execute a simulation.

## Usage

``` r
hammel_model()
```

## Value

A named list of class `hammel_model` containing `population` (initially
`NULL`) and a named list of `processes`. Once set, `population` is a
`hammel_population` object whose `end_time` is the model's starting
time.

## Examples

``` r
persons <- data.frame(person_id = 1:2, sex = c("m", "f"),
  birth_time = c(-25, -23), death_time = NA_real_,
  mother_id = NA_integer_, father_id = NA_integer_)
model <- hammel_model() |>
  set_population(hammel_population(persons, end_time = 0)) |>
  set_mortality(data.frame(sex = c("m", "f"), age_start = 0,
    age_end = Inf, rate = c(0.01, 0.008)))
```
