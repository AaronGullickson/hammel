# Create an empty Hammel model specification

This function creates the empty `hammel_model` object which can then be
built out by adding populations and processes.

## Usage

``` r
hammel_model()
```

## Value

A named list of class `hammel_model` containing `population` (initially
`NULL`) and a named list of `processes`.

## Details

The `hammel_model` is a container for a single `population` and a named
list of `processes`. These can all be added by the various setters
through a pipable interface or by feeding the model into
[`set_population()`](https://aarongullickson.github.io/hammel/reference/set_population.md)
or
[`set_process()`](https://aarongullickson.github.io/hammel/reference/set_process.md).

## Examples

``` r
model <- hammel_model() |>
  set_population(generate_uniform_population(1000)) |>
  set_mortality(usa_2000_mortality) |>
  set_fertility(usa_2000_fertility)
```
