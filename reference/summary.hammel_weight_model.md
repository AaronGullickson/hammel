# Summarize a weight model

Summary method for class `"hammel_weight_model"`

## Usage

``` r
# S3 method for class 'hammel_weight_model'
summary(object, ...)

# S3 method for class 'summary.hammel_weight_model'
print(x, ...)
```

## Arguments

- object:

  A `hammel_weight_model` object returned by
  [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_weight_model` object returned by
  `summary.hammel_weight_model()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_weight_model` containing:

- custom:

  Whether the model uses a custom function rather than a formula.

- formula:

  The weight formula, or `NULL` for a custom function.

- coefficients:

  The named log hazard ratios of a formula model.

[`print()`](https://rdrr.io/r/base/print.html) also shows each
coefficient's hazard ratio, `exp(coefficient)`.

## See also

[`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)

## Examples

``` r
partnered_mortality <- weight_model(~partnered, c(partnered = log(0.9)))
summary(partnered_mortality)
#> 
#> ── Weight model ──
#> 
#> • Formula: ~partnered
#> • Coefficients: partnered = -0.1054 (hazard ratio 0.9)
```
