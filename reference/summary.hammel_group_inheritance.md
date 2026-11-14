# Summarize a group inheritance model

Summary method for class `"hammel_group_inheritance"`

## Usage

``` r
# S3 method for class 'hammel_group_inheritance'
summary(object, ...)

# S3 method for class 'summary.hammel_group_inheritance'
print(x, ...)
```

## Arguments

- object:

  A `hammel_group_inheritance` object returned by
  [`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md).

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_group_inheritance` object returned by
  `summary.hammel_group_inheritance()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_group_inheritance` containing:

- method:

  The inheritance method.

- order:

  For `"hierarchy"`, the group labels, highest-ranked first; otherwise
  `NULL`.

- passing:

  For `"hierarchy"`, the probability of reversing the order; otherwise
  `NULL`.

- fixed_value:

  For `"fixed"`, the category every newborn receives; otherwise `NULL`.

## See also

[`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md),
[`group()`](https://aarongullickson.github.io/hammel/reference/group.md)

## Examples

``` r
summary(group_inheritance("hierarchy", order = c("B", "W"), passing = 0.05))
#> 
#> ── Group inheritance model ──
#> 
#> • Method: hierarchy
#> • Order: B > W
#> • Passing probability: 0.05
```
