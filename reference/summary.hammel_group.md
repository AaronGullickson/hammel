# Summarize a group dimension

Summary method for class `"hammel_group"`

## Usage

``` r
# S3 method for class 'hammel_group'
summary(object, categories = NULL, ...)

# S3 method for class 'summary.hammel_group'
print(x, ...)
```

## Arguments

- object:

  A `hammel_group` object returned by
  [`group()`](https://aarongullickson.github.io/hammel/reference/group.md).

- categories:

  Optional categories to report instead of the declared ones, such as
  the categories a population resolves for the dimension. `NULL`
  (default) uses the declared categories.

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_group` object returned by `summary.hammel_group()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_group` containing:

- categories:

  The dimension's categories, or `NULL` when they are inferred from the
  population.

- ancestry:

  Whether ancestry shares are tracked.

- inheritance:

  A summary of the newborn inheritance model: a
  `summary.hammel_group_inheritance`, or a list of `periods` and
  `models` for a calendar schedule.

## See also

[`group()`](https://aarongullickson.github.io/hammel/reference/group.md),
[`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md)

## Examples

``` r
ethnicity <- group(
  categories = c("W", "B", "O"),
  ancestry = TRUE,
  inheritance = group_inheritance("ancestry")
)
summary(ethnicity)
#> 
#> ── Group dimension ──
#> 
#> • Categories: W, B, O
#> • Ancestry tracked: yes
#> • Inheritance: ancestry
```
