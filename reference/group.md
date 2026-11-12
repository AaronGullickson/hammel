# Describe and declare a group dimension

The `group()` function describes a single group dimension, including
valid categories, whether to track ancestry, and inheritance rules. The
`set_groups()` function declares all group dimensions for a
`hammel_model`.

For the `set_groups` function, each named argument should identify a
group name and provide the returned value from a `group()` function. See
examples for syntax.

## Usage

``` r
group(categories = NULL, ancestry = FALSE, inheritance = group_inheritance())

set_groups(model, ...)
```

## Arguments

- categories:

  Optional complete, ordered vector of the dimension's categories
  (character labels or finite numeric codes). `NULL` (the default)
  either keeps the categories the population already declares, if any,
  or they are inferred from the population's values in first-occurrence
  order, plus transition destinations and categories with ancestry
  shares. When given, every value in the population and every transition
  destination must be listed; categories can be declared even if
  starting population contains no members of that category.

- ancestry:

  Logical, default `FALSE`. Track ancestry shares for this dimension in
  the simulation output. The returned population then has an `ancestry`
  table indicating the share of ancestry from each category of the
  group. Required for use of the `ancestry` inheritance method.

- inheritance:

  Either a single character string, a
  [`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md)
  object, or a calendar schedule tibble with `time_start`, `time_end`
  and a list-column `model` providing the group inheritance model for a
  time period. Determines how group membership is assigned to newborns.
  Defaults to `random_parent`. Character sting entries will use built-in
  methods provided gy the `group_inheritance` model.

- model:

  A `hammel_model`, used by `set_groups()`.

- ...:

  Named `group()` objects, one per dimension, named by the persons
  column. No arguments removes the declaration.

## Value

`group()` returns a `hammel_group` object describing one dimension.
`set_groups()` returns an updated `hammel_model`.

## See also

[`hammel_population()`](https://aarongullickson.github.io/hammel/reference/hammel_population.md),
[`generate_uniform_population()`](https://aarongullickson.github.io/hammel/reference/generate_uniform_population.md),
[`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md),
[`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)

## Examples

``` r
# create two groups - ethnicity and education
model <- hammel_model() |>
  set_groups(
    ethnicity = group(
      categories = c("W", "B", "O"),
      ancestry = TRUE,
      inheritance = group_inheritance("ancestry")
    ),
    education = group(
      categories = c("LHS", "HS", "SC", "C"),
      ancestry = FALSE,
      inheritance = group_inheritance("fixed", fixed_value = "LHS")
    )
  )
```
