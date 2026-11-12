# Specify newborn group identity inheritance

A model specification indicating how the group category of newborns
should be assigned for a particular group dimension. Used in
[`group()`](https://aarongullickson.github.io/hammel/reference/group.md).

## Usage

``` r
group_inheritance(
  method = "random_parent",
  fun = NULL,
  order = NULL,
  passing = 0,
  fixed_value = NULL
)
```

## Arguments

- method:

  One of `"mother"`, `"father"`, `"random_parent"`, `"same_sex"`,
  `"opposite_sex"`, `"ancestry"`, `"fixed"`, `"hierarchy"`, or
  `"custom"`. See details for further information.

- fun:

  For `"custom"`, an R function defined by the user. See details for
  variables exposed for this function.

- order:

  For `"hierarchy"`, unique group labels, highest-ranked first. Chooses
  between parental group membership, not ancestry.

- passing:

  For `"hierarchy"`, probability in `[0,1]` of reversing the order for a
  birth. Defaults to zero.

- fixed_value:

  For `"fixed"`, one category label or finite numeric code in the
  dimension's categories. Every newborn receives this category.

## Value

An inspectable `hammel_group_inheritance` object.

## Details

Whenever groups are in use, an inheritance rule is required to determine
how the group membership of newborns is determined. The function
provides several out of the box methods, as well as the ability for
users to define a custom method:

- `random_parent`: The default method. The newborn is assigned the group
  membership of a randomly selected parent with equal probability.

- `mother` - The newborn is assigned the group membership of the mother.

- `father` - The newborn is assigned the group membership of the father.

- `same_sex` - The newborn is assigned the group membership of their
  same-sex parent.

- `different_sex` - The newborn is assigned the group membership of
  their opposite-sex parent.

- `ancestry` - The newborn is assigned a group membership randomly,
  weighted by their ancestry share from each group. `ancestry = TRUE`
  must be set for the group to use this method.

- `fixed` - All newborns are assigned a fixed group membership given by
  the `fixed_value` argument.

- `hierarchy` - The newborn is assigned the group membership from the
  parent whose membership is higher on the `order` argument. If the
  `passing` argument is provided, it defines the probability that this
  ordering is inverted.

- `custom` - The newborn is assigned a group membership based on a
  custom R function provided by the user with the `fun` argument.
  Simulations with custom R functions are slower to process.

For methods that require no additional arguments (`random_parent`,
`mother`, `father`, `same_sex`, `different_sex`, `ancestry`), the name
of the method may be passed to arguments requiring an inheritance model
like
[`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
instead of the model itself.

For the custom function, the following variables are exposed and can be
called from the function declaration:

- `mother_group`: The mother's group membership.

- `father_group`: The father's group membership. NA if unknown.

- `mother_ancestry`: Numeric vector of ancestry share from each group
  for mother, as well as an "unknown" value for unknown ancestry.

- `father_ancestry`: Numeric vector of ancestry share from each group
  for father, as well as an "unknown" value for unknown ancestry.

- `child_ancestry`: Numeric vector of ancestry share from each group for
  child, as well as an "unknown" value for unknown ancestry.

- `child_sex`: Either "f" or "m".

- `parents_in_union`: Logical indicating whether the parents were in a
  union, including on-demand union formation.

## Examples

``` r
group_inheritance("fixed", fixed_value = "none")
#> $method
#> [1] "fixed"
#> 
#> $fun
#> NULL
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> $fixed_value
#> [1] "none"
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
group_inheritance("mother")
#> $method
#> [1] "mother"
#> 
#> $fun
#> NULL
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> $fixed_value
#> NULL
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
group_inheritance("hierarchy", order = c("B", "W"), passing = 0.05)
#> $method
#> [1] "hierarchy"
#> 
#> $fun
#> NULL
#> 
#> $order
#> [1] "B" "W"
#> 
#> $passing
#> [1] 0.05
#> 
#> $fixed_value
#> NULL
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
group_inheritance("custom",
  fun = function(mother_group, father_group, ...) {
    if (is.na(father_group) || mother_group == father_group) {
      mother_group
    } else {
      "mixed"
    }
  }
)
#> $method
#> [1] "custom"
#> 
#> $fun
#> function (mother_group, father_group, ...) 
#> {
#>     if (is.na(father_group) || mother_group == father_group) {
#>         mother_group
#>     }
#>     else {
#>         "mixed"
#>     }
#> }
#> <environment: 0x5bafa677fc40>
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> $fixed_value
#> NULL
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
```
