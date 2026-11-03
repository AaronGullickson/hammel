# Specify newborn group identity inheritance

Assign one group independently of genealogical ancestry shares. Attach
this object with
[`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
or
[`set_fertility()`](https://aarongullickson.github.io/hammel/reference/set_process.md).
The default selects either parent with equal probability, with maternal
fallback. Existing people are not reassigned, and ancestry is unchanged
by assignment. Fertility also accepts a schedule tibble with
`time_start`, `time_end`, and a `model` list column containing complete
inheritance objects. Half-open absolute-year intervals select the model
at birth; all scheduled possible identities enter one shared registry
before execution. See
[`fertility()`](https://aarongullickson.github.io/hammel/reference/mortality.md).

## Usage

``` r
create_group_inheritance(
  method = "random_parent",
  fun = NULL,
  groups = NULL,
  order = NULL,
  passing = 0
)
```

## Arguments

- method:

  One of `"mother"`, `"father"`, `"random_parent"`, `"same_sex"`,
  `"opposite_sex"`, `"ancestry"`, `"hierarchy"`, or `"custom"`. Random
  parent chooses each with probability 0.5. Same/opposite sex chooses
  the parent whose sex matches/differs from the child's. Ancestry
  samples known child shares, conditioning out unknown ancestry. Missing
  information required by a built-in rule falls back to the mother,
  including missing fathers, entirely unknown ancestry, and missing
  hierarchy ranks.

- fun:

  For `"custom"`, an R function receiving named arguments
  `mother_group`, `father_group`, `mother_ancestry`, `father_ancestry`,
  `child_ancestry`, `child_sex` (`"f"`/`"m"`), and `parents_in_union`
  (logical union status at birth, including on-demand formation).
  Nonmarital births and unknown fathers give FALSE. Older callbacks
  without this new argument or `...` remain supported. Use `...` for
  unused inputs. Ancestry vectors have `group_<label>`/`unknown` names.
  A missing father has NA group and NULL ancestry. Return exactly one
  nonmissing identity in the possible group set. Ordinary R random draws
  are supported; functions must not mutate simulation/global state.
  Invalid returns and errors fail.

- groups:

  Optional unique group labels extending the inferred possible
  identities (starting groups plus ancestry groups with positive
  shares). Declare any additional custom outputs here. Grouped rates and
  categorical matching levels must cover the complete set before
  simulation.

- order:

  For `"hierarchy"`, unique group labels, highest-ranked first. Choose
  between parental identities, not ancestry. Complex orders use `fun`.

- passing:

  For `"hierarchy"`, probability in `[0,1]` of reversing the order for a
  birth. Defaults to zero. Equal parental identities and missing
  ranks/parents retain maternal assignment without a passing draw.

## Value

An inspectable `hammel_group_inheritance` object.

## Examples

``` r
create_group_inheritance("mother")
#> $method
#> [1] "mother"
#> 
#> $fun
#> NULL
#> 
#> $groups
#> NULL
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
create_group_inheritance("hierarchy", order = c("A", "B"), passing = 0.05)
#> $method
#> [1] "hierarchy"
#> 
#> $fun
#> NULL
#> 
#> $groups
#> NULL
#> 
#> $order
#> [1] "A" "B"
#> 
#> $passing
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
create_group_inheritance("custom", groups = "mixed",
  fun = function(mother_group, father_group, ...) {
    if (is.na(father_group) || mother_group == father_group) mother_group else "mixed"
  })
#> $method
#> [1] "custom"
#> 
#> $fun
#> function (mother_group, father_group, ...) 
#> {
#>     if (is.na(father_group) || mother_group == father_group) 
#>         mother_group
#>     else "mixed"
#> }
#> <environment: 0xaaabcc727068>
#> 
#> $groups
#> [1] "mixed"
#> 
#> $order
#> NULL
#> 
#> $passing
#> [1] 0
#> 
#> attr(,"class")
#> [1] "hammel_group_inheritance"
```
