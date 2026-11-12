# Define who moves together in household transitions

Define a set of rules for determining household membership for
transitions.

## Usage

``` r
household(age_independence = 18)
```

## Arguments

- age_independence:

  Age in years at which a child stops being a dependent: below it, an
  unpartnered child with a living parent is a dependent. Defaults to 18.

## Value

A `hammel_household`, to pass as `household` to
[`set_group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md).
Each transition process may use its own rule.

## Details

A household rule describes a household for
[`set_group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md).
The simulation does not track households independently, but determines
who belongs to a household when a household would move and assembles the
moving party.

A person's household is their current partner (married or cohabiting)
and the dependents attached to either of them. A dependent is alive,
unpartnered, younger than `age_independence`, and attached to their
mother while she is alive, otherwise to their father while he is alive.
People with neither parent alive are independent. A person under
`age_independence` who leaves a union is a dependent again if a parent
is alive.

## Examples

``` r
household(age_independence = 18)
#> $age_independence
#> [1] 18
#> 
#> attr(,"class")
#> [1] "hammel_household"
```
