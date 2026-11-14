# Summarize a group transition process

Summary method for class `"hammel_group_transition"`

## Usage

``` r
# S3 method for class 'hammel_group_transition'
summary(object, categories = NULL, ...)

# S3 method for class 'summary.hammel_group_transition'
print(x, ...)
```

## Arguments

- object:

  A `hammel_group_transition` object returned by
  [`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md).

- categories:

  Optional character vector of the dimension's categories, used to count
  possible transitions. `NULL` (default) uses the categories that appear
  in the `origin` or `destination` columns.

- ...:

  Additional arguments, currently unused.

- x:

  A `summary.hammel_group_transition` object returned by
  `summary.hammel_group_transition()`.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.hammel_group_transition` containing:

- group:

  The name of the group dimension whose categories change, or `NULL`
  when the process applies to the model's only dimension.

- categories:

  The categories used to define possible transitions.

- transitions:

  A tibble of the `origin` and `destination` pairs that have a positive
  rate in at least one band or stratum.

- missing_transitions:

  A tibble of the possible `origin` and `destination` pairs between
  distinct categories that have no positive rate in any band or stratum.

- rate_strata:

  A list describing the transition rate schedule apart from its `origin`
  and `destination` columns, containing `stratified_by` (stratum column
  names), `strata` (a tibble of observed stratum combinations), `bands`
  (a tibble of exact band combinations), and `bands_by_stratum` (whether
  band combinations differ by stratum). Shared bands are stored once;
  otherwise `bands` includes the stratum columns. Band endpoints define
  half-open intervals `[start, end)`.

- weights:

  A summary of the weight model used for selection within rate cells: a
  `summary.hammel_weight_model`, a list of `periods` and `models` for a
  calendar schedule, or `NULL` when selection is unweighted.

- age_independence:

  The household age of independence, in years, or `NULL` when
  transitions move individuals rather than households.

## See also

[`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md),
[`household()`](https://aarongullickson.github.io/hammel/reference/household.md)

## Examples

``` r
rates <- tibble::tribble(
  ~age_start, ~age_end, ~origin, ~destination, ~rate,
  0,          Inf,      "A",     "B",          0.05,
  0,          Inf,      "B",     "A",          0.02
)
process <- group_transition(rates, household = household())
summary(process)
#> 
#> ── Group transition process ──
#> 
#> • Group dimension: the model's only dimension
#> • Transitions (2 of 2 possible): A → B, B → A
#> • Missing transitions: none
#> • Stratified by: none
#> • Age bands (years): [0, Inf)
#> • Household moves: yes, age of independence 18 years
#> • Weight model: none
```
