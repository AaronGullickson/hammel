# Create a Hammel population

Collect person and union histories in an inspectable population object.
This constructor does not run a simulation. It validates consistent
tables against the population schema. Custom columns are retained. An
optional `group` column contains nonmissing character labels or finite
numeric codes (character labels must be nonempty). Newborns inherit
their group from either parent with equal probability by default, with
maternal fallback, or use fertility's configured inheritance model.

## Usage

``` r
hammel_population(
  persons,
  unions = NULL,
  end_time = NULL,
  ancestry = NULL,
  transitions = NULL
)
```

## Arguments

- persons:

  A data frame with one row per person, including deceased persons.
  Columns are preserved and the table is returned as a tibble.

- unions:

  A data frame with one row per union episode, or `NULL` for no unions.
  `wife_id` references the female partner and `husband_id` the male
  partner; these roles do not imply legal marriage. An empty data frame
  is normalized to the union schema. Optional `union_type` is cohab or
  marriage. A conversion closes cohab with `end_reason = "conversion"`
  and links the same-time marriage through `previous_cohab_id`.

- end_time:

  Optional finite simulation endpoint in years. Its calendar origin is
  supplied by the user. `NULL` means that the endpoint is unknown.

- ancestry:

  Optional data frame with exactly one row per person, keyed by
  `person_id`, known shares in `group_<label>` columns, and an `unknown`
  share. Numeric shares must be finite, between zero and one, and sum to
  one per row (tolerance `1e-8`). Supplied shares are preserved, with
  rows aligned to persons. Without it, grouped persons start with 100%
  ancestry in their recorded group, regardless of parent histories.
  Without groups, no ancestry table is created automatically. Newborns
  average parental shares; a missing parent contributes its half to
  `unknown`. Identity remains separate.

- transitions:

  Optional data frame of individual group changes with `person_id`,
  `time`, `group_origin`, and `group_dest`. Histories must chain to
  current person groups and lie within life dates and the known
  endpoint. Grouped populations always include this tibble, even when
  empty.

## Value

A named list of class `hammel_population` with components `persons`,
`unions`, and `end_time`, plus an `ancestry` tibble when supplied or
when persons has groups, and a `transitions` tibble when groups are
present. Ancestry is retained across restarts and deaths.

## Examples

``` r
sim <- hammel_population(
  persons = data.frame(person_id = 1:2, sex = c("m", "f"),
    birth_time = c(-25, -23), death_time = NA_real_,
    mother_id = NA_integer_, father_id = NA_integer_),
  end_time = 0
)
sim$persons
#> # A tibble: 2 × 6
#>   person_id sex   birth_time death_time mother_id father_id
#>       <int> <chr>      <dbl>      <dbl>     <int>     <int>
#> 1         1 m            -25         NA        NA        NA
#> 2         2 f            -23         NA        NA        NA
sim$unions
#> # A tibble: 0 × 6
#> # ℹ 6 variables: union_id <int>, wife_id <int>, husband_id <int>,
#> #   start_time <dbl>, end_time <dbl>, end_reason <chr>
```
