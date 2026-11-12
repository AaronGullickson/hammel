# Create a Hammel population

Create a `hammel_population` object to serve as the starting population
for a simulation.

## Usage

``` r
hammel_population(
  persons,
  unions = NULL,
  end_time,
  ancestry = NULL,
  transitions = NULL,
  groups = NULL
)
```

## Arguments

- persons:

  A data frame with one row per person, including deceased persons. Must
  include a unique `person_id`, nonmissing `sex`, nonmissing
  `birth_time`, `death_time`, `mother_id`, and `father_id`.

- unions:

  A data frame with one row per union episode, or `NULL` for no unions.
  `wife_id` references the female partner and `husband_id` the male
  partner. Optional `union_type` is cohab or marriage.

- end_time:

  Required finite numeric reference time in years. Its calendar origin
  is supplied by the user. This supplies the next simulation's starting
  time. Recorded events must not occur after this time.

- ancestry:

  Optional ancestry shares for dimensions declared with
  `group(ancestry = TRUE)`: a named list with one data frame per such
  group dimension (a bare data frame when only one dimension tracks
  ancestry), each with exactly one row per person, keyed by `person_id`,
  one column of shares per category (named by the category) and an
  `unknown` share. Shares must be finite, between zero and one, and sum
  to one per row (tolerance `1e-8`).

- transitions:

  Optional data frame of individual category changes with `person_id`,
  `time`, `group` (the dimension; may be omitted when one is declared),
  `origin`, and `destination`. Each person's history in each dimension
  must chain to their current category and lie within life dates and the
  known endpoint. Grouped populations always include this tibble, even
  when empty, with character `group`, `origin` and `destination`.

- groups:

  Optional group declaration: a character vector of dimension names
  (each with
  [`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
  defaults) or a named list of
  [`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
  objects. Persons columns that are not declared are ordinary custom
  columns. Models usually declare groups with
  [`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md).

## Value

A named list of class `hammel_population` with components `persons`,
`unions`, and `end_time`, plus an `ancestry` tibble when supplied or
when persons has groups, and a `transitions` tibble when groups are
present.

## Details

In most cases, it will be more common to create a `hammel_population`
object by using
[`generate_uniform_population()`](https://aarongullickson.github.io/hammel/reference/generate_uniform_population.md)
or to receive a `hammel_population` object as the output of
[`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md).

The `hammel_population` object must always contain a `persons` data
frame and an `end_time` value. The `end_time` value defines the starting
point for the simulation and no dates in the `hammel_population` may
exceed the `end_time` value.

## Examples

``` r
population <- hammel_population(
  persons = data.frame(
    person_id = 1:2, sex = c("m", "f"), birth_time = c(-25, -23),
    death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_
  ),
  end_time = 0
)
```
