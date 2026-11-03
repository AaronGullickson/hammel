# Validate a Hammel population

Check the population structure, population schemas, IDs, relationship
references, and date consistency. This function does not normalize or
change the object. Invalid objects raise an error; valid objects are
returned invisibly.

## Usage

``` r
validate_population(results)
```

## Arguments

- results:

  A `hammel_population` object.

## Value

The supplied object, invisibly and unchanged.

## Examples

``` r
persons <- data.frame(person_id = 1L, sex = "f", birth_time = -20,
  death_time = NA_real_, mother_id = NA_integer_, father_id = NA_integer_)
results <- hammel_population(persons, end_time = 0)
validate_population(results)
```
