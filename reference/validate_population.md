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
