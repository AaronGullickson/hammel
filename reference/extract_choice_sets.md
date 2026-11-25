# Extract matching choice sets from simulation results

Build wife-side conditional-logit choice sets for unions formed in a
simulation run, with each candidate man's log matching weight as an
offset.

## Usage

``` r
extract_choice_sets(
  population,
  unions = 1000,
  alternatives = 50,
  union_types = c("marriage", "cohabitation")
)
```

## Arguments

- population:

  A `hammel_population` returned by
  [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)
  with `union_weight_sample` above 0.

- unions:

  The number of recorded unions to sample, or all of them if fewer.

- alternatives:

  The number of other eligible men sampled for each union, or all of
  them if fewer.

- union_types:

  The union types to sample from: `"marriage"`, `"cohabitation"` or
  both.

## Value

A tibble with one row per candidate: `choice_id` (the union ID),
`union_type`, `time`, `chosen`, `wife_id`, `man_id`, `birth_time_w`,
`birth_time_h`, each matching dimension's category at `time` with `_w`
and `_h` suffixes, rounded ancestry shares named
`<dimension>_<category>_w`/`_h` when matching used ancestry, and
`log_weight`.

## Details

The run must keep matching weights, with the `union_weight_sample`
argument in
[`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)
set to a value above zero. Each choice set contains the actual recorded
union (`chosen = TRUE`) and a specified number of counterfactual unions
(`chosen = FALSE`) created from a uninform sample of the other men
eligible when it formed, who were alive, at least the union process's
`min_male_age`, had arrived if immigrants, and were in no existing
union. Known-kin exclusions and hard group restrictions (`NA` in
`group_exogamy`) are not applied.

Each candidate's `log_weight` is the calibrated weight of his matching
type (categories and rounded ancestry profile) at that search, 0 when
his model has no calibrated weights. With it as an offset, a conditional
logit of `chosen` on the matching model's terms estimates its
coefficients: for example,
`survival::clogit(chosen ~ age_diff + I(age_diff^2) + offset(log_weight) + strata(choice_id), data = sets)`.
