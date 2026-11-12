# Package index

## Model setup and simulation

Construct a model, attach its population and processes, and run a
simulation.

- [`hammel_model()`](https://aarongullickson.github.io/hammel/reference/hammel_model.md)
  : Create an empty Hammel model specification
- [`set_population()`](https://aarongullickson.github.io/hammel/reference/set_population.md)
  : Set the starting population of a Hammel model
- [`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
  [`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md)
  : Describe and declare a group dimension
- [`set_process()`](https://aarongullickson.github.io/hammel/reference/set_process.md)
  : Set a demographic process specification
- [`run_simulation()`](https://aarongullickson.github.io/hammel/reference/run_simulation.md)
  : Run a demographic simulation

## Population preparation

Create and validate starting populations and describe their group
dimensions.

- [`hammel_population()`](https://aarongullickson.github.io/hammel/reference/hammel_population.md)
  : Create a Hammel population
- [`generate_uniform_population()`](https://aarongullickson.github.io/hammel/reference/generate_uniform_population.md)
  : Generate a population with uniformly distributed ages
- [`validate_population()`](https://aarongullickson.github.io/hammel/reference/validate_population.md)
  : Validate a Hammel population
- [`group()`](https://aarongullickson.github.io/hammel/reference/group.md)
  [`set_groups()`](https://aarongullickson.github.io/hammel/reference/group.md)
  : Describe and declare a group dimension

## Rate preparation

Import demographic rates and prepare proportional rate strata.

- [`get_mortality_hmd()`](https://aarongullickson.github.io/hammel/reference/get_mortality_hmd.md)
  : Download prepared mortality rates from the HMD
- [`get_fertility_hfd()`](https://aarongullickson.github.io/hammel/reference/get_fertility_hfd.md)
  : Download prepared fertility rates from the HFD
- [`stratify_rates()`](https://aarongullickson.github.io/hammel/reference/stratify_rates.md)
  : Add proportional strata to a rate table

## Demographic processes

Specify mortality, fertility, unions, and group transitions.

- [`mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
  [`set_mortality()`](https://aarongullickson.github.io/hammel/reference/mortality.md)
  : Specify an age-specific mortality process
- [`fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md)
  [`set_fertility()`](https://aarongullickson.github.io/hammel/reference/fertility.md)
  : Specify an age-specific fertility process
- [`union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md)
  [`set_union_formation()`](https://aarongullickson.github.io/hammel/reference/union_formation.md)
  : Specify a union-formation process
- [`union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md)
  [`set_union_dissolution()`](https://aarongullickson.github.io/hammel/reference/union_dissolution.md)
  : Specify duration-dependent union dissolution
- [`group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)
  [`set_group_transition()`](https://aarongullickson.github.io/hammel/reference/group_transition.md)
  : Specify individual group transitions
- [`household()`](https://aarongullickson.github.io/hammel/reference/household.md)
  : Define who moves together in household transitions

## Matching, inheritance, and weight models

Configure partner selection, newborn group inheritance, and selection
within rate cells.

- [`matching_model()`](https://aarongullickson.github.io/hammel/reference/matching_model.md)
  : Specify a conditional-logit matching model
- [`age_gap_matching_model()`](https://aarongullickson.github.io/hammel/reference/age_gap_matching_model.md)
  : Create a quadratic matching model from interpretable age-gap
  preferences
- [`group_inheritance()`](https://aarongullickson.github.io/hammel/reference/group_inheritance.md)
  : Specify newborn group identity inheritance
- [`weight_model()`](https://aarongullickson.github.io/hammel/reference/weight_model.md)
  : Create a weight model for selection within rate cells

## Demographic summaries

Calculate population sizes, observed event rates, and synthetic
demographic measures.

- [`calculate_population_size()`](https://aarongullickson.github.io/hammel/reference/calculate_population_size.md)
  : Calculate the living population at time points
- [`calculate_mortality_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_fertility_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_union_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_marriage_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_dissolution_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_cohab_convert_rate()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  [`calculate_transition_rates()`](https://aarongullickson.github.io/hammel/reference/calculate_rates.md)
  : Calculate observed demographic event rates
- [`calculate_life_table()`](https://aarongullickson.github.io/hammel/reference/calculate_life_table.md)
  : Construct period life tables from age-specific mortality rates
- [`calculate_life_expectancy()`](https://aarongullickson.github.io/hammel/reference/calculate_life_expectancy.md)
  : Calculate life expectancy at given ages
- [`calculate_tfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md)
  [`calculate_patfr()`](https://aarongullickson.github.io/hammel/reference/calculate_tfr.md)
  : Calculate total fertility and parity-adjusted total fertility rates

## Family and union histories

Extract kin relationships and enrich recorded person, birth, and union
histories.

- [`extract_kin()`](https://aarongullickson.github.io/hammel/reference/extract_kin.md)
  : Extract known kin within generation-distance bounds
- [`count_kin()`](https://aarongullickson.github.io/hammel/reference/count_kin.md)
  : Calculate the average number of living kin over time
- [`augment_persons()`](https://aarongullickson.github.io/hammel/reference/augment_persons.md)
  : Extract persons data from simulation results with additional detail
- [`augment_unions()`](https://aarongullickson.github.io/hammel/reference/augment_unions.md)
  : Extract union data from simulation results with additional detail
- [`extract_unpartnered_births()`](https://aarongullickson.github.io/hammel/reference/extract_unpartnered_births.md)
  : Extract births with known parents to unpartnered mothers

## Example data

Packaged demographic inputs for examples and exploration.

- [`usa_2000_mortality`](https://aarongullickson.github.io/hammel/reference/usa_2000.md)
  [`usa_2000_fertility`](https://aarongullickson.github.io/hammel/reference/usa_2000.md)
  [`usa_2000_fertility_union`](https://aarongullickson.github.io/hammel/reference/usa_2000.md)
  [`usa_2000_marriage`](https://aarongullickson.github.io/hammel/reference/usa_2000.md)
  [`usa_2000_migration`](https://aarongullickson.github.io/hammel/reference/usa_2000.md)
  : United States example rates, around 2000
- [`england_wales_1992_onward`](https://aarongullickson.github.io/hammel/reference/england_wales_1992_onward.md)
  : England and Wales inputs for a second demographic transition example
