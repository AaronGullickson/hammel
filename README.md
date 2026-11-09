# Hammel <img src="man/figures/logo.png" align="right" width="120" alt="Hammel hex sticker showing a branching kinship network" />

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/AaronGullickson/hammel/graph/badge.svg)](https://app.codecov.io/gh/AaronGullickson/hammel)
[![R-CMD-check](https://github.com/AaronGullickson/hammel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AaronGullickson/hammel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Hammel is an R package for doing demographic microsimulation. Hammel is a spiritual successor to the venerable SOCSIM program but provides a modern interface, enhanced features, greater user customization, and faster computation. Hammel models demographic processes in continuous time and features:

- Stratification of rate schedules by time, group membership, union type, sex, and parity (for fertility).
- Union formation processes that including union dissolution, cohabitation and remarriage, and customizable rules for partner matching.
- Multiple user-defined group identifiers (e.g. race, education, region) complete with age or duration-specific transition rates, customizable partner matching rules across groups, customizable inheritance rules for newborns, and ancestry tracking.
- Helper functions to extract important details from simulations such as age-specific rates, population size, and kinship relations.
- Helper functions to build starting populations and to stratify rate schedules.
- A pipable interface for building and modifying simulation models.

## Installation

Hammel is currently only available through GitHub. You can install the latest development version with `pak` or `devtools`:

``` r
# install.packages("pak")
pak::pak("AaronGullickson/hammel")

# install.packages("devtools")
devtools::install_github("AaronGullickson/hammel")
```

## Example

Several detailed examples of Hammel simulations are available on the [documentation website](https://hammel.github.io). Here I provide a simple toy example.

To use Hammel, you first set up a simulation model, which can be constructed using a pipable interface. The primary inputs for Hammel are tables of age-specific rate schedules.

```r
data(usa_2000_mortality)
data(usa_2000_fertility)

model <- hammel_model() |>
  set_fertility(usa_2000_fertility) |>
  set_mortality(usa_2000_mortality) |>
  set_union_formation("on_demand")
```

All models must have fertility and mortality rates and some way to deal with union formation. Instead of supplying marriage rates, I am using an "on demand" model for union formation that partners off unpartnered women when they have a birth, but `set_union_formation` allows customization based on marriages rates, dissolution rates, and cohabitation.

Once I have the model, I can add a starting population and simulate. In many cases, you would "burn in" a simulation to get a stable population distribution for the starting population, but I won't do that here.

```r
pop_start <- generate_uniform_population(10000, end_time = 0)
model <- model |> 
  set_population(pop_start)
pop_end <- run_simulation(model, 100, seed = 42)
```

Now, I can compute summary statistics.

```r
calculate_population_size(pop_end, seq(0, 100, 5))
```

