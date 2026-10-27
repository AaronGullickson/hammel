
# hammel

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/AaronGullickson/hammel/graph/badge.svg)](https://app.codecov.io/gh/AaronGullickson/hammel)
[![R-CMD-check](https://github.com/AaronGullickson/hammel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AaronGullickson/hammel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Hammel is a standalone demographic microsimulator with an R interface and a
C++17 engine. It models continuous-time mortality, fertility, and union
formation, with group-specific rates, partner matching, ancestry, and newborn
identity inheritance. Results expose population and union histories as tibbles.

Milestones 1–3 are complete by user acceptance: mortality/fertility, union
formation, and groups. See the [documentation index](docs/README.md),
[scientific contract](docs/scientific-contract.md), and
[validation record](docs/testing.md) for implemented behavior and limitations.
Future mechanisms are listed in the [roadmap](docs/roadmap.md).

## Installation

You can install the development version of hammel from [GitHub](https://github.com/AaronGullickson/hammel) with:

``` r
# install.packages("pak")
pak::pak("AaronGullickson/hammel")
```

## Example

The [union-formation example](inst/examples/milestone2.R) demonstrates model
construction and simulation. The user-authored
[milestone 3 report](validation/milestone3.qmd) assesses group-specific rates,
matching, ancestry, inheritance, and performance. Rendering the report installs
the package, retrieves HMD rates, and runs substantial simulations; ordinary
package tests use offline inputs unless live database tests are explicitly enabled.

