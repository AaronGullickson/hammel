# hammel

Hammel is a standalone demographic microsimulator with an R interface
and a C++17 engine. It models continuous-time mortality, fertility, and
union formation, dissolution, remarriage, and optional
cohabitation/conversion, with group-specific rates, individual group
transitions, partner matching, ancestry, and newborn identity
inheritance. Results expose person, union, and transition histories as
tibbles.

Milestones 1–6 are complete by user acceptance: mortality/fertility,
union formation, groups, dissolution/remarriage/cohabitation,
parity-specific fertility with persistent fertility heterogeneity, and
individual group transitions (D116). Milestone 7, time-varying inputs,
is next. See the [documentation
index](https://aarongullickson.github.io/hammel/design/README.md),
[scientific
contract](https://aarongullickson.github.io/hammel/design/scientific-contract.md),
and [validation
record](https://aarongullickson.github.io/hammel/design/testing.md) for
implemented behavior and limitations. Future mechanisms are listed in
the
[roadmap](https://aarongullickson.github.io/hammel/design/roadmap.md).

## Installation

You can install the development version of hammel from
[GitHub](https://github.com/AaronGullickson/hammel) with:

\
`# install.packages("pak")`\
`pak``::`[`pak`](https://pak.r-lib.org/reference/pak.html)`(``"AaronGullickson/hammel"``)`

## Example

The [union-formation
example](https://aarongullickson.github.io/hammel/inst/examples/milestone2.R)
demonstrates model construction and simulation. The user-authored
[milestone 3
report](https://aarongullickson.github.io/hammel/validation/milestone3.qmd)
assesses group-specific rates, matching, ancestry, inheritance, and
performance. Rendering the report installs the package, retrieves HMD
rates, and runs substantial simulations; ordinary package tests use
offline inputs unless live database tests are explicitly enabled.

The [milestone 4
report](https://aarongullickson.github.io/hammel/validation/milestone4.qmd)
and its [rendered
assessment](https://aarongullickson.github.io/hammel/validation/milestone4.md)
cover dissolution, age/duration remarriage, cohabitation, conversion,
separate matching roles, and performance. The user accepted milestone 4
on October 29, 2026 (D093). Rendering retrieves HMD mortality, runs 100
independently initialized replications per scenario, and executes the
offline test suite. See [validation
contents](https://aarongullickson.github.io/hammel/validation/README.md)
for retained reports, inputs, and supporting diagnostic evidence.

The accepted [milestone 5
report](https://aarongullickson.github.io/hammel/validation/milestone5.qmd)
and [rendered
assessment](https://aarongullickson.github.io/hammel/validation/milestone5.md)
cover parity-specific fertility, cell-accounted birth spacing, restart
and union-type rate checks, and lognormal fertility weights with
maternal inheritance (D105). Completion retains the documented
sparse-cell/unrealized-intensity and reproductive-selection limitations.
Rendering installs the package, retrieves HFD/HMD inputs, runs
replications and the offline tests; it is not an ordinary package check.

The accepted [milestone 6
report](https://aarongullickson.github.io/hammel/validation/milestone6.qmd)
and [rendered
assessment](https://aarongullickson.github.io/hammel/validation/milestone6.md)
cover age/duration transitions, sex/union-type strata, demographic rate
and formation-group recovery, population shares and performance (D116,
October 30, 2026). The final render was active at acceptance; preserve
its intermediates until it finishes. Household moves are deferred to the
regions milestone (D114).
