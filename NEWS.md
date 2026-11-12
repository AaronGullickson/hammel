# hammel 0.1.0

First public release.

Hammel is a continuous-time demographic microsimulate inspired by SOCSIM. It uses an R interface to allow users to construct, run, and, analyze microdemographic simulations. The current release features:

- Stratification of rate schedules by time, group membership, union type, sex, and parity (for fertility).
- Union formation processes that including union dissolution, cohabitation and remarriage, and customizable rules for partner matching.
- Multiple user-defined group identifiers (e.g. race, education, region) complete with age or duration-specific transition rates, customizable partner matching rules across groups, customizable inheritance rules for newborns, and ancestry tracking.
- Helper functions to extract important details from simulations such as age-specific rates, population size, and kinship relations.
- Helper functions to build starting populations and to stratify rate schedules.
- A pipable interface for building and modifying simulation models.

Hammel is under active development and we appreciate your feedback!
