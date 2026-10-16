# Hammel project brief

Updated October 15, 2026. The first mortality/fertility simulator is implemented.
The user accepted milestone 1 as complete after rendering the validation report
on October 15, 2026. See D029 and [testing](testing.md) for the evidence and limits.

## Purpose

Hammel is a standalone demographic microsimulator with a modern R package
interface, informed by SOCSIM's approach and experience. It is a substantial
rewrite, not a SOCSIM wrapper, rsocsim companion package, or upstream contribution.
Legacy file compatibility must not determine its architecture.

The intended workflow combines pipeable model construction, inspectable R
populations and rate tables, separate simulation execution, and eventual
customizable demographic mechanisms. The user owns demographic requirements
and scientific interpretation; agents investigate, implement, document, and
validate software.

## Milestone 1 scope and status

The accepted initial scope is fertility and mortality schedules constant over
calendar time, but potentially varying by age (D003/D019).

| Component | Current state |
| --- | --- |
| Population/model construction | Implemented; public tables are tibbles |
| Mortality | Piecewise constant age/sex hazards; continuous sampling |
| Fertility | Population ASFR inputs; midpoint spacing compensation (D028) |
| Birth spacing | Fertility-process setting, default 0.75 years |
| Parent assignment | Uniformly sampled living male; unknown father if none exists |
| Relationships | Parent histories retained; no new unions formed |
| Analysis | Exposure/counts, ASFR/ASMR, TFR, life tables, population sizes |
| Restarts | Output can initialize another run; not an exact checkpoint |
| Validation | Milestone 1 accepted by the user October 15, 2026; completed report in `validation/milestone1.qmd` |

D019 replaced the original basic-partnership requirement with random male
parent assignment for the first executable milestone. D024 replaced monthly
execution with continuous event times. D028 replaced monthly and exact
age-varying fertility compensation with one adjusted rate per original age
band, evaluated at its midpoint. Those older descriptions are historical.

## Milestone 1 completion

The completed [validation report](../validation/milestone1.qmd) documents the
user's assessment of replicated fertility/mortality measures, population growth
and stable age structure against a two-sex Leslie projection, restart behavior,
and performance. Its conclusions and the user's explicit acceptance establish
milestone completion; no additional numerical acceptance thresholds are inferred.

Sparse oldest-age rates, initial fertility transients, and fertility compensation
after rate changes remain limitations. The restart issue is deferred beyond this
milestone. A controlled SOCSIM comparison has not been performed. Completion
does not claim accuracy for every schedule or authorize new mechanisms.
Union formation, calendar-varying rates, and initialization correction remain
outside milestone 1's scope. Milestone 2 now targets union formation (D030), as specified below.

## Milestone 2 scope

The user selected union formation for milestone 2: female search initiation
separate from conditional-logit matching, birth-triggered and scheduled search
modes, initial age-difference matching, kinship screening, and partnered/unpartnered
fertility in scheduled mode. Only death ends unions in this milestone. Parity-specific fertility is deferred
to a later milestone.
Candidate sampling defaults to 1000 men before pair screening (zero means the
whole pool). Optional waiting, detailed eligibility, and fertility-transition
rules remain open. See [milestone 2 notes](union-formation-models.md)
for scope, unresolved decisions, and validation. The on-demand C++ path now forms unions at births with formula-based
conditional-logit selection and handles widowhood and restart histories.
Minimum male age screening (default 15 years) and known-kin exclusion within
two generations of each partner are implemented. `extract_kin()` shares the
native kinship module. Exclusive ongoing unions are enforced before sampling; female rate-based searches and status-specific fertility schedules are implemented.
Overall fertility is also supported in rate mode with a warning; optional waiting remains pending.

## Boundaries and future work

C++17/Rcpp is selected and implemented for the first engine (D019); tibble is
used at the R boundary (D023). Package metadata and licensing are not yet ready
for release. A universal simulation framework and legacy compatibility are not
requirements. Union formation and matching are now in milestone 2 development. Custom
process execution, non-death union dissolution, and richer initialization remain
future work, not automatically authorized features.

The user has postponed finding an African historical alternative to USA/2000.
No replacement case or new rate-source integration is currently selected.

## Provenance and related documents

The project direction originated in the October 11, 2026 conversation,
summarized in `/home/aarong/Work/Hammel-project-notes.md`. Current decisions are
recorded in the [decision log](decisions.md). See [design](design.md) for the
implemented workflow and [scientific contract](scientific-contract.md) for
units, timing, rate interpretation, and limitations.
