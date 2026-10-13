# Hammel project brief

Updated October 13, 2026. The first mortality/fertility simulator is implemented.
Milestone 1 remains in demographic validation; scientific acceptance thresholds
and the completed validation report are not yet established.

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
| Validation | Software/analytical checks pass; broader demographic assessment ongoing |

D019 replaced the original basic-partnership requirement with random male
parent assignment for the first executable milestone. D024 replaced monthly
execution with continuous event times. D028 replaced monthly and exact
age-varying fertility compensation with one adjusted rate per original age
band, evaluated at its midpoint. Those older descriptions are historical.

## Remaining milestone work

1. Complete the [Quarto validation report](../validation/milestone1.qmd), using
   saved inputs, recorded seeds/settings, and identifiable code versions.
2. Quantify replicated ASFR, ASMR, TFR, and life-expectancy agreement, with
   sampling uncertainty, exposure coverage, and comparable age bands.
3. Assess higher-fertility inputs, five-year fertility bands, and alternative
   birth intervals; define acceptable approximation errors with the user.
4. Record which criteria are met and which limitations remain.

Reported visual agreement in USA/2000 replications is encouraging evidence,
not a completed scientific assessment. Union formation, calendar-varying rate
support, and initialization correction are not requirements for closing this
current milestone.

## Boundaries and future work

C++17/Rcpp is selected and implemented for the first engine (D019); tibble is
used at the R boundary (D023). Package metadata and licensing are not yet ready
for release. A universal simulation framework and legacy compatibility are not
requirements. Custom process execution, union formation/dissolution, realistic
matching, and richer initialization remain future work, not automatically
authorized features.

The user has postponed finding an African historical alternative to USA/2000.
No replacement case or new rate-source integration is currently selected.

## Provenance and related documents

The project direction originated in the October 11, 2026 conversation,
summarized in `/home/aarong/Work/Hammel-project-notes.md`. Current decisions are
recorded in the [decision log](decisions.md). See [design](design.md) for the
implemented workflow and [scientific contract](scientific-contract.md) for
units, timing, rate interpretation, and limitations.
