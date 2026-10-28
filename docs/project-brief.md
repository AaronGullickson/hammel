# Hammel project brief

Updated October 27, 2026. The user accepted milestone 1 (mortality/fertility)
on October 15, 2026 (D029), milestone 2 (union formation) on October 19,
2026 (D050), and milestone 3 (groups) on October 27, 2026 (D077).
See [testing](testing.md) for evidence and limits.

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

## Accepted milestone 1 scope and status

The accepted initial scope was fertility and mortality schedules constant over
calendar time, but potentially varying by age (D003/D019). The table records
that milestone at acceptance, not the current mandatory-union engine (D032).

| Component | State at milestone 1 acceptance |
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
Partners are now drawn exactly from the whole eligible pool, without a
candidate sample (D072). Eligibility and basic fertility transitions are implemented; the
failed-match birth cancellation is implemented (default TRUE); transition-aware fertility
compensation remains deferred (D045), including after milestone acceptance (D050). The poor-pool wait alternative is deferred to a later milestone (D044). See [milestone 2 notes](union-formation-models.md)
for implemented rules, deferred decisions, and validation. The on-demand C++ path now forms unions at births with formula-based
conditional-logit selection and handles widowhood and restart histories.
Minimum male age screening (default 15 years) and known-kin exclusion within
two generations of each partner are implemented. `extract_kin()` shares the
native kinship module. Exclusive ongoing unions are enforced; female rate-based searches and status-specific fertility schedules are implemented.
Overall fertility is also supported in rate mode with a warning. Previous-father
continuity is deferred to a later milestone. The user-authored [milestone 2 report](../validation/milestone2.qmd) now includes
USA circa 2000–2001 demographic comparisons, matching-parameter recovery, and
R formula-fallback performance comparisons (D042). The user accepted milestone 2
as complete on October 19, 2026 (D050), retaining the documented limitations.
Matching defaults to
`create_age_gap_matching_model(0, 5)`, and optional union diagnostics record
the applied log male-group weights on new union records (D070/D072).

## Milestone 3 completion and future scope

The user accepted milestone 3 as complete on October 27, 2026 (D077), after
reviewing and revising the [validation report](../validation/milestone3.qmd).
Implemented features include group-specific schedules and proportional hazards,
group-stratified observations, separate ancestry proportions, symmetric pair
preferences, calibrated weighted choice with delta target adjustments, exact
full-pool partner selection, and built-in/custom newborn identity assignment.
Random-parent inheritance is the default, with maternal fallback (D076).

The report assesses rate and preference recovery, directional pairing under
different demographic conditions, inheritance and custom mixed identities,
sensitivity to starting populations, and performance. Scientific acceptance
belongs to the user; completion does not establish universal accuracy or remove
sparse-cell, calibration, or fertility-transition limitations. Weighted matching
does not model marriage squeezes; fuller two-sex models remain future work.

The [roadmap](roadmap.md) records proposed dissolution/remarriage, parity,
time-varying inputs, group transitions, summaries, and migration work. No future
feature is authorized by milestone completion. Fixed newborn assignment was
explicitly deferred to the group-transitions discussion; custom inheritance
already supports assignment of mixed-parent births to a specified category.

## Milestone 4 in progress

The first component is optional duration-dependent dissolution (D078), with
annual hazards and one female-carried clock per union, group schedules and
multipliers by her group, and birth before dissolution on exact same-person
ties. See [union dissolution](union-dissolution.md). Separate remarriage rates
and cohabitation remain deferred; no milestone 4 acceptance is implied.

## Boundaries and future work

C++17/Rcpp is selected and implemented for the first engine (D019); tibble is
used at the R boundary (D023). Package metadata and licensing are not yet ready
for release. A universal simulation framework and legacy compatibility are not
requirements. Union formation and matching are implemented and milestone 2 is complete. Custom
process execution, separate remarriage schedules, cohabitation and richer initialization remain
future work, not automatically authorized features.

The user has postponed finding an African historical alternative to USA/2000.
No replacement case or new rate-source integration is currently selected.

## Provenance and related documents

The project direction originated in the October 11, 2026 conversation,
summarized in `/home/aarong/Work/Hammel-project-notes.md`. Current decisions are
recorded in the [decision log](decisions.md). See [design](design.md) for the
implemented workflow and [scientific contract](scientific-contract.md) for
units, timing, rate interpretation, and limitations.
