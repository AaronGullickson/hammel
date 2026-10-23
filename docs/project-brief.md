# Hammel project brief

Updated October 19, 2026. The user accepted milestone 1 (mortality/fertility)
on October 15, 2026 (D029), and milestone 2 (union formation) on October 19,
2026 (D050). See [testing](testing.md) for evidence and limits.

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
eligible sample counts and selected utility ranks.

## Beyond milestone 2

The user selected groups as milestone 3's theme: intergroup matching, baseline
rate multipliers, flexible child group identity and possible proportional group
ancestry. Other future priorities are analysis helpers, input-rate preparation,
fertility heterogeneity, inherited/life-course characteristics (such as education),
and within-run immigration with generational status. Initial optional group schedules, maternal inheritance, and weighted baseline
sampling are specified in [groups](groups.md) (D052). Newborn identity now has
built-in/custom inheritance rules (D058), including declared mixed labels and
independently configured schedules; see [group inheritance](group-inheritance.md).
Multiple simultaneous identities and other richer demographic rules
remain open. See the [roadmap](roadmap.md); milestone 2 completion does
not authorize unspecified future demographic behavior.

## Boundaries and future work

C++17/Rcpp is selected and implemented for the first engine (D019); tibble is
used at the R boundary (D023). Package metadata and licensing are not yet ready
for release. A universal simulation framework and legacy compatibility are not
requirements. Union formation and matching are implemented and milestone 2 is complete. Custom
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
