# Hammel project brief

Updated November 6, 2026. The user accepted milestone 1 (mortality/fertility)
on October 15, 2026 (D029), milestone 2 (union formation) on October 19,
2026 (D050), milestone 3 (groups) on October 27, 2026 (D077), and milestone 4
(dissolution/remarriage/cohabitation) on October 29, 2026 (D093). Milestone 5
(parity-specific fertility, cell-accounted spacing and fertility heterogeneity)
was accepted October 29, 2026 (D105). The user accepted individual group transitions, milestone 6, on October 30,
2026 (D116). Milestone 7 calendar-varying rates (D118) and cohabitation
probabilities (D122), plus matching/inheritance model schedules (D123), are
accepted complete on November 2, 2026 (D127), based on the user-authored
[validation report](../validation/milestone7.qmd) and explicit user acceptance.
Milestone 8, cell scheduling for every rate process (D128) and weight models
(D129), was accepted complete on November 5, 2026 (D130).
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
| Fertility | Population ASFR inputs; midpoint compensation under spacing (D028, later superseded by D101) |
| Birth spacing | Fertility-process setting, default 0.75 years |
| Parent assignment | Uniformly sampled living male; unknown father if none exists |
| Relationships | Parent histories retained; no new unions formed |
| Analysis | Exposure/counts, ASFR/ASMR, TFR, life tables, population sizes |
| Restarts | Output can initialize another run; not an exact checkpoint |
| Validation | Milestone 1 accepted by the user October 15, 2026; completed report in `validation/milestone1.qmd` |

D019 replaced the original basic-partnership requirement with random male
parent assignment for the first executable milestone. D024 replaced monthly
execution with continuous event times. D028 replaced monthly and exact
age-varying fertility compensation with midpoint compensation, which D101
replaced with cell-accounted fertility. Those older descriptions are historical.

## Milestone 1 completion

The completed [validation report](../validation/milestone1.qmd) documents the
user's assessment of replicated fertility/mortality measures, population growth
and stable age structure against a two-sex Leslie projection, restart behavior,
and performance. Its conclusions and the user's explicit acceptance establish
milestone completion; no additional numerical acceptance thresholds are inferred.

At milestone 1 acceptance, sparse oldest-age rates, initial fertility transients
and compensation after rate changes remained limitations. Milestone 5 later
replaced that correction with cell accounting and revisited the restart case.
No broad SOCSIM equivalence comparison was part of milestone 1; the later
parity/spacing reference experiment addresses a narrower question. Completion
does not claim accuracy for every schedule or authorize new mechanisms.
Union formation, calendar-varying rates, and initialization correction remain
outside milestone 1's scope. Milestone 2 now targets union formation (D030), as specified below.

## Milestone 2 scope

The user selected union formation for milestone 2: female search initiation
separate from conditional-logit matching, birth-triggered and scheduled search
modes, initial age-difference matching, kinship screening, and partnered/unpartnered
fertility in scheduled mode. Only death ended unions in the milestone 2 acceptance scope. Parity-specific
fertility was subsequently added in milestone 5.
Partners are now drawn exactly from the whole eligible pool, without a
candidate sample (D072). Eligibility and basic fertility transitions are implemented; the
failed-match birth cancellation is implemented (default TRUE); transition-aware fertility
compensation was deferred at milestone 2 acceptance (D045/D050), then replaced
by cell accounting in milestone 5 (D101). The poor-pool wait alternative is deferred to a later milestone (D044). See [milestone 2 notes](union-formation-models.md)
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

The [roadmap](roadmap.md) records multiple group columns as the current
milestone 9 (D131/D134). Under D136, once it is completed and accepted, work
turns to the first version bump and public release. External immigration and
two-sex marriage models are post-release projects. The release version,
criteria and scope of proposed summaries/helpers remain unselected; the
two-sex model and API also remain unselected (D117). No future feature or
publication is authorized by milestone completion. Fixed newborn assignment was
selected November 7, 2026 (D142), using an explicit category per dimension; see
[group inheritance](group-inheritance.md).

## Milestone 4 complete

The user accepted milestone 4 on October 29, 2026 (D093), based on the
[final validation report](../validation/milestone4.qmd). Its scope includes
optional duration-dependent dissolution (D078), with
annual hazards per union (originally one female-carried clock per union; cells since D128), symmetric group-pair
schedules/multipliers (D079), and birth before dissolution on exact same-person
ties. See [union dissolution](union-dissolution.md). Search schedules optionally stratify by `previous_union`
(D080) and select age or duration since the latest ending per category in one
rate table (D083); see [remarriage](remarriage.md). Optional rate-mode cohabitation adds type selection, separate matching,
type-specific rates and same-partner conversion (D085/D086); see
[cohabitation](cohabitation.md). Separate demographic, partnership-entry,
marriage-entry, dissolution and conversion observations use their documented
risk denominators. `get_unpartnered_births()` supports father-matching assessment.
Acceptance retains the documented fertility-transition, sparse-pool,
initialization, and calibration limitations; it does not authorize later work.

## Milestones 5 and 6 complete; time-varying inputs next

The user accepted milestone 5 on October 29, 2026 (D105), based on the
[report](../validation/milestone5.qmd) and its rendered assessment. Its final
render was running at closeout and is preserved. The milestone implements
conditional parity schedules and HFD preparation, history-derived actual parity,
PATFR/TFR dispatch, parity-specific observations, cell-accounted fertility and
persistent lognormal fertility weights with maternal inheritance (D094–D104).
The report revisits baseline, restart and partnered/unpartnered fertility,
then assesses parity rates, heterogeneity and performance.

Cell accounting is the only birth-spacing mechanism. Cells with no eligible
women count unrealized births (D128); adding them to an analysis is an accounting
diagnostic, not actual births. Sparse and structurally infeasible young/high-
parity cells can still fall below input rates. Retained histories determine
parity and spacing; missing prehistory is not imputed. Reproductive selection
can change inherited weight distributions, and restarts are not exact
checkpoints. Completion retains these limitations rather than asserting
universal rate recovery.

The user accepted milestone 6 individual transitions on October 30 (D116),
following implementation and report review; see [the transition contract](group-transitions.md). The
process uses annual hazards, common age or duration bands, directed pair
columns and optional sex/union-type strata. Group histories are retained and
affected risks, cells and matching state update at moves. Ancestry stays
unchanged. Household moves remain deferred (D114). Multiple group columns are
now milestone 9 (D121, D131), potentially covering the former standalone
regions/internal-migration work; schema and household rules remain unselected. The final milestone 6 render was active at
acceptance; preserve the report and intermediates. Acceptance retains history,
finite-cell and non-checkpoint-restart limitations. The user subsequently
approved calendar-varying rate semantics and authorized their implementation
under D118; see [the rate contract](time-varying-rates.md). Whole-model matching and
inheritance schedules are implemented under D123; see
[the model contract](time-varying-models.md). The user accepted milestone 7
as complete on November 2, 2026 (D127). Acceptance retains the report and
current-contract limitations.

## Milestone 8 complete

The user accepted milestone 8 on November 5, 2026 (D130), based on the
[validation report](../validation/milestone8.qmd) and its rendered assessment.
Every rate process (mortality, fertility, rate-mode search, dissolution,
conversion and group transitions) is scheduled by cells of at-risk units with
one clock per process, cells stratifying only as finely as the input rates
(D128). Weight models, as log-linear formulas or custom functions, optionally
on a calendar schedule, choose who in a cell experiences an event without
changing cell intensities (D129). Unrealized births are counted, calendar rate
changes update cells in place where the layout allows, and matching
calibration takes intensities from the cells. See
[cell scheduling and weight models](binned-group-hazards.md). Acceptance
retains the documented limitations: weights are normalized within cells, so
population-level effects depend on cell size and composition, and no
eligibility rules beyond fertility spacing or user-defined characteristics are
implemented.

## Boundaries and future work

C++17/Rcpp is selected and implemented for the first engine (D019); tibble is
used at the R boundary (D023). Package metadata and licensing are not yet ready
for release. A universal simulation framework and legacy compatibility are not
requirements. Union formation and matching are implemented and milestone 2 is complete. Custom
process execution and richer initialization remain
future work, not automatically authorized features.

The user has postponed finding an African historical alternative to USA/2000.
No replacement case or new rate-source integration is currently selected.

## Provenance and related documents

The project direction originated in the October 11, 2026 conversation,
summarized in `/home/aarong/Work/Hammel-project-notes.md`. Current decisions are
recorded in the [decision log](decisions.md). See [specification](specification.md) for the
implemented workflow and [scientific contract](scientific-contract.md) for
units, timing, rate interpretation, and limitations.
