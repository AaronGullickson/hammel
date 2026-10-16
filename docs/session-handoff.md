# Current project handoff

Updated October 15, 2026. Milestone 1 is complete by explicit user acceptance
following the rendered validation report (D029). This is a dated status record,
not authorization for additional features. Start with the [documentation index](README.md).

## Implemented foundation

- R model/population/process construction; all public tables are tibbles.
- Continuous C++17/Rcpp mortality/fertility scheduling with R RNG and call-local
  state; midpoint compensation uses ordinary Band and shared waiting_time (D028).
- Fertility-process `birth_interval`, default 0.75 years, accepting fractional
  durations. Recorded histories enforce spacing on restart.
- Finite fertility inputs only; HFD open endpoints are omitted. Mortality
  requires complete age/sex coverage ending in an open interval.
- Conditional-logit on-demand union formation, partner-based father assignment,
  minimum male age (default 15), and known-kin exclusion (default two generations
  on each side). Failed matches cancel births by default (D043); exclusive
  ongoing unions are enforced before sampling (D037). Shared native kinship also backs R `extract_kin()`.
- Output reuse with resampled candidates, not exact checkpoints.
- Lexis exposure/counts, TFR, life tables, and population sizes; testthat suite
  and offline analytical checks pass.

See [design](design.md), [scientific contract](scientific-contract.md), and
[testing](testing.md) for details and evidence.

## Completed validation and ownership

Aaron Gullickson completed [milestone1.qmd](../validation/milestone1.qmd), rendered
`validation/milestone1.html`, and accepted milestone 1 on October 15, 2026.
The report includes demographic coverage/deviation summaries, a two-sex Leslie
projection, stable growth and age structure comparisons, extended runs,
fertility-rate restart assessments, and performance/test results.

The agent reviewed the source and identified the local rendered artifact;
it did not rerun simulations or independently verify all scientific conclusions.
[Testing](testing.md) records the artifact hash, source commit, and remaining
archive/provenance limits. Preserve the user's analysis and source comments.

## Next work

Milestone 2 union-formation scope is selected (D030); see
[the specification notes](union-formation-models.md). Detailed modeling choices
remain unresolved. Process construction and a conditional-logit on-demand C++
path are implemented (D031), including current-partner caches, union histories,
widowhood, and restart reconstruction. D032 now requires an explicit union
process for every simulation; partner state is initialized unconditionally. Minimum male age screening is implemented (default 15 years, user-selectable);
known-kin exclusion is also implemented (default two generations on each side),
with shared C++ traversal and R extraction. Male exclusivity is now enforced (D037);
and rate-based initiation is implemented as described below. Matching
uses the user formula, with native linear/quadratic evaluation and an R fallback
for other transformations. The report proposes a
GitHub issue for fertility compensation after changing rates at restart; no
issue was created as part of recording completion. Follow the user's next request.

## Deferred work

- Benchmark R formula-fallback matching (for example a cubic age-gap term)
  against the native linear/quadratic path later; explicitly deferred by the
  user October 15, 2026; now scheduled for the final milestone validation
  report under D042. See [testing](testing.md#deferred-matching-performance-comparison).

- Initial fertility calibration: retain simple initialization and assess burn-in.
- Fertility compensation across status/rate transitions: required remaining
  milestone 2 work (D042); no corrective algorithm selected.
- African historical benchmark/new rate-source preparation: postponed.
- Non-death union dissolution remains beyond milestone 2. Union formation and
  matching are now active specification work, with detailed choices still open.
- Calendar-varying schedules, arbitrary custom processes, and exact checkpoints.

Earlier inspection/run records are available in Git history. Do not use a historical
pause request, deprecated API, or superseded timing rule to govern current work.


### Exclusive ongoing unions (October 15, 2026; D037)

The available male pool now contains only living, unpartnered men, before
candidate sampling. Formation removes a man; widowhood returns a surviving
husband. Both sexes have a current-partner cache, reconstructed on restart;
conflicting active input unions are rejected. Failed searches still preserve
births with unknown fathers. Prior concurrent-union outputs need rerunning to
apply this policy. Polygamy/polyamory is deferred.


### Rate-based searches (D038; October 15, 2026)

`initiation = "rate"` now schedules female unpartnered search attempts as
competing events. Failed attempts draw a new waiting time. Both initiation
modes share search/formation code. Fertility follows the current status schedule
or an overall schedule for all women (D039); widowhood invalidates pending
events and restores search risk. Existing
fertility compensation is retained without a population-ASFR recovery claim.
Status-specific fertility is implemented (D039); refusal/wait remains future work.


### Requested coherent empirical example

England and Wales 2015 is now the selected empirical example. `validation/data`
contains female marriage rates, partnered/unpartnered and overall fertility,
tidy ONS Table 6 for 1999-2015, the source workbook, provenance and extraction
script. Published age groups are retained. Legal marital status includes
cohabiting women among unmarried women; non-death dissolution is not simulated.


### Status-specific fertility (D039; October 15, 2026)

Fertility accepts optional `union_status` with both `partnered` and `unpartnered`.
No column means an overall schedule. Rate mode applies overall fertility to
both statuses with one warning per call; on-demand rejects conditional tables.
Rate-mode unpartnered births select fathers without unions under D040. Transitions
switch schedules with global maternal spacing retained. Compensation is computed
separately per status; transition compensation is required remaining work (D042). Full suite:
2269 passing expectations, no failures/warnings, one intentional live-data skip.

### Nonmarital father selection (D040; October 15, 2026)

Rate-mode nonmarital births select a father from the existing available male
pool using the same sampling, eligibility and formula model as partners.
Selection creates parentage only and does not reserve the man or form a union.
No qualifying sample leaves father unknown. Paternal links enter native kinship
immediately. Separate models are explicitly deferred; preferential reuse of a
previous father at another birth/union is deferred to a later milestone (D042). The pinned ordinary
SOCSIM paths have no such preferential reuse; its optional random-father fallback
also differs by allowing partnered men and not using spouse scores. Full suite:
2313 passing expectations, no failures/warnings, one intentional live-data skip.


### Optional father-matching model (D041; October 15, 2026)

`union_formation()` and `set_union_formation()` now accept `father_matching = NULL`
as a trailing argument. NULL uses the current general `matching` specification;
an explicit model affects only unpartnered rate-mode father selection. General
matching still governs all union searches and on-demand initiation. Both models
are validated/prepared independently; candidate sampling and eligibility remain
shared. This supersedes D040's separate-model deferral. Previous-father
preference is deferred to a later milestone (D042).


### Revised completion queue (D042; October 15, 2026)

Required before milestone 2 completion: fertility compensation and demographic
validation/report. The cancellation toggle is implemented under D043 below.
Specify the compensation method before implementing dependent behavior. The wait alternative is under
consideration and may be postponed. Previous-father continuity is deferred to
a later milestone. Test R formula-fallback performance in the final validation
report. No new simulation behavior or numerical acceptance criteria were adopted
by this scope update.


### Birth cancellation (D043; October 15, 2026)

`cancel_birth_if_no_match = TRUE` now defaults on in both union constructors.
Failed on-demand searches and nonmarital father selections cancel the attempted
birth without allocating IDs, recording births/kinship, drawing sex, or updating
last_birth. Draw new waiting times immediately; no new spacing interval starts.
FALSE explicitly preserves unknown-father births. Initial missing parent IDs
remain unchanged. Old process objects lacking the required logical toggle must
be reconstructed with the current constructor. Full suite: 2403 passing
expectations, no failures/warnings, one intentional live-data skip.
