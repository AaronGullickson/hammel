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
  on each side). Failed searches retain births with unknown father; exclusive
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
  user October 15, 2026. See [testing](testing.md#deferred-matching-performance-comparison).

- Initial fertility calibration: retain simple initialization and assess burn-in.
- Fertility transition adjustment after replacing rates between runs: document
  the mismatch, with no correction implemented.
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

Find a future validation/example population-period with union-entry (or legal
marriage) rates plus partnered/married and unpartnered/unmarried fertility rates.
Five-year age bands are acceptable. Prefer a common source and consistent
status/exposure definitions; do not combine unrelated schedules and present
them as a coherent empirical case. No dataset has yet been selected.


### Status-specific fertility (D039; October 15, 2026)

Fertility accepts optional `union_status` with both `partnered` and `unpartnered`.
No column means an overall schedule. Rate mode applies overall fertility to
both statuses with one warning per call; on-demand rejects conditional tables.
Rate-mode unpartnered births select fathers without unions under D040. Transitions
switch schedules with global maternal spacing retained. Compensation is computed
separately per status; history-aware compensation remains deferred. Full suite:
2269 passing expectations, no failures/warnings, one intentional live-data skip.

### Nonmarital father selection (D040; October 15, 2026)

Rate-mode nonmarital births select a father from the existing available male
pool using the same sampling, eligibility and formula model as partners.
Selection creates parentage only and does not reserve the man or form a union.
No qualifying sample leaves father unknown. Paternal links enter native kinship
immediately. Separate models are explicitly deferred; preferential reuse of a
previous father at another birth/union remains unresolved. The pinned ordinary
SOCSIM paths have no such preferential reuse; its optional random-father fallback
also differs by allowing partnered men and not using spouse scores. Full suite:
2313 passing expectations, no failures/warnings, one intentional live-data skip.
