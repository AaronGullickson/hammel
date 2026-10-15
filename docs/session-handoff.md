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
- Random living-male parent assignment; unknown father if none exists. No
  new unions, minimum paternal age, or kinship restrictions.
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
process for every simulation; partner state is initialized unconditionally. Eligibility is deferred, and rate-based initiation remains rejected. Matching
uses the user formula, with native linear/quadratic evaluation and an R fallback
for other transformations. The report proposes a
GitHub issue for fertility compensation after changing rates at restart; no
issue was created as part of recording completion. Follow the user's next request.

## Deferred work

- Initial fertility calibration: retain simple initialization and assess burn-in.
- Fertility transition adjustment after replacing rates between runs: document
  the mismatch, with no correction implemented.
- African historical benchmark/new rate-source preparation: postponed.
- Non-death union dissolution remains beyond milestone 2. Union formation and
  matching are now active specification work, with detailed choices still open.
- Calendar-varying schedules, arbitrary custom processes, and exact checkpoints.

Earlier inspection/run records are available in Git history. Do not use a historical
pause request, deprecated API, or superseded timing rule to govern current work.
