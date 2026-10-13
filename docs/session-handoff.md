# Current project handoff

Updated October 13, 2026. This is a concise status record, not a stop instruction
or authorization for additional features. Follow the user's current request.
Start with the [documentation index](README.md).

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

## Active validation work

The user is developing [milestone1.qmd](../validation/milestone1.qmd) from the
[replication example](../inst/examples/milestone1_bootstrap.R). The template
renders successfully with execution disabled; it is not a completed analysis.
The example remains USA/2000, with 100,000 initial persons and 100 runs.
No alternative data source/country is currently selected.

The user reported decadal TFR means around 2.05–2.06 in 100 runs, and first-year
TFR around 2.12 settling to roughly 2.05–2.06 in later years. Those reported
summaries preceded the final midpoint rewrite and are not current-version
verification. More recently, the user reported close ASFR/ASMR agreement with
input curves, with smoother ASMRs under 20-year pooling. The plots and
replication outputs were not independently rerun or archived by the agent.

Life-expectancy summaries use 100+ pooling to avoid sparse older cells. Direct
ASMR comparisons use single-year bands through 110 with an appended open tail.
These are distinct comparisons; the input life table retains its original
older-age bands. Missing ASMR cells and zero terminal estimates can cause NA
or infinite life expectancy. Report coverage rather than silently ignoring it.

## Next work and ownership

The remaining milestone work is to quantify agreement and uncertainty, broaden
rate/spacing cases, agree on acceptance criteria, and complete the reproducible
report. The user will write its analysis and interpretation; do not fill in
scientific conclusions or start a new replication study without a request.
Preserve the user's edits and explanatory code comments.

## Deferred work

- Initial fertility calibration: retain simple initialization and assess burn-in.
- Fertility transition adjustment after replacing rates between runs: document
  the mismatch, with no correction implemented.
- African historical benchmark/new rate-source preparation: postponed.
- Union formation and matching: [planning notes](union-formation-models.md),
  with no algorithm accepted and no implementation authorization.
- Calendar-varying schedules, arbitrary custom processes, and exact checkpoints.

Earlier inspection/run records are available in Git history. Do not use a historical
pause request, deprecated API, or superseded timing rule to govern current work.
