# Current project handoff

Updated October 19, 2026. This is a dated status record, not authorization for
additional features. Start with the [documentation index](README.md); the
current user request controls work.

## Implemented state

- Call-local C++17/Rcpp continuous event scheduling, R RNG isolation, tibble
  interfaces, and recorded population/union histories.
- Mortality and fertility with shared hazard sampling, per-schedule midpoint
  compensation, and actual-birth spacing preserved across restarts.
- A mandatory union process with on-demand and female rate-based initiation,
  conditional-logit matching, sampling before age/kin screening, exclusive
  ongoing unions, widowhood, and event invalidation after status changes.
- Overall fertility or, in rate mode, partnered/unpartnered schedules. Nonmarital
  births select fathers without creating unions; `father_matching` can override
  the general matching model. Failed matches cancel births by default (D043).
- Default matching through `create_age_gap_matching_model(0, 5)`; optional
  `.n_eligible` and `.utility_rank` diagnostics on new union rows.
- Shared native kinship screening and `extract_kin()`. Lexis summaries optionally
  include overall/partnered/unpartnered exposure, births, deaths and formation
  rates, using `rate_fertility`, `rate_mortality` and `rate_marriage`.

See [design](design.md), [scientific contract](scientific-contract.md), and
[union formation](union-formation-models.md) for current rules and limitations.

## Validation and ownership

Milestone 1 was accepted by the user on October 15, 2026 (D029).
`validation/milestone1.qmd` is the accepted historical assessment. Reproduce it
at its recorded commit: current simulation requires a union process that was
not part of that assessment. [Testing](testing.md) retains the acceptance
artifact provenance and archive limitations.

The user is writing `validation/milestone2.qmd`. It now includes demographic
comparisons for both search modes, union diagnostics, age-gap comparisons,
conditional-logit parameter recovery, and performance comparisons including
R formula fallback. The current report uses USA circa 2000–2001 mortality,
fertility and marriage inputs. England and Wales 2015 preparation remains
available in `validation/data`, but is not the current report's empirical case.
The report's conclusions section is empty and no milestone 2 acceptance is
recorded. Preserve the user's report and explanatory source comments.

The October 19 documentation review inspected report source and implementation;
it did not rerender the reports, download rates, independently reproduce the
scientific findings, or set acceptance thresholds. Earlier implementation
verification counts and timings remain in Git history.

## Remaining and deferred work

The user confirmed on October 19 that milestone 2 is nearing completion. The
remaining sequence is to compile the report, obtain agent feedback, and then
write the conclusions before final user acceptance. Parameter recovery and
R formula-fallback timing are already present. Report conclusions and scientific
interpretation belong to the user.

Transition-aware fertility compensation (D045), poor-pool waiting (D044), and
previous-father continuity (D042) are deferred. The current adjustment does not
establish calibrated conditional ASFRs through status transitions or rate
changes between runs. Retain this limitation when interpreting report results.

The [roadmap](roadmap.md) records groups as the milestone 3 theme and proposed
non-death dissolution for milestone 4. Detailed group rules, parity conditioning,
calendar-varying schedules, arbitrary custom mechanisms, initialization
correction, exact checkpoints, and other future work are not implicitly
authorized. An African historical benchmark remains postponed.
