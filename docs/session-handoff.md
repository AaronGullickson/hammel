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

The user accepted milestone 2 as complete on October 19, 2026 (D050), after
review of `validation/milestone2.qmd` and its rendered HTML and user revisions.
The report includes demographic comparisons for both search modes, union
diagnostics, age-gap comparisons, conditional-logit parameter recovery, and
performance comparisons including R formula fallback. It uses USA circa
2000–2001 mortality, fertility and marriage inputs. England and Wales 2015
preparation remains available as an additional example. Preserve the user's
report and explanatory source comments.

Report review inspected source, rendered tables and figures; the acceptance
documentation update did not rerender reports, download rates, or independently
rerun simulations or tests. No new scientific acceptance thresholds are set.
The report and documentation are to be committed together; this record does
not assign a final artifact hash or acceptance commit.

## Deferred work and next scope

Milestones 1 and 2 are complete by user acceptance. Report findings and
scientific interpretation belong to the user. Completion does not imply every
review suggestion was implemented or remove documented limitations.

Transition-aware fertility compensation (D045), poor-pool waiting (D044), and
previous-father continuity (D042) are deferred. The current adjustment does not
establish calibrated conditional ASFRs through status transitions or rate
changes between runs. Retain this limitation when interpreting report results.

The [roadmap](roadmap.md) records groups as the milestone 3 theme and proposed
non-death dissolution for milestone 4. Detailed group rules, parity conditioning,
calendar-varying schedules, arbitrary custom mechanisms, initialization
correction, exact checkpoints, and other future work are not implicitly
authorized. An African historical benchmark remains postponed.
