# Current project handoff

Updated October 29, 2026. This is a dated status record, not authorization for
additional features. The current user request controls work. Start with the
[documentation index](README.md) and [scientific contract](scientific-contract.md).

## Accepted milestones

The user accepted milestones 1–4 on October 15, October 19, October 27 and
October 29, 2026 (D029, D050, D077, D093). Preserve all four report sources,
rendered assessments, rate provenance and the author's explanatory comments.
Milestone 1 requires its recorded historical commit because current execution
requires a union process. The final milestone 4 report was rendering when the
user requested documentation closeout; do not interrupt that render or delete
its temporary files. No new scientific analysis is implied by closeout.

## Implemented state

- Independent C++17/Rcpp continuous-event simulations, isolated R RNG policy,
  tibble population/union histories, shared known-kin screening and extraction.
- Mortality and fertility with annual hazards and midpoint birth-spacing
  compensation. Actual birth spacing persists through union transitions.
- On-demand or female rate-based partnership searches, exact full-pool matching,
  group schedules/multipliers, calibrated group weights and optional delta.
  Current female hazards are summed by group/cohort; a common boundary queue
  handles age, duration, spacing and cohabitation probability changes.
- Duration-specific dissolution with one female-carried clock, symmetric group
  pairs/multipliers, optional cohab/marriage schedules and retained end reasons.
- Searches by previous_union none/dissolved/widowed, requiring all categories
  when stratified. None uses age; each other category may use age or duration
  since the latest actual ending. State is reconstructed from complete histories.
- Rate-mode cohabitation: cohab_prob by woman's age and optional group/previous
  category; optional separate cohab_matching; age/duration union_convert_rate
  optionally by woman's group. Conversion creates same-partner linked episodes
  with previous_cohab_id and resets dissolution duration. Exact ties prioritize
  death, birth, dissolution, conversion, then search.
- Fertility input union_type accepts unpartnered/partnered or
  unpartnered/cohab/marriage; hazard-list names match supplied types. Father
  matching at unpartnered births can use a separate model. Failed matches cancel
  births by default. Cohabitors remain unavailable for new matching.
- Cohabitation endings recover earlier legal status by default; the configured
  toggle can count them as dissolved/widowed. Latest simulation metadata carries
  the policy for observation. Changing policy reclassifies supplied histories.
- Newborn identities and ancestry are separate; random_parent inheritance is
  default with maternal fallback. Custom inheritance is supported.

## Observation APIs

- calculate_fertility_rates() and calculate_mortality_rates(): separate process
  summaries; by_union_type = TRUE returns unpartnered plus discovered types,
  without overall rows. Fertility accepts parity_breaks for current-parity ASFRs.
- calculate_union_rates(): partnership entries over unpartnered years, excluding
  conversion. Optional by_union_type repeats shared unpartnered exposure;
  by_previous_union and optional dur_breaks retain mixed age/duration observations.
- calculate_marriage_rates(): direct marriages plus conversions over nonmarried
  (unpartnered plus cohabiting) years; optional previous categories, age only.
- calculate_dissolution_rates(): dissolutions per union-year, optionally by
  symmetric group_pair and union_type.
- calculate_cohab_convert_rate(): conversions per cohabitation union-year,
  either woman's age or cohabitation duration, optionally woman's group.
- get_unpartnered_births(): known-parent births outside partnerships, enriched
  with parent ages/groups. Preserve repeated pairs; use to assess father choice.

## Validation and cleanup

The accepted [milestone 4 report](../validation/milestone4.qmd) assesses rate
recovery, separate matching contexts and performance over independently seeded
baselines. Its final render includes the test suite. The earlier reviewed HTML
had 321 passing tests and one opt-in live-database skip; final-render results
remain the author's record. No new tests/downloads were run for docs cleanup.

Reports remain at validation/ top level; the 14 supplementary benchmark/recovery
files moved to [validation/milestone4/](../validation/milestone4/README.md), with
links and script paths updated. Rate/source inputs remain in validation/data/.
The author deleted milestone4_tests.R; the final report supersedes it. Quarto
intermediates are ignored and disposable only after rendering finishes.

## Retained limitations and future scope

Transition-aware fertility compensation (D045), poor-pool waiting (D044), former-
partner/father continuity (D042), sparse exposure, initialization transients and
approximate calibration ages remain limitations. Pooled cross-tab log odds are
not universal exact constraints; conditional-choice recovery is separate evidence.
Static SOCSIM audits are not executable reference comparisons.

The roadmap proposes parity (milestone 5), calendar-time schedules (milestone 6),
and later group transitions, summaries, migration and characteristics. The user
has discussed adding calendar bounds alongside age/duration, but this does not
authorize implementation or resolve gap/coverage and fertility-compensation rules.
Fixed newborn assignment remains deferred to group transitions. Marriage-order
effects, richer couple covariates, general two-sex models, API restructuring,
arbitrary callbacks and exact checkpoints remain future scope. Acceptance does
not authorize those mechanisms. An African historical benchmark remains postponed.


## Milestone 5 current work

The user confirmed the milestone 4 PR was merged and authorized parity work.
D094 selects half-open parity_start/parity_end, including Inf. HFD conditional
mi rates map birth-order 5+ (m5px) to current parity 4+. D095 retains 12- as
[12,13) and 55+ as [55,Inf), supports upper-open fertility ages with interior-
limit spacing adjustment, and warns when TFR drops nonzero open-age rates.
D096 zero-fills unavailable HFD conditional cells with policy/count provenance.
D097 adds calculate_patfr(), with automatic calculate_tfr() delegation and all
existing strata, missing included cells yielding NA and upper-age truncation.

D098/D099 implement actual parity counts from recorded maternal births only.
No optional persons parity column overrides history; deceased children and
unknown fathers count. Initialize random baselines at zero if there are no
recorded births. Increment only on successful births; reconstruct at restart.
Complete parity partitions are required for every group/union-type schedule.
Parity zero is unadjusted; higher parity uses the preceding actual parity's raw
age schedule for the selected SOCSIM-style compensation. Entry to the open
highest category differs from later births within it. Refresh father-matching
fertility contributions and scheduling after births, with no parity matching bins.

**Important limitation:** the offline initialized-cohort check at
[validation/milestone5/](../validation/milestone5/README.md) finds mean parity-1
and 2+ rates .2665 and .3335 with .75-year spacing against inputs .3 and .4.
Zero-spacing controls recover approximately .2992 and .4035. This is the selected
approximation's limitation, not a scientific acceptance claim. Do not silently
change to population-ratio compensation. D100 adds `calculate_fertility_rates(parity_breaks = ...)` and separates
`calculate_mortality_rates()`. Fertility heterogeneity is not implemented. Preserve the user's
validation/milestone5.qmd, which is outside these implementation changes.


## Current investigation for collaboration (October 29)

The user confirmed persistent high-parity/young-age shortfalls after 50 years
of burn-in, then close schedule recovery with `birth_interval = 0`. The user
requested an executable SOCSIM toy comparison and documentation for collaboration
with Claude. See [the investigation](parity-spacing-investigation.md) and
[reference experiment](../validation/milestone5/socsim/README.md). No new
compensation method is accepted. Investigate a shared eligibility framework
before choosing between pre-calibration and runtime accounting. Preserve the
user-authored milestone 5 report and do not implement a replacement on the basis
of the agent's proposals alone.

Later on October 29 the user agreed to pursue a general runtime mechanism
instead of a parity-only correction and requested its specification:
[cell-accounted fertility](fertility-cell-accounting.md). Births are generated
per fertility cell at m_c * N_c and assigned to an eligible woman, with
reserved individual weights for future heterogeneity. Its six open decisions
need the user's answers before any engine change; D028/D099 remain current.
The user then agreed to all six decisions (D101): drop births in cells with no
eligible women and report them; individual weights set selection probability
within a cell rather than multiplying rates. Next step is the opt-in engine
implementation sharing the matching-calibration tracker; D028/D099 remain the
current behavior until it is validated.

Implementation status (October 29): cell accounting is implemented as the
opt-in `spacing_method = "cell"` with tests, internal scan checks, offline
validation and a benchmark. After validation the user selected 1-year age cells
split by actual parity (revised decision 5); young-age parity rates now recover
except structurally infeasible cells (for example parity 5+ at age 18), whose
intensity is reported as unrealized. Midpoint remains the default and its
histories are unchanged. Making cell mode the default needs the user's
decision and re-validation of earlier milestones.
