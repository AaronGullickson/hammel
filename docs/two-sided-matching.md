# Experimental fixed-pool two-sided matching

Authorized for implementation and investigation October 22, 2026 on
`experiment/two-sided-sequential-matching`, branched from `milestone3-groups`
at `4d03f55`. This experiment is not accepted demographic validation or a change
to the default matching mechanism. Preserve the user's milestone 3 report.

## Scope and inputs

`union_formation()` / `set_union_formation()` accept
`matching_method = "female_choice"` (existing default) or `"two_sided"`.
Older process objects lacking this field retain female choice. The method
applies through the shared selection routine to on-demand union formation,
rate-based union searches, and rate-mode nonmarital father selection. An explicit
`father_matching` still supplies the latter's preferences and hard restrictions.
Existing unions retain their partners. Initiation clocks, annual input hazards,
birth cancellation, spacing, and event invalidation remain unchanged.

All scores use the existing matching model: continuous male age minus female
age, fixed female/male group orientation, and additive directed/symmetric
pair effects. No conversion of utility coefficients into equilibrium margins
is performed. Finite utilities give conditional probabilities proportional to
`exp(utility - max(utility))` among the fixed allowed alternatives.

## Sampling and eligibility

1. Sample men from the existing living, unpartnered male pool before screening,
   using the existing `sample_size` (zero means the full pool).
2. Apply the existing minimum male age, ego kinship, and applicable hard group
   exclusions. If no men remain, use the existing failure behavior.
3. The alternate-woman pool is all other living, unpartnered women whose age
   is at least the same `min_male_age` cutoff. Do not weight by fertility,
   spacing, or search hazards. The actual ego woman retains her existing event
   eligibility; this does not impose a new minimum age on her.
4. With M eligible men, sample M - 1 alternate women without replacement.
   If there are too few women, uniformly reduce the eligible men to one plus
   the number of available alternate women, then sample that many women.
   Ego plus alternates and the retained men form a square fixed market.
5. Enforce hard group exclusions for every hypothetical pair. Do not apply
   kinship screening to alternate women. All retained men were screened for
   ego, so her row always has an allowed alternative.

Alternate women may themselves have children or prior ended unions. Their
only pool conditions are living, unpartnered, and the common age cutoff.
The pool updates at actual births, deaths, formation, and widowhood. Ages are
in continuous years; exactly the cutoff is eligible. Under-age women become
available at their cutoff birthday, without requiring a demographic event at
that birthday. Pool state belongs to each simulation call and is reconstructed
on restart.

## Fixed-pool probability calculation

The user authorized replacing sequential depletion on October 22, 2026.
Keep the sampled market fixed. Calculate each allowed pair's utility U[i,j],
with women as rows and men as columns. Hard exclusions have zero weight.
Normalize ego's row to obtain a[j] = P(ego chooses man j). Normalize each
column across women to obtain b[j] = P(man j chooses ego). Draw once with
probabilities (a[j] + b[j]) / sum(a + b).

This exactly integrates repeated rounds that flip a fair sex coin, uniformly
select an initiator, and draw a conditional-logit partner, stopping when ego
is involved. Non-ego pairs are discarded without removing anyone. Equal pool
sizes make the common initiator factor cancel. A person with no allowed
partners produces no pair in a hypothetical round; they remain in the pool.
Ego's eligible row ensures termination with probability one. No hypothetical
pairs are actually generated, and only one uniform partner draw is consumed.

Scores retain female-first group orientation and male-minus-female age gaps.
For custom R formulas, evaluate each woman's full allowed male vector once;
use those same fixed utilities for both normalizations. Data-dependent
transforms therefore operate within each woman's allowed row, never on a
changing vector or separately re-evaluated column. Female-only score terms
remain relevant to column probabilities. Row and column sums use stable
log-scale accumulation; no full matrix allocation is necessary.

No hypothetical outcome changes unions, parent links, availability, spacing,
or event schedules. Only ego's returned match changes actual state through
the existing event. Nonmarital father selection still creates no union.
Default female choice retains its existing seeded histories and RNG consumption;
the experimental method changes histories relative to the prior prototype.

Worked two-by-two example: pair weights [[1,3],[2,4]] give ego's row
probabilities [1/4,3/4] and men's probabilities of choosing ego [1/3,3/7].
Adding and normalizing gives partner probabilities [49/148,99/148].

## Observations and examples

Existing `.n_eligible` reports ego's screened men before any reduction for a
female shortage. `.utility_rank` ranks the returned man among that original
eligible set using ego's original utilities. Diagnostics do not change draws.
No temporary market histories are added to public population tables.

- With one eligible man and no alternate women, ego selects that man.
- With two men and one alternate woman, the combined probabilities reflect
  both choosing directions without generating hypothetical pairs.
- With three eligible men but only one alternate woman available, randomly
  retain two men and use the resulting two-by-two market.
- With a B ego, only B men eligible under strict endogamy, and an A alternate,
  the alternate has zero weight in those columns. Ego can still select a B man.

## Validation and benchmarking

Software checks cover analytical two-by-two probabilities, male-choice score
orientation, empty/single/shortage markets, no hypothetical state changes,
age eligibility and ageing-in, widowhood and restart pool reconstruction,
kinship for ego, group restrictions for all pairs, father overrides, native/R
score agreement, diagnostic invariance, input validation, reproducibility and
independent runs. Run the ordinary suite through `devtools::test()`.

Keep scientific assessment separate: use milestone 3's group pair tables,
age gaps, parameter-recovery analysis, and inheritance trajectories to evaluate
whether input coefficients retain useful benchmarks. No aggregate odds-ratio
recovery or inheritance balance is assumed. Benchmark both fixed-market
selection and complete simulations, with sample sizes, formulas, seeds,
population sizes, event counts, and timings recorded. No live rate downloads
are needed for software checks or synthetic benchmarks.

Run `validation/two_sided_matching_benchmark.R` from the repository root for
synthetic optimized-build timings. The measured October 22 results and raw
tables are recorded in
[the timing assessment](../validation/two_sided_matching_benchmark.md).

Score work is quadratic in market size, with linear scratch storage and one
R evaluation per nonempty female row for custom formulas. There is no iterative
matching or early termination. Historical sequential timings in the timing
assessment describe the superseded prototype; compare against fresh measurements.

## Native bulk calculation contract

Native linear/quadratic age-gap models may calculate whole matrix columns in
reusable C++ scratch arrays, with group offsets and hard exclusions applied
in female-first orientation. This is a software optimization of D062, not a
change to its probabilities, square sampling, eligibility, or single draw.
No new demographic rule or dependency is introduced. The R formula path
retains one evaluation per female row, including contextual transforms.

Validation compares the bulk native kernel against the original row-streamed
calculation and explicitly materialized matrices, including extreme scores,
restrictions, directed effects, and RNG state. Column summation may differ
at floating-point roundoff; probability agreement and sampled fixed-seed
histories are checked. Benchmarks use clean optimized builds and report
market size, full-simulation event counts, and workload limitations.

The measured optimization and explicit vectorized R comparison are recorded
in [native performance](../validation/native_two_sided_performance.md).
