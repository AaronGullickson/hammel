# Experimental sequential two-sided matching

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
`exp(utility - max(utility))` among the remaining allowed alternatives.

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
   Ego plus alternates and the retained men form a square initial market.
5. Enforce hard group exclusions for every hypothetical pair. Do not apply
   kinship screening to alternate women. All retained men were screened for
   ego, so her row always has an allowed alternative while she is present.

Alternate women may themselves have children or prior ended unions. Their
only pool conditions are living, unpartnered, and the common age cutoff.
The pool updates at actual births, deaths, formation, and widowhood. Ages are
in continuous years; exactly the cutoff is eligible. Under-age women become
available at their cutoff birthday, without requiring a demographic event at
that birthday. Pool state belongs to each simulation call and is reconstructed
on restart.

## Temporary pairing

At each round:

1. Draw a fair independent coin for female or male choice.
2. Uniformly select one remaining person of that sex.
3. Draw a partner using conditional-logit probabilities over their remaining
   allowed alternatives. When a man chooses, score each candidate woman with
   female-first group orientation and male-minus-female age gaps.
4. If an alternate woman has no allowed men remaining, remove only her from
   the temporary market and continue; she remains actually unpartnered.
5. If ego is in the selected pair, stop and return the selected man. Otherwise
   remove both participants from the temporary market and continue.

There are no actual unions, parent links, removals, spacing changes, or event
rescheduling for hypothetical pairs. Only ego's returned match is used by the
existing union/father-selection event. Nonmarital father selection still does
not form a union or remove the father from actual availability.

If one man and ego remain, the same round/draw procedure is used; there is no
special deterministic last-pair shortcut. All sampling uses R's RNG within the
existing seed-isolation contract. Experimental histories need not match the
female-choice method. Explicit/default female choice must preserve existing
seeded histories and RNG consumption.

Common native linear/quadratic age scores and group offsets can be evaluated
on demand for the selected row/column without allocating a full matrix.
For a general R formula, pass only the selected person's remaining allowed
row/column alternatives to the R evaluator, once per choice. No full-market
information is exposed to the formula. Data-dependent transforms therefore
operate on that round's choice vector, just as in ordinary ego-only selection. Female-only terms can affect
male-initiated choices and must not be removed as row constants; the formula
intercept still cancels within every choice. Custom native/fallback parity
is tested for pointwise-equivalent formulas.

## Observations and examples

Existing `.n_eligible` reports ego's screened men before any reduction for a
female shortage. `.utility_rank` ranks the returned man among that original
eligible set using ego's original utilities. Diagnostics do not change draws.
No temporary market histories are added to public population tables.

- With one eligible man and no alternate women, ego selects that man.
- With two men and one alternate woman, either sex may choose. A hypothetical
  pair removes only its local participants; the remaining actual persons and
  their event schedules are unchanged.
- With three eligible men but only one alternate woman available, randomly
  retain two men and use the resulting two-by-two market.
- With a B ego, only B men eligible under strict endogamy, and an A alternate,
  the alternate may leave unmatched. Ego can still select a B man.

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

Worst-case score work is quadratic in market size; early termination at ego's
match reduces work in individual searches. Native score evaluation uses linear
scratch storage, as does the R formula fallback (with an R evaluation per
round that has allowed alternatives). Full-pool
sampling can be expensive. Measure overhead rather than infer it solely from
complexity.
