# Experimental weighted-choice matching

Specified with the user and authorized for implementation October 23, 2026 on
`experiment/two-sided-sequential-matching` (D066). This describes implemented
behavior of `matching_method = "weighted_choice"`. Scientific assessment is
pending; no acceptance thresholds are set. Additive two-sided matching remains
available and remains the default for new specifications (D064).

## Problem

The user requires partner matching that (1) preserves age-specific female
union timing from the inputs, (2) recovers the input matching preferences, and
(3) yields equal female-A/male-B and female-B/male-A union counts when
preferences are symmetric and each group's sexes are balanced.

| Method | Preferences recovered | Directional balance |
| --- | --- | --- |
| Female choice | Yes | No: 75/25 groups, A::B = -3 give AB/BA = 0.377 |
| Additive two-sided (D062) | No: same case gives cross-product -5.77, not -6 | Approximately |
| Sinkhorn balancing of each temporary market | No: age and group effects attenuated | No: overcorrected |

The female-choice imbalance arises because nothing specifies how often each
type of man forms a union. Women's union timing is an input; men's is whatever
female choice leaves. With 75/25 groups, 1.6% of A women but 13% of B women
choose out-group men, so B men form about 8% fewer unions than B women.
Weighted choice adds the missing rule for men without altering female choice's
preference structure.

The 75/25 figures are exact group-only calculations. The Sinkhorn and sampling
findings come from an exploratory static-snapshot analysis on October 23, 2026
(synthetic mortality fixture; scripts in a session scratch directory, not
archived). The user clarified that `sample_size` is a speed shortcut, not an
encounter model; at 250 men, sampling did not materially attenuate
female-choice recovery in that analysis.

## Selection rule

Union initiation, timing, candidate sampling, minimum male age, kinship
screening, hard group exclusions, failure, cancellation, and spacing rules are
unchanged. Only the probability of selecting among ego's screened candidates
changes. For ego woman i with screened candidate set C,

    P(j | i, C) = b[g(j)] exp(U[i,j]) / sum over k in C of b[g(k)] exp(U[i,k])

U is the existing matching utility (formula plus pair effects). b[g] > 0 is a
weight for the man's group, recalculated at each search from the current
population state. Only ratios matter; the first group with a positive target
has weight 1. With one group, or no recorded groups, b = 1: the method then
reproduces explicit female choice exactly, including random-number use.

Because b depends only on the man's group, it multiplies whole columns of the
expected pairing table. It changes how often each male group forms unions, not
the association between partners' characteristics. The age model, including
`create_age_gap_matching_model()`, keeps its current meaning.

## Calibrating the weights (continuous time)

At each search time t, with no execution step and no random draws:

1. **Female intensity.** The calibration needs how many women of each group
   are about to need a partner now. For each living, unpartnered woman w at or
   above `min_male_age`, lambda[w] is her current chance per year of
   triggering a partner search:
   - union searches in rate mode: her search hazard at her current age, after
     group multipliers;
   - on-demand union formation: her adjusted unpartnered fertility hazard,
     because a birth triggers the search;
   - nonmarital father selection in rate mode: the same fertility hazard.
   For fertility-driven intensities, a woman within her birth-spacing interval
   contributes zero. Failed searches and cancelled births need no special
   treatment: intensity counts attempts, whatever their outcome.
2. **Male target.** Let N_F[g] and N_M[g] count living, unpartnered women and
   men of group g at or above `min_male_age`, and r[g] be group g's summed
   female intensity divided by N_F[g]. The mirrored share is
   M[g] = N_M[g] r[g] / sum over h of N_M[h] r[h]: each unpartnered man has the
   per-capita intensity of an unpartnered woman of his group. If a group has
   no counted women, its men use the population-wide female per-capita
   intensity. F[g] is group g's share of female intensity among groups with
   counted men. The target is
   tau[g] = F[g] + kappa (M[g] - F[g]),
   where kappa is the matching model's `imbalance_absorption` (below). Women's
   timing still determines the total number of unions; the target only
   distributes them across male groups.
3. **Expected male flows.** For woman w, S[w,g] sums exp(U[w,j]) over counted
   men j of group g, with hard-excluded pairs contributing zero. Her expected
   choice share is pi[w,g] = b[g] S[w,g] / sum over h of b[h] S[w,h], and the
   expected male flow share is phi[g] = sum over w of lambda[w] pi[w,g] divided
   by the summed intensity. Candidate subsampling and kinship screening are
   ignored here.
4. **Solve** phi(b) = tau for the free log weights by Newton iterations on a
   convex objective, warm-started from the previous search's weights, until
   every share is within 1e-10. Expected flow shares rise monotonically with a
   group's own weight, so a solution is unique when it exists.

Crude per-capita counts were selected by the user. Age standardization was not:
it would require a male age-specific intensity (a woman's rate at the same or
a shifted age), reintroducing the male age assumption and preferred-gap
question that group-only weights avoid.

### Imbalance absorption

`imbalance_absorption` (kappa) belongs to each matching model, like
`sample_size`, defaulting to 1 (user, October 23, 2026). It sets the share of a
group's sex imbalance among unpartnered people absorbed through cross-group
unions. At 1, men receive the full mirrored target. At 0, each male group's
share of unions equals its female group's share, so the excess sex waits
longer. With symmetric preferences and balanced sex composition, F = M and
every value gives equal directional counts. With no group preferences, values
below 1 still disfavor men of an over-represented group. Other matching
methods ignore the argument; older matching objects without it use 1.

### Nonmarital fathers

Rate-mode father selection has its own weights, calibrated from fertility
intensity with the applicable father utilities (`father_matching`, or
`matching` when NULL) and that model's `imbalance_absorption` (user, October
23, 2026). Union searches keep weights calibrated from search hazards. The
male pool, target rule, and solver are shared. On-demand mode has one
calibration, driven by fertility.

### Numerical representation

The engine keeps counts of unpartnered men and women by group and 0.25-year
birth-time bin, and of women past their birth-spacing interval. Counts change
only at events (birth, death, formation, widowhood) and when a queued spacing
interval ends; restarts rebuild them from histories. At a search, each bin's
age is the current time minus its centre birth time, which determines the age
cutoff and hazards.

Pair utilities depend only on the lag between bin centres, which equals male
age minus female age. A table over every lag in the run's birth-time range is
built once per matching model, including one evaluation of a custom R formula
on that grid. The formula must therefore give finite values at every such lag.
Bin width is a numerical tolerance, not a demographic parameter. Candidate
selection itself uses exact ages.

### When no solution exists

A male group may have a positive target but no feasible demand, for example
under `prevent_exogamy` when no compatible women have positive intensity. Its
log weight is set to 50 and the search is counted as approximate; Newton
non-convergence within 100 iterations is counted likewise. A group with zero
target receives log weight -1000. If no women have positive intensity, all
weights are 1. As selected by the user, `simulate()` warns once, reporting the
number of searches that used approximate weights.

## Worked examples

Exact group-only expectations with equal per-woman intensities, zero
within-group utility, and absorption 1 unless stated:

| Case | b_B | AA | AB | BA | BB | Log cross-product | AB/BA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 75/25, A::B = -3, female choice | 1 | .73776 | .01224 | .03249 | .21751 | -6 | 0.377 |
| 75/25, A::B = -3 | 1.683 | .72962 | .02038 | .02038 | .22962 | -6 | 1.000 |
| 50/50, A::B = -3 | 1.000 | .47629 | .02371 | .02371 | .47629 | -6 | 1.000 |
| 75/25, A:B = -1, B:A = -2 | 0.924 | .67368 | .07632 | .07632 | .17368 | -3 | 1.000 |
| 75/25, A::B = -1.5 | 1.523 | .67368 | .07632 | .07632 | .17368 | -3 | 1.000 |
| Women 75/25, men 70/30, A::B = -3 | 3.875 | .69272 | .05728 | .00728 | .24272 | -6 | 7.87 |
| Same, absorption 0.5 | 2.400 | | .0365 | .0115 | | -6 | 3.17 |
| Same, absorption 0 | 1.309 | | .0204 | .0204 | | -6 | 1.00 |

Consequences:

- **Directed effects in two groups.** A:B = -1 with B:A = -2 produces exactly
  the same expected table as A::B = -1.5. With both sexes' group flows fixed,
  a two-group table has one free quantity, its cross-product. Directed effects
  change which weights are needed, not the flows. With more groups, only
  combinations not absorbed by the group flows affect pairing.
- **Sensitivity to composition.** AB - BA equals the difference between the
  male and female group-B shares of unions. When exogamy is rare, small
  differences in the sexes' unpartnered group shares move AB/BA strongly.
  Directional balance should be assessed with those pool compositions.
- **Pooled cross-products.** These examples have no ages. With age preferences
  and group-specific age structures, a cross-tabulation pooled over ages need
  not equal the summed group utilities exactly, under any method.

## Diagnostics and interpreting estimates

`.n_eligible` and `.utility_rank` keep their meanings: rank uses U, not the
weighted score. Weights are not recorded and consume no random draws.

A conditional logit fitted to wife-side choice sets with the input terms only
estimates separate directed group terms offset by log(b) differences, while
their sum and the age terms are unaffected. Adding a male-group main effect
absorbs the weights. In the static analysis, calibrated weights left the age
term (-0.0278 against -0.0277) and symmetric group term (-2.87 to -2.91,
against -2.92 to -2.96 for female choice) essentially unchanged.

## Validation

`tests/testthat/test-weighted-choice.R` covers the analytical examples above
(including absorption 0, 0.5 and 1, and balanced symmetry at every absorption),
binned calibration against an exact individual-level R reference for native
and custom formulas in both intensity types, three-group convergence, zero and
unattainable targets, exact equality with female choice for one or no groups
in both initiation modes, reproducibility, restarts, separate father models,
diagnostics without changed histories, the single warning, and input
validation. Timing has not yet been measured.

Scientific assessment (user): the milestone 3 group-matching scenarios, with
directional counts, unpartnered pool composition by sex and group over time,
cross-tabulations, conditional logit with and without partner-group effects,
and realized age gaps.
