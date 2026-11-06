# Weighted-choice matching

Weighted choice is the sole matching mechanism (D075, October 26, 2026).
There is no `matching_method` argument. Both partner and nonmarital father
selection use calibrated male-group weights; separate father models retain
independent preferences, restrictions and delta adjustments. Group pair effects
use a single colon (`A:B`) for symmetric preferences. This describes implemented
behavior, not scientific acceptance thresholds.

## Problem

The user requires partner matching that (1) preserves age-specific female
union timing from the inputs, (2) recovers the input matching preferences, and
(3) yields equal female-A/male-B and female-B/male-A union counts when
preferences are symmetric, delta is zero, and each group's sexes are balanced.

| Method | Preferences recovered | Directional balance |
| --- | --- | --- |
| Female choice | Yes | No: 75/25 groups, A:B = -3 give AB/BA = 0.377 |
| Additive two-sided (D062; removed in D067) | No: same case gives cross-product -5.77, not -6 | Approximately |
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

Union initiation, timing, minimum male age, kinship screening, hard group
exclusions, failure, cancellation, and spacing rules are unchanged. Only the
probability of selecting among ego's eligible men changes. For ego woman i
with eligible pool C,

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
   - union searches in rate mode: her search hazard at her current age and, when supplied,
     `previous_union` category, after group multipliers;
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
   intensity. The target is
   tau[g] = M[g] exp(delta[g]) / sum over h of M[h] exp(delta[h]),
   where `delta` is the matching model's named vector of target adjustments (below). Women's
   timing still determines the total number of unions; the target only
   distributes them across male groups.
3. **Expected male flows.** For woman w, S[w,g] sums exp(U[w,j]) over counted
   men j of group g, with hard-excluded pairs contributing zero. Her expected
   choice share is pi[w,g] = b[g] S[w,g] / sum over h of b[h] S[w,h], and the
   expected male flow share is phi[g] = sum over w of lambda[w] pi[w,g] divided
   by the summed intensity. Kinship screening is ignored here.
4. **Solve** phi(b) = tau for the free log weights by Newton iterations on a
   convex objective, warm-started from the previous search's weights, until
   every share is within 1e-10. Expected flow shares rise monotonically with a
   group's own weight, so a solution is unique when it exists.

The engine now maintains actual current hazard contributions per unpartnered
woman (D082). These are summed into group/birth-cohort cells before solving;
previous union and other rate predictors do not require separate matching
rows. Separate marriage and fertility-driven totals share a boundary-update
queue. Female counts and mirrored male targets keep their existing definition;
there are no male previous-category targets. See
[current-hazard contributions](calibration-intensity.md) for timing and state
updates and [remarriage](remarriage.md) for schedule selection.

### Objective and uniqueness

With Lambda the summed intensity and the first group's log weight fixed at 0,
the engine minimizes

    L(log b) = (1/Lambda) sum over w of lambda[w] log sum over h of b[h] S[w,h]
               - sum over g of tau[g] log b[g]

(`WeightCalibrator::solve` in `src/simulate.cpp`). Its gradient with respect
to log b[g] is phi[g] - tau[g], so the stationary point is exactly the target
condition phi(b) = tau. The first term is an intensity-weighted sum of
log-sum-exp functions of log b, which is convex, and the second is linear, so
the minimum is unique when it exists. Newton steps use the analytical gradient
and Hessian.

### Why the weights leave the association unchanged

Because b depends only on the man's group, it multiplies whole columns of the
expected pairing table. In any two-by-two cross-product of female groups g, h
and male groups k, l, each woman's normalizing denominator and each column
weight appear once in the numerator and once in the denominator. For the
group-only case below,

    (AB x BA) / (AA x BB) = exp(U[A,B] + U[B,A] - U[A,A] - U[B,B])

whatever b and the group sizes are. The weights therefore change how often
each male group forms unions, not the association between partners'
characteristics, and wife-side conditional-logit estimates of the input terms
use calibrated log group weights as offsets when reconstructing choice sets. When unweighted flows already meet their targets, b = 1; one or no
recorded groups always use unit weights.

Crude per-capita counts were selected by the user. Age standardization was not:
it would require a male age-specific intensity (a woman's rate at the same or
a shifted age), reintroducing the male age assumption and preferred-gap
question that group-only weights avoid.

### Male target adjustments

`delta` is an optional named finite numeric vector on each matching model.
Omitted groups have delta 0. Targets are proportional to mirrored male flows
multiplied by exp(delta); normalization preserves a total target of one.
Only differences matter: adding the same constant to all group adjustments
changes nothing. Positive values increase a group's male target relative to
baseline and other groups; negative values decrease it. This is not an
individual utility term, a pair-specific ratio, or an adjustment to women's
search rates. Named groups must belong to
the simulation's identity registry, including declared future identities.
With no identities, nonempty delta is an error. Zero baseline targets remain
zero; no available men or no positive search rates are not repaired by delta.
Strict endogamy and other restrictions may make targets unattainable; existing
approximate-weight diagnostics and warnings remain applicable.

For a static 75/25 population of each sex with symmetric exogamy -3,
delta A = 0 and delta B = 0.12274682836 produce expected shares AA =
0.71443059512, AB = 0.03556940488, BA = 0.01185646829, BB = 0.23814353171.
B men receive 0.27371293659 of unions; AB/BA = 3 and the log cross-product
remains -6. This is an analytical illustration, not simulation evidence.
Changing composition or preferences changes the ratio associated with delta.
Delta can also offset composition imbalance. In a static group-only case with
women 75/25 and men 70/30, delta_A = 0 and delta_B = log(7/9), approximately
-0.2513144283, shift B men's target from 0.30 to 0.25. Symmetric preferences
then give AB = BA. This matches female union shares at that composition;
a fixed delta need not retain symmetry as composition or search rates change.
The engine evaluates the normalized targets in log space to avoid overflow.

### Nonmarital fathers

Rate-mode father selection has its own weights, calibrated from fertility
intensity with the applicable father utilities (`father_matching`, or
`matching` when NULL) and that model's `delta` (D074). Both use calibrated
weights. Union searches keep weights calibrated from search hazards. The
male pool, target rule, and solver are shared. On-demand mode has one
calibration, driven by fertility.

### Numerical representation

The engine keeps counts of unpartnered men and women by group and 0.25-year
birth-time bin. It maintains each woman's actual current search and eligible
fertility hazards, summed into those cells. A boundary queue updates hazards at
age-band starts/ends and birth-spacing releases; pool and state changes replace
or remove contributions and pending entries. There is at most one pending
entry per counted woman. Updates consume no randomness and do not reschedule
individual demographic events. Restarts rebuild the state from histories.
At each search, bin-centre age still determines whole-bin cutoff eligibility
for both sexes. Female hazard contributions use actual ages, including within
narrow bands and gaps, replacing the previous full-bin hazard averaging (D082).

Pair utilities depend only on the lag between bin centres, which equals male
age minus female age. A table over every lag in the run's birth-time range is
built once per matching model, including one evaluation of a custom R formula
on that grid. The formula must therefore give finite values at every such lag.
Bin width is a numerical tolerance, not a demographic parameter. Candidate
selection itself uses exact ages. Remaining binning approximations are the
whole-bin age cutoff and pair utilities at bin-centre ages (gap error up to
about 0.25 years, which matters only for very sharp age preferences). In the
earlier age-band-averaged implementation, with one-year rate bands and a 5-year half-weight distance, default bins gave B
log weights within 0.0001 to 0.0003 of an exact person-level calculation;
0.01-year bins agreed to within 0.00001 (October 23, 2026).

The pool sums S are maintained incrementally for every female bin, group pair
and matching model: a man joining or leaving the eligible pool updates them,
and a male birth cohort is added once when its bin reaches the age cutoff. A
search then reads its rows directly instead of summing over all male bins.
Entries below 1e-6 are recomputed exactly at the search, and all sums are
rebuilt from the counts after every 100,000 updates, bounding accumulated
roundoff. This is a software optimization (October 23, 2026): probabilities
are unchanged up to floating-point roundoff, which can in principle alter
fixed-seed histories relative to the first implementation.

### Strict endogamy

When a matching model prohibits all exogamy (`group_exogamy = NA` for its
dimension, D135; formerly `prevent_exogamy = TRUE`), every woman sees only men
of her own group, so a group weight scales all her candidates equally and
cannot change her choice. Calibration is skipped for that model, giving
female choice exactly, with no warnings (D070). Partial restrictions (listed
pairs) still calibrate.

### When no solution exists

A male group may have a positive target but no feasible demand, for example
under prohibitions when no compatible women have positive intensity. Its
log weight is set to 50 and the search is counted as approximate; Newton
non-convergence within 100 iterations is counted likewise. After solving,
expected flows are recomputed in log space with every woman who has an
allowed man, including those whose only allowed men carry near-zero weights;
any flow missing its target by more than 1e-8 also marks the search
approximate (D071). Women with no allowed man are left out, since no weights
can match them. A group with zero
target receives log weight -1000. If no women have positive intensity, all
weights are 1. As selected by the user, `simulate()` warns once, reporting the
number of searches that used approximate weights.

## Worked examples

Unweighted female choice below is a mathematical comparison, not an API option.

Exact group-only expectations with equal per-woman intensities, zero
within-group utility, and zero delta unless stated:

| Case | b_B | AA | AB | BA | BB | Log cross-product | AB/BA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 75/25, A:B = -3, female choice | 1 | .73776 | .01224 | .03249 | .21751 | -6 | 0.377 |
| 75/25, A:B = -3 | 1.683 | .72962 | .02038 | .02038 | .22962 | -6 | 1.000 |
| 50/50, A:B = -3 | 1.000 | .47629 | .02371 | .02371 | .47629 | -6 | 1.000 |
| 75/25, A:B = -1.5 | 1.523 | .67368 | .07632 | .07632 | .17368 | -3 | 1.000 |
| 75/25, A:B = -3, delta_B = .122747 | 3.000 | .71443 | .03557 | .01186 | .23814 | -6 | 3.000 |
| Women 75/25, men 70/30, delta_B = log(7/9), A:B = -3 | 1.309 | .72962 | .02038 | .02038 | .22962 | -6 | 1.000 |
| Women 75/25, men 70/30, A:B = -3 | 3.875 | .69272 | .05728 | .00728 | .24272 | -6 | 7.87 |

These are analytical expectations, not simulation output. With equal
per-woman intensities and no ages, every woman of a group has the same choice
probabilities, so S[w,g] reduces to group g's male share times exp(U). With
b_A = 1 and b_B = b, female shares F_A, F_B, male shares M_A, M_B and pair
utilities U:

    pi_AB = b M_B exp(U[A,B]) / (M_A exp(U[A,A]) + b M_B exp(U[A,B]))
    pi_BA = M_A exp(U[B,A]) / (M_A exp(U[B,A]) + b M_B exp(U[B,B]))
    AA = F_A (1 - pi_AB)    AB = F_A pi_AB
    BA = F_B pi_BA          BB = F_B (1 - pi_BA)

Female choice sets b = 1. Weighted choice solves AB + BB = tau_B, with
tau_B = M_B exp(delta_B) / (M_A exp(delta_A) + M_B exp(delta_B)). For 75/25 groups and A:B = -3: at b = 1,
AB = 0.01224 and BA = 0.03249, so group-B men receive AB + BB = 0.2298 of
unions while group-B women form 0.25; scarce B women often choose A men, few
A women choose B men. Solving for tau_B = 0.25 gives b = 1.683, which raises
AB to 0.02038 and lowers BA to the same value. Equal groups give b = 1 by
symmetry. With women 75/25 and men 70/30, zero delta targets B men at
0.30. With balanced 75/25 composition and delta_B = 0.12274682836,
the B target is 0.27371293659 and AB/BA = 3 (see target adjustments above).
`validation/wc_table.R` reproduces these cases with a one-dimensional root
solve. Tests check default weights, signed delta, a common shift, extreme
adjustments, and the analytical 3:1 case.

Consequences:

- **Formula group effects.** A `male_group` main effect is absorbed by
  calibrated weights and does not change group flows; interactions can be
  partly absorbed. Such formulas warn once at construction. Use `delta` to
  adjust male-group targets. Pair effects are symmetric, specified once per
  unordered pair with `A:B`; directional utility entries and `::` are errors.
- **Sensitivity to composition.** AB - BA equals the difference between the
  male and female group-B shares of unions. When exogamy is rare, small
  differences in the sexes' unpartnered group shares move AB/BA strongly.
  Directional balance should be assessed with those pool compositions.
- **Directional tuning.** Delta adjusts total male-group targets and can
  tune directional ratios or offset sex imbalance. It does not directly set
  pair-specific ratios, which change as population composition changes.
- **Pooled cross-products.** These examples have no ages. With age preferences
  and group-specific age structures, a cross-tabulation pooled over ages need
  not equal the summed group utilities exactly, under any method.

## Diagnostics and interpreting estimates

Weights consume no random draws. In a grouped population, `diagnostics = TRUE` records `.log_weight_<group>` on each new union: the log group weights applied
at the search that formed it (October 23, 2026, D070). One group is the
reference with log weight 0; only differences between groups matter, and
searches where calibration cannot affect choice record 0 for every group. Earlier unions and new groups'
columns for earlier rows are NA. In a wife-side choice set, using the column
for each alternative's group as an offset gives the exact benchmark, in place
of approximating the weights with a male-group main effect. Nonmarital father
selections form no union and are not recorded.

A conditional logit fitted to reconstructed wife-side choice sets should
include the recorded log group weights as offsets. A male-group main effect
can approximate varying calibration weights but cannot reproduce each search
exactly. In the static analysis, calibrated weights left the age
term (-0.0278 against -0.0277) and symmetric group term (-2.87 to -2.91,
against -2.92 to -2.96 for female choice) essentially unchanged.

Choice sets reconstructed from the wife's side match the simulated choice
process and are the appropriate recovery benchmark. Husband-side sets treat
each husband as choosing among unpartnered women, but men do not choose in
this model: which women form unions, and when, follows women's hazards, which
those sets read as male age preferences. Weights apply only to men's groups,
so men's ages at union formation remain an outcome, not an input.

## User-reported assessment (October 23, 2026)

The user ran these replicated comparisons in the milestone 3 report's uneven
two-group scenario (on-demand formation, `create_age_gap_matching_model(0, 5)`
with A:B = -3, 100 years). The 250-man runs of both methods and the
1,000-man female-choice run used the same 20 seeds; the 1,000-man
weighted-choice runs used 40 replications. Values are means across replications,
with standard errors from between-replication variation. These are the user's
results and assessment; they were not independently rerun for this document.

| Method, sample size | Pooled cross-product (input -6) | AB/BA |
| --- | --- | --- |
| Female choice, 250 | -5.87 (0.017) | 0.596 (0.008) |
| Weighted choice, 250 | -5.88 (0.015) | 1.01 (0.039) |
| Female choice, 1,000 | -5.98 (0.017) | 0.615 (0.008) |
| Weighted choice, 1,000 | -5.95 (0.0095) | 0.983 (0.023) |
| Weighted choice, exact draws (D072) | -5.96 (0.0101) | 0.984 (0.021) |

Wife-side conditional logit (inputs: linear 0, quadratic -0.0277, symmetric
group -3):

| Method, sample size, terms | Linear | Quadratic | ABsym | Male group B |
| --- | --- | --- | --- | --- |
| Female choice, 250 | -0.0030 (0.0027) | -0.0275 (0.0005) | -3.03 (0.060) | |
| Weighted choice, 250 | 0.0023 (0.0028) | -0.0262 (0.0005) | -3.02 (0.056) | |
| Weighted choice, 250, + male group | -0.0025 (0.0039) | -0.0265 (0.0006) | -3.06 (0.064) | 0.322 (0.043) |
| Weighted choice, 1,000, + male group | 0.0013 (0.0019) | -0.0277 (0.0004) | -3.02 (0.037) | 0.311 (0.038) |
| Weighted choice, 1,000 | 0.0005 (0.0019) | -0.0274 (0.0004) | -3.03 (0.039) | |
| Same, new choice-set draws | 0.0011 (0.0020) | -0.0283 (0.0003) | -3.00 (0.037) | |
| Weighted choice, exact draws (D072) | -0.0002 (0.0019) | -0.0278 (0.0003) | -3.06 (0.034) | |

With weighted choice at 250 men, husband-side sets alone gave a linear term
of -0.0383 (0.0030), quadratic -0.0263 (0.0005), and directed terms AB -3.25
and BA -2.57; combining husband- and wife-side sets gave linear -0.0233
(0.0022), quadratic -0.0267 (0.0003), and ABsym -3.01 (0.031).

Rows labelled 250 or 1,000 used candidate samples of that many men. After exact
full-pool draws replaced sampling (D072), the user reran the scenario with 40
replications (rows labelled exact draws): results match the 1,000-man runs,
with direction balanced and the input preferences recovered in wife-side sets.
In the user's setup a run took about 1.1 s, down from about 2.2 s with
1,000-man samples.

Findings from these runs and their discussion:

- Weighted choice removes female choice's directional imbalance (AB/BA about
  1 rather than 0.6) and recovers the input preferences in wife-side sets.
- At 250 men, both methods' pooled cross-products fell about 0.12 short of -6,
  and weighted choice's quadratic term was slightly flat. At 1,000 men both
  largely disappeared. Each woman's sample contains few men near her age, and
  varying group mixes across samples attenuate pooled associations; candidate
  sampling is a speed shortcut, so the larger sample is closer to the intended
  model. A small residual pooled shortfall remains under both methods.
- At 1,000 men the male-group term is not needed for recovery; its estimate,
  about 0.31, reflects the B men's average log weight.
- Husband-side sets are not a valid benchmark (see above).
- After the incremental pool sums, the user's rerun with an installed build
  reproduced the cross-products and AB/BA exactly; the clogit differed only
  through new choice-set draws. In the user's setup, weighted choice at 1,000
  men took 6.4 s per run before that optimization and 2.18 s after (5.5 s at
  250 before); build settings may differ from the timings above.

## Validation

`tests/testthat/test-weighted-choice.R` covers the analytical examples above
(including default targets, signed delta adjustments and the analytical 3:1 case),
binned calibration against an exact individual-level R reference for native
and custom formulas in both intensity types, three-group convergence, zero and
unattainable targets, a male group with no counted women, identities created
during a run (custom AB inheritance) and restarts with them, strict endogamy
with unit weights, recorded log weights across restarts, one/no-group
unit-weight shortcuts in both initiation modes, reproducibility, restarts, separate father models,
diagnostics without changed histories, the single warning, and input
validation, and incremental pool sums against a full rebuild after thousands
of pool changes and cutoff crossings.

Timing with candidate sampling, before D072, October 23, 2026 (optimized build; 20,000 starting persons in two
groups, on-demand, 50 years, one seed; synthetic test mortality): female
choice 0.13 s and 0.29 s at sample sizes 250 and 1,000; weighted choice 0.28 s
and 0.44 s; additive two-sided (since removed) 1.05 s at 250. Before incremental sums,
weighted choice took 1.1 s at 250. These are single-scenario measurements, not
performance guarantees.

With exact full-pool draws (D072; same scenario, optimized build): female
choice 0.10 s and weighted choice 0.25 s at 20,000 persons over 50 years;
0.39 s and 1.12 s at 50,000 persons over 100 years; weighted choice with
three groups 0.40 s. Single-seed measurements.

## Accepted milestone 3 assessment

The user accepted milestone 3 on October 27, 2026 (D077), after review and
revision of the [validation report](../validation/milestone3.qmd). It assesses
uneven and even three-group populations, pooled cross-products, directional
union ratios with delta, wife-side conditional-logit recovery using recorded
log-weight offsets, marriage rates, and two-group fertility/mortality/age-gap
scenarios. These are the user's scientific findings, not independent reruns
or guarantees of exact pooled odds ratios or fixed directional balance.

The report retains limitations for small or sex-imbalanced markets and the
absence of marriage-squeeze feedback. Starting-population sensitivity in the
inheritance experiments did not establish a calibration bug or a systematic
3:1 group-size effect. Full two-sex models remain deferred.
