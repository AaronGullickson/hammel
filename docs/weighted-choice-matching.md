# Weighted-choice matching

Specified with the user and implemented October 23, 2026 on
`experiment/two-sided-sequential-matching` (D066). The user made
`matching_method = "weighted_choice"` the default and removed additive
two-sided matching (D067); partners are now drawn exactly from the whole
eligible pool (D072). The
method is an argument of the matching-model constructors (D068), so a separate
father model can use its own; `"female_choice"` remains available. This describes
implemented behavior; no acceptance thresholds are set.

## Problem

The user requires partner matching that (1) preserves age-specific female
union timing from the inputs, (2) recovers the input matching preferences, and
(3) yields equal female-A/male-B and female-B/male-A union counts when
preferences are symmetric and each group's sexes are balanced.

| Method | Preferences recovered | Directional balance |
| --- | --- | --- |
| Female choice | Yes | No: 75/25 groups, A::B = -3 give AB/BA = 0.377 |
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
   by the summed intensity. Kinship screening is ignored here.
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
`group_exogamy`, defaulting to 1 (user, October 23, 2026). It sets the share of a
group's sex imbalance among unpartnered people absorbed through cross-group
unions. At 1, men receive the full mirrored target. At 0, each male group's
share of unions equals its female group's share, so the excess sex waits
longer. With symmetric preferences and balanced sex composition, F = M and
every value gives equal directional counts. With no group preferences, values
below 1 still disfavor men of an over-represented group. Female-choice models
ignore the argument.

### Nonmarital fathers

Rate-mode father selection has its own weights, calibrated from fertility
intensity with the applicable father utilities (`father_matching`, or
`matching` when NULL) and that model's `imbalance_absorption` (user, October
23, 2026). Each model applies its own `matching_method`: for example, weighted
unions with female-choice fathers, or the reverse (D068). Union searches keep weights calibrated from search hazards. The
male pool, target rule, and solver are shared. On-demand mode has one
calibration, driven by fertility.

### Numerical representation

The engine keeps counts of unpartnered men and women by group and 0.25-year
birth-time bin, and of women past their birth-spacing interval. Counts change
only at events (birth, death, formation, widowhood) and when a queued spacing
interval ends; restarts rebuild them from histories. At a search, each bin's
centre age (current time minus centre birth time) determines whether the whole
bin is at or above the age cutoff, for both sexes. A bin's intensity uses the
average hazard over its full quarter-year age span, so rate bands narrower than
a bin or edges inside it are not missed (D071); the only remaining
approximation is how births are spread within a bin.

Pair utilities depend only on the lag between bin centres, which equals male
age minus female age. A table over every lag in the run's birth-time range is
built once per matching model, including one evaluation of a custom R formula
on that grid. The formula must therefore give finite values at every such lag.
Bin width is a numerical tolerance, not a demographic parameter. Candidate
selection itself uses exact ages. Remaining binning approximations are the
whole-bin age cutoff and pair utilities at bin-centre ages (gap error up to
about 0.25 years, which matters only for very sharp age preferences). With
one-year rate bands and a 5-year half-weight distance, default bins gave B
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

When a matching model sets `prevent_exogamy = TRUE`, every woman sees only men
of her own group, so a group weight scales all her candidates equally and
cannot change her choice. Calibration is skipped for that model, giving
female choice exactly, with no warnings (D070). Partial restrictions (listed
pairs) still calibrate.

### When no solution exists

A male group may have a positive target but no feasible demand, for example
under `prevent_exogamy` when no compatible women have positive intensity. Its
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

Because of the first two points below, constructing a weighted-choice matching
model with directed pair effects (two directions differing, including a single
direction) or `male_group` formula terms warns once (D069).

- **Directed effects in two groups.** A:B = -1 with B:A = -2 produces exactly
  the same expected table as A::B = -1.5. With both sexes' group flows fixed,
  a two-group table has one free quantity, its cross-product. Directed effects
  change which weights are needed, not the flows. With more groups, only
  combinations not absorbed by the group flows affect pairing. A
  `male_group` main effect in the formula has no effect at all, since the
  group weights absorb it. Under female choice directed effects do move group
  flows, but they mix with female choice's directional bias when group sizes
  differ (for example, AB/BA = 0.52 with symmetric -1.5 effects at 75/25).
- **Sensitivity to composition.** AB - BA equals the difference between the
  male and female group-B shares of unions. When exogamy is rare, small
  differences in the sexes' unpartnered group shares move AB/BA strongly.
  Directional balance should be assessed with those pool compositions.
- **No user-specified directional ratios.** Because the weights fix both
  sexes' group flows, a gender-asymmetric intermarriage pattern (for example
  3:1 between the two cross-group directions) cannot be specified with
  balanced sexes. The user declined adding directional targets (D073); see
  [two-sex union models](two-sex-union-models.md) for the analysis and
  future alternatives.
- **Pooled cross-products.** These examples have no ages. With age preferences
  and group-specific age structures, a cross-tabulation pooled over ages need
  not equal the summed group utilities exactly, under any method.

## Diagnostics and interpreting estimates

Weights consume no random draws. When the general matching model uses
weighted choice in a grouped population, `diagnostics = TRUE` records `.log_weight_<group>` on each new union: the log group weights applied
at the search that formed it (October 23, 2026, D070). One group is the
reference with log weight 0; only differences between groups matter, and
unweighted searches record 0 for every group. Earlier unions and new groups'
columns for earlier rows are NA. In a wife-side choice set, using the column
for each alternative's group as an offset gives the exact benchmark, in place
of approximating the weights with a male-group main effect. Nonmarital father
selections form no union and are not recorded.

A conditional logit fitted to wife-side choice sets with the input terms only
estimates separate directed group terms offset by log(b) differences, while
their sum and the age terms are unaffected. Adding a male-group main effect
absorbs the weights. In the static analysis, calibrated weights left the age
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
with A::B = -3, 100 years). The 250-man runs of both methods and the
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
(including absorption 0, 0.5 and 1, and balanced symmetry at every absorption),
binned calibration against an exact individual-level R reference for native
and custom formulas in both intensity types, three-group convergence, zero and
unattainable targets, a male group with no counted women, identities created
during a run (custom AB inheritance) and restarts with them, strict endogamy
reproducing female choice, recorded log weights across restarts, exact equality with female choice for one or no groups
in both initiation modes, reproducibility, restarts, separate father models,
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

Scientific assessment (user): the milestone 3 group-matching scenarios, with
directional counts, unpartnered pool composition by sex and group over time,
cross-tabulations, conditional logit with and without partner-group effects,
and realized age gaps.
