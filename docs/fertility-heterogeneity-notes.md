# Fertility weights: distribution and inheritance design notes

October 29, 2026. Discussion notes requested by the user, not an accepted
scientific specification or authorization to implement a distribution. The
[accepted cell-accounting contract](fertility-cell-accounting.md) reserves
individual selection weights; they are currently all one. Heterogeneity
parameters, defaults, initialization, persistence and inheritance remain to be
selected. These notes distinguish a candidate mathematical construction from
its unresolved behavior under reproductive selection.

## Separate dispersion from inheritance

The user wants control of fertility heterogeneity independently of the degree
of maternal inheritance. The same dispersion conventions should apply to
initial women and newborn daughters. SOCSIM does not provide that separation:
its fresh distribution is fixed, arithmetic blending shrinks variance, and
betaT transforms both mean and dispersion. See the [SOCSIM audit](socsim-fertility-heterogeneity-audit.md).

Under Hammel's cell semantics, births occur at intensity m_c*N_c and an eligible
woman is chosen in proportion to her weight. Multiplying all eligible weights
in a cell by the same positive constant does not change selection. Mean-one
weights provide an interpretable reference scale; they do not enforce the cell
birth rate, which is controlled by the intensity itself. They are not direct
multipliers of a woman's unadjusted schedule hazard.

## Candidate distribution: mean-one lognormal

The agent recommended lognormal weights as a starting proposal, not a selected
default. They are positive and allow continuous control of dispersion. Gamma
weights are another candidate, with a less extreme upper tail than lognormal
weights at the same positive coefficient of variation (CV). Neither has a
finite upper bound. Potential concentration of births among a few eligible
women needs empirical assessment.

Parameterize by desired CV c >= 0. With Z standard normal, set

    sigma^2 = log(1 + c^2)
    W = exp(-sigma^2/2 + sigma*Z).

Then E(W)=1, Var(W)=c^2 and SD(W)=c. Thus c also equals SD on this mean-one
scale. At c=0 all weights are one; there is no meaningful empirical inheritance
correlation when everyone has the same weight. No positive numerical default
CV has been chosen. Sample means and variances fluctuate around these
population moments; this is not exact per-sample normalization.

The [NIST lognormal formulas](https://www.itl.nist.gov/div898/handbook/apr/section1/apr164.htm)
support the mean/variance parameterization. Selection of lognormal rather than
gamma or another distribution remains a demographic/modeling decision.

## Candidate inheritance: preserve variance on a latent normal scale

For positive sigma, recover a mother's standardized latent value as

    Z_m = (log(W_m) + sigma^2/2) / sigma.

Draw an independent epsilon ~ Normal(0,1), and for 0 <= alpha <= 1 set

    Z_d = alpha*Z_m + sqrt(1-alpha^2)*epsilon
    W_d = exp(-sigma^2/2 + sigma*Z_d).

If Z_m ~ Normal(0,1), Z_d is also Normal(0,1). Consequently daughters have
the same complete lognormal distribution, including mean and variance.
The variance contributions are alpha^2 and 1-alpha^2; ordinary arithmetic
weights alpha and 1-alpha would instead shrink variance at intermediate alpha.

| alpha | Interpretation under the stated assumptions |
| --- | --- |
| 0 | Fresh independent daughter weight |
| .5 | Latent/log-weight mother–daughter correlation .5 |
| 1 | Exact maternal weight inheritance |

This alpha is a correlation on the latent normal, equivalently log-weight,
scale. It is not generally the correlation of raw weights. Under the same
standard-normal maternal-distribution assumption and sigma>0, raw-weight
correlation is

    rho_W = [exp(alpha*sigma^2)-1] / [exp(sigma^2)-1].

If users instead want alpha to denote raw-weight correlation, the latent
coefficient would need to be transformed accordingly. The public parameter's
interpretation has not been selected. The formula also assumes the same sigma
for mother and daughter; changes in dispersion or arbitrary supplied initial
weights require a separate contract.

## Crucial qualification: mothers are selected

The user explicitly set reproductive selection aside while discussing the
variance-preserving formula. It must not be forgotten when specifying the
simulation. A transition preserving the lognormal distribution for a randomly
sampled mother does not automatically preserve it among daughters born to a
fertility-selected set of mothers.

More generally, if the actual reproducing mothers' latent values have mean mu_m
and variance v_m, the proposed transition yields

    E(Z_d) = alpha*mu_m
    Var(Z_d) = alpha^2*v_m + (1-alpha^2).

These equal zero and one only under the corresponding maternal moments (or
alpha=0); normality is additionally needed for the full lognormal claim. If
maternal latent values are non-normal, their latent moments alone do not
establish the daughters' weight variance. The transformation therefore does
not justify promising a common fixed weight variance across all generations,
parities or union types under reproductive selection.

It remains unresolved whether a future model should preserve a fixed marginal
birth-cohort distribution, permit selection-driven distribution changes, or
use another explicit calibration. Within-cell birth-rate calibration does not
settle that question: it fixes cell intensity but still selects mothers by
weight. No corrective centering, ranking or renormalization is authorized here.

## Decisions and validation needed before implementation

- Select distribution, dispersion parameter and any defaults; decide whether
  users may supply individual weights and whether zero weights are supported.
- Define alpha's scale and whether inherited variation is maternal only.
- Specify assignment/storage for founders, newborns and restarts, including
  changed dispersion, missing parents and user-supplied weights.
- Decide what distributional changes from reproductive selection are intended.
- Implement weighted eligible-member selection consistently with the existing
  cell contract; uniform selection is the current implementation.
- Validate independent draws' moments, alpha endpoints, latent variance and
  correlations under the theoretical assumptions, and actual simulated cohort
  distributions under selection. Check conditional rates, spacing, group/type
  transitions and sparse-cell/unrealized-intensity behavior separately.

This note records the discussion; it does not add engine behavior, adopt
SOCSIM's inheritance transformation, or mark milestone 5 complete.
