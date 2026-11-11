# Fertility weights: distribution and inheritance design notes

**Status (D176, November 11, 2026):** Hetfert is not available to users.
`het_fert()` is no longer exported and `fertility()`/`set_fertility()` no
longer accept `hetfert`. The engine path below is retained and tested through
the internal `enable_het_fert(model, het_fert(cv, alpha))`, pending a general
quantitative-characteristics feature. The rules below describe that internal
path.

October 29, 2026. The user selected mean-one lognormal fertility weights and
maternal inheritance on the normal/log-weight scale, exposed as
`set_fertility(..., hetfert = het_fert(cv, alpha))`. This specifies the accepted
behavior below; implementation and validation accompany it. The user accepted
milestone 5 under D105 on October 29, 2026, based on
[the validation report](../validation/milestone5.qmd). The reproductive-selection
qualification remains part of the contract: no per-generation distribution correction is used.

## Execution, inputs and observations

`het_fert(cv, alpha)` requires finite scalar cv >= 0 and alpha in [0,1]. Neither
parameter has an implicit default. No `hetfert` (the default) disables heterogeneity and
retains uniform mother selection. Model construction consumes no random draws.
When enabled, simulation uses the following rules:

- If persons lacks `hetfert`, draw mean-one lognormal weights for every recorded
  woman, including deceased women, under the simulation RNG/seed. Male rows get
  NA. This also applies to a zero-year run. CV zero draws no randomness and
  gives weight one.
- If the column exists, retain supplied female weights rather than redrawing or
  reconstructing them from ancestry. Female entries must be finite and strictly
  positive; an incomplete column is an error. Existing male entries are ignored
  and preserved. Female zero is not an infertility flag. New male rows get NA.
- Daughters inherit from their actual mother's stored weight at birth using the
  formula below. Alpha zero gives a fresh draw; alpha one copies exactly,
  without an innovation draw. Weights persist across all age/parity/group/type
  cells and are retained on restart. Supplied weights need not follow the model
  distribution; selecting a new CV on restart does not rescale existing weights.
- Use the equivalent log-space formula
  `log(W_d) = mu + alpha*(log(W_m)-mu) + sigma*sqrt(1-alpha^2)*epsilon`,
  where mu=-sigma^2/2. This avoids division by sigma. At CV zero, daughters
  of ordinary unit-weight mothers stay at one; for supplied nonunit mothers it
  gives `W_d = W_m^alpha`, with no random innovation.
- When enabled, output persons always includes numeric `hetfert`, including
  founder/descendant females and zero-year or empty runs. When disabled, output
  omits `hetfert` even if the starting population supplied it (D104). The input
  population is unchanged. Re-enabling from a disabled run draws fresh values
  for all recorded women, rather than retaining a partially populated column.
- Cell intensity remains input rate times member count. Eligible mothers are
  selected proportional to weight; equivalent individual hazard is
  `cell_intensity * weight / eligible_weight_sum`. Matching calibration uses
  these same individual hazards. Spacing, failed-attempt rules, death priority
  and unrealized-intensity reporting are unchanged.
- Generated weights that underflow to zero or overflow to a nonfinite value
  cause an explicit error rather than clipping. Extremely large supplied totals
  must likewise be representable. No weights are silently capped or normalized.

For example, two eligible women with weights 1 and 3 in the same cell have
selection probabilities 1/4 and 3/4. Their cell still generates attempts at
`rate * 2`. A birth removes the mother's weight from eligible totals until her
spacing release, while she continues contributing to total member count in her
current parity cell. Cancelled attempts do not remove eligibility or weight.

Validation must cover supplied/drawn/restarted weights, RNG isolation and the
NULL path, CV zero, alpha endpoints, intermediate-alpha innovation moments,
weighted draws and removal/re-entry, group/type/parity transitions, and matching
calibration's equality with an independent scan of individual weighted hazards.

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

## Mean-one lognormal weights

The user selected lognormal weights; gamma remains an unimplemented alternative. They are positive and allow continuous control of dispersion. Gamma
weights are another candidate, with a less extreme upper tail than lognormal
weights at the same positive coefficient of variation (CV). Neither has a
finite upper bound. Potential concentration of births among a few eligible
women needs empirical assessment.

Parameterize by desired CV c >= 0. With Z standard normal, set

    sigma^2 = log(1 + c^2)
    W = exp(-sigma^2/2 + sigma*Z).

Then E(W)=1, Var(W)=c^2 and SD(W)=c. Thus c also equals SD on this mean-one
scale. At c=0 all weights are one; there is no meaningful empirical inheritance
correlation when everyone has the same weight. CV must be supplied explicitly; no positive numerical default
is supplied. Sample means and variances fluctuate around these
population moments; this is not exact per-sample normalization.

The [NIST lognormal formulas](https://www.itl.nist.gov/div898/handbook/apr/section1/apr164.htm)
support the mean/variance parameterization. The implementation uses the user-selected lognormal construction.

## Maternal inheritance on a latent normal scale

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
coefficient would need to be transformed accordingly. The user selected alpha on the latent/log-weight scale. The formula also assumes the same sigma
for mother and daughter; existing weights are retained when dispersion changes, and subsequent
daughters use the new specification in the log-space formula above.

## Crucial qualification: mothers are selected

The user explicitly set reproductive selection aside while discussing the
variance-preserving formula. It must not be forgotten when specifying the
simulation. A transition preserving the lognormal distribution for a randomly
sampled mother does not automatically preserve it among daughters born to a
fertility-selected set of mothers.

More generally, if the actual reproducing mothers' latent values have mean mu_m
and variance v_m, the selected transition yields

    E(Z_d) = alpha*mu_m
    Var(Z_d) = alpha^2*v_m + (1-alpha^2).

These equal zero and one only under the corresponding maternal moments (or
alpha=0); normality is additionally needed for the full lognormal claim. If
maternal latent values are non-normal, their latent moments alone do not
establish the daughters' weight variance. The transformation therefore does
not justify promising a common fixed weight variance across all generations,
parities or union types under reproductive selection.

The selected transition permits selection-driven distribution changes.
Within-cell birth-rate calibration fixes cell intensity but still selects
mothers by weight. There is no corrective centering, ranking or renormalization;
preserving a fixed birth-cohort distribution would require a separate decision.

## Remaining scope

Gamma weights, permanent zero-weight infertility, cohort-distribution
normalization, paternal inheritance, heterogeneity distributions varying across
strata and runtime weight changes are not included. None is implied by the
selected maternal lognormal model. Validation of weights and conditional rates
does not by itself establish universal demographic plausibility. Milestone
acceptance rests on the user's report assessment (D105), with these limitations
retained.
