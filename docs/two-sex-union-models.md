# Two-sex union formation: future directions

Recorded October 24, 2026 from a discussion with the user. **Status: research
notes and proposals only.** Nothing here is implemented, specified for
implementation, or an accepted modeling decision, except the user decisions
recorded in D073. The current union process remains as described in
[union-formation models](union-formation-models.md) and
[weighted-choice matching](weighted-choice-matching.md).

Calculations below are exact group-only (no age) expectations derived during
the discussion unless stated otherwise; they were not produced by simulation.

## Starting question: directional asymmetry in intermarriage

The user's example: in the US, black-man/white-woman unions outnumber
white-man/black-woman unions roughly 3 to 1. A single symmetric odds ratio
can describe the association, but not this gender asymmetry. Could users
supply a directional ratio that weighted choice preserves when each group's
sexes are balanced?

### Why weighted choice cannot produce it

With two groups, AB - BA equals the male group-B share of unions minus the
female group-B share (A:B means female A, male B). With balanced sexes the
weighted-choice target sets those shares equal, so AB = BA however the
preferences are written; directed pair effects are absorbed by the weights
(see the weighted-choice consequences). A 3:1 ratio therefore means men of
group B form more unions than women of group B. Because women's timing is an
input, those extra unions are taken from men of group A. The asymmetry is
about how often each male group forms unions, not about who prefers whom.

### Two mechanisms considered

1. **Male union propensity.** A multiplier c[g] on each group's male
   per-capita intensity in the target, applied before imbalance absorption so
   kappa handles only sex imbalance. Simple and estimable from male/female
   union rates by group, but the ratio is an outcome that moves with group
   sizes and preference strength. The c_B/c_A giving exactly AB/BA = 3 with
   balanced sexes:

   | Groups (A/B) | Log cross-product | c_B/c_A |
   | --- | ---: | ---: |
   | 50/50 | -6 | 1.115 |
   | 75/25 | -6 | 1.131 |
   | 90/10 | -6 | 1.184 |
   | 75/25 | -3 | 1.554 |

2. **Directional ratio target.** The user supplies `c("A:B" = 3)`. At each
   search, solve for the male group shares giving expected AB/BA = R at the
   current female intensities and utilities, then add kappa (M - F) for sex
   imbalance. With two groups the ratio is monotone in log b_B, so the solve
   is unique. R = 1 reproduces the current method. With K groups the weights
   have K - 1 degrees of freedom, so at most K - 1 non-cyclic pair ratios can
   be set; cyclic directional variation comes from directed preferences,
   which do affect pairing with three or more groups.

### Feedback through the male pools

Women still form unions at their input hazards, so women's unpartnered pools
do not respond. Group-B men partner faster than group-B women and group-A men
slower than group-A women, so the unpartnered pools become sex-imbalanced
within groups (B short of men, A with extra men). With kappa > 0 the
imbalance correction then pulls the realized ratio back toward 1; at kappa = 0
the ratio holds and the cost falls only on group-A men's waiting times. Fewer
group-B men also raises group-B women's chance of choosing group-A men, which
a per-search ratio solve would correct and a fixed propensity would not.

The user noted research suggesting this asymmetry contributes to a marriage
squeeze for black women. The current model cannot represent that squeeze:
female union timing never responds to male availability (D073 declines both
mechanisms above).

## Two classes of union-formation models

The user's framing:

- **Timing independent of matching (current).** The union event (a female
  search, or an unpartnered birth) occurs independently of partner matching.
  This preserves desired input rates (female union rates, or overall
  fertility) but prevents group and age composition from affecting event
  timing, and does not treat the sexes symmetrically.
- **Full two-sex models.** Union timing responds to both sexes' availability
  and composition, allowing marriage-squeeze effects, at the cost of not
  necessarily reproducing desired input rates.

Users may want one or the other for specific purposes.

## Possible extensions

1. **Participation (outside-option) function.** Keep the current structure,
   but make the probability that a woman's search produces a union depend on
   market conditions, for example an alternative "remain unpartnered" with
   utility V0, so P(union) = sum b exp(U) / (exp(V0) + sum b exp(U)). Input
   hazards become reference rates, reproduced only at a reference pool. This
   is the most direct extension of the current code and overlaps the deferred
   poor-pool waiting item (D044).
2. **Schoen harmonic-mean model**, possibly modified along the lines of
   Dagsvik's model (below).
3. **Choo–Siow transferable-utility model**, as a later possibility.

## Schoen's harmonic-mean model

For female type i and male type j with unpartnered counts F_i and M_j:

    Y_ij = alpha_ij F_i M_j / (F_i + M_j)

(Schoen 1981; Schoen 1988). alpha_ij is the force of attraction.

### Properties relevant to Hammel

- **Directional asymmetry is native.** alpha is indexed by (female type, male
  type). With balanced sexes the two cross-group harmonic means are equal, so
  Y_AB / Y_BA = alpha_AB / alpha_BA exactly; with imbalance the ratio moves
  through the harmonic means.
- **Marriage squeeze is native.** A shortage of group-B men reduces
  H(F_B, M_B) and H(F_A, M_B), lowering group-B women's union rates.
- **Continuous time.** Schoen used the form as an instantaneous rate. A type-i
  woman's hazard toward type j is alpha_ij M_j / (F_i + M_j); each specific
  pair's hazard is alpha_ij / (F_i + M_j).

### Implementation sketch

**Types and cells.** A type is defined for one sex: group by birth-time bin
(as in weighted choice), plus any other attribute affecting pairing. A cell is
a pair of types. Birth-time bins never change membership as people age; age
enters through alpha evaluated at the bins' current ages. Quarter-year bins
over a 50-year range with two groups give about 200 x 200 x 4 = 160,000
cells.

**Market clock.** Women's and men's hazards describe the same event, so
simulating both sides' clocks would double-count. Unions leave the per-person
event queue: the market holds one pending next-union event whose rate is the
sum of all cell rates R_ij = alpha_ij F_i M_j / (F_i + M_j), kept in a sum tree.

**Rate changes.** Unlike the current engine, union times cannot be drawn once
and kept: every change to F_i or M_j (union, death of an unpartnered person,
widowhood, dissolution, reaching the minimum age) changes all cell rates in
that row and column. The update touches one row and one column, then redraws
the single market time; no individual clocks are redrawn, and mortality,
fertility within unions and other per-person events keep their precomputed
times. Cost scales with stock-changing events times a row-plus-column update,
comparable to weighted choice's incremental pool sums. Not benchmarked.

**Rates drifting between events.** alpha changes as people age (and later with
period inputs) even without events. Options: thinning against an upper bound
on the total rate (exact, continuous; costly if the bound is loose);
recomputing rates at bin boundaries (a timing approximation); or fixed-interval
stock snapshots (reintroduces discrete execution, would require revisiting
D024; not recommended). Thinning with a bin-boundary fallback was suggested.

**Pair selection.** When the clock fires: choose cell (i, j) with probability
R_ij / R_total (one sum-tree walk), then one woman uniformly from type i and
one man uniformly from type j. Uniform draws are exact because all pairs in a
cell share the same hazard, so each pair has probability 1/(F_i M_j). Keep
membership arrays per type with swap-removal for O(1) draws. Uniform selection
requires every pairing-relevant characteristic in the type definition and
alpha; coarse types with within-cell weighting would break the equal-pair
hazard and make cell rates inconsistent.

**Barred pairs (kinship).** Two choices with different demography: a null
event (reject, form no union, redraw the clock) is exact thinning if barred
pairs have zero hazard, so their unions do not happen; redrawing within the
cell keeps the cell rate and redirects those unions. Hard group exclusions set
alpha_ij = 0.

**Inputs.** Users would supply a parametric alpha, not cell values, for
example

    log alpha = a_F(female age) + a_M(male age) + age-gap terms + group effects

Compared with current inputs, a_F replaces the female search schedule
(reinterpreted), a_M is new, and the gap and group terms (including directed
effects) carry over. Group multipliers become group main effects. Separable
schedules with a gap term are a restriction; it can be loosened with splines,
group-specific age schedules, or gap terms varying by group pair or age.

Two input routes:

1. **Estimation.** From new unions in a period, cross-classified by both
   partners' ages and groups, plus unpartnered stocks by sex, age and group:
   a Poisson regression with log alpha as the linear predictor and
   log(F_i M_j / (F_i + M_j)) as offset. Estimate at single-year ages with a
   parametric form rather than from coarse tables (see aggregation below).
2. **Calibration.** Female and male union rates by age (and group) for a
   reference population, plus gap and group preferences; solve for a_F and a_M.
   The two sexes' rates must imply the same total unions, so they can be
   matched only up to that consistency.

### Built-in homogamy and aggregation

The denominator F_i + M_j is local to the cell. Each pair type behaves as a
separate market, so with flat alpha and unequal group sizes the harmonic mean
produces endogamy. Example, 90/10 groups, balanced sexes, alpha = 1:

| Pair | F_i + M_j | Per-pair hazard |
| --- | ---: | ---: |
| A-A | 1.8 | 0.56 |
| A-B, B-A | 1.0 | 1.00 |
| B-B | 0.2 | 5.00 |

Union table: AA 0.45, AB 0.09, BA 0.09, BB 0.05 (random mating would give BB
0.01). Log odds ratio 1.02. Group B forms unions at 1.4 per capita against 0.6
for group A. In general the odds ratio is (a + b)^2 / 4ab >= 1, with equality
only for equal groups.

Splitting a homogeneous population into two meaningless labels doubles total
unions at a 50/50 split (0.5 to 1.0) and gives a log odds ratio of 1.02 at a
90/10 split. Type definitions therefore change behavior. Fitted alpha
reproduces the data at the estimation composition, but alpha is not a pure
preference and projections are sensitive to composition and type definition.

Bin width: for fixed alpha, halving bin width roughly doubles total unions, so
alpha must be defined per unit of age density (scaled by bin width) for bin
width to remain a numerical tolerance; this converges to Schoen's
continuous-age form. Alpha estimated from coarse published tables does not
transfer exactly to fine ages.

### Translating odds ratios loosely

The user wants a simple way to translate published log odds ratios into
comparable Schoen parameters, not to force the model to reproduce odds ratios.
For two female and two male types, the F M products cancel in the flow-table
cross-product:

    log OR_unions = log OR_alpha
                    + log[(F_i + M_j')(F_i' + M_j) / ((F_i + M_j)(F_i' + M_j'))]

A woman-side conditional logit on available partners recovers the same
quantity. With balanced sexes and unpartnered group shares a and b, the
correction is log((a + b)^2 / 4ab): 0 at 50/50, 0.288 at 75/25, 1.022 at 90/10.
So a crude translation is alpha log OR = input log OR - correction, valid when
the shares are of the unpartnered at-risk population, the odds ratio comes from
new unions rather than current marriages, and there is no age structure.

With age, the correction applies to each female-age/male-age pair. With the
same group shares and age distributions in both groups and unpartnered
densities f and m at the paired ages:

    log[1 + ((a - b)^2 / ab) fm / (f + m)^2]

fm / (f + m)^2 <= 1/4, with equality at f = m; it falls toward 0 as one side
dominates (at a 3:1 unpartnered sex ratio the 90/10 correction is about 0.85).
The crude correction is therefore an upper bound.

A possible helper (illustrative name `convert_log_odds_to_attraction()`) would
take log odds ratios in `group_exogamy` pair syntax, an optional directional
ratio, and either group shares or a population; apply the correction pair by
pair against within-group pairing; carry a directional ratio over unchanged
under balanced sexes (corrected by the ratio of the cross-group harmonic means
otherwise); and return a tibble of pair effects with the correction shown. Its
documentation would state that it ignores age structure and gives starting
values, not reproduced odds ratios. Alpha stays fixed during simulation, so
realized odds ratios move with composition.

### Known weaknesses

- No spillovers: a shortage of one male type does not shift demand toward
  other types (the Choo–Siow critique).
- No accounting constraints: implied rates can exceed what is plausible far
  from the estimation data.
- The built-in homogamy and type-definition sensitivity above.

## Related models

**Choo–Siow (2006).** An econometric implementation of Becker's
transferable-utility matching, positioned as an answer to demography's two-sex
problem; the companion Demography paper ("Estimating a marriage matching
model with spillover effects") criticizes mean functions for lacking a
behavioral basis, spillovers and accounting constraints. Matching function
mu_ij = exp(Pi_ij / 2) sqrt(mu_i0 mu_0j) with population constraints for each
type, which link all cells. The square-root terms cancel in the odds ratio, so
flat gains give no association. Splitting a type still raises total unions
(at a 50/50 split with gain k, from k/(1+k) to 2k/(1+2k)), likely because
taste shocks are drawn per partner type.

**Dagsvik (2000); Cheng, Dagsvik, Han & Jia (2026).** A stochastic two-sided
matching model from random utility and stable matching, with shocks per
potential partner. Matching probabilities phi_i(j) = (M_j / F) c_ij phi_i(0)
phi_j(0), with single probabilities solving a joint system; c_ij = a_ij b_ij
is joint attraction, not identified separately by side, but directional
c_AB != c_BA is allowed. Cheng et al. compare c_ij with Schoen's alpha_ij
algebraically (their Section 5.4): c_ij ~ Y_ij F / (F'_i M'_j) uses end-of-period
singles and a global scale F, while alpha_ij uses mid-period stocks and the
local scale F_i + M_j. They do not run both models on the same data. The
global denominator is what removes the built-in homogamy and type-splitting
sensitivity, and the joint system supplies spillovers. Cheng et al. estimate
from ACS unions in the past 12 months plus singles, using 18 types (2
education x 3 race x 3 age); their age brackets (women under 24, 24–37, 38+;
men under 26, 26–41, 42+) cannot identify realistic age gaps. A parametric
attraction surface over fine ages avoids that. The model describes a market
clearing each period, so it needs a continuous-time analogue for Hammel.

| Property | Schoen | Choo–Siow | Dagsvik |
| --- | --- | --- | --- |
| Endogamy from group size with flat preferences | Yes | No | No |
| Splitting a type changes total unions | Yes | Yes (modestly) | No |
| Spillovers between types | No | Yes | Yes |
| Accounting constraints | No | Yes | Yes |
| Directional asymmetry | Yes | Yes | Yes |

The table's Schoen and Choo–Siow entries follow from the formulas; the
type-splitting and spillover interpretations in the econometric literature
(for example Menzel 2015; Galichon and Salanié) were not checked.

## Interface considerations

`union_formation()` mixes timing (`initiation`, `rates`, `group_prop_hazard`),
partner selection (`matching`), pair eligibility (`min_male_age`,
`kinship_generations`), unpartnered-birth handling (`father_matching`,
`cancel_birth_if_no_match`) and diagnostics, with validation and the engine
branching on `initiation`. A two-sex model needs eligibility, unpartnered-birth
handling and diagnostics, but neither initiation nor a matching model.

A sketched alternative: `set_union_formation()` takes a union-model object from
class-specific constructors (for example on-demand, rate search, and later
harmonic mean), each with its own validation, requirements on other processes
and engine preparation, plus a shared eligibility object; internally, a C++
interface implemented by each built-in model. These names and shapes are
placeholders.

The user identified unresolved connections:

- Cohabitation (possibly milestone 4) would add union types and transitions.
- Nonmarital father selection: setting harmonic-mean alphas for unions and
  conditional-logit odds ratios for unpartnered fathers would be inconsistent,
  so father selection, cohabitation and unpartnered partnering must be designed
  together with any two-sex model.

The user deferred API restructuring (D073): the current API is unchanged until
a single working implementation incorporates all union-formation components.

## Open questions

- Whether alpha or the odds ratio should be the invariant under composition
  change (they predict differently when group sizes change).
- Treatment of barred pairs (null event or within-cell redraw).
- Exact handling of alpha drift between events (thinning bounds, refresh).
- Spillovers: whether to extend Schoen's form, adopt a Dagsvik-type global
  denominator, or use Choo–Siow.
- Reference conditions for estimated or calibrated inputs.
- How unions, cohabitation and nonmarital fathers share one model.

## References

- Schoen, R. 1981. The harmonic mean as the basis of a realistic two-sex
  marriage model. *Demography* 18:201–216.
  <https://link.springer.com/article/10.2307/2061093>
- Schoen, R. 1988. *Modeling Multigroup Populations*. Plenum.
- Choo, E., and A. Siow. 2006. Who marries whom and why. *Journal of Political
  Economy* 114(1):175–201.
- Choo, E., and A. Siow. 2006. Estimating a marriage matching model with
  spillover effects. *Demography* 43(3).
  <https://link.springer.com/article/10.1353/dem.2006.0023>
- Dagsvik, J. K. 2000. Aggregation in matching markets. *International
  Economic Review* 41(1).
- Cheng, Y., J. K. Dagsvik, X. Han, and Z. Jia. 2026. Force of attraction and
  partner availability in the U.S. marriage market: a two-sided matching
  model. *Sociological Science* 13:178–213. DOI 10.15195/v13.a8.
  <https://sociologicalscience.com/articles-v13-8-178/>
