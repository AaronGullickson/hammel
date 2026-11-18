# Two-sex union formation: future directions

Recorded October 24, 2026 from a discussion with the user. **Status: research
notes and proposals only.** Nothing here is implemented, specified for
implementation, or an accepted modeling decision, except the user decisions
recorded in D073 and D074. D074 authorizes the male-target adjustment
described below; D075 subsequently removes unweighted matching and directed
pair effects. The user assigned a two-sex marriage model to
tentative milestone 11 on October 30, 2026 (D117), subsequently moved to
tentative milestone 12 on October 31 (D121) and then to
[post-release](roadmap.md#post-release-two-sex-marriage-model) (D136). This
sets a planning priority, not a model or API; the fuller two-sex proposals
remain deferred. A second discussion on November 12, 2026 revised the
implementation sketch for the cell-scheduled engine (D128) and compared
denominator choices and spillover mechanisms; see
[denominators and spillover](#denominators-and-spillover-november-12-2026).
The current union process remains as described in
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

### Why zero-delta weighted choice cannot produce it

With two groups, AB - BA equals the male group-B share of unions minus the
female group-B share (A:B means female A, male B). With balanced sexes and zero delta the
weighted-choice target sets those shares equal, so AB = BA however the
preferences are written; historically, directed pair effects were absorbed by the weights
(see the weighted-choice consequences). A 3:1 ratio therefore means men of
group B form more unions than women of group B. Because women's timing is an
input, those extra unions are taken from men of group A. The asymmetry is
about how often each male group forms unions, not about who prefers whom.

### Two mechanisms considered

1. **Male union propensity.** A multiplier c[g] on each group's mirrored
   male target. D074 now implements c[g] = exp(delta[g]), followed by
   normalization, replacing imbalance absorption. The directional
   ratio is an outcome that moves with group
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
   current female intensities and utilities, with a separate rule needed for
   how sex imbalance affects the target. With two groups the ratio is monotone in log b_B, so the solve
   is unique. R = 1 reproduces the current method. With K groups the weights
   have K - 1 degrees of freedom, so at most K - 1 non-cyclic pair ratios can
   be set. The historical directed-preference mechanism could also create
   cyclic directional variation with three or more groups; D075 removes that
   pair-effect option.

### Feedback through the male pools

Women still form unions at their input hazards, so women's unpartnered pools
do not respond. Group-B men partner faster than group-B women and group-A men
slower than group-A women, so the unpartnered pools become sex-imbalanced
within groups (B short of men, A with extra men). The mirrored baseline responds to this imbalance and can pull the realized
ratio back toward 1 despite a fixed delta. Fewer
group-B men also raises group-B women's chance of choosing group-A men, which
a per-search ratio solve would correct and a fixed propensity would not.

The user noted research suggesting this asymmetry contributes to a marriage
squeeze for black women. The current model cannot represent that squeeze:
female union timing never responds to male availability (full two-sex feedback remains deferred; D074 authorizes target adjustments).

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

### Marriage squeeze (October 25, 2026)

In the current model a squeeze never changes whether or when a woman forms a
union, only whom she partners with: her search fails only when no eligible man
exists. A female outside option (extension 1) would let women respond, but not
symmetrically with men. A woman's union rate would be capped at her input
search hazard and would respond through a choice parameter (V0), while men's
rates respond mechanically and without a cap through being chosen more or
less often, so the response depends on which sex is scarce. The user decided
to defer marriage-squeeze modeling to a full two-sex choice model rather than
approximate it with a one-sided outside option. At the time, Schoen's
harmonic mean alone seemed insufficient: it lowers union rates under scarcity
but, at fixed stocks, does not redirect squeezed women toward other ages or
groups, and it is a rate function rather than a choice model. Choo–Siow and
Dagsvik are two-sex choice models with redirection through equilibrium
accounting; Dagsvik avoids the type-splitting sensitivity noted below. Given
the close algebraic relationship between Dagsvik's attraction parameter and
Schoen's force of attraction, the user hopes to find a form between the two.
The November 12 analysis
([denominators and spillover](#denominators-and-spillover-november-12-2026))
found that the Choo–Siow and Dagsvik spillovers act only through who remains
single, which unpartnered stocks reproduce in continuous time; whether that
passive redirection is sufficient remains open.

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

The sketch was revised on November 12, 2026 for the cell-scheduled engine
(D128); the October 24 version assumed per-person event queues.

**Types and cells.** A type is defined for one sex: current age band by group
category (by every declared dimension that affects pairing), plus any other
attribute affecting pairing, such as previous-union category for remarriage.
A cell is a pair of types. As in other D128 processes, people move between
types at age-band edges, and alpha is constant within a band. One-year bands
over ages 15–80 with four group combinations give about 260 types per sex and
68,000 cells. Union type (cohabitation or marriage) could be chosen after the
couple, as `cohab_prob` does now, or added as a cell dimension.

**Market clock.** Women's and men's hazards describe the same event, so
simulating both sides' clocks would double-count. The market is one process
replacing rate-mode female search, with one exponential clock whose rate is
the sum of all cell rates R_ij = alpha_ij F_i M_j / (F_i + M_j), kept in a sum
tree. It would take search's place in the tie order. Cells with F_i = 0 or
M_j = 0 have zero rate, so there are no failed searches.

**Rate changes.** Unlike other D128 processes, a cell's rate is not its own
member count times a rate: F_i appears in every cell of row i and M_j in every
cell of column j. Every change to F_i or M_j (union, death of an unpartnered
person, widowhood, dissolution, reaching the minimum age, crossing an age
band, group transition, immigration) recomputes that row or column, then
redraws the market time. Rates are piecewise constant between such changes,
so redrawing the memoryless clock is exact, as for other processes; age-band
and calendar-band edges are handled by the existing refresh and boundary
queue, so the October 24 thinning options are no longer needed. With 260
types per sex each change updates about 260 cells. Annual band crossings by
every unpartnered person will probably dominate the cost. Not benchmarked.

**Weights (D129).** Selection stays exact if a pair's hazard factors as
alpha_ij w_woman w_man; the cell rate then uses the summed weights on each
side and each partner is drawn by weight within the type. Non-factoring pair
weights break exact selection.

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
The correction below arises from Schoen's local denominator; with pooled or
global denominators it vanishes for group contrasts (see
[odds ratios and alpha](#odds-ratios-and-alpha)).
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

- No spillovers at fixed stocks: a shortage of one male type does not raise
  women's hazards toward other types (the Choo–Siow critique). In continuous
  time, squeezed women remain unpartnered and so form more unions with other
  types through the stocks; the November 12 analysis suggests this is the
  same mechanism as the Choo–Siow and Dagsvik spillovers.
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
sensitivity, and the joint system supplies spillovers through end-of-period
singles (see [denominators and spillover](#denominators-and-spillover-november-12-2026)).
Cheng et al. estimate
from ACS unions in the past 12 months plus singles, using 18 types (2
education x 3 race x 3 age); their age brackets (women under 24, 24–37, 38+;
men under 26, 26–41, 42+) cannot identify realistic age gaps. A parametric
attraction surface over fine ages avoids that. The model describes a market
clearing each period, so it needs a continuous-time analogue for Hammel.

| Property | Schoen | Choo–Siow | Dagsvik |
| --- | --- | --- | --- |
| Endogamy from group size with flat preferences | Yes | No | No |
| Splitting a type changes total unions | Yes | Yes (modestly) | No |
| Spillovers between types | Through stocks only | Through singles only | Through singles only |
| Accounting constraints | No | Yes | Yes |
| Directional asymmetry | Yes | Yes | Yes |

The table's Schoen and Choo–Siow entries follow from the formulas; the
type-splitting and spillover interpretations in the econometric literature
(for example Menzel 2015; Galichon and Salanié) were not checked. The
spillover row follows the November 12 analysis below.

## Denominators and spillover (November 12, 2026)

Discussion notes and the agent's own algebra, checked where stated against
Dagsvik, Flaatten and Brunborg (1998) and Cheng et al. (2026). Nothing here is
selected.

### Mass-action forms differ only in the denominator

Schoen's model is one member of a family Y_ij = alpha_ij F_i M_j / D, where
Y_ij is the flow of new unions per year between unpartnered type-i women and
type-j men, F_i and M_j are their unpartnered counts, and alpha_ij (per year)
is the force of attraction. A type-i woman's hazard toward type j is
alpha_ij M_j / D, and each specific pair's hazard is alpha_ij / D. With
D = F_i + M_j, alpha is roughly the hazard of a member of the scarce side when
the other side is plentiful, and alpha / 2 each person's hazard when the sexes
are balanced. The meaning and scale of alpha depend on D, so denominator
choices should be compared at the same reference rates, not the same alpha.

| D | Group-size endogamy with flat alpha | Splitting a type changes unions | Who dilutes rates | Update per stock change |
| --- | --- | --- | --- | --- |
| Local, F_i + M_j (Schoen) | Yes | Yes | No one else | Row i and column j |
| Pooled over groups within ages, F_a + M_b | No | Not across groups; age bins need density scaling | Others of the same ages | Cells with female age a or male age b |
| Global, F + M | No | No | Everyone eligible, for example older widows | Numerators only; divide the total by one scalar |

All are homogeneous of degree one: doubling every count doubles the flows, so
per-capita rates do not depend on simulated population size. The pooled form
is the agent's proposal, not from the literature. All three use the same market process and engine work;
only the cell-rate formula and update scope differ. Under global D, cell
selection uses the numerators alone.

### Odds ratios and alpha

For female types i, i' and male types j, j', every F and M appears once in
the numerator and once in the denominator of the flow-table cross-product, so

    log OR(flows) = log OR(alpha)
                    + log[D_ij' D_i'j / (D_ij D_i'j')]

where log OR(alpha) = log[alpha_ij alpha_i'j' / (alpha_ij' alpha_i'j)]. The
denominator term depends on D:

| D | Denominator term | Result |
| --- | --- | --- |
| Local, F_i + M_j | Nonzero, composition-dependent | The correction in [translating odds ratios](#translating-odds-ratios-loosely) |
| Pooled, F_a + M_b | Zero for group contrasts at the same ages (one shared D) | Group log OR equals log OR(alpha) exactly at every age pair |
| Global | Always zero | Exact for any contrast |

Under pooled D, contrasts between age pairs involve different D and keep a
correction; published odds ratios are usually group contrasts. If
log alpha = a_F(age) + a_M(age) + gap terms + group effects with group
effects constant across ages, a group interaction in alpha is the log odds
ratio of new unions in every age cell, whatever the sex ratio or group sizes.
Published log odds ratios could then be entered directly, much as
`group_exogamy` works now, and the illustrative conversion helper would be
unnecessary. A Poisson log-linear model of new-union counts with offset
log F_i + log M_j estimates alpha's group interactions; a D constant across
groups is absorbed by the intercept and age terms.

Caveats:

- The odds ratio must describe new unions among the unpartnered. Odds ratios
  from current marriages also reflect durations, dissolution and mortality.
- The equality holds within each age pair. A table aggregated over ages has
  the same odds ratio only if the groups' unpartnered age distributions are
  proportional; otherwise aggregation changes it, as in any non-collapsible
  table. Aggregation over calendar time with changing composition is
  similar.
- One odds ratio fixes only the association. Group main effects (union
  levels) and directional asymmetry need other information, such as group
  union rates or a directional ratio.
- The equality holds in expectation at each instant; simulated tables vary,
  especially for small groups.

### Spillovers in Dagsvik and Choo–Siow act through singles

Dagsvik, Flaatten and Brunborg (1998, eq. 2.5–2.10) give the number of
marriages as X_ij = c_ij X_i0 X_j0, where X_i0 and X_j0 are the women and men
of each type remaining single at the end of the period; Cheng et al. (eq. 3–4)
add the scale 1/F. The odds that a type-i woman matches a type-j man rather
than remaining single are therefore c_ij times the type-j men remaining
single (over F), with no term for the availability of other male types. A
shortage of type-j men raises matches with other types only because more
type-i women remain single (independence of irrelevant alternatives). The
same holds for Choo–Siow: mu_ik / mu_i0 = exp(Pi_ik / 2) sqrt(mu_0k / mu_i0)
does not involve type j (from the standard formula; the paper was not
rechecked). Dagsvik et al. note that McFarland's axiom that X_ij decreases
with other types' counts can fail in a special case, though not at their
Norwegian estimates.

In continuous time, a squeezed woman who does not form a union remains in
the unpartnered stock, so her type's flows toward other types rise without
any change in hazards. These period models' spillovers therefore appear to
correspond to stock dynamics, not to an active substitution that mass-action
forms lack. The relevant distinction is whether a shortage lowers a person's
total union hazard or redirects it:

| Model | Hazard toward other types when preferred partners are scarce | Total union hazard |
| --- | --- | --- |
| Current weighted choice (fixed search rate) | Rises; all of the rate is redirected | Unchanged; no squeeze |
| Outside option, lambda alpha_k M_k / (V0 + sum alpha M) | Rises partly, depending on V0 | Falls partly |
| Mass action; Dagsvik and Choo–Siow at fixed singles | Unchanged; redirection through stocks only | Falls |

A two-sided outside-option form is the natural intermediate, with the
difficulty noted for extension 1.

### Continuous-time limits

Dagsvik's model is not the integral of continuous mass action over a period.
With one homogeneous type and balanced sexes, continuous flow k F M gives
period unions kT F_start F_end, while Dagsvik gives c F_end F_end. They agree
only as the period shortens, where Dagsvik becomes Y_ij = c_ij F_i M_j / F,
matching Cheng et al.'s eq. 32. So mass action with a global scale is the
short-period limit of Dagsvik, and its behavioral foundation carries over only
in that sense.

Dagsvik's scale is F, the number of women (from the normalization kappa /
sqrt(F) -> 1), not F + M. In continuous time with D = F, a man's hazard
depends only on women's type shares, not their scarcity, so when women are
scarce total flow grows with the number of men; Dagsvik's period accounting
prevents this, but a continuous version would not. D = F + M is the agent's
symmetric modification. The choice of denominator is therefore not settled by
the literature.

The short-period limit of Choo–Siow is a geometric mean, Y_ij = gamma
sqrt(F_i M_j) over unpartnered stocks (the agent's derivation, not checked
against the literature). It is not mass action: a woman's hazard,
gamma sqrt(M_j / F_i), falls as more women of her own type compete. It is a
candidate for comparison alongside the mass-action forms.

### Inputs

Most users have female, and sometimes male, union rates by age rather than
alpha. One possible bridge is to calibrate alpha so the model reproduces input
rates at a reference population, such as the initial population, and let
composition change rates afterward. This would reproduce inputs at the start
while allowing squeezes later, combining the two model classes distinguished
above. Not designed.

### Suggested exploration

A deterministic stock-flow model in R, outside the package, with two groups
and a few age bands, comparing local, pooled and global denominators, the
geometric mean and perhaps a two-sided outside option on balanced sexes, a
male shortage in one group, and splitting a type into meaningless labels.
This would show squeeze size, redirection and odds-ratio behavior before any
specification. Not started.

## Cohabitation and nonmarital fathers (November 13, 2026)

Discussion notes for a two-sex market model. The cohabitation and father
directions are the user's tentative preferences, not accepted decisions.

### Cohabitation: partnership attraction and a couple-based share

Two ways to bring union type into the market:

- **A. Type drawn after the couple.** The market forms partnerships with one
  alpha; after the couple is chosen, cohabitation is drawn with probability
  p, as `cohab_prob` does now.
- **B. Type as a cell dimension.** Separate surfaces, with
  Y^marr_ij = alpha^marr_ij F_i M_j / D and
  Y^cohab_ij = alpha^cohab_ij F_i M_j / D.

If p can depend on both partners, the two are equivalent:
alpha^cohab = alpha p and alpha^marr = alpha (1 - p). **Tentative preferred
choice (user, November 13, 2026):** partnership alpha plus a couple-based
cohabitation share p, because it offers the most flexibility. A p that depends
only on the woman's age, group and previous union reproduces the current
`cohab_prob` semantics, and a couple-based p replaces `cohab_matching`. The
odds-ratio relationship above holds separately within each type.

Because p is drawn after a specific couple is chosen, it may use any
individual characteristic of either partner, such as the man's previous-union
category, parity or a shared child, without breaking the equal-pair hazard
within a cell; partnership alpha must use only type-defining characteristics.
A possible interface is a `cohab_prob`-style table (woman's age, group,
previous union, optional calendar bands under D122) for the woman-only case,
and a logistic model over both partners' variables, similar to D129 weight
formulas, for the couple-based case. Tables with both partners' dimensions
would grow quickly. Not designed.

**Caveat.** The parameterizations can represent the same models but hold
different quantities fixed when a user changes one input. Under alpha plus p,
raising the cohabitation share lowers direct marriage one-for-one and leaves
the total partnership rate unchanged, consistent with calendar-varying
`cohab_prob` today (D122). Under separate surfaces, raising cohabitation
attraction adds unions without reducing marriage. A user who wants
cohabitation to grow on top of marriage under alpha plus p must also raise
alpha. This resembles the open question of whether alpha or the odds ratio
is invariant under composition change, and must be documented clearly.

Conversion and dissolution act within existing couples and do not use the
market; dissolution returns both partners to it.

### Nonmarital fathers

Currently an unpartnered woman's births occur at her input fertility rate, and
the father is drawn from unpartnered men by conditional logit
(`father_matching`), without forming a union. Combining that with market
alpha for unions would give two inconsistent preference structures. Options:

- **F1. One-sex timing, father from the market row.** The birth occurs at the
  woman's input rate; the father's type is drawn with probability
  proportional to Y_ij over her row of the market, then a man uniformly within
  type j. One preference structure, fertility inputs and cell accounting
  (D101) preserved, and cheap given row access in the market sum tree. A
  father-specific difference (for example, younger men as fathers than as
  husbands) would be an optional multiplier on alpha, replacing
  `father_matching`. A shortage of men changes who the father is, not whether
  the birth occurs; births still fail only without eligible men. Squeezed
  women remain unpartnered longer and have births at unpartnered rates, so a
  squeeze shifts births between union and nonmarital fertility by the rate
  difference (an interpretation, not a validated result).
- **F2. Nonmarital conception as a market event,** at rate
  beta_ij F_i M_j / D, so nonmarital fertility responds to male availability.
  It would sit outside the fertility cells and need reconciliation with
  parity, spacing and cell accounting (D101), and would not preserve
  unpartnered fertility inputs.
- **F3. Fathers drawn also from partnered men.** A separate extension, not
  needed for consistency.

**Tentative preferred choice (user, November 13, 2026):** F1. The goal is to
keep father selection simple and preserve unpartnered fertility input rates.
This is not an accepted decision while the two-sex model remains unselected.
The user also wants an equivalent of `father_matching`: an optional
father-specific multiplier on alpha, so that the partners chosen as fathers
can differ from those chosen for unions (for example in age gap). Because
father selection is a draw conditional on a birth rather than a market
clock, the multiplier could use individual characteristics of either partner
beyond the type definition, as the cohabitation share can; drawing within a
type would then be weighted rather than uniform. Form and interface not
designed. Behavior with no eligible man (currently cancellation by default,
or an unknown father) is not yet discussed.

Previous-father continuity (D042, deferred) would raise alpha for specific
co-parent pairs, breaking the equal-pair hazard within a cell. It would need
its own process on linked co-parent pairs outside the market, like
conversion. Deferred.

## Interface considerations

`union_formation()` mixes timing (`initiation`, `rates`),
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

- Rate-mode cohabitation and conversion are implemented in milestone 4 (D086,
  D093); their separate matching roles would need integration with any new
  two-sex formation model.
- Nonmarital father selection: setting harmonic-mean alphas for unions and
  conditional-logit odds ratios for unpartnered fathers would be inconsistent,
  so father selection, cohabitation and unpartnered partnering must be designed
  together with any two-sex model.

See [cohabitation and nonmarital fathers](#cohabitation-and-nonmarital-fathers-november-13-2026)
for the November 13 discussion of both connections.

The user deferred API restructuring (D073): the current API is unchanged until
a single working implementation incorporates all union-formation components.

## Open questions

- Whether alpha or the odds ratio should be the invariant under composition
  change (they predict differently when group sizes change).
- Treatment of barred pairs (null event or within-cell redraw).
- The denominator: local (Schoen), pooled over groups within ages, global
  F + M, or Dagsvik's F; or a non-mass-action form such as the Choo–Siow
  geometric mean.
- Whether redirection through unpartnered stocks is sufficient, or a
  two-sided outside option is needed for partial active substitution.
- Reference conditions for estimated or calibrated inputs, including
  calibrating alpha to input rates at a reference population.
- How unions, cohabitation and nonmarital fathers share one model.
  Tentatively, cohabitation uses partnership alpha with a couple-based share,
  and nonmarital fathers are drawn from the woman's market row at input
  fertility rates (F1).

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
- Dagsvik, J. K., A. S. Flaatten, and H. Brunborg. 2001. A behavioral two-sex
  marriage model. *Mathematical Population Studies* 9(2):97–121. Discussion
  Paper 238 version (1998) consulted:
  <https://ssb.no/forskning/discussion-papers/a-behavioral-two-sex-marriage-model>
- Cheng, Y., J. K. Dagsvik, X. Han, and Z. Jia. 2026. Force of attraction and
  partner availability in the U.S. marriage market: a two-sided matching
  model. *Sociological Science* 13:178–213. DOI 10.15195/v13.a8.
  <https://sociologicalscience.com/articles-v13-8-178/>
