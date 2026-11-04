# Milestone 8: binned-group hazards and formula weights

Recorded October 31, 2026 (D121). **Status: step 1 (cell scheduling for every
rate process) is implemented on the `milestone8-binned-intensity` branch and
not yet accepted; step 2 (weights and eligibility) is not implemented.**
Milestone 7 was accepted complete with its previous engine (D127, November 2,
2026). Milestone 8 extends the fertility cell approach to all existing rate
processes and will allow custom formula weights on those processes. This
document records the implemented step 1 design, the discussion, selected
principles, proposed diagnostics and questions still needing a contract. It
does not select a public API or authorize an unspecified implementation.

The engine uses [cell-accounted fertility](fertility-cell-accounting.md) with
[persistent fertility weights](fertility-heterogeneity-notes.md), and the same
cell scheduling for the other five processes (see
[step 1 implementation](#step-1-implementation)). Calendar-varying rates remain
governed by [the milestone 7 contract](time-varying-rates.md). Other design
documents that describe individual waiting-time clocks, such as
[waiting-time sampling](waiting-time-sampling.md), predate step 1 and still
need revision.

## Phasing

User direction (November 4, 2026):

1. **Step 1 (implemented):** move every rate process to cell scheduling, then
   measure performance against the previous engine. Fertility keeps its spacing
   eligibility and persistent `hetfert` weights unchanged. The other five
   processes use weights of one and have no eligibility rules. Fertility's cell
   machinery, including its handling of members crossing age boundaries, is
   reused and generalized.
2. **Step 2:** add weighting and eligibility mechanisms. The weight and
   eligibility material below is notes for this step.

One scheduling approach serves all rate processes. The engine will not keep
individual clocks for some processes alongside cells for others. The aim is a
simpler code base. If step 1 performance is unacceptable, the approach may be
abandoned rather than combined with individual clocks.

## Step 1 implementation

This section describes the engine on the milestone 8 branch (November 4,
2026). Choices marked *user direction* were selected by the user during step 1;
the rest are implementation choices within that direction.

### Processes, units and cells

Each rate process keeps cells of at-risk units. A unit belongs to at most one
cell of a process at a time. A cell is one band of a stratum's schedule,
measured on the unit's own clock from its origin.

| Process | Unit | Strata | Clock origin |
| --- | --- | --- | --- |
| Mortality | Living person | Group x sex | Birth |
| Fertility | Living woman (members); eligible women receive births | Group x union type x actual parity (when parity-specific) | Birth |
| Union search (rate mode) | Living unpartnered woman | Group, or group x previous-union category | Birth, or most recent union end for duration schedules |
| Dissolution | Ongoing union, by its history index | Unordered group pair x union type | Union start |
| Conversion | Ongoing cohabitation, by its history index | Wife's group | Union start (duration schedules) or wife's birth (age schedules) |
| Group transition | Living person in the origin group | One process per outgoing-schedule slot k: the origin group's k-th schedule (destination, sex, union type) | Birth, or group entry for duration schedules |

- **Cell width (user direction):** fertility splits bands into cells of at
  most one year and by actual parity, as in D101, because spacing is an
  eligibility rule. The other five processes use their input bands as cells.
  Zero-rate bands are gaps with no cells.
- **Symmetric pairs:** dissolution schedules are defined by unordered group
  pairs, so both orders of a pair (A:B and B:A) share cells.
- **Transitions:** one process per outgoing-schedule slot keeps each origin
  group's competing destinations as separate processes, as the previous
  engine did.

### One clock per process (user direction)

A process does not keep one clock per cell. Its cells are independent Poisson
streams with intensity m_c * N_c, and their sum is one stream at the total
intensity Lambda = sum over cells of m_c * N_c, in which each event belongs to
cell c with probability m_c * N_c / Lambda. Each process therefore keeps:

- an `IntensityTree`, a sum tree whose leaf c holds the cell's current
  intensity. Sums are recomputed from children along each updated path, so they
  do not drift;
- a `ProcessClock`, one Exp(1) threshold spent at the current total intensity.
  Before any intensity change the elapsed interval consumes threshold; the due
  time is then `last + threshold / Lambda`. Membership and rate changes
  therefore recompute the due time without random draws.

When the clock fires, the engine draws one uniform to choose the cell in
proportion to its intensity, and one uniform to choose the unit within it
(uniformly, or by `hetfert` weight in fertility). It then draws a new threshold.

Fertility cells with members but no eligible woman contribute zero to the tree.
Their intensity keeps accumulating as expected unrealized births, exactly as in
D101 (*user direction*: smoothed unrealized births are retained rather than
counting discarded events).

### Boundaries and membership updates

- Each process has a `BoundaryQueue` with at most one pending boundary per
  unit: its next band edge. Boundaries are bookkeeping, consume no randomness,
  and are processed in time order before any event at or after them.
- The queue is a bucket queue (buckets of 1/365 year over a two-year ring, plus
  an overflow heap). Buckets are sorted when reached, so pops follow exactly
  (time, unit) order.
- `refresh_person()` recomputes a person's membership in every process from
  current state. Every event refreshes each person whose state it changes.
  `refresh_union()` does the same for a union's dissolution and conversion
  cells when a union starts or ends and when a partner changes group.

### Event selection and priority

Each loop iteration processes the earliest boundary if it is not later than
the earliest process clock. Otherwise the earliest clock fires. Exact ties
follow the source order death, birth, dissolution, conversion, search, then
transition slots. Calendar changes still precede same-time boundaries and
events. A rate change rebuilds every process from the active schedules and
re-enters every person and union. New thresholds are valid after survival to
the boundary, since the clocks are memoryless.

### Random-number use

At the start, persons are refreshed in index order. Each process's clock draws
its first threshold when its total intensity first becomes positive, so
processes with no positive cells draw nothing. Each event draws the cell and
unit uniforms and the next threshold, followed by the event's own draws (for
example, union type, partner choice and newborn sex). Fixed-seed histories
differ from the previous engine, and models that stratify equal rates
differently (for example, one shared search schedule against three equal
previous-union schedules) create different cells and so differ draw by draw,
while agreeing in distribution.

### Removed with step 1

The per-person event queue and its event versions, `schedule()`, the stored
per-union dissolution and conversion clocks, the special handling of ties
between cell clocks and person events, the individual waiting-time sampler
`waiting_time()` and its internal test entry point `hazard_wait_cpp()`.

The matching-calibration intensity tracker (`FemaleIntensityTracker`) is
retained. It keeps its own per-woman boundaries for weighted matching. Feeding
it from the search cells is a possible later simplification.

### Validation and performance (November 4, 2026)

These are software checks, not user-selected acceptance thresholds.

- **Equivalence with the previous engine:** a model with every process (two
  groups, mortality, fertility, rate-mode search, cohabitation and conversion,
  typed duration dissolution, two-way transitions) was run for ten seeds on
  each engine. All 23 summary measures agree within replication noise (largest
  |z| = 2.05), and observed dissolution, conversion and transition rates match
  their inputs.
- **Fertility:** the milestone 5 cell-spacing recovery check agrees with the
  accepted results across its 29 comparisons (largest |z| = 2.4), including
  unrealized births.
- **Mortality:** observed single-year death rates match Gompertz inputs
  (chi-square over about 155 age-by-sex cells, p = 0.30).
- **Tests:** the suite passes, also with internal consistency checks that
  rescan every process's cells, and its tree, at each event. Tests that
  constructed exact ties from the previous draw order were rewritten for the new
  order; tests asserting draw-for-draw identity between differently stratified
  models now compare distributions.

Timings on the development machine (seconds, optimized builds):

| Model | Previous engine | Step 1 | SOCSIM 1.9.18 |
| --- | --- | --- | --- |
| Full model above, 50K people, 100 years | 7.7 | 9.2 | not run |
| Mortality and fertility, 50K people, 100 years | 0.51 | 1.29 | 2.00 |
| Mortality and fertility, 100K people, 200 years, replacement fertility | 3.86 | 7.2 | 13.9 |
| Demographic transition article, main run (yearly rates, on-demand; timed by the user) | 42 | 32 | not run |
| Second demographic transition article, three main runs (rate mode) | 4.6 | 6.3 | not run |

The SOCSIM runs used comparable inputs with random fathers instead of matching.
The second demographic transition timings are the means of two runs of the
article's code, extracted with `knitr::purl`.
Per-boundary updates are the main added cost; calendar rate changes and
fertility became cheaper.

## Purpose and selected direction

Separate two aspects of each demographic rate process:

1. The input schedule determines the aggregate event intensity within a cell.
2. Individual or union weights determine which eligible member experiences an
   event, preserving relative risks within that cell.

This generalizes the current fertility mechanism. Users can supply detailed
rate strata, use formulas to distribute risk according to characteristics, or
combine the two. Formula weighting should preserve the input cell intensity,
rather than let the cell's covariate composition multiply its aggregate rate.

The intended scope is mortality, fertility, rate-based union-search attempts,
union dissolution, cohabitation-to-marriage conversion and directed individual
group transitions. These are the six existing inputs in [the rate census](rate-inputs.md).
On-demand partner matching remains triggered by a birth attempt; it is not a
seventh input-rate process. Conditional partner choice and group inheritance
are separate mechanisms, not replaced by event-weight formulas.

A binned group here means a demographic rate cell. It need not correspond to
one named population `group`, nor require a grouped population.

## Cell boundaries selected by the user

User direction (November 4, 2026), revising the original single-year rule:
cells are split to single years only for a process with eligibility rules,
which is currently fertility alone (spacing; D101 decision 5). The other five
processes use their input bands as cells. Without eligibility rules every
member's hazard equals the cell rate whatever the cell's width, so splitting
cannot change results and only adds boundary updates. In the second
demographic transition article, single-year splitting roughly doubled the cost
of a rate-mode model with five-year bands. Weights will not depend on age
within a band. Step 2 must still consider weights on traits correlated with
age inside a wide band, which shift events toward the ages where the trait is
common while keeping the band total.

The rules below describe single-year splitting where it applies:

- Use single-year attained-age or elapsed-duration bands, even when input bands
  are wider. The process keeps its existing choice of clock; this does not add
  a joint age-by-duration surface.
- Preserve every other stratum in the input schedule: sex, group, union type,
  previous-union category, directed origin/destination pair or symmetric group
  pair, as applicable to that process.
- Retain original fractional boundaries. Cells span at most one year: input
  `[15,20)` gives `[15,16)`, ..., `[19,20)`; `[15.5,17)` gives `[15.5,16)` and
  `[16,17)`. Duration bands split the same way, measured from each member's
  own duration origin.
- Retain fertility's existing actual-parity cells when parity-specific schedules
  are supplied. An input band covering several parities applies to each
  contained actual parity, as in D101.
- Calendar periods select active schedules and trigger updates, as they do now.
  They are not a persistent additional cross-classification of members by
  historical period. Past events and exposure retain their recorded dates.

An input age band of five years crossed with two sexes and two input groups
therefore yields 20 potential cell definitions, before any additional
input strata. Formula predictors do not automatically create further cells.

For upper-open input ages in fertility, single-year cells are created up to the
oldest attainable age in the run and the rest stays one open cell. Fertility's
small-rate open-band tolerance (an open band stays one cell when spacing could
affect only a negligible share) is retained. The other processes keep each
upper-open band as one cell.

## Aggregate intensity and relative individual risk

For cell c, let:

- m_c be its effective input rate, after applicable input-rate multipliers;
- N_c be its current number of reference exposure units;
- E_c be the members currently eligible to receive the event;
- w_i be a member's relative selection weight.

The proposed shared mathematical contract is:

    Lambda_c = m_c * N_c
    P(i | event in c) = w_i / sum(w_j for j in E_c)
    h_i = Lambda_c * w_i / sum(w_j for j in E_c), for i in E_c

Ineligible members have zero event hazard. With a positive eligible weight sum,
summing individual hazards gives Lambda_c. Input-rate recovery means this
conditional aggregate intensity, with observed counts fluctuating stochastically
around their exposure-based expectations. It does not promise exact event quotas
or exact observed rates in every finite cell or time interval.

A process clock generates an event, a cell is chosen in proportion to its
intensity, and selection assigns the recipient; see
[step 1 implementation](#step-1-implementation) for the clocks, data
structures, competing processes and same-time priority. This remains
continuous-time execution, not annual simulation steps.

With all eligible weights one and N_c equal to eligible membership, each member
has hazard m_c. Aggregate scheduling can then reproduce the current independent
hazard model in distribution, provided boundaries and competing events are
handled correctly. Fixed-seed histories need not match after the algorithm changes.

### Worked mortality example

Two people share a cell with mortality rate 0.01 per person-year and weights
1 and 3. The cell intensity is 0.02 deaths per year. Selection probabilities
are 1/4 and 3/4, giving hazards 0.005 and 0.015. The hazard ratio is three;
the average hazard is the supplied 0.01.

After a death, membership and normalization change. If the higher-weight member
dies, the remaining member has hazard 0.01 while alone in the cell, rather than
retaining 0.005. This composition dependence is part of preserving the cell target.

## Interpretation of formulas

A log-linear form was discussed:

    w_i = exp(x_i * beta)

It gives h_i / h_j = exp((x_i - x_j) * beta) for eligible members of the same
cell at the same instant. Formula syntax, coefficient storage, supported terms
and whether arbitrary weight expressions will be allowed are not selected.
The log-linear example is not an implemented API or a demographic coefficient default.

The input rate is a cell-average target, not a reference individual's hazard.
Normalization can change with time and composition, including selective deaths
or transitions. A cell that loses high-weight members continues to experience
its prescribed aggregate intensity. Formula effects are conditional within-cell
relative hazards; comparisons across cells also involve their schedules and
normalization factors.

For comparison, ordinary unnormalized proportional hazards would use
h_i = m_c * w_i and total intensity m_c * sum(w_i). That alternative permits
composition to alter aggregate rates. It was discussed to distinguish computation
from demographic interpretation; it is not the selected target-preserving
weight mechanism for milestone 8.

User direction (November 4, 2026): the first milestone 8 version offers
normalized weights only. An unnormalized option may be considered as a later
feature, but it needs its own contract first. A consequence to record as a
limitation: normalized weights cannot produce frailty selection. Under
normalization, a cell keeps its input rate however its members' weights shift.
As high-weight members die, the survivors' hazards rise to compensate. Old-age
mortality deceleration or cross-overs arising from selection therefore cannot
come from weights. The current engine has no frailty mechanism either, so this
is not a regression.

### Cell composition and population-level effects

Within a cell, the ratio of two eligible members' hazards equals their weight
ratio exactly. In a cell where a share p has a binary trait with weight 1.1 and
the rest weight 1, members with the trait have hazard 1.1 m_c / (1 + 0.1 p) and
the others m_c / (1 + 0.1 p). The mix changes the levels, not the ratio.

A cell whose members all share one value has no contrast: every member gets
m_c whatever the weight. A singleton cell is the extreme case. Such cells pull
population-level estimates of the effect toward no effect. A model fitted to
simulated output therefore recovers the weight's effect only to the extent that
cells are large and mixed in the weighted trait. The same model can show
smaller effects in a small population than in a large one. Persistent fertility
weights are already subject to this. Coarser cells for weighted processes are a
possible later remedy if the effect proves to matter in practice; they are not
selected.

Existing proportional input-rate multipliers act on m_c and therefore can change
cell targets. They are conceptually distinct from normalized formula weights.
Their API relationship to the new weight interface remains to be specified.

## Two ways to express stratification

Stratifying a rate table creates separate aggregate targets. Adding a formula
predictor redistributes events within those targets. These approaches need not
give the same model, even if they use the same characteristic.

A positive factor common to every eligible member of a cell cancels during
normalization. In log-linear weights, a main effect of an existing input-rate
stratum therefore has no additional effect. This was explicitly recognized by
the user as a potential specification mistake worth warning about.

| Log-weight formula term | If input rates already stratify by group |
| --- | --- |
| `group` | Main effect cancels within each cell |
| `education` | Redistributes risk where education varies within the cell |
| `group + education` | Group main effect cancels; education can remain effective |
| `group:education` | Can specify different education effects within each group |

### Proposed warning behavior

The user requested consideration of a redundancy warning. The following policy
was discussed, but its exact implementation and diagnostic rules remain proposals:

- Warn once for a user-specified term whose cancellation is structurally provable.
- Identify the ineffective term, not the whole formula, and allow other terms to
  retain their effect. Execution can continue.
- Do not warn for the ordinary intercept, whose cancellation is expected.
- Do not infer an error merely because a realized small cell happens to contain
  identical covariate values.
- Do not reject or warn about every reference to a stratification variable:
  interactions can remain effective. Arbitrary custom expressions require
  algebraic care; a constant additive contribution to raw weights need not
  cancel like a common multiplicative factor.

Illustrative wording: "The main effect of group has no effect on relative
hazards because the input rates already stratify by group. Weights are
normalized separately within each cell."

## Process-specific exposure and eligibility

A shared scheduler must retain each process's demographic denominator and unit.
The fertility spacing eligibility rule must not silently become a universal
eligibility rule for other mechanisms.

### Proposed general eligibility mechanism

Status: proposal awaiting the user's acceptance. The user wants eligibility as a
general mechanism, so that rules can be added later without redesigning the
scheduler. Milestone 8 adds no new eligibility rules; fertility spacing remains
the only one. A possible future rule is a grace period after divorce or
widowhood before a new union search.

Proposed rules:

- **Rates are per member, not per eligible member.** Lambda_c = m_c * N_c counts
  all reference members, including ineligible ones, as fertility does now
  (D101). Eligible members' hazards rise to carry the cell's full intensity.
  This keeps output rates equal to input rates when inputs are observed rates
  for the whole population at risk. A future rule may choose per-eligible-member
  rates only through its own decision.
- **A weight of zero means ineligible.** Weights and eligibility are one
  mechanism.
- **A cell with no eligible member produces no events,** and its intensity is
  reported as unrealized, as with fertility.
- **Mortality has no eligibility rules,** and any mortality weights must be
  strictly positive, so that no cell can make its members immortal.
- **Model checks should warn when a rule and a schedule conflict,** for example
  when a grace period covers a whole duration band whose input rate is nonzero.

Worked example (illustrative rule, not implemented): the age-60 widowed
search rate is 0.10 per year. The cell has 100 widows, 20 widowed less than a
year ago and in a one-year grace period. The cell produces 100 x 0.10 = 10
searches per year, shared by the 80 eligible widows at 0.125 each. Observed
output is 10 searches per 100 widow-years, matching the input. Per-eligible-
member rates would instead give 8 per year and an observed rate of 0.08.

Duration-specific search schedules can already express a delay after a union
ends: a schedule with `[0,1): 0` gives no searches in the first year
([remarriage](remarriage.md)). A grace-period rule would add something new only
for age-based schedules.

| Process | Reference exposure or distinction to preserve |
| --- | --- |
| Mortality | Living person-years in the specified age/sex/group cell |
| Fertility | All female cell members contribute exposure, including spacing-excluded women; only eligible women receive births |
| Union search | Unpartnered female exposure; rates describe attempts, not successful formations |
| Dissolution | One union exposure unit per active episode; weights may involve both partners |
| Conversion | Current cohabitation exposure, with the documented woman's age or episode-duration clock |
| Group transition | Origin-group exposure for each directed destination risk; destination denominators can share the same members |

Dissolution and conversion units are unions, so future weight predictors need
not be restricted to the woman's characteristics. Couple covariate names and
evaluation context remain to be designed.

No feasible recipient or no positive eligible weight sum can make a target
unattainable. The zero-weight proposal above awaits acceptance; invalid-weight
policy, reporting and other process-specific failure rules need specification. Existing fertility unrealized-intensity and cancelled-
birth rules remain current behavior; they are not automatically generalized or
removed. Failed matching can still reduce union formation or realized births,
even when their attempt intensities are preserved.

## User-defined characteristics

The user identified formula weights as a route to using future individual
characteristics such as education, income or health without supplying every
cross-classified rate schedule. Rate formulas determine how current values affect
relative event risk; characteristic mechanisms determine the values themselves.

This does not implement characteristic initialization, newborn assignment,
inheritance, longitudinal updates or missing-value handling. Those contracts
remain separate. Formula evaluation must not permit arbitrary callbacks to mutate
engine state or bypass event/cell invalidation rules. Multiple independent group
columns are separately assigned to tentative milestone 10 (D121).

## Scheduling, implementation and validation work still required

Changing weights alone does not change Lambda_c when eligibility and reference
membership are unchanged. This could permit evaluation of current covariates at
recipient selection instead of redrawing individual clocks. It is a computational
opportunity, not a selected evaluation policy or demonstrated speed improvement.
Time-dependent formulas and efficient weighted selection need investigation.

Before step 2 implementation, specify:

- Formula API, weight scale and validation, supported person/partner/history
  covariates, persistent heterogeneity interactions and handling of missing values.
- For each process, exact membership, exposure, eligibility, zero-weight behavior,
  affected cells after events and duration-origin preservation/reset rules.
- How weighted selection and eligibility fit the step 1 clocks, trees and
  updates, including matching-calibration integration.
- Characteristic/formula update timing and safe evaluation context, including
  continuously changing age/duration predictors and new people.
- Diagnostic provenance, observation semantics and any reporting of infeasible
  target intensity beyond the existing fertility contract.

Validation should include analytical uniform-weight equivalence in distribution,
within-cell hazard-ratio recovery, how much population-level effects shrink
across population sizes, cell intensity/exposure accounting, selective
removal and changing normalization, formula redundancy/interactions, changing
rates and covariates, duration boundaries, competing destinations, joint union
state changes, infeasible risks, reproducibility and independence between runs.
Use replications for comparisons after changed random-number consumption.
Benchmark performance without assuming cell scheduling is always faster.
These are proposed verification tasks, not user-selected scientific thresholds.

No milestone 7 result or acceptance claim depends on this engine rewrite.
Its current validation report and artifacts must be preserved.
