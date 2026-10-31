# Milestone 8: binned-group hazards and formula weights

Recorded October 31, 2026 (D121). **Status: user-selected future engine direction;
not implemented.** Finish milestone 7 with its current engine. Milestone 8 will
extend the fertility cell approach to all existing rate processes and allow
custom formula weights on those processes. This document records the discussion,
selected principles, proposed diagnostics and questions still needing a contract.
It does not select a public API or authorize an unspecified implementation.

The current engine still uses [cell-accounted fertility](fertility-cell-accounting.md)
with [persistent fertility weights](fertility-heterogeneity-notes.md), alongside
individual/union hazard scheduling for other processes. Calendar-varying rates
remain governed by [the milestone 7 contract](time-varying-rates.md).

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

For upper-open input ages/durations, implementation must represent the attainable
one-year cells without allocating infinitely many cells. The treatment of the
existing fertility small-rate open-band optimization needs explicit review;
it is not automatically a general exception to the selected cell definition.

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

A cell clock generates an event and selection assigns its recipient. The
precise clock implementation, threshold retention and data structures remain to
be designed. It must handle competing events, changes in membership, eligibility
and rates, and preserve the engine's documented state-consistency and invalidation
rules. Calendar changes precede same-time events; the current demographic
priority is death, birth, dissolution, conversion, search, then group transition.
How cell-level selection preserves these priorities needs explicit specification. This remains continuous-time execution, not annual simulation steps.

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
The fertility spacing compensation rule must not silently become a universal
eligibility rule for other mechanisms.

| Process | Reference exposure or distinction to preserve |
| --- | --- |
| Mortality | Living person-years in the specified age/sex/group cell |
| Fertility | All female cell members contribute exposure, including spacing-excluded women; only eligible women receive births |
| Union search | Unpartnered female exposure; rates describe attempts, not successful formations |
| Dissolution | One union exposure unit per active episode; weights may involve both partners |
| Conversion | Current cohabitation exposure, with the documented woman's age or episode-duration clock |
| Group transition | Origin-group exposure for each directed destination risk; destination denominators can share the same members |

Dissolution's current female clock carrier is an implementation detail, not a
reason to restrict future weight predictors to the woman's characteristics.
Couple covariate names and evaluation context remain to be designed.

No feasible recipient or no positive eligible weight sum can make a target
unattainable. Zero/invalid weight policy, reporting and process-specific failure
rules need specification. Existing fertility unrealized-intensity and cancelled-
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

Before implementation, specify:

- Formula API, weight scale and validation, supported person/partner/history
  covariates, persistent heterogeneity interactions and handling of missing values.
- For each process, exact membership, exposure, eligibility, zero-weight behavior,
  affected cells after events and duration-origin preservation/reset rules.
- Aggregate clocks, boundary ordering, destination competition, updates after
  births/deaths/unions/moves, and matching-calibration integration.
- Characteristic/formula update timing and safe evaluation context, including
  continuously changing age/duration predictors and new people.
- Diagnostic provenance, observation semantics and any reporting of infeasible
  target intensity beyond the existing fertility contract.

Validation should include analytical uniform-weight equivalence in distribution,
within-cell hazard-ratio recovery, cell intensity/exposure accounting, selective
removal and changing normalization, formula redundancy/interactions, changing
rates and covariates, duration boundaries, competing destinations, joint union
state changes, infeasible risks, reproducibility and independence between runs.
Use replications for comparisons after changed random-number consumption.
Benchmark performance without assuming cell scheduling is always faster.
These are proposed verification tasks, not user-selected scientific thresholds.

No milestone 7 result or acceptance claim depends on this engine rewrite.
Its current validation report and artifacts must be preserved.
