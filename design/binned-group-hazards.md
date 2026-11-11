# Milestone 8: cell scheduling and weight models

Updated November 5, 2026. **Milestone 8 is accepted complete (D130).** Every
rate process is scheduled by cells with one clock per process (D128), and
weight models choose who in a cell experiences an event (D129). This is the
current specification; the [scientific contract](scientific-contract.md)
summarizes it, and the user's [validation report](../validation/milestone8.qmd)
is the acceptance evidence. The direction was selected on October 31, 2026
(D121); the earlier discussion and the step-by-step development record are in
Git history.

The work was done in two steps at the user's direction. Step 1 moved all six
rate processes to cell scheduling with weights of one, and was benchmarked
against the previous engine and SOCSIM before continuing. Step 2 added weight
models. One scheduling approach serves all rate processes: the engine keeps no
individual waiting-time clocks, so as to simplify the code base. Fertility
keeps its spacing eligibility (D101) and persistent `hetfert` weights (D102).
Calendar-varying rates remain governed by [their contract](time-varying-rates.md).

## Cell scheduling (D128)

Choices marked *user direction* were selected by the user during the
milestone; the rest are implementation choices within that direction.

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

*Unrealized births are counted (user direction, November 4, 2026, revising
D101's reporting and an earlier step 1 choice to keep expected values).* Every
fertility cell contributes rate * members, including cells whose members are
all within the birth interval. When the clock draws a cell with no eligible
woman, no child is born and one unrealized birth is recorded for its one-year
period from the run start, group, union type, parity band and age band. In a
cell shared by groups, the draw that would pick a mother picks the group in
proportion to its members. Realized births are unchanged in distribution;
`simulation$unrealized_births` now holds counts rather than expected values,
`years_without_eligible` is removed, and records no longer split at fertility
rate changes (the per-epoch tables are gone). The accepted milestone 5 report
used the earlier expected values. Simpler code and output motivated the
change; it also removes the unrealized accounting that complicated in-place
rate updates.

### Boundaries and membership updates

- Each process has a `BoundaryQueue` with at most one pending boundary per
  unit: its next band edge. Boundaries are bookkeeping, consume no randomness,
  and are processed in time order before any event at or after them.
- The queue is a bucket queue (buckets of 1/365 year over a 16-year ring, plus
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
events.

*In-place rate updates (user direction, November 4, 2026).* At a calendar rate
change, a process whose new schedules keep its cells (same strata and nonzero
band edges; for fertility also the same union types, group stratification and
age splitting, and for transitions the same destination, sex and union type
per slot) changes its cells' rates in place: the clock spends its threshold up
to the change, the bins take the new rates, and the clock continues, with no
draw and no one moving. This is exact, because a Poisson process whose rate
changes at a known time continues at the new rate. Otherwise that part is
rebuilt from the active schedules and every living person and ongoing union
re-enters; the new thresholds are valid after survival to the boundary. With
yearly rates this cut the demographic transition article's main run from 30.9
to 25.3 seconds and the second demographic transition's changing-rate run
from 2.66 to 1.45 seconds (now faster than the previous engine's 1.82).

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

### Removed from the previous engine

The per-person event queue and its event versions, `schedule()`, the stored
per-union dissolution and conversion clocks, the special handling of ties
between cell clocks and person events, the individual waiting-time sampler
`waiting_time()` and its internal test entry point `hazard_wait_cpp()`.

*Calibration tracker replaced (user direction, November 4, 2026, for
simplicity).* The former `FemaleIntensityTracker` kept a second, per-woman
record of search intensity with its own boundary queue, duplicating the search
cells, and an unused fertility channel. Matching calibration now takes all
search intensities from the search cells (SearchCalibration, with weights of
one when search is unweighted) and fertility intensities from the fertility
cells. A small `ProbabilityTracker` remains only for each unpartnered woman's
cohabitation probability, whose age bands can differ from the search clock;
it exists only with rate-mode search and cohabitation. Calibrated runs gave
identical histories before and after the change, and run time was unchanged.
`inspect_intensity_tracker_cpp()` and its tests were removed; a
simulation-level test exercises every calibration path under
`hammel.check_cells`.

### Step 1 validation and performance (November 4, 2026)

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

## Weight models (D129)

Specified November 4, 2026 from decisions the user made in discussion.

### Scope

Weights apply to mortality, fertility, rate-mode union search,
dissolution, cohabitation-to-marriage conversion and group transitions.
Eligibility rules are out of scope: every weight is strictly positive.

*Calendar schedules (user direction, November 4, 2026).* Any `weights` or
`weights_conversion` argument may be a tibble with `time_start`, `time_end`
and a list column `model` of weight models, the structure of the matching and
inheritance model schedules (D123), validated the same way: half-open
intervals that do not overlap and cover the run, each entry a complete model.
A schedule's models must all be formulas or all custom functions, because the
evaluation path is set per process; a formula with zero coefficients turns
weighting off for an interval. At a boundary the engine switches the active
model and recomputes every stored weight. Weights never change cell
intensities, so event clocks are kept and no random draws are used: a
schedule of identical models reproduces the static model exactly.

### Weights and selection

For a unit i in cell c, the weight model gives a log weight eta_i and the
weight w_i = exp(eta_i). The cell's intensity is unchanged:

    Gamma_c = m_c * N_c
    P(i | event in c) = w_i / sum(w_j for eligible j in c)

Weights decide who in the cell experiences the event; they never change a
cell's intensity, the process clock or which cell is chosen. exp(beta) for a
term is the hazard ratio between two members of the same cell that differ by
one unit in it. Population-level ratios can be smaller where cells are small
or homogeneous (see [cell composition](#cell-composition-and-population-level-effects)).
In fertility only eligible women are weighted and chosen, as now.

### API

    weight_model(formula, coefficients)
    weight_model(fun = f)

The formula form mirrors `matching_model()`: a one-sided formula, a
finite named numeric vector of coefficients.
Categorical variables use treatment coding with fixed levels, and coefficient
names follow the model-matrix column names (for example `raceB`,
`union_typecohab`). Group variables are named by declared dimension (D134,
stage 5), and a group-valued variable's levels are its dimension's categories
(declared, or in order of first appearance), so variables of one dimension
share its reference: `~ race + partner_race` with categories A and B has
coefficients `raceB` and `partner_raceB`. Names that could be group variables
are resolved, and their coefficients checked, when the model meets the
declaration at simulation; other names are checked when the model is set. Other categorical variables have
engine-fixed levels, with the first listed in [Variables](#variables) as the
reference. The intercept has no effect and is ignored. Non-finite
coefficients are errors.

Formula terms are main effects of supported variables and their interactions
(`:` or `*`). Transformations such as `I()` or `log()` are not supported; a
custom function covers anything else.

Weight models attach to processes through new arguments:

| Process | Argument |
| --- | --- |
| Mortality | `set_mortality(model, rates, weights = NULL)` |
| Fertility | `set_fertility(..., weights = NULL)` |
| Dissolution | `set_union_dissolution(..., weights = NULL)` |
| Conversion | `set_union_formation(..., weights_conversion = NULL)` |
| Group transitions | `set_group_transition(model, rates, weights = NULL)` |
| Union search (rate mode) | `set_union_formation(..., weights = NULL)` |

`NULL` keeps unweighted selection.

### Variables

Person processes (mortality, fertility, transitions) can use:

| Variable | Values |
| --- | --- |
| `sex` | `"f"`, `"m"` |
| `<dim>` (e.g. `race`) | the person's category in a declared dimension |
| `union_type` | `"unpartnered"`, `"cohab"`, `"marriage"` |
| `partnered` | 1 if in a union, else 0 |
| `parity` | number of births (women; 0 for men) |
| `previous_union` | `"none"`, `"separated"`, `"widowed"` |
| `partner_<dim>` | partner's category; NA when unpartnered |
| `ancestry_<dim>_<category>`, `ancestry_<dim>_unknown` | ancestry shares of a tracked dimension, each in [0, 1] |
| `destination` | transitions only: the destination of the risk in the transition's dimension (formerly `group_dest`) |

Union processes (dissolution, conversion) use partner variables with the
suffix `_w` (wife) or `_h` (husband): `<dim>_w`, `<dim>_h`, `parity_w`,
`previous_union_w`, `previous_union_h`, `ancestry_<dim>_<category>_w` and
`ancestry_<dim>_<category>_h` (births are counted for mothers only, so there is no
`parity_h`), plus the couple terms:

| Variable | Values |
| --- | --- |
| `same_<dim>` | 1 if both partners are in the same category of the dimension, else 0 |
| `age_difference` | husband's age minus wife's age, fixed for the union |
| `union_type` | `"cohab"` or `"marriage"` |

All variables change only at events, so a unit's weight is constant between
its refreshes. Continuously changing quantities (exact age, durations) are
excluded, as are user-supplied population columns, which belong to a separate
characteristics feature. The ancestry shares sum to 1; a formula using all of
them is collinear, which is harmless with fixed coefficients.

Parity is tracked whenever a weight model uses it, as well as for
parity-specific fertility.

### Cell stratification

*User direction (November 4, 2026):* cells stratify only as finely as the input
rate schedules. A process whose active rates have no group-specific schedules
pools all groups into one set of cells, so group terms in its weight model take
effect. Tables from `stratify_rates()` are group-specific even when every
factor is one (D132). Pooled fertility cells still report unrealized births by
group: each group's share is the cell rate times its members, so the reported
values match per-group cells. A woman's matching-calibration bin uses her own
group. Group transitions always stratify by origin group, and dissolution by
group pair only when its rates do. Pooling changes cell structure, and so
fixed-seed histories, for grouped models with shared rates; results are
equivalent in distribution.

### Missing values

- A formula term whose variable is NA (for example `partner_group` for an
  unpartnered person) contributes zero to eta. Other terms keep their effect.
- An NA (or NaN) returned by a custom function gives eta = 0 for that unit.
- Either case warns once per run, naming the process and variable or function.
- An infinite value from a custom function is an error.

### Custom functions

`fun(members, bin)` is called with:

- `members`: a tibble with one row per eligible unit in the chosen cell. It
  holds `person_id` (or `union_id`, `wife_id`, `husband_id` for union
  processes) and every variable listed above, under the same names.
- `bin`: a list describing the cell: `process`, its strata (for example sex,
  group, union type, parity range, previous-union category, group pair or
  destination) and its band, as `age_start`/`age_end` or `dur_start`/`dur_end`,
  plus `time_start`/`time_end` of the active rate period (NA for static rates).

No exact times are passed: a function can vary weights by band but not by
exact age or duration within a band.

It must return a numeric vector of log weights, one per row. A wrong length,
a non-numeric result or an infinite value is an error; NA follows the rule
above. The function cannot change engine state. It shares the engine's random
number stream, as custom inheritance does, so runs remain reproducible from
the seed.

### Evaluation

- **Formulas: stored at refresh** (*user direction, November 4, 2026*,
  replacing selection-time evaluation after timing). Formula variables change
  only at events, so each unit's weight is computed when it is refreshed and
  stored. A chosen cell picks its member by rejection sampling: a uniform
  candidate is accepted with probability w / bound, where the cell keeps an
  upper bound on its members' weights (raised as weights arrive, reset when
  the cell empties, and reset to the exact largest weight after 32
  consecutive rejections, so a bound left by a departed extreme-weight member
  cannot stall selection). This is exact for any bound at or above the largest
  weight; a candidate at the bound is accepted without a further draw, so
  equal weights keep the unweighted pick and its random-number use. No
  structure is updated when units change cells.
- **Custom functions: at selection.** They are called once per event for the
  chosen cell's eligible members, because they can depend on the band.
- **Search with matching calibration:** a searching woman's hazard is
  rate_c * N_c * w_i / W_c, so her calibration contribution depends on her
  whole cell. Search weights are then computed when each woman is refreshed and
  stored. Each search cell keeps, per calibration bin, its member count and
  the sums of w (1 - p) and w p, where p is her current cohabitation
  probability (reported by the probability tracker, which also handles
  probability-band boundaries). A bin's marriage and cohabitation intensities
  are the sum over cells of rate_c * N_c / W_c times those sums. Cells whose
  contents changed are recomputed before the next search solve, and the two
  channels are assigned to the calibrator, as fertility's channel already is.
  Searchers are chosen by the same rejection sampling. Without search weights
  the same sums are kept with weights of one. A custom search function is called at
  every refresh of every searching woman, with a warning. An internal check
  (`hammel.check_cells`) compares the maintained totals with an explicit scan
  before each solve, recomputing each woman's probability from her schedule.
- **Fertility with matching calibration:** when weighted matching calibration
  is active (more than one group and group terms in matching), calibration
  needs the weighted intensity of unpartnered women by calibration bin at all
  times. Fertility weights are then computed when each woman is refreshed and
  kept in the cell's existing weight tree and calibration sums, as `hetfert`
  is now. A custom function is then called at every refresh of every woman;
  simulation warns that this can be slow.
- **`hetfert`** (internal only, D176): combines additively on the log scale,
  eta_i = log(hetfert_i) + x_i * beta, so accepted milestone 5 behavior is
  unchanged when no fertility weight model is given. Replacing `hetfert` by a
  formula term may follow once persistent, inherited characteristics exist.
- **Equal weights reproduce unweighted runs:** a formula whose coefficients
  are all zero is skipped. When every eligible member of the chosen cell has
  the same weight, the member is chosen with the unweighted rule and the same
  uniform draw, so a custom function returning a constant also reproduces the
  unweighted run.

### Cost (November 4, 2026)

Optimized builds, development machine:

| Model | Without weights | With formula weights |
| --- | --- | --- |
| Rate-mode second-transition model, two groups, 150K people, 50 years; search weights under calibration | 7.91 s | 8.53 s (+8%) |
| Same, 20K people; custom search function | 1.01 s | 2.73 s |
| Mortality and fertility, 100K people, 100 years; mortality weights | 9.25 s | 9.07 s (no measurable cost) |

Models without weights are unaffected. Earlier designs cost more: scanning the
whole chosen cell at every event (+29% for the mortality model, +78% for
search) and per-cell sum trees, which must be updated at every birthday move
(+19% for mortality). Search's remaining cost is mostly the calibration sums.

### Diagnostics

Simulation warns once for each formula term that cancels within every cell
of its process, because the process's rates already stratify by its
variables. For example, `group` in mortality, `union_type` in dissolution, or
`parity` in parity-specific fertility. The warning names the term. Interactions
with a non-cancelling variable do not warn. Custom functions are not
inspected.

### Worked examples

*Mortality.* A cell (men aged 70 in group A) has 100 members, 60 of them
partnered, and death rate m. The formula `~ partnered` with coefficient
log(0.8) gives weights 0.8 and 1. The cell's intensity stays 100 m. A
partnered man has hazard 100 m x 0.8 / (60 x 0.8 + 40) = 0.909 m, and an
unpartnered man 1.136 m. Their ratio is 0.8, and the cell's average hazard is
still m.

*Dissolution.* `~ age_difference` with coefficient log(1.05) makes a union
with a five-year age gap 1.05^5 = 1.28 times as likely to be the one dissolved
as an otherwise similar union with no gap in the same cell. Terms such as
`same_group` cancel when dissolution rates stratify by group pair, because
every union in a pair cell shares the same value.

*Fertility with `hetfert`.* Fertility rates are not stratified by union type,
and the weight model is `~ partnered` with coefficient log(1.2). A partnered,
eligible woman with hetfert 1.5 has eta = log(1.5) + log(1.2) and weight 1.8,
against weight 1.2 for a partnered woman with hetfert 1 and 1 for an
unpartnered woman with hetfert 1 in the same cell. With union-type-specific
fertility rates, `partnered` would cancel and trigger the warning.

### Validation

The user's [milestone 8 report](../validation/milestone8.qmd) assesses rate
recovery and within-cell hazard ratios for every process, population-level
attenuation and performance. The test suite checks every cell's input rate,
exact reproduction of unweighted runs with zero coefficients or constant
custom weights, the cancellation warnings, NA handling, the custom-function
contract, weight schedules, rejection sampling and RNG reproducibility, and
compares maintained cell and calibration totals with explicit scans under
`hammel.check_cells`. These are software checks, not user-selected
acceptance thresholds.

## Aggregate intensity and relative individual risk

For cell c, let:

- m_c be its effective input rate, after applicable input-rate multipliers;
- N_c be its current number of reference exposure units;
- E_c be the members currently eligible to receive the event;
- w_i be a member's relative selection weight.

The contract is:

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
[cell scheduling](#cell-scheduling-d128) for the clocks, data structures,
competing processes and same-time priority. Execution is continuous, not in
annual steps.

With all eligible weights one and N_c equal to eligible membership, each member
has hazard m_c, which reproduces the previous independent-hazard engine in
distribution; fixed-seed histories differ.

### Worked mortality example

Two people share a cell with mortality rate 0.01 per person-year and weights
1 and 3. The cell intensity is 0.02 deaths per year. Selection probabilities
are 1/4 and 3/4, giving hazards 0.005 and 0.015. The hazard ratio is three;
the average hazard is the supplied 0.01.

After a death, membership and normalization change. If the higher-weight member
dies, the remaining member has hazard 0.01 while alone in the cell, rather than
retaining 0.005. This composition dependence is part of preserving the cell target.

## Interpretation of formulas

Formula weights are log-linear, w_i = exp(x_i * beta), so
h_i / h_j = exp((x_i - x_j) * beta) for eligible members of the same cell at
the same instant. Coefficients are user inputs; there are no defaults.

The input rate is a cell-average target, not a reference individual's hazard.
Normalization can change with time and composition, including selective deaths
or transitions. A cell that loses high-weight members continues to experience
its prescribed aggregate intensity. Formula effects are conditional within-cell
relative hazards; comparisons across cells also involve their schedules and
normalization factors.

For comparison, ordinary unnormalized proportional hazards would use
h_i = m_c * w_i and total intensity m_c * sum(w_i), letting composition alter
aggregate rates.

User direction (November 4, 2026): milestone 8 offers normalized weights only. An unnormalized option may be considered as a later
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
selected. Composition change over calendar time has the same effect on
estimates pooled across periods; estimating within strata (for example age by
sex by year) and excluding an open top age group recovers the within-cell
ratio, as the validation report shows.

Proportional stratification (`stratify_rates()`, D132) acts on m_c and so
changes cell targets; it is distinct from normalized weights.

## Two ways to express stratification

Stratifying a rate table creates separate aggregate targets. Adding a formula
predictor redistributes events within those targets. These approaches need not
give the same model, even if they use the same characteristic.

A positive factor common to every eligible member of a cell cancels during
normalization. In log-linear weights, a main effect of an existing input-rate
stratum therefore has no additional effect. This was explicitly recognized by
the user as a potential specification mistake worth warning about.

Here `education` stands for a hypothetical future characteristic:

| Log-weight formula term | If input rates already stratify by group |
| --- | --- |
| `group` | Main effect cancels within each cell |
| `education` | Redistributes risk where education varies within the cell |
| `group + education` | Group main effect cancels; education can remain effective |
| `group:education` | Can specify different education effects within each group |

Simulation warns about such terms (see [Diagnostics](#diagnostics)).

## Process-specific exposure and eligibility

A shared scheduler must retain each process's demographic denominator and unit.
The fertility spacing eligibility rule must not silently become a universal
eligibility rule for other mechanisms.

| Process | Reference exposure or distinction to preserve |
| --- | --- |
| Mortality | Living person-years in the specified age/sex/group cell |
| Fertility | All female cell members contribute exposure, including spacing-excluded women; only eligible women receive births |
| Union search | Unpartnered female exposure; rates describe attempts, not successful formations |
| Dissolution | One union exposure unit per active episode; weights may involve both partners |
| Conversion | Current cohabitation exposure, with the documented woman's age or episode-duration clock |
| Group transition | Origin-group exposure for each directed destination risk; destination denominators can share the same members |

## Deferred issues

These are recorded for future work; none is selected or implemented.

- **Eligibility rules beyond fertility spacing** (user direction, November 4,
  2026: deferred until a concrete rule is needed, such as a grace period after
  divorce or widowhood before a new search). The proposed direction keeps rates
  per member, as fertility does (D101): eligible members carry the cell's full
  intensity, and a cell with no eligible member would report unrealized events.
  Duration-specific search schedules can already express a delay after a
  union ends with a zero-rate first band ([remarriage](remarriage.md)).
- **`NULL` entries in weight schedules,** to turn weighting off for an
  interval, only if offered in parallel for matching, inheritance and weight
  schedules (roadmap, milestone 11). A zero-coefficient formula does this now.
- **Unnormalized weights,** which would permit frailty selection, need their own
  contract.
- **User-defined characteristics** (education, income, health) as weight
  variables need initialization, inheritance, updating and missing-value
  contracts; user-supplied population columns are not weight variables.
  Multiple group columns are milestone 9 (D121, D131).
- **Weights that vary with exact age or duration within a band** are not
  offered; a custom function can vary weights by band only.
