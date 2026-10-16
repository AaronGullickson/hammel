# Milestone 2: union formation

Updated October 15, 2026 from the user's milestone 2 request.
**Status: milestone scope accepted (D030); detailed demographic rules below
remain proposals or open decisions. Process construction and formula scoring
and conditional-logit on-demand execution with minimum male age and known-kin
screening are implemented; exclusive ongoing unions are enforced; rate-based female search initiation is implemented.**
Every simulation requires a union-formation process. Current-partner lookup
and union histories are common infrastructure for all initiation modes.

## Terminology

Use **union**, **union formation**, and **partnered/unpartnered** in code and
documentation. The user sometimes says “marriage” conversationally, but the
model is deliberately agnostic about legal marriage, long-term cohabitation,
and other union types. Milestone 2 does not introduce a legal-marital-status
classification. At the user's request, union tables nevertheless use explicit
`wife_id` (female) and `husband_id` (male) roles so sex is clear. Conditional fertility refers to whether a woman is in an
ongoing union, subject to the detailed exposure rules still to be specified.

## Accepted scope

Separate a woman's decision to search from the process that selects a partner.
Both initiation approaches must share eligibility screening and matching:

1. **On demand:** no independent union-formation waiting time. A birth to an
   unpartnered woman triggers a search and possible union at the birth date
   (the user's “shotgun wedding” approach).
2. **Scheduled search:** add a union-search event and waiting time analogous
   to birth and death events in the continuous-time engine.

Partner selection uses conditional-logit choice probabilities over eligible men.
Users must be able to specify a flexible matching model; milestone 2 considers
age difference only, while leaving room for later characteristics. The user selected configurable sampling from the available male pool before
pair-specific screening, following the inspected SOCSIM ordering.
A wait/no-partner alternative for a poor pool is under consideration, not yet
selected. Eligibility screening must support kin exclusions, including incest
and cousin relationships; exact prohibited relationships remain to be specified.
SOCSIM is an investigation reference, not an adopted matching algorithm.

Scheduled-search mode must support partnered and unpartnered fertility inputs.
Parity-specific fertility is deferred to a later milestone; milestone 2 does
not add parity conditioning. Union dissolution is excluded except when a
partner dies (widowhood). The user
has not yet decided overlapping-union rules, remarriage eligibility, or search
behavior following widowhood. The construction interface below adds no major dependency.

## Implemented construction interface

`create_matching_model(formula, coefficients)` stores a one-sided R formula
and finite named coefficients matching its non-intercept model-matrix columns.
Only `age_difference` is available initially, defined as male minus female age
in continuous years. Ordinary R formula transformations are supported. An
intercept is ignored because it cancels from partner-only choice probabilities;
no outside/wait alternative is implemented. Formula functions must be
deterministic and must not mutate state. Data-dependent transformations require
care: their basis can depend on the candidate pool. There is no fitted-model
estimation interface yet.

`union_formation(initiation, matching, sample_size = 1000, rates = NULL, min_male_age = 15, kinship_generations = 2)` constructs
a process; `set_union_formation(model, ...)` attaches it. Initiation is explicitly
`"on_demand"` or `"rate"`; initiation and coefficients have no defaults.
The user selected a default sample size of 1000. `sample_size = 0` means use the full available pool; a positive
integer requests up to that many men without replacement BEFORE pair-specific
eligibility screening. Screening can leave no candidates. Rate-mode tables
contain female search-attempt hazards by age (`age_start`, `age_end`, `rate`),
in years and per person-year, with no calendar or sex columns; on-demand mode
rejects rates. This construction layer retains schedules without assigning
unresolved execution eligibility, retry, or gap behavior. Tables are tibbles
and rate provenance is preserved. Edited specifications are checked on attachment.

Worked scoring example: with `~ age_difference` and coefficient `log(2)`,
three eligible candidates whose gaps are 0, 1, and 2 have relative weights
1, 2, and 4 and probabilities 1/7, 2/7, and 4/7. With utilities 1000 and 2000,
subtracting the maximum before exponentiation avoids overflow. An empty pool
has no partner probabilities; its demographic consequences remain undecided.
These are controlled mathematical cases, not suggested demographic defaults.

The internal scoring helper is tested against these analytical cases. The
engine supports the conditional-logit on-demand increment specified below;
Rate-based initiation uses the execution increment below.

## On-demand conditional-logit execution specification

At each birth, use the woman's existing living male partner when she has an
ongoing union. Otherwise draw `min(sample_size, living unpartnered men)` distinct men
uniformly without replacement (`0` means all living unpartnered men), then exclude sampled
men below `min_male_age` or within the configured kinship bound without refilling. For every remaining candidate calculate `age_difference = male age - female age` in continuous
years. The shared search time cancels, so this equals the woman's birth time
minus the man's birth time. Evaluate the user's matching formula and named
coefficients to obtain utility V_j, then choose one candidate with probability

    P(j) = exp(V_j - max(V)) / sum_k exp(V_k - max(V)).

Create the union at the birth time and use the selected man as father. A
nonempty screened pool always produces a partner; there is no wait alternative. Empty
pools preserve the birth with unknown father and no union. A future choice to
cancel the birth is deferred. All-zero utilities give uniform probabilities.
Nonfinite formula values/utilities are errors, not silently excluded candidates.
Formula evaluation cannot directly mutate engine population state. Candidate
sets are not retained.

Worked example requested by the user: formula
`~ age_difference + I(age_difference^2)` with coefficients 0.072 and -0.014
has utility `0.072*d - 0.014*d^2`. Its maximum is at
`d = -0.072 / (2 * -0.014) = 2.5714286` years. Among equally represented
candidate ages, this is the preferred gap; the resulting mean gap depends on
candidate availability. With candidate gaps 0, 3, and 10, probabilities are
proportional to `exp(c(0, 0.09, -0.68))`. The coefficients are test inputs,
not demographic defaults.

The temporary uniform-choice increment is superseded by this formula choice.
Common linear and raw quadratic age-difference terms can be evaluated directly
in C++ with the same formula semantics; other supported R formula transforms
use the existing model-matrix evaluator on the sampled candidate data. There
is no age discretization, cached age-grid approximation, or new dependency.

Minimum male age and known kinship (default two generations on each side)
are checked after sampling. Partnered men are excluded BEFORE sampling: the
available male vector contains living, unpartnered men only. Formation removes
the selected man immediately; widowhood returns a surviving husband immediately.
A woman retains her current partner until death. Input must supply active
female/male unions with living partners and at most one ongoing union per
person; conflicting active unions are rejected. Ended histories are retained.
Death clears both partners' current-partner state; a woman's next birth can
search again. A deceased husband never returns to the available pool.

Worked availability cases: with 1000 partnered men and one unpartnered man,
sample size 1 always samples the unpartnered man before age/kin screening.
If he forms a union, a second woman cannot select him. If his wife dies,
he returns to the pool and may form a later union; if he dies, he stays absent.
If no unpartnered men remain, a birth proceeds with unknown father as before.
Polygamy/polyamory is deferred and requires an explicit future model choice.

Union formation/widowhood do not change
fertility or mortality schedules in this increment, so a surviving partner's
pending event stays valid. Rate-based initiation is specified below.

Validation covers the analytical peak/probabilities, native versus R formula
utilities, sampled choice frequencies, uniform and singleton/empty pools,
extreme utilities and invalid formula output, reproducibility/independence,
partner retention, union accounting, widowhood and restarts. Test thresholds
for stochastic checks are sampling-error bounds, not demographic acceptance
criteria. No scientific validity is claimed for the remaining unrestricted relationships.

## Minimum male age: eligibility increment

The user selected a configurable `min_male_age` in continuous years, default
15. It must be a finite, nonnegative numeric scalar; fractional thresholds
are allowed and zero removes this age exclusion. A sampled man passes when
`search_time - birth_time >= min_male_age`. Check at the search, without
rounding ages. No separate female minimum age is imposed; initiation eligibility
is governed by fertility or, once implemented, union-search rates.

Sample living unpartnered men first at the configured sample size, then remove sampled
men below the minimum, then evaluate matching utilities and normalize
probabilities over survivors. Do not refill the sample or search elsewhere
when sampled candidates fail. Existing ongoing unions are not retroactively
ended by this search restriction. Other exclusions remain pending.

Worked cases: with a threshold of 15, ages 14.999 and 15 are respectively
excluded and included. With threshold 16.5, age 16.25 is excluded and age 16.5
included. A boy initially below the threshold can become eligible later in the
same continuous run. With sample size 1, drawing a boy can leave an empty pool
even if an older man exists outside the sample. At the user's direction, a pool emptied by screening preserves the birth
with unknown father and no new union, just like an empty living-male pool.

Validate construction, attachment, and simulation of edited specifications;
exercise exact/fractional boundaries, aging into eligibility, sampling before
screening, zero threshold, existing partner retention, and no added female-age
restriction. Eligibility changes matching alternatives; coefficients themselves
remain unchanged.

## Initiation, matching, and state consistency

The initiation process determines when a search occurs. The matching process
receives the woman, current time, and eligible candidate information and returns
a supported selection or failure/wait outcome. Matching specifications must not
mutate population state directly. The engine owns union creation, history,
partner availability, parent assignment, and dependent event invalidation.

Public times remain years and hazards per person-year. A scheduled search hazard
is not automatically a realized union-formation hazard: failures and optional waiting
can make realized union formation lower than search initiation. Decide which
quantity the input specifies before defining exposure and validation.

Birth, death, union formation, and widowhood can change eligibility and conditional
fertility. Their timing and invalidation rules need explicit worked examples
before implementation. In on-demand mode, determine whether a union created at
a birth changes that birth's classification and father assignment. There is no
conception process, so no implicit gestational timing should be introduced.

## Matching probabilities and candidate pools

For a woman i, let V_ij be the specified utility of candidate j. Conditional
logit selects among the eligible candidate set C_i with

    P(j | i, C_i) = exp(V_ij) / sum_{k in C_i} exp(V_ik).

Evaluate this with a numerically stable softmax. These are relative partner
choice probabilities, not independent binary union-formation probabilities. This
milestone uses age difference, but its sign and unit are accepted as male age minus female age in continuous
years. Users specify transformations and coefficients through a formula. For example, a signed linear age-gap
term implies a directional preference; a quadratic or absolute-deviation term
can represent a preferred gap. No demographic coefficients are defaulted.

A common additive utility constant cancels from partner-only probabilities.
Without an outside alternative, the model chooses someone whenever the pool is
nonempty, even if all candidates have low utilities. A wait alternative would
require an explicit utility scale relative to partners and defined retry timing.
Whether it belongs in on-demand mode, scheduled mode, or both remains open.

Candidate subsampling changes the choice set and generally changes the induced
matching distribution. With a wait alternative it can also change the chance of
forming any union. Specify sampling, any correction, and the intended target
before implementing it as an efficiency option. Begin by investigating full-pool
costs; this is not a decision to adopt a particular sampling scheme.

## Rationale and matching-model validation

The user prefers probabilistic conditional-logit matching because it gives
less-preferred but eligible candidates nonzero probabilities and accommodates
future characteristics through model terms. For example, rare group exogamy
could later be represented by a finite utility coefficient, rather than an
extra acceptance mechanism around a maximum-score rule. Group effects are
future extensions, not milestone 2 matching inputs. A hard eligibility ban
would still prevent exogamy regardless of utility; distinguish exclusions from
preferences. Coefficients govern relative choice probabilities conditional on
the available pool, not fixed realized exogamy proportions.

The second rationale is testability: fit the specified choice model back to
simulated matching decisions and assess approximate recovery of its parameters.
At the user's direction, do not retain candidate sets. Instead, reconstruct
approximate alternatives from eligible men in time slices using recorded
population and union histories. Time-slice width, reference instant, eligibility
reconstruction, and any sampling scheme remain to be specified for analysis.

This is approximate model recovery rather than the exact likelihood of the
simulation's choices. Within-slice births, deaths, aging, and union formation or
widowhood can change availability. If simulation searches use sampled pools,
the reconstructed eligible population will also differ from the actual sampled
alternatives. An optional wait alternative would require an additional analysis
design because formed unions alone omit decisions to wait or failed searches.
These limits do not authorize candidate-set retention or imply extra search logs
are accepted. Assess sensitivity to time-slice width when developing validation.

Combine this approximate recovery analysis with controlled known-utility choice
experiments that check probabilities without retaining candidate sets in
simulation outputs. Model form, sample sizes, and scientific acceptance criteria
remain to be agreed. Neither demographic realism nor parameter recovery has
yet been demonstrated.

## Eligibility and ancestry

The [shared kinship module](kinship.md) now supports configurable generation
bounds and `extract_kin()` using the same recorded-parent traversal. The
default two-generation bound screens demonstrated relationships only, including
half kin and paternal links when a mother is unknown. Ancestry cycles are
rejected by simulation and extraction; missing parent links remain unknown.

Screen sampled candidates before scoring or selection without refilling. Male
minimum age is implemented and there is no extra female age restriction. Exclusive ongoing unions are required before sampling. Expanded
social or categorical kin rules remain unresolved. The generation bound also limits direct
ancestry: two does not automatically exclude a great-grandparent. Unknown
ancestry must not be represented as verified absence of kinship. More extensive
parent-age/sex validation remains separate from cycle-safe traversal.

The [SOCSIM audit plan](socsim-audit.md) identifies the pinned reference and
provenance requirements. The [static union-formation trace](socsim-union-formation-audit.md) records
queue behavior, matching, and kin-screening caveats in the pinned source.
These have not been verified with executable cases; do not claim complete parity.

## Fertility and widowhood

On-demand mode retains birth-triggered initiation. When no sampled candidate
passes screening, the accepted behavior preserves the birth with unknown father
and no new union. Requiring a match before a
birth would change the supplied fertility target through partner availability.

Partnered and unpartnered fertility rates use woman-years in the current union
status; the implemented schedule and transition rules are specified under D039 below. D028's population-ASFR midpoint compensation cannot simply be
assumed correct across status transitions; specify the target and spacing
compensation before implementing conditional fertility.

A partner's death ends the ongoing union. Specify the survivor's next search
eligibility and timing, fertility changes, and subsequent father assignment.
No non-death dissolution process is included in milestone 2.

## Open decisions before implementation

| Decision | Consequence |
| --- | --- |
| Failed on-demand search: birth proceeds or depends on a match | Preserves or changes fertility targets |
| One ongoing union per person; remarriage after widowhood | Defines availability and subsequent searches |
| Initiation-rate denominator and failed-search retries | Separates attempted from realized formation |
| Matching coefficients and transformations | User-supplied formula; male-minus-female continuous-year gap |
| Candidate sample size | Defaults to 1000; zero means full available pool, sampled before screening |
| Wait alternative and its scope | Allows refusal; needs utility and retry semantics |
| Further eligibility | Exclusive ongoing unions, minimum male age and two-sided known-kin distance are implemented; further policies remain open |
| Conditional fertility inputs and spacing after status changes | Defines rates, exposure, and event resampling |
| Birth-time union classification and father assignment | Defines observed histories and parentage |
| Equal-time priorities and effects on both partners | Preserves state consistency |

## Specification and validation work

Once dependent decisions are settled, write the executable scientific contract,
worked examples, and observation definitions before engine changes. Proposed
validation includes:

- Uniform matching at equal utilities and known probability ratios at unequal
  utilities; stable behavior for large utilities; fixed-seed reproducibility.
- Empty/single-candidate pools, exact age boundaries, prohibited kin pairs,
  incomplete ancestry, availability, and any selected sampling/wait rules.
- On-demand births with successful and failed searches, including classification
  and father assignment under the agreed rules.
- Analytical search waiting times where every search succeeds; attempted versus
  realized rates where searches can fail.
- Union history accounting, widowhood, and event invalidation for both partners;
  no non-death dissolution or prohibited overlap under the selected rules.
- Status-specific fertility exposure, transitions, birth spacing, and accounting;
  unchanged milestone 1 behavior when union formation is not enabled.

These are proposed validation cases, not completed findings or numerical
scientific acceptance thresholds. A narrow SOCSIM source audit can inform the
remaining choices without adopting legacy architecture or copying its code.

## Rate-mode execution increment (October 15, 2026)

The user requested female-only age-specific search events competing with death.
Only living, unpartnered women face the supplied search hazard (per woman-year);
never-partnered and widowed women use the same schedule. Age-band gaps and
outside support have zero search hazard, using the shared `waiting_time()`.
Men remain available for selection but never initiate a rate search.
Search execution uses the same sampling, screening, utility and union-creation
code as on-demand execution, without creating a child at union formation.

In rate mode, both statuses can have birth events according to the fertility
schedule selected under D039. Formation and widowhood switch schedules; a
husband's death also restores search risk immediately. Birth spacing is preserved.
Death wins exact same-person event ties. Cross-person ties retain queue order.
Pending events invalidated by a status change must never execute. Surviving
women are rescheduled from the transition time using their attained age.
Midpoint fertility adjustment operates independently within each schedule,
without a claim of calibrated conditional ASFRs across status transitions.

Worked cases: zero search hazards and no initial union yield no unions; births then follow
the unpartnered schedule. With explicit zero unpartnered rates, there are no births. With zero mortality and a constant search hazard
h and sufficient men, the fraction forming a first union by t is 1-exp(-h*t).
A woman whose husband dies before her pending birth cannot execute that birth;
birth risk now follows the unpartnered schedule. A search
outside fertility ages can still form a union, but cannot produce births.

After a failed search, draw another waiting time from the same age-specific
schedule starting at the attempt time, as selected by the user. Waiting/refusal for poor nonempty pools is
still deferred until after this rate increment.

## Fertility schedules by union status (D039)

An optional character `union_status` column in the fertility table must contain
both `partnered` and `unpartnered`, with no missing or other values. Without
that column the table is one overall schedule. Mixing overall and conditional
rows or providing only one status is an error; explicitly provide zero rates
for a status intended to have no births. Validate finite nonoverlapping age
bands separately within each status; their age boundaries may differ. Gaps
and ages outside each schedule have zero fertility. Rates are births per
woman-year in the indicated status; an overall schedule applies to all women.

On-demand initiation requires an overall schedule and rejects conditional
schedules at simulation validation, including zero-duration runs. Rate initiation
accepts either form. With overall rates it emits one warning per simulate call
that the same schedule applies to both statuses. This replaces D038's temporary
suppression of all unpartnered births in rate mode.

In rate mode an unpartnered woman faces birth, search and death simultaneously;
a birth does not initiate a union but selects a father under D040 below. Partnered
births use the current husband. Formation and widowhood switch the applicable
fertility schedule and reschedule from the transition time. Birth spacing is
shared across statuses and survives transitions/restarts. Failed searches retain
unpartnered fertility. Death wins exact ties; birth precedes search on an exact
same-person birth/search tie. Cross-person ties retain queue insertion order.

Apply existing midpoint compensation independently within each schedule; do not
integrate across statuses or claim calibrated status-specific ASFR recovery
around transitions. History-aware compensation remains deferred.

Worked cases: overall rate 0.1 applies to both statuses with a rate-mode warning.
Unpartnered rate zero and partnered rate 0.2 reproduces the former gate explicitly.
If union formation occurs within a birth-spacing interval, no birth can occur
before that interval ends, regardless of the new schedule. A widow switches to
unpartnered fertility, so a previously scheduled partnered birth is invalidated.
Validate these cases alongside different age bands, rejection of invalid status
labels/incomplete tables, independent compensation and RNG isolation.

## Nonmarital father selection (D040)

For a birth to an unpartnered woman in rate mode, select a father from the same
living, currently unpartnered male vector used for union searches. Use the union
process's sample size, sample-before-screening rule, minimum male age and known-kin
exclusion. Use `father_matching` when provided, otherwise the general
conditional-logit matching formula/coefficients. Selection is at
birth time, before the newborn enters the population. Empty or entirely rejected
samples preserve the birth with an unknown father and no refill.

Assign parentage without forming a union, updating current partners, removing
the man from availability, or changing his event schedule. Record paternal links
in the shared kinship graph immediately. Subsequent nonmarital births/searches
use the same ordinary selection procedure, with no preferential reuse of an
earlier child's father. Such preferential reuse is an unresolved extension;
there is no new persistence parameter. Father selection and partner selection can use separate matching models
under D041 below.

Worked cases: with one eligible unpartnered man, multiple nonmarital births can
all name him while union history stays empty. If that man forms a union with
another woman, he is unavailable for further nonmarital father selection. With
no eligible men, births continue with unknown fathers. Selection can repeat the
same man by chance without implying a continuing partnership. Paternal half
siblings share the recorded father in kinship extraction/screening.

## Separate nonmarital father matching (D041)

`union_formation()` and `set_union_formation()` accept `father_matching = NULL`.
A non-NULL value must be a `create_matching_model()` specification with the same
formula contract as general `matching`. Only rate-mode births to unpartnered
women use this model. NULL resolves to the current general model at execution,
not a copy frozen at process construction. Existing objects without the field
also resolve to the general model.

Union searches, including on-demand searches, always use `matching`; partnered
births still use the current partner. On-demand mode does not use the optional
father model. Validate any supplied model at construction, attachment and
simulation, including zero-duration runs. Candidate pool, sampling, age/kin
eligibility, parentage-only state changes and unknown-father fallback are shared.
Both matching models support native linear/quadratic evaluation and R formula
fallback. No new previous-father preference is introduced.

Worked case: with equally available men at age gaps 0 and 10, general utility
`100*d` selects the older man for unions, while father utility `-100*d` selects
the younger man for nonmarital fatherhood. Setting father_matching to NULL
restores the general utility and must reproduce an explicitly identical model
under the same seed. Separate coefficients affect selection, not search hazards
or fertility schedules.
