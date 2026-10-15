# Milestone 2: union formation

Updated October 15, 2026 from the user's milestone 2 request.
**Status: milestone scope accepted (D030); detailed demographic rules below
remain proposals or open decisions. Process construction and formula scoring
and conditional-logit on-demand execution with minimum male age screening
are implemented; other pair restrictions and rate-based initiation remain pending.**
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

`union_formation(initiation, matching, sample_size = 1000, rates = NULL, min_male_age = 15)` constructs
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
`simulate()` rejects rate-based initiation.

## On-demand conditional-logit execution specification

At each birth, use the woman's existing living male partner when she has an
ongoing union. Otherwise draw `min(sample_size, living men)` distinct men
uniformly without replacement (`0` means all living men), then exclude sampled
men below `min_male_age` without refilling. For every remaining candidate calculate `age_difference = male age - female age` in continuous
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

Minimum male age is checked after sampling, as specified below. Other pair
restrictions remain deferred: partnered men and relatives are not excluded. Men may have concurrent unions; this is not an
accepted final eligibility policy. A woman retains her current partner until
death. Input must supply unambiguous active female/male unions with living
partners and at most one ongoing union per woman. Ended histories are retained.
A partner's death ends affected unions and clears current-partner state; a
woman's next birth can search again. Union formation/widowhood do not change
fertility or mortality schedules in this increment, so a surviving partner's
pending event stays valid. Scheduled initiation remains rejected.

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

Sample living men first at the configured sample size, then remove sampled
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

Screen candidates before scoring or selection. Specify minimum ages, availability,
existing-union exclusions, remarriage, and exact kin prohibitions. Distinguish
parent/child, siblings and half-siblings, other ancestor/descendant relationships,
aunt/uncle relationships, and cousin degrees rather than relying on “incest” as
a complete rule. Decide how unknown parents and incomplete baseline ancestry
limit screening. Unknown ancestry must not be represented as verified absence
of kinship. Structural validation of input histories and ancestry cycles also
needs a defined policy.

The [SOCSIM audit plan](socsim-audit.md) identifies the pinned reference and
provenance requirements. The [static union-formation trace](socsim-union-formation-audit.md) records
queue behavior, matching, and kin-screening caveats in the pinned source.
These have not been verified with executable cases; do not claim complete parity.

## Fertility and widowhood

On-demand mode is intended to retain birth-triggered initiation; the outcome
when no partner is found still needs a decision. Requiring a match before a
birth would change the supplied fertility target through partner availability.

For scheduled mode, partnered and unpartnered fertility rates require exposure denominators
based on union status. Decide whether these are separate status-specific ASFRs,
how gaps/zero rates work, and what happens to an already scheduled birth when
status changes. D028's population-ASFR midpoint compensation cannot simply be
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
| Exact age/kin exclusions and unknown ancestry policy | Defines candidate eligibility |
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
