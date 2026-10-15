# Milestone 2: union formation

Updated October 15, 2026 from the user's milestone 2 request.
**Status: milestone scope accepted (D030); detailed demographic rules below
remain proposals or open decisions. No union-formation process is implemented.**
The current engine still assigns fathers without forming new unions.

## Terminology

Use **union**, **union formation**, and **partnered/unpartnered** in code and
documentation. The user sometimes says “marriage” conversationally, but the
model is deliberately agnostic about legal marriage, long-term cohabitation,
and other union types. Milestone 2 does not introduce a legal-marital-status
classification. Conditional fertility refers to whether a woman is in an
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
age difference only, while leaving room for later characteristics. Whether to
use every eligible man or a smaller candidate sample remains undecided.
A wait/no-partner alternative for a poor pool is under consideration, not yet
selected. Eligibility screening must support kin exclusions, including incest
and cousin relationships; exact prohibited relationships remain to be specified.
SOCSIM is an investigation reference, not an adopted matching algorithm.

Scheduled-search mode must support partnered and unpartnered fertility inputs.
Parity-specific fertility is deferred to a later milestone; milestone 2 does
not add parity conditioning. Union dissolution is excluded except when a
partner dies (widowhood). The user
has not yet decided overlapping-union rules, remarriage eligibility, or search
behavior following widowhood. No new union API or package dependency is selected.

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
milestone uses age difference, but its definition, functional form, and
coefficients still need user decisions. For example, a signed linear age-gap
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
| Age gap sign, transformations, coefficients | Defines the matching model without demographic defaults |
| Full pool versus subsampling | Changes cost and potentially selection probabilities |
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
