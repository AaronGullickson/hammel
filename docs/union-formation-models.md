# Union-formation planning notes

Updated October 13, 2026 from the user's notes and subsequent discussion.
**Status: proposals, not an accepted demographic algorithm or implementation
instruction.** Union formation is deferred from the initial milestone scope.
The current simulator assigns fathers without forming unions.

## Fertility and union formation interact

If population ASFRs are applied only to partnered women, overall births can
fall below the input target. Unpartnered women still belong in the population
ASFR denominator; marital fertility instead uses exposure within unions.
The interpretation of the supplied rate, not merely the event rate, matters.

The user's “driver's seat” distinction asks which outcomes should be prescribed:
aggregate fertility or union formation and conditional fertility. These inputs
must be compatible. Preserving total fertility may require calibration or
simplifying unions; realistic union dynamics may make total fertility an outcome.
Neither priority is selected here.

## Separate initiation from partner choice

The user proposes distinguishing the decision to seek a union from the choice
of a partner. Both formation approaches below could share eligibility and
matching machinery. Partner supply and attractiveness may also affect the
initiation decision—a marriage squeeze need not be ignored permanently.

### Fertility-triggered formation

An unpartnered woman scheduled to give birth triggers a partner search and
possible union at that event. No separate nuptiality schedule is required.
This deliberately links union timing to births rather than realistic independent
marital histories. Widowhood or dissolution could permit subsequent searches.

Preserving fertility depends on what happens when no acceptable partner exists:
does the birth proceed, proceed without a union, or wait? No fallback is decided.
A union created at the birth date also requires an explicit marital/non-marital
birth classification. There is no conception process in the current engine.

### Formation as a competing event

An age-specific initiation hazard schedules a search independently of birth
and death. The user prefers considering female-initiated searches among
available men over a two-sided waiting queue.

A search-attempt hazard differs from an observed formation rate if searches can
fail. Restricting population fertility to partnered women can lower total ASFR;
using conditional marital/non-marital fertility requires corresponding exposure
and validation. Formation-rate units, eligibility, remarriage, retries, and
failed-search behavior still need decisions.

The original notes describe SOCSIM's one- and two-queue approaches and express
a preference against the legacy two-queue approach. Those specific matching
claims are investigation leads, not findings established by the current narrow
[SOCSIM audits](socsim-audit.md). Verify source/version before relying on them.

## Partner eligibility and probabilistic selection

The user proposes searching all available men or a sample, first excluding
ineligible pairs such as prohibited kin relationships, then scoring candidates.
A flexible user-specified score would support stochastic selection instead of
always choosing the highest-scoring partner. Exact exclusions, age rules,
availability, and candidate sampling remain unspecified.

A conditional-logit choice among candidate utilities s_j would use

    P(j) = exp(s_j) / sum_k exp(s_k).

These are relative choice probabilities, not separately calculated binary
marriage probabilities. A “wait/no union” alternative could make the chance of
formation depend on both partner supply and attractiveness. Threshold-based
waiting was also suggested in the original notes; neither mechanism is chosen.
Sampling a subset of candidates can change choice and waiting probabilities,
so its interpretation would need a separate specification.

## Decisions needed before implementation

| Topic | Consequence |
| --- | --- |
| Population versus conditional fertility | Determines targets and exposure denominators |
| Initiation versus realized formation rates | Failed searches change observed formation |
| Failed-search fallback and retry timing | Can alter births and union histories |
| Pair eligibility and overlap rules | Controls kinship, ages, and partnership availability |
| Choice score, waiting alternative, candidate sampling | Determines matching and supply effects |
| Same-time union/birth classification | Determines marital/non-marital birth observations |
| Dissolution, widowhood, and remarriage | Changes eligibility and future searches |

Forming/ending a union affects both partners, so the engine must consistently
invalidate or reschedule dependent candidates. Future user callbacks should
return supported decisions/scores rather than directly mutate relationship state.
These are software contract requirements; they do not settle the demographic
choices above. Write worked examples and validation expectations before coding.
