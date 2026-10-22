# Milestone 2: union formation

Updated October 19, 2026. Implemented milestone 2 specification; the user
accepted the milestone as complete on this date (D050).
**Status: core milestone implementation is in place: both search modes,
conditional-logit matching, age/kin/exclusivity restrictions, status-specific
fertility, widowhood, and separate optional nonmarital father matching. The
failed-match birth cancellation is implemented (default TRUE); transition-aware fertility
compensation (D045) and poor-pool waiting (D044) remain deferred.
The user accepted milestone 2 as complete on October 19, 2026 (D050).**
Every simulation requires a union-formation process. Current-partner lookup
and union histories are common infrastructure for all initiation modes.

An optional sequential two-sided selection experiment was added October 22,
2026 (D061), with `matching_method = "two_sided"`. The specification below
continues to describe default female choice; see
[two-sided matching](two-sided-matching.md) for the experiment's temporary
markets and unchanged event-initiation rules. Its demographic benchmarks
remain under investigation.

## Terminology

Use **union**, **union formation**, and **partnered/unpartnered** in code and
documentation. The user sometimes says “marriage” conversationally, but the
model is deliberately agnostic about legal marriage, long-term cohabitation,
and other union types. Milestone 2 does not introduce a legal-marital-status
classification. At the user's request, union tables nevertheless use explicit
`wife_id` (female) and `husband_id` (male) roles so sex is clear. Conditional fertility refers to whether a woman is in an
ongoing union; input rates use woman-years in the indicated status.

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
A wait/no-partner alternative for a poor pool is deferred to a later milestone
(D044); its utility, scope and retry details have not been selected. Eligibility uses the
accepted age, generation-distance kin and exclusivity restrictions below.
SOCSIM is an investigation reference, not an adopted matching algorithm.

Scheduled-search mode must support partnered and unpartnered fertility inputs.
Parity-specific fertility is deferred to a later milestone; milestone 2 does
not add parity conditioning. Union dissolution is excluded except when a
partner dies (widowhood). Exclusive ongoing unions, return of widowed men to
availability and renewed female rate searches after widowhood are implemented. The construction interface below adds no major dependency.

## Implemented construction interface

`create_matching_model(formula, coefficients, group_exogamy = NULL,
group_levels = NULL)` stores a one-sided R formula and finite named coefficients
matching its non-intercept model-matrix columns. Predictors include
`age_difference` (male minus female age in continuous years), `female_group`,
`male_group`, and numeric `same_group`. Categorical group terms require explicit
fixed `group_levels`. Optional directed/symmetric pair effects add to the formula
score. See [group matching](group-matching.md) for syntax, factor contrasts, and
worked probabilities. Ordinary R formula transformations are supported. An
intercept is ignored because it cancels from partner-only choice probabilities;
no outside/wait alternative is implemented. Formula functions must be
deterministic and must not mutate state. Data-dependent transformations require
care: their basis can depend on the candidate pool. There is no fitted-model
estimation interface yet.

```r
union_formation(
  initiation, matching = create_age_gap_matching_model(0, 5),
  sample_size = 1000, rates = NULL, min_male_age = 15,
  kinship_generations = 2, father_matching = NULL,
  cancel_birth_if_no_match = TRUE, diagnostics = FALSE
)
```

This constructs
a process; `set_union_formation(model, ...)` attaches it. Initiation is explicitly
`"on_demand"` or `"rate"`; initiation has no default. Omitting `matching` uses
`create_age_gap_matching_model(0, 5)`: a preferred male-minus-female age gap of
zero years and a symmetric half-weight distance of five years. Explicit matching
models override this default in both the constructor and setter. With
`father_matching = NULL`, nonmarital father selection uses this same model.
The user selected a default sample size of 1000. `sample_size = 0` means use the full available pool; a positive
integer requests up to that many men without replacement BEFORE pair-specific
eligibility screening. Screening can leave no candidates. Rate-mode tables
contain female search-attempt hazards by age (`age_start`, `age_end`, `rate`),
in years and per person-year, with no calendar or sex columns; on-demand mode
rejects rates. Search hazards apply to unpartnered women; gaps have zero hazard and failed
searches draw a new waiting time. Tables are tibbles
and rate provenance is preserved. Edited specifications are checked on attachment.

Worked scoring example: with `~ age_difference` and coefficient `log(2)`,
three eligible candidates whose gaps are 0, 1, and 2 have relative weights
1, 2, and 4 and probabilities 1/7, 2/7, and 4/7. With utilities 1000 and 2000,
subtracting the maximum before exponentiation avoids overflow. An empty pool
has no partner probabilities: follow the configured birth-cancellation policy
or reschedule a failed rate search, as appropriate.
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
pools cancel the birth by default, or preserve it with unknown father when
`cancel_birth_if_no_match = FALSE` (D043). All-zero utilities give uniform probabilities.
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
If no unpartnered men remain, cancellation follows D043; explicitly disabling
cancellation permits a birth with unknown father.
Polygamy/polyamory is deferred and requires an explicit future model choice.

In on-demand mode, union formation/widowhood do not change fertility or
mortality schedules, so a surviving partner's pending event stays valid.
Rate mode invalidates and resamples women's events as specified below.

Validation covers the analytical peak/probabilities, native versus R formula
utilities, sampled choice frequencies, uniform and singleton/empty pools,
extreme utilities and invalid formula output, reproducibility/independence,
partner retention, union accounting, widowhood and restarts. Test thresholds
for stochastic checks are sampling-error bounds, not demographic acceptance
criteria; broader demographic assessment is recorded separately in the report.

## Minimum male age: eligibility increment

The user selected a configurable `min_male_age` in continuous years, default
15. It must be a finite, nonnegative numeric scalar; fractional thresholds
are allowed and zero removes this age exclusion. A sampled man passes when
`search_time - birth_time >= min_male_age`. Check at the search, without
rounding ages. No separate female minimum age is imposed; initiation eligibility
is governed by fertility or union-search rates, according to initiation mode.

Sample living unpartnered men first at the configured sample size, then remove sampled
men below the minimum, then evaluate matching utilities and normalize
probabilities over survivors. Do not refill the sample or search elsewhere
when sampled candidates fail. Existing ongoing unions are not retroactively
ended by this search restriction. Known-kin screening follows the age check.

Worked cases: with a threshold of 15, ages 14.999 and 15 are respectively
excluded and included. With threshold 16.5, age 16.25 is excluded and age 16.5
included. A boy initially below the threshold can become eligible later in the
same continuous run. With sample size 1, drawing a boy can leave an empty pool
even if an older man exists outside the sample. A pool emptied by screening follows the configured birth-cancellation policy
(D043), just like an empty available-male pool.

Validate construction, attachment, and simulation of edited specifications;
exercise exact/fractional boundaries, aging into eligibility, sampling before
screening, zero threshold, existing partner retention, and no added female-age
restriction. Eligibility changes matching alternatives; coefficients themselves
remain unchanged.

## Group eligibility restrictions

`prevent_exogamy` belongs to each matching model and defaults to FALSE. TRUE
requires same-group candidates; a character vector such as `c("A::B", "A:C")`
prohibits selected symmetric/directed cross-group combinations. Sample before
these restrictions, and do not refill exclusions. The selected father/partner
model supplies the restrictions, finite score effects cannot override them,
and an empty sample follows existing failure/cancellation rules. Existing
partners are retained. See [group matching](group-matching.md).

## Initiation, matching, and state consistency

The initiation process determines when a search occurs. Native matching scores
the screened age differences and selects a partner, or reports failure for an
empty set. Formula callbacks receive age differences and, when requested,
fixed group-identity predictors; they must not mutate state. Add pair offsets
only after the formula score is computed. The engine owns union creation, history, partner availability, parent
assignment, and dependent event invalidation. No wait outcome is implemented.

Public times remain years and hazards per person-year. Scheduled input hazards
specify search attempts per unpartnered woman-year, not guaranteed formations:
failed searches can reduce realized formation rates. In on-demand mode a
successful search forms the union at birth, assigns its husband as father, and
classifies that birth as partnered in union-status summaries (D046). There is
no conception process or implicit gestational timing. Rate-mode transitions
and tie priorities are specified below and in the scientific contract.

## Matching probabilities and candidate pools

For a woman i, let V_ij be the specified utility of candidate j. Conditional
logit selects among the eligible candidate set C_i with

    P(j | i, C_i) = exp(V_ij) / sum_{k in C_i} exp(V_ik).

Evaluate this with a numerically stable softmax. These are relative partner
choice probabilities, not independent binary union-formation probabilities. This
milestone uses age difference, but its sign and unit are accepted as male age minus female age in continuous
years. Users specify transformations and coefficients through a formula. For example, a signed linear age-gap
term implies a directional preference; a quadratic or absolute-deviation term
can represent a preferred gap. The accepted default uses a zero-year preferred
gap and five-year half-weight distance; explicit specifications override it.

A common additive utility constant cancels from partner-only probabilities.
Without an outside alternative, the model chooses someone whenever the pool is
nonempty, even if all candidates have low utilities. A wait alternative would
require an explicit utility scale relative to partners and defined retry timing.
It is deferred beyond milestone 2 (D044); its scope and semantics remain open.

Candidate subsampling changes the choice set and can change the induced
matching distribution. The accepted implementation samples up to 1000 living,
unpartnered men by default before screening, without refilling; zero requests
the full pool. No sampling correction or outside alternative is implemented.

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
At the user's direction, candidate sets are not retained. The
[user-authored milestone 2 report](../validation/milestone2.qmd) now reconstructs
alternatives at each sampled union's exact formation time: 5,000 on-demand
unions, each combined with 50 alternative living, unpartnered men aged at least
15. It fits linear and quadratic age-gap coefficients and reports recovery of
the generating values 0.072 and -0.014. This replaces the earlier proposed
coarse time-slice analysis; it does not introduce a package estimation API.

The reconstructed alternatives are not the actual simulation samples. The
report's helper checks age, survival and partnership status but does not apply
pair-specific kin screening or reconstruct ordering of coincident events.
Keep these limits distinct from controlled known-utility software checks.
The report's findings are the user's analysis, not an independent replication.
The user accepted milestone 2 on October 19, 2026 (D050), retaining these
limitations. No new scientific thresholds are imposed.

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
passes screening, the birth is cancelled by default (D043). Explicitly disabling
cancellation preserves the birth with unknown father and no new union. Requiring
a match allows partner availability to reduce realized fertility.

Partnered and unpartnered fertility rates use woman-years in the current union
status; the implemented schedule and transition rules are specified under D039
below. Midpoint compensation is computed independently per status. It does not
account for prior-status rate history; correction is deferred under D045.

A partner's death ends the ongoing union. In rate mode, a surviving woman
resumes search risk and switches fertility status with pending events
invalidated; in on-demand mode, her next birth attempt can trigger a new search.
Surviving men return to the available pool. No non-death dissolution process
is included in milestone 2.

## Deferred decisions and validation status

Poor-pool waiting, previous-father continuity, transition-aware fertility
compensation, parity conditioning and additional matching characteristics are
future work. Their detailed rules remain open; see the [roadmap](roadmap.md).
Current initiation, eligibility, retries, birth cancellation, status fertility,
event priorities and birth-time classification are implemented, not open choices.

Software checks cover analytical choice probabilities, numerical stability,
empty/singleton samples, age/kin/exclusivity restrictions, search hazards,
widowhood and invalidation, birth spacing, restart histories, diagnostics and
reproducibility. [Testing](testing.md) maps these to the current suite.
Demographic assessment, parameter recovery and timing comparisons appear in
`validation/milestone2.qmd`. The user accepted milestone 2 as complete on
October 19, 2026 (D050), based on the report and review, with deferred
limitations retained. Passing software tests alone does not establish
demographic validity.

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
deferred to a later milestone (D044).

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
around transitions. Transition-aware compensation is deferred under D045;
the report examines case-specific discrepancies, but no general bias bound or
corrective algorithm has been established.

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
samples cancel the birth by default, or preserve it with unknown father when
cancellation is explicitly disabled. Do not refill.

Assign parentage without forming a union, updating current partners, removing
the man from availability, or changing his event schedule. Record paternal links
in the shared kinship graph immediately. Subsequent nonmarital births/searches
use the same ordinary selection procedure, with no preferential reuse of an
earlier child's father. Such preferential reuse is deferred to a later milestone (D042);
there is no new persistence parameter. Father selection and partner selection can use separate matching models
under D041 below.

Worked cases: with one eligible unpartnered man, multiple nonmarital births can
all name him while union history stays empty. If that man forms a union with
another woman, he is unavailable for further nonmarital father selection. With
no eligible men, births are cancelled unless unknown fathers are explicitly allowed. Selection can repeat the
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
eligibility, parentage-only state changes and failed-match cancellation settings
are shared.
Both matching models support native linear/quadratic evaluation and R formula
fallback. Previous-father preference is deferred to a later milestone (D042).

Worked case: with equally available men at age gaps 0 and 10, general utility
`100*d` selects the older man for unions, while father utility `-100*d` selects
the younger man for nonmarital fatherhood. Setting father_matching to NULL
restores the general utility and must reproduce an explicitly identical model
under the same seed. Separate coefficients affect selection, not search hazards
or fertility schedules.

## Failed-match birth cancellation (D043)

`union_formation()` and `set_union_formation()` accept the scalar logical
`cancel_birth_if_no_match`, default TRUE. Reject missing, nonlogical or nonscalar
values at construction, attachment and simulation validation. In on-demand mode
this applies when the birth-triggered union search fails; in rate mode it applies
when an unpartnered birth's father selection fails. An empty available pool or
sample emptied by age/kin screening both count as failure. Do not refill samples.

TRUE prevents creating the child when no father is selected. Do not allocate a
person ID, append kinship links, draw newborn sex, increment births, or create a
union for that failed attempt. FALSE explicitly retains the existing unknown-
father birth behavior. Partnered births and successful matches are unaffected.
Input persons with missing historical father IDs remain unchanged; this controls
new births only. Search events without births retain their existing retry policy.
Cancellation allows partner scarcity to depress realized fertility below inputs;
fertility compensation must not silently restore those suppressed births.

After cancellation redraw competing event waiting times from the failed attempt,
using attained age and current status. Do not change last_birth or start a new
birth-spacing interval; any spacing from an actual prior birth remains valid.
A man becoming eligible after several cancellations can therefore be selected
at the next sampled birth attempt, without an artificial spacing delay.

## Optional union diagnostics

`diagnostics = TRUE` on the union process records `.n_eligible` and
`.utility_rank` for each newly formed union in either initiation mode. The count
is the sampled candidate set remaining after all eligibility screening, not the
whole male pool. Rank is one plus the number of eligible candidates with strictly
greater deterministic utility; 1 is best and ties share the minimum rank. For
utilities 3, 3, 1, choosing either first candidate records rank 1 of 3; choosing
the last records rank 3 of 3. A singleton records count 1 and rank 1.

Diagnostics default off, consume no random draws, and do not change matching or
event timing. Failed searches and nonmarital father selections create no union
row. Existing recorded values survive restarts; missing historical diagnostics
are NA, including when enabled for a zero-duration run. With diagnostics off,
existing columns survive and newly formed unions have NA in those columns.
Rank requires one linear pass over already computed utilities, without sorting
or evaluating the formula again; that pass is skipped when disabled.
## Age-gap matching helper

`create_age_gap_matching_model(preferred_age_gap, half_weight_distance)` converts
two explicit inputs in years into the existing formula-based matching model.
Age gap is male age minus female age; the preferred gap may be negative or zero.
The half-weight distance must be strictly positive. There are no defaults.
With preferred gap `d` and distance `h`, utility is equivalent, up to an ignored
constant, to `-log(2) * ((age_difference - d) / h)^2`. The quadratic coefficient
is `-log(2) / h^2`; the linear coefficient is `-2 * d * quadratic`.
For example, `d = 3` and `h = 5` give equal half-peak weights at gaps -2 and 8.
Decline is symmetric. These are relative candidate selection weights, not
guaranteed realized age-gap frequencies or a halving of utility itself.
The helper returns an ordinary matching model and retains the native quadratic
execution path. Inputs must be finite numeric scalars and yield representable
finite coefficients with strictly negative quadratic curvature.
