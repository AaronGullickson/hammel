# Newborn group identity

Requested and specified by the user October 20, 2026. This assigns one identity
at birth, separately from deterministic ancestry averaging. Multiple simultaneous
identities remain deferred.

Milestone 3 was accepted by the user October 27, 2026 (D077). The
[validation report](../validation/milestone3.qmd) records the assessment;
acceptance retains the documented limitations.

## Interface and rules

`create_group_inheritance(method = "random_parent", fun = NULL, groups = NULL,
order = NULL, passing = 0)` constructs an inspectable inheritance object.
`fertility(..., inheritance = create_group_inheritance())` and its setter attach
it. Unspecified inheritance, including older specifications without this field,
uses random-parent assignment with maternal fallback when the father is unknown.
Configuration is validated
at construction, attachment, and execution.

For rules requiring no additional parameters, `inheritance` also accepts a
method name: `"mother"`, `"father"`, `"random_parent"`, `"same_sex"`,
`"opposite_sex"`, or `"ancestry"`. For example,
`set_fertility(model, rates, inheritance = "father")` is equivalent to using
`create_group_inheritance("father")`. Stored process specifications contain the
full model object. Hierarchy and custom rules still require the constructor
to supply `order` or `fun`; their bare names produce the corresponding error.

| Method | Identity assignment |
| --- | --- |
| `mother` | Mother's group |
| `father` | Father's group |
| `random_parent` | Select a parent with equal probability |
| `same_sex` | Mother for daughters, father for sons |
| `opposite_sex` | Father for daughters, mother for sons |
| `ancestry` | Sample using the child's known ancestry shares |
| `hierarchy` | Choose the higher-ranked parental identity in `order` |
| `custom` | Call `fun` and validate its returned identity |

Missing information required by a built-in rule falls back to the mother.
Missing fathers therefore never create an unknown identity. Ancestry sampling
conditions on known ancestry: divide known shares by `1 - unknown`. Entirely
unknown ancestry falls back to the mother. A hierarchy compares identities,
not ancestry. Missing parental ranks fall back to the mother. Its strict order
lists highest first; ties, cycles, and pair-specific exceptions use a custom
function. `passing` is a probability in `[0,1]`, default 0: independently reverse
the order for that birth. For A above B and A/B parents, passing 0.05 yields
A with probability 0.95 and B with probability 0.05. Reversal has no effect
when identities are equal or a required parent/rank is missing.

Deterministic rules, identical-parent choices, and fallback cases consume no
identity draw. Random-parent choice draws only for different known parental
identities; ancestry draws only when more than one known component is positive;
passing at 0 or 1 is deterministic. Grouped models use the new rule; without
population identity data the default random-parent rule and explicit maternal
inheritance are valid and have no effect on identity or random-number consumption.

## Possible identities and coverage

Infer possible identities from all starting identities plus reference ancestry
groups with any positive initial share. Explicit `groups` extends this set with
identities a custom function may introduce; it never removes inferred groups.
All-zero ancestry columns cannot introduce a group through parental averaging.
Keep this registry fixed during a run. Validate grouped mortality/fertility/
search coverage and categorical matching levels against the full registry,
including identities not yet present. Shared schedules cover everyone and group
rate multipliers use 1.0 when omitted. Newborn identity selects all subsequent
demographic schedules, matching effects and eligibility rules immediately.

Declared identities do not add ancestry components. A new `mixed` identity may
have ancestry expressed entirely in A and B. Preserve existing group-column
types when the extended labels can be represented in that type; otherwise
promote identity labels to character. Restarts infer the registry again from
their population and configured rule; callers reusing custom rules retain the
explicit declaration of possible outputs.

## Custom functions and timing

`fun` receives named arguments `mother_group`, `father_group`, `mother_ancestry`,
`father_ancestry`, `child_ancestry`, `child_sex` (`"f"` or `"m"`), and
`parents_in_union` (one nonmissing logical). The latter is TRUE when the parents
are in an ongoing union at birth, including a union formed on demand for that
birth. It is FALSE for nonmarital births, including those with a known father,
and births with an unknown father. It describes birth-time status, not status
at conception. Older callbacks without this argument or `...` remain supported.
Ancestry
inputs are the rule's own dimension's shares, numeric vectors named like its
public ancestry columns (`A`, `B`, `unknown`). Missing fathers have an NA group
and NULL ancestry.
Functions may use `...` for unused arguments and ordinary R RNG draws. They
receive copied inputs, not mutable population/queue state, and must not mutate
simulation/global state. Custom functions must handle missing inputs as needed,
using the mother's group as the fallback. Errors, missing/invalid returns,
or identities outside the fixed possible set fail explicitly.

Sex is drawn and ancestry calculated before identity assignment. The female
birth probability is configured through fertility's `ffab` (default 0.4886);
the same-sex/opposite-sex rules and custom callbacks use that realized sex.
Each custom call returns exactly one nonmissing character label or finite
numeric code in the registry. Synchronize R and native RNG state around calls so ordinary R
draws participate in the same stream and explicit seeds restore caller state,
including errors. Assignment finishes before scheduling the newborn. Ancestry
shares are unaffected by identity assignment. Existing people are not reassigned.

## Worked cases and validation

- Default rule, mother A/father B: A or B with probability 0.5 each. Unknown
  father gives A without an identity draw; equal parental groups need no draw.
- Mother A/father B: father method gives B; same-sex gives A for daughters and
  B for sons; opposite-sex reverses those outcomes. Missing father gives A.
- Child A ancestry 0.35, B 0.15, unknown 0.5: ancestry probabilities are 0.7/0.3.
- Parents both identify as A but have B ancestry: ancestry sampling can assign B.
- Hierarchy A/B with A first and passing 1 always gives B; missing rank gives
  the maternal identity.
- A custom function may return `mixed` only if inferred or declared in `groups`.
  All grouped schedules and categorical matching must cover it before execution.

Tests cover each rule and fallback, analytical sampling probabilities, passing,
custom multinomial-style assignment, new identities and rate/matching coverage,
native/R RNG agreement, errors, empty populations, restarts, independent runs,
default/explicit random-parent equivalence and explicit maternal histories, and conserved ancestry shares.

## Directional matching and sex-specific rules

Matching always uses calibrated male-group weights (D075), with directional
adjustments specified by delta (D074). The obsolete female-choice inheritance
warning has been removed. The rules `mother`, `father`, `same_sex`, and
`opposite_sex` depend on parental sex, so directional pairing differences can
still affect identity composition. Inspect those differences alongside group
preferences, delta, sex-specific pool composition and inheritance rules;
weighted matching does not guarantee identical maternal and paternal outcomes.

## User-reported inheritance assessment (October 27, 2026)

The user completed the inheritance-methods section of the milestone 3 report.
Replications sharing a baseline showed sensitivity to that baseline: directional
union imbalances and corresponding cross-group birth imbalances varied with
the baseline seed, including changes in direction. These observations do not
establish a calibration bug or a systematic effect of the A:B size ratio.

With independently generated baselines per replication and maximum starting
age reduced from 100 to 50, the user reported similar ending group and ancestry
shares for mother, father, random-parent, same-sex, and opposite-sex inheritance.
Both changes were made together, so their contributions were not isolated.
Between-replication Monte Carlo standard errors in this design include baseline
generation variability; fixed-baseline errors condition on that baseline.

The user selected random-parent inheritance as the default with maternal
fallback (D076). This is a default choice, not a guarantee of directional
balance or absence of sensitivity to starting populations. See the
[user-authored report](../validation/milestone3.qmd) for the assessment.


## Calendar-varying inheritance (D123)

Fertility's `inheritance` argument also accepts a tibble with `time_start`,
`time_end` and list-column `model`, containing complete inheritance objects.
Select the active rule at birth, without reassigning existing identities or
changing ancestry. All declared possible identities from all intervals enter
one shared simulation registry, including future custom outputs. Other grouped
processes and matching models must cover that registry upfront. Direct models
and existing static method-name aliases remain supported. See
[calendar-varying models](time-varying-models.md) for interval validation,
coverage, callback behavior, worked examples and dynamic matching fallbacks.
