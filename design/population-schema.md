# Person and union tables

All public tables are returned as tibbles. Base data frames remain valid input;
constructors normalize the table class while preserving column contents and
order. Both empty and populated `persons` and `unions` tables follow this rule.

Updated October 19, 2026. The schema is implemented under D012/D015/D016/D023.
Normalization and structural validation are shared by population construction
and attachment. See the [scientific contract](scientific-contract.md) for
event and restart semantics.

## Persons

One row per person, including deceased persons.

| Column | Meaning |
| --- | --- |
| `person_id` | Unique positive whole-number person ID |
| `sex` | Character labels: `"m"` = male, `"f"` = female |
| `birth_time` | Birth time in years |
| `death_time` | Death time in years; `NA` for a living person at the result endpoint |
| `mother_id` | Mother's person ID; `NA` if unknown |
| `father_id` | Father's person ID; `NA` if unknown |
| `group` (optional) | Nonmissing character label or finite numeric code; newborn assignment uses fertility's inheritance model (random parent by default, with maternal fallback) |

Omit marital status and latest-union pointers. Derive relationship state from `unions`; the engine may cache derived state for performance without requiring another user-facing column. Unknown earlier histories must not be silently interpreted as known lifetime singlehood. Optional persistent fertility selection weights use `hetfert` under D102;
sibling and last-born-child pointers are omitted. Group labels must be nonempty; their input type is preserved. See [groups](groups.md) for schedule lookup and inheritance rules.

## Uniform-age population generation

`generate_uniform_population()` draws independent continuous ages uniformly
between `min_age` (default zero) and `max_age` (default 70), in years (D168).
Both bounds must be finite nonnegative numeric scalars, with
`min_age <= max_age`. Birth time is `end_time - age`; everyone is alive at the
reference time, with unknown parents and no unions. Sex and groups retain
their independently specified sampling distributions. This is an initial
population, without any claim that its age distribution is stable.

For example, `min_age = 18, max_age = 65, end_time = 2000` produces ages in
18–65 years, with birth times in 1935–1982. Equal bounds of 25 produce birth
times exactly 1975. Empty populations retain the typed empty schema and still
validate the bounds. Invalid or reversed bounds fail before drawing random
numbers. With `min_age = 0`, existing calls retain their output and RNG state
for the same seed when other arguments are unchanged. The new argument follows
`end_time`, before `female_probability` and `groups`; positional calls using
those later arguments must account for the new order.

## Ancestry

An optional `ancestry` component is a named list with one tibble per dimension
declared with `ancestry = TRUE`. Each has one row per person, including deceased
people, keyed by `person_id`. Known shares are in columns named by category; `unknown`
records unknown ancestry. Shares are finite and in `[0,1]`, summing to one per
person within `1e-8`. Supplied rows are aligned to persons; reference groups may
differ from current identity groups. Without supplied shares, grouped populations
initialize everyone at 100% of their recorded group. Without groups or supplied
shares, no table is created. Newborns average parental shares; a missing parent
contributes half to unknown. See [ancestry](ancestry.md) for examples and restarts.

## Unions

One row per union episode; repeated unions have distinct IDs.

| Column | Meaning |
| --- | --- |
| `union_id` | Unique positive whole-number union ID |
| `wife_id` | Female partner's person ID (`sex = "f"`) |
| `husband_id` | Male partner's person ID (`sex = "m"`) |
| `start_time` | Union formation time in years |
| `end_time` | Union ending time in years; `NA` for an ongoing union |
| `end_reason` | `"death_husband"`, `"death_wife"`, or `"dissolution"`; `NA` for an ongoing union |

At the user's direction, partner roles are explicit: `wife_id` must reference
a female person and `husband_id` a male person, for ended and ongoing unions.
These column names do not imply legal marriage; the relationship remains a
union. The former unordered partner columns are no longer accepted as the
required schema. Existing tables must be mapped by sex, not blindly renamed
by column position. For example, an old row with male partner 1 and female
partner 2 must place partner 2 in `wife_id` and partner 1 in `husband_id`. `unions = NULL` on person-table input denotes no supplied unions, as already accepted. NULL unions and zero-row, zero-column data frames are normalized to a typed empty union table. Other inputs must supply every schema column; extra columns and data-frame subclasses are preserved.

Optional union diagnostics `.log_weight_<group>` are populated for new unions
when the process sets `diagnostics = TRUE` and weighted choice runs in a
grouped population; they are not required core columns. Historical unknowns
remain NA. See
[diagnostic semantics](union-formation-models.md#optional-union-diagnostics).

## Time representation

Dates are numeric years, including fractions, with a user-supplied common
origin. Calendar years such as 2000 and relative time such as zero are both
valid; dates before the origin may be negative. Do not round dates to integer
months or birthdays. Execution is continuous under D024, and midpoint fertility
compensation under D028 has no monthly grid.

`end_time` is required and must be a single finite numeric reference time in
years (D166). Omission, `NULL`, missing/nonfinite values, complex numbers,
vectors and arrays are invalid. The next simulation starts there; after a run,
it records the new endpoint. Tables must not contain recorded events later
than this reference time. Negative and fractional reference times are valid
when consistent with the histories.

For example, a population with births at -25 and -23 and `end_time = 0` can
start a simulation at zero. With `end_time = 2000`, a 50-year run ends at 2050.
Supplying `NULL` or omitting `end_time` fails at construction. Direct validation
and model attachment reject a population whose endpoint was edited to `NULL`.

Recorded lifetimes and analysis periods are half-open. Exposure is exact overlap
in person-years, split by age and period boundaries. Population counts include
births at an instant and exclude deaths at that instant. See
[demographic analysis](demographic-analysis.md) for observation details.

## Implemented validation scope

Require the listed columns and their types, finite dates (except allowed NA), unique positive whole-number IDs, sex labels, and references to existing persons. Parents cannot reference themselves. Union partners must be distinct, with female wife and male husband references. Death cannot precede birth; union ending cannot precede formation, and formation cannot precede a partner's birth. End reasons must match whether a union has ended. Death reasons identify the
partner whose executed death event ended the union; if both die at the same
time, the first processed death closes it under the existing event ordering. When a reference time is supplied, recorded events cannot lie after it. Equal event dates are permitted; monthly alignment is not required by table validation.

This is structural validation, not a full demographic eligibility contract.
Simulation additionally rejects conflicting ongoing unions and ongoing unions
with dead partners; it enforces the implemented search eligibility rules.
Broader validation of historical overlaps, parent ages/sexes, and ended-union
consistency is not provided by the structural validator. Ancestry cycles are rejected by the shared kinship module when simulating or extracting kin; the structural population validator itself does not traverse ancestry. No missing core fields are silently filled in. Unknown parents use NA, not SOCSIM's zero sentinel.

Normalization and validation are separate operations. `hammel_population()` normalizes absent unions before checking the constructed object. `validate_population()` verifies a `hammel_population` without mutation and returns it invisibly. `set_population()` accepts only a population object, validates it, normalizes
its person/union tables to tibbles, and stores it as `model$population` while
preserving extra components. Population construction/validation lives in
`R/hammel_population.R`, alongside the population setter (D165).


## Typed union episodes (D086)

Cohabitation-enabled rate runs add union_type (cohab/marriage) and
previous_cohab_id. Conversion closes cohab with end_reason conversion and
starts marriage with the same partners at the same time. The new row points
to the cohab ID; links must be unique and consistent. Both episodes remain in
history, with continuous partnered exposure. Untyped mode retains its six core
union columns. Enabling cohabitation treats legacy untyped history as marriage.


## Optional fertility weight column

With `hetfert = het_fert(cv, alpha)` on fertility (D102), `persons$hetfert` stores
persistent dimensionless positive female selection weights. The simulator draws
values if the entire column is absent; a supplied column requires complete,
finite, strictly positive female values. Existing values are retained on
restart. Newly drawn male entries are NA, and existing male values are unused.
When the feature is disabled, output omits the column even if supplied (D104).
Re-enabling from that output draws fresh weights for every woman. See
[heterogeneity](fertility-heterogeneity-notes.md) for assignment and inheritance.

## Group-transition histories

Grouped populations always include `transitions`, a tibble with columns
`person_id`, `time` (years), `group` (the dimension), `origin`, and
`destination`. It contains one row per individual move, preserves all earlier
rows across restarts, and has the same schema when empty. `group`, `origin`
and `destination` are character. Histories are validated against life dates,
the known endpoint, and chronological chaining to the current category in each
dimension. `hammel_population(..., transitions = ...)` accepts
base data frames. See [group transitions](group-transitions.md) for event and
observation semantics. Ancestry proportions are unaffected by moves.
