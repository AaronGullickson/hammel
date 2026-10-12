# Person and union tables

All public tables are returned as tibbles. Base data frames remain valid input;
constructors normalize the table class while preserving column contents and
order. Both empty and populated `persons` and `unions` tables follow this rule.

Status: initial schema direction accepted October 12, 2026 (D012). Shared structural validation and empty-union normalization are implemented in both population entry points. Based on SOCSIM's person and marriage histories, with engine pointers and group omitted.

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

Omit marital status and latest-union pointers. Derive relationship state from `unions`; the engine may cache derived state for performance without requiring another user-facing column. Unknown earlier histories must not be silently interpreted as known lifetime singlehood. Group, fertility multipliers, sibling pointers, and last-born-child pointers are also omitted for now.

## Unions

One row per union episode; repeated unions have distinct IDs.

| Column | Meaning |
| --- | --- |
| `union_id` | Unique positive whole-number union ID |
| `partner1_id` | First partner's person ID |
| `partner2_id` | Second partner's person ID |
| `start_time` | Union formation time in years |
| `end_time` | Union ending time in years; `NA` for an ongoing union |
| `end_reason` | `"death"` or `"dissolution"`; `NA` for an ongoing union |

Partner columns do not encode wife/husband roles. `unions = NULL` on person-table input denotes no supplied unions, as already accepted. NULL unions and zero-row, zero-column data frames are normalized to a typed empty union table. Other inputs must supply every schema column; extra columns and data-frame subclasses are preserved.

## Time representation

User-facing table dates are numeric years, including fractional years, using a common origin. For example, monthly dates relative to year 2000 can be represented as `2000 + month_index / 12`. Relative years are also possible; the calendar origin convention and metadata still need specification. Times before the origin may be negative.

Year-based storage does not change D004 or D008: the first engine still uses monthly time and truncates candidate waiting times to integer months. Use generic integer step indices if needed: an event at index `k` is exposed as `origin + k * dt` years, with configured `dt = 1/12` for the initial milestone. Month assumptions must be confined to explicit time configuration, not table schemas or demographic algorithms. Do not round to a fixed number of decimal places: one month is generally not exactly representable as a finite decimal, and ordinary floating-point arithmetic also introduces small errors. Boundary comparisons and input alignment need a documented numerical policy before implementation.

Exposure must be recorded in person-years, but year-based dates alone do not define exposure accounting. Event placement within a month, observation boundaries, partial intervals, and age-band splitting require their own rules.

Population reference-time metadata uses `end_time` in years; the model takes its starting time from that value. The former month-named arguments have been removed rather than reinterpreted. NULL continues to represent an unknown reference time.

## Implemented validation scope

Require the listed columns and their types, finite dates (except allowed NA), unique positive whole-number IDs, sex labels, and references to existing persons. Parents cannot reference themselves. Union partners must be distinct. Death cannot precede birth; union ending cannot precede formation, and formation cannot precede a partner's birth. End reasons must match whether a union has ended. When a reference time is supplied, recorded events cannot lie after it. Equal event dates are permitted; monthly alignment is not required by table validation.

This is structural validation, not a full demographic eligibility contract. Parent-age and sex constraints, ancestry cycles, overlapping unions, death-related union consistency, and process-specific eligibility still require scientific specification. No missing core fields are silently filled in. Unknown parents use NA, not SOCSIM's zero sentinel.

Normalization and validation are separate operations. `hammel_population()` normalizes absent unions before checking the constructed object. `validate_population()` verifies a `hammel_population` without mutation and returns it invisibly. `set_population()` accepts only a population object, validates it, and stores it unchanged as `model$population`. These functions and helpers live in `R/hammel_population.R`.
