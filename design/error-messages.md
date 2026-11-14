# Error and warning messages

Status: style, architecture, condition classes and native message tags
selected by the user on November 14, 2026 (D183); migration in progress under
issue #25. Until the migration finishes, unconverted call sites still use base
`stop()` and `warning()` with `call. = FALSE`.

This guide covers every error and warning Hammel raises, in R and from the
native engine. Most are input validation errors, but the same rules apply to
failures detected during a simulation, database download failures, internal
consistency failures and all warnings. It follows the
[tidyverse style guide for error messages](https://style.tidyverse.org/errors.html),
adapted to Hammel's terms and units.

## Goals

A message should let a user fix the problem without reading Hammel's source:

1. Say what is wrong, in the user's terms (arguments, columns, categories,
   years), not the engine's (cells, banks, routes, codes).
2. Show the offending value: the missing column, the bad category, the time
   of the coverage gap.
3. Suggest a fix when one is known.
4. Name the function the user called.

## Anatomy

Messages are built with `cli` and raised with `rlang` through Hammel's helpers
(see [Raising conditions](#raising-conditions)). A message has a header and
optional bullets:

```
Error in `set_mortality()`:
! `rates` is missing required columns.
✖ Missing: `age_end` and `rate`.
ℹ Mortality rates need `sex`, `age_start`, `age_end` and `rate` columns.
```

- **Header (`!`).** One sentence stating the problem or requirement. Prefer
  "must" for a known requirement (`` `years` must be a single nonnegative
  number. ``) and "can't" when describing an impossible request (`` Can't find
  group dimension `region`. ``).
- **Problem bullets (`✖`, name `"x"`).** What is actually wrong, with values.
  One fact per bullet.
- **Information bullets (`ℹ`, name `"i"`).** Context or the fix: what is
  expected, which function or argument to use, where to read more.
- Bullets are optional. A header alone is fine when it already shows the
  value and the fix is obvious.

## Wording rules

- Use plain language and complete sentences ending with a period, in the
  header and in bullets.
- Refer to things with cli inline markup so they are styled consistently:
  `{.arg rates}` for arguments, `{.field age_end}` for columns,
  `{.val {value}}` for values, `{.fn set_groups}` for functions,
  `{.cls hammel_model}` for classes.
- Interpolate values rather than pasting them into the format string, so
  braces in user data print literally.
- Report units as the public API does: times and durations in years,
  exposure in person-years, rates per person-year.
- Show at most about five offending values, then a count of the rest
  (`and 12 more`). cli truncation of interpolated vectors does this.
- Do not blame the user ("you supplied the wrong ...") or apologize. State
  the requirement and the actual value.
- Do not mention internal objects, engine codes, file names, or decision
  numbers. Mention a design concept only when the user can act on it, for
  example a model schedule's calendar coverage.
- Keep the wording of a message identical wherever the same check is made.
  Shared checks live in one helper for that reason.

## Naming the function

Errors name the exported function the user called, not the internal
validator that detected the problem.

- Internal helpers that can raise take a `call` argument defaulting to
  `rlang::caller_env()` and pass it on unchanged to further helpers and to
  `abort_hammel()`.
- An exported function that validates directly needs no extra code: the
  helper's default `caller_env()` is the exported function's frame.
- A setter that delegates to its constructor (`set_mortality()` calling
  `mortality()`) should name the setter. The constructor and setter both call
  an internal builder that takes `call`, so each names itself.
- Validation repeated inside `run_simulation()` (for example re-attaching an
  edited process) names `run_simulation()`.
- Native engine errors name the R function that made the native call (see
  [Native engine messages](#native-engine-messages)).

## Condition classes

Every condition Hammel raises has a specific class and a base class:
`c("hammel_error_<kind>", "hammel_error")` for errors and
`c("hammel_warning_<kind>", "hammel_warning")` for warnings. The classes let
users handle one kind of failure in code, let tests check the kind of
failure instead of matching message text, and can carry data about the
failure as fields of the condition object.

Classes are public API: once released, a class name and the meaning of its
fields should not change without a decision entry.

**Error classes:**

| Class | Raised when | Typical fields |
| --- | --- | --- |
| `hammel_error_argument` | A scalar argument or option has the wrong type, length or value | `argument` |
| `hammel_error_table` | A data-frame input is not a data frame, lacks or duplicates columns, or has invalid column values | `argument`, `columns` |
| `hammel_error_rates` | A rate or probability schedule has invalid bands, strata, overlaps or gaps within a period | `argument` |
| `hammel_error_coverage` | A rate or model schedule does not cover the requested simulation time | `argument`, `time` |
| `hammel_error_groups` | Group declarations, dimensions or categories are missing, unknown or inconsistent | `group`, `categories` |
| `hammel_error_population` | Person, union, transition or ancestry histories are invalid or inconsistent | `table` |
| `hammel_error_specification` | A model or component model (matching, weight, inheritance) is incomplete or inconsistent | `argument` |
| `hammel_error_database` | An HMD/HFD request fails: credentials, network, or unavailable data | `database` |
| `hammel_error_runtime` | A user-supplied input fails during a simulation, such as a custom function returning invalid values | `process` |
| `hammel_error_internal` | An internal consistency check fails; this indicates a bug in Hammel | none |

**Warning classes:**

| Class | Raised when |
| --- | --- |
| `hammel_warning_no_effect` | An input is accepted but has no or only partial effect (an unused weight term, a cohabitation setting without conversion) |
| `hammel_warning_estimate` | A computed summary is biased or incomplete for the stated reason (TFR omitting an open age band) |
| `hammel_warning_simulation` | A simulation completed but could not fully meet a target (matching calibration shortfalls) |
| `hammel_warning_performance` | A setting is valid but likely to be slow |
| `hammel_warning_deprecated` | A retired option was supplied and ignored |

Each class gets at least one test. A condition fits one specific class; when
two seem to apply, use the one that tells the user what to change.

## Raising conditions

Shared helpers live in a dedicated R script, `R/errors.R`, authorized by the
user on November 14, 2026 (D183):

- `abort_hammel(message, kind, ..., call)` wraps `cli::cli_abort()`. `kind`
  is the short class name, such as `"table"`; the helper adds the
  `hammel_error_<kind>` and `hammel_error` classes, passes `call`, and stores
  named data in `...` as condition fields. `message` is a cli character
  vector: the header followed by named `"x"` and `"i"` bullets, interpolated
  in the caller's environment.
- `warn_hammel(message, kind, ...)` does the same for warnings with
  `cli::cli_warn()`. Warnings do not name a calling function.
- A native-call wrapper translates engine errors and warnings (below).

The general checks in `R/checks.R` (`check_columns()`, `check_numeric_column()`,
`check_time()`) and domain `check_*`/`validate_*` helpers call these helpers,
so every caller gets the same class and wording. Call sites do not call
`stop()`, `warning()`, `cli_abort()` or `cli_warn()` directly once converted.

## Native engine messages

The engine's C++ errors and warnings are translated at the R boundary rather
than formatted in C++.

- A user-caused native message begins with a tag in square brackets, for
  example `"[weight_function] A custom fertility weight function returned 3
  values for 5 members."`. The R wrapper maps the tag to a condition class
  and an R-side message with bullets and a suggested fix. Adding a tag changes
  no simulation result.
- An untagged native error is an internal consistency failure. The wrapper
  raises `hammel_error_internal`, keeps the original text, and adds a bullet
  saying this is likely a bug, with the issue tracker URL from `DESCRIPTION`.
- Native warnings are re-raised through `warn_hammel()` with their class.
- Where a check made during the run could instead be made in R before the
  engine starts, moving it is a behavior change and is listed for the user
  rather than made silently.

**Tags.** These cover the user-caused native messages known when D183 was
made. A native message found during conversion that fits none of them gets a
new tag only with the user's agreement:

| Tag | Class | Example |
| --- | --- | --- |
| `weight_function` | `hammel_error_runtime` | Custom weight function returns the wrong length or an infinite value |
| `weight_value` | `hammel_error_runtime` | A computed weight is not finite and positive |
| `matching_utility` | `hammel_error_runtime` | A matching formula yields non-finite utilities |
| `schedule_coverage` | `hammel_error_coverage` | No fertility schedule covers a woman's parity; no weight model covers the run start |
| `limits` | `hammel_error_runtime` | The population, union history or kin counts exceed supported sizes |
| `weight_missing` | `hammel_warning_no_effect` | Missing weight values treated as a zero log-weight contribution |
| `slow_custom_weights` | `hammel_warning_performance` | Custom weights evaluated at every calibration refresh |

Native changes follow the C++ verification policy (D175): hashed native
results must be identical before and after with `-ffp-contract=off`.

## Testing

- Check the class and, where relevant, a field, rather than message text:
  `expect_error(mortality(rates), class = "hammel_error_table")`.
- Check how messages look with `expect_snapshot(error = TRUE)` for a
  representative message of each class. Snapshot changes are reviewed like
  code.
- Matching on text remains acceptable for a key value that must appear (a
  column or category name).

## Examples

Current messages and their proposed replacements.

**Missing columns** (`check_columns()`):

```
Before: `rates` is missing columns: age_end, rate.

After:  Error in `set_mortality()`:
        ! `rates` is missing required columns.
        ✖ Missing: `age_end` and `rate`.
```

**Calendar coverage** (`check_calendar_coverage()`):

```
Before: mortality rates have a calendar coverage gap at time 2010.

After:  Error in `run_simulation()`:
        ! Mortality rates don't cover the whole simulation.
        ✖ No rates apply at time 2010.
        ℹ The run covers 2000 to 2030. Add rate rows with `time_start` and
          `time_end` covering 2010, or shorten `years`.
```

**Setting that requires another setting** (`union_formation()`):

```
Before: Search `weights` require rate-based union formation.

After:  Error in `set_union_formation()`:
        ! `weights` can only be used with rate-based union formation.
        ✖ `initiation` is "on_demand".
        ℹ Use `initiation = "rate"` with search `rates`, or remove `weights`.
```

**Custom weight function during a run** (native):

```
Before: A custom fertility weight function returned 3 values for 5 members.

After:  Error in `run_simulation()`:
        ! A custom fertility weight function returned the wrong number of
          values.
        ✖ It returned 3 values for a cell with 5 members.
        ℹ The function must return one log weight per row of `members`.
```

**Internal failure** (native, untagged):

```
Before: Invalid matching group index.

After:  Error in `run_simulation()`:
        ! Internal error: Invalid matching group index.
        ℹ This is likely a bug in Hammel. Please report it at
          <https://github.com/AaronGullickson/hammel/issues>.
```

The suggested fixes above are drafts; the wording of each message is settled
as its file is converted.

## Migration plan

1. Add `rlang` to `Imports`, the helper script, and tests for the helpers.
2. Convert the general checks in `R/checks.R` and add `call` passing.
3. Convert the remaining call sites in reviewable batches, updating tests in
   the same change: process constructors and rate utilities; population and
   groups; matching; analysis (`calculate_*`); database import.
4. Translate native messages: add tags, the wrapper and the internal-error
   fallback, verified per D175.
5. Update the R source organization guide and the decision log when the
   migration is complete.
