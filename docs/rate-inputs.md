# Census of existing rate inputs

Updated October 30, 2026. This document describes implemented inputs through
accepted milestone 6, extended by milestone 7 calendar bands under D118.
All six annual-rate tables now optionally take `time_start`/`time_end`; see
[calendar-varying rates](time-varying-rates.md). `cohab_prob` also accepts
calendar bands under D122. Matching/inheritance arguments accept separate
[whole-model schedules](time-varying-models.md) under D123.

## Rates, bands and strata

There are six rate-based inputs. Tables use a numeric `rate` column. Public
ages and durations are years, with half-open bands `[start, end)`. Starts are
finite and nonnegative; ends exceed starts and may be `Inf`. Rates are finite
and nonnegative. Bands cannot overlap within a schedule's strata. Base data
frames are accepted and stored as tibbles.

Optional strata can be combined subject to the restrictions below. Additional
columns do not automatically create supported rate strata.

| Input and interface | Supported bands | Required strata or risk population | Optional strata |
| --- | --- | --- | --- |
| Mortality: `set_mortality(model, rates)` | Age only | `sex`: both `m` and `f` | Individual `group` |
| Fertility: `set_fertility(model, rates)` | Age only | Women; no `sex` column | Mother's `group`; `union_type`; current-parity bands `parity_start`/`parity_end` |
| Union-search attempts: `set_union_formation(model, "rate", rates = rates)` | Age, or duration since the most recent partnership ended | Currently unpartnered women | Woman's `group`; `previous_union`: `none`, `dissolved`, `widowed` |
| Union dissolution: `set_union_dissolution(model, rates)` | Duration since the current union episode began | One risk per union, carried by the woman | Symmetric `group_pair`; `union_type`: `cohab`/`marriage` |
| Cohabitation-to-marriage conversion: `union_convert_rate` on union formation | Woman's age or duration since the current cohabitation began | Currently cohabiting women | Woman's `group` |
| Individual group transitions: `set_group_transition(model, rates)` | Age or duration since latest group entry; age without earlier transition history | Directed `group_origin` and `group_dest` | `sex`; current `union_type` |

Age bands use `age_start`/`age_end`; duration bands use
`dur_start`/`dur_end`. No existing input supports a joint age-by-duration
surface. Fertility can cross age bands with parity bands.

## Restrictions and interpretation

### Mortality

Both sex schedules are required within each supplied group. Simulation requires
complete age coverage; uncovered ages and gaps are not treated as zero hazard.
There are no union-type, parity or duration strata.

See [the scientific contract](scientific-contract.md) and
[groups](groups.md).

### Fertility

Without `union_type`, rates apply to all women. With that column, each group
must supply either `{unpartnered, partnered}` or
`{unpartnered, cohab, marriage}`. Shared partnered rates apply to both partnership
types. Explicit type-specific rates require cohabitation mode. On-demand union
formation requires overall fertility without union-type strata; rate-based
formation also accepts overall fertility, with a warning.

Optional parity partitions must cover `[0, Inf)` without gaps or overlaps in
every group/union-type stratum. They describe current maternal parity, derived
from recorded births, rather than the birth order of the next child. A pooled
parity band applies its rate to each contained actual parity. Age gaps have
zero fertility intensity.

Rates are births per woman-year in the indicated cell, including women
temporarily ineligible under birth spacing. Cell accounting assigns births
to eligible women; intensity in cells with no eligible woman is reported as
unrealized. Persistent fertility weights affect assignment within cells.

See [parity fertility](parity-fertility.md) and
[cell-accounted fertility](fertility-cell-accounting.md).

### Union-search attempts

These are annual female search hazards, not guaranteed union-formation rates.
Failed matching leaves the woman unpartnered. Men have no search clock.
Without `previous_union`, a shared age schedule applies to all unpartnered
women. Supplying that column requires all three categories within each group.

`none` must use age. `dissolved` and `widowed` may each independently use age
or duration, but each category must use one clock consistently across all its
bands and groups. Each row supplies exactly one boundary pair, with unused
age/duration boundaries `NA`. The age columns remain part of the table schema
even when some categories use duration. Gaps have zero hazard.

Search duration runs from the most recent partnership ending, not simulation
start or first union ending. The `cohabitation_counts_as_previous_union`
setting controls whether ended cohabitations determine the previous category;
every actual partnership ending resets the duration origin even when a prior
legal-marriage category is restored. There is no input `union_type` split of
search hazards; the separate `cohab_prob` schedule selects intended type.

See [remarriage schedules](remarriage.md) and [cohabitation](cohabitation.md).

### Union dissolution

`group_pair` is symmetric: `A:B` and `B:A` select the same rate. Pair-stratified
tables must cover all unordered pairs in the possible identity registry.
There are no directional wife-group/husband-group strata. If `union_type` is
supplied, both `cohab` and `marriage` schedules are required within each pair
stratum; omission shares rates across types. Duration gaps have zero hazard.

Rates are per union-year, with one clock per episode. Conversion closes the
cohabitation and starts a marriage episode with dissolution duration zero.
Age and sex are not supported rate strata.

See [union dissolution](union-dissolution.md).

### Cohabitation-to-marriage conversion

Use one age or duration clock throughout the table, including across groups.
Gaps have zero hazard. There are no sex, previous-union, union-type or
group-pair strata. Conversion requires cohabitation enabled in rate-based
union formation. Rates are per cohabiting woman-year; conversion retains the
same partners and requires no new matching.

See [cohabitation](cohabitation.md).

### Individual group transitions

One age or duration clock applies to the whole table. Self-transitions are
invalid. Optional `sex` uses `m`/`f`; optional `union_type` uses
unpartnered/partnered or unpartnered/cohab/marriage, without mixing `partnered`
with explicit partnered types. Omitted columns apply universally. Missing
directed pairs, categories and gaps have zero hazard.

Duration resets at each group entry and uses birth as its origin without
earlier transition history. Union changes do not reset it. Rates are per
person-year in the origin group and applicable strata.

See [group transitions](group-transitions.md).

## Related probability schedule

| Input | Supported bands | Optional strata | Interpretation |
| --- | --- | --- | --- |
| `cohab_prob` on union formation | Woman's age only | Woman's `group`; `previous_union` | Probability that a search attempt seeks cohabitation rather than marriage |

This table uses `prob` in `[0, 1]`, not an annual rate. It requires rate-based
union formation. Bands cannot overlap within active strata; age gaps have
probability zero. Optional calendar bands follow the rate-table coverage rules
(D122), and select probability at search execution without relabeling existing
unions or changing the overall search hazard. If previous
categories are supplied, all three are required within each group. Type is
drawn before matching. The overall search intensity is split into intended
cohabitation and marriage searches by this probability.

See [cohabitation](cohabitation.md).

## Proportional rate multipliers

`group_prop_hazard` supplies static, dimensionless multipliers on mortality,
fertility, union-search and dissolution rates. Mortality, fertility and search
use individual group labels; dissolution uses symmetric group-pair labels.
Omitted factors use one. Fertility also permits separate group multiplier
vectors by the supplied union-type categories. Conversion and group transitions
have no corresponding multiplier argument.

These vectors have no age or duration bands of their own; they scale the
selected rate schedule. They are distinct from matching preferences,
inheritance models and persistent individual fertility weights.

See [groups](groups.md), [cohabitation](cohabitation.md), and
[union dissolution](union-dissolution.md).
