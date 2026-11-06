# Demographic analysis

Updated November 5, 2026. All outputs are tibbles; base data frames are valid
inputs. Dates and age boundaries use years, and exposure uses person-years.
See the [scientific contract](scientific-contract.md) for model assumptions.

## Exposure, counts, and rates

For a population simulated from 2000 through 2100:

```r
fertility_summary <- result |>
  calculate_fertility_rates(
    time_breaks = seq(2020, 2100, by = 20),
    age_breaks = 0:100
  )
```

`calculate_mortality_rates()` takes the same arguments and returns deaths and
`rate_mortality` for both sexes; fertility returns female rows only.
Neither function includes the other process's counts or rates.

Both functions append `Inf` to finite final age breaks. Thus `0:100` ends with
`[100, Inf)` and `0:110` ends with `[110, Inf)`; those ages are pooled, not
removed. Existing `Inf` breaks are retained.

| Column | Meaning |
| --- | --- |
| `time_start`, `time_end` | Half-open observation period |
| `sex`, `age_start`, `age_end` | Sex and half-open age band |
| `person_years` | Exact recorded-lifetime exposure in the cell |
| `deaths`, `rate_mortality` | Deaths and deaths/person-years |
| `births`, `rate_fertility` | Maternal births and births/all female person-years |

Fertility returns female rows only. Female exposure includes the time
excluded by birth spacing. Zero events with positive exposure give a zero rate;
zero exposure gives `NA`. Keep those situations distinct.

`by = "group"` returns a separate grid for every group in the person
history, preserving character labels or numeric codes. The default pools
groups. It requires a `group` column; no pooled group row is added. Exposure
and deaths follow each person's group; births follow the mother's group, even
when the child's group differs. It can be combined with `by = "union_type"`:
each spouse's formations and exposure follow their own group. Complete histories
are retained during extraction, so cross-group family links remain available.

```r
group_counts <- calculate_fertility_rates(
  result, seq(2020, 2100, 20), 0:100,
  by = c("group", "union_type")
)
group_tfr <- calculate_tfr(group_counts)
group_mortality <- calculate_mortality_rates(
  result, seq(2020, 2100, 20), 0:100, by = "group"
)
group_life_tables <- calculate_life_table(group_mortality)
```

Groups appear in first-occurrence order, with all period/age/sex cells retained,
including zero exposure. Empty grouped histories produce typed empty tables.
See [groups](groups.md) for exact accounting examples and observation rules.

Births with unknown mothers cannot be assigned maternal age or group. Their period
counts are in `attr(summary, "unassigned_births")`, also a tibble. Initial
persons born within the observation window can appear there if parents are
unknown; a baseline is not complete historical birth follow-up. The attribute
remains one pooled count per period when split by group; the child's group
does not establish an unknown mother's group.

## Stratifying results with `by` (D133)

Every rate calculator and `calculate_population_size()` take one argument,
`by`, a character vector of categorical variables to split results over, for
example `by = c("group", "union_type")`. `NULL` (the default) pools them. Bands
on continuous scales keep their own arguments (`time_breaks`, `age_breaks`,
`dur_breaks`, and `parity_breaks` in fertility). `by` replaced the
`by_group`, `by_union_type`, `by_previous_union` and `by_sex` arguments on
November 5, 2026.

| Function | `by` values | Always split by |
| --- | --- | --- |
| `calculate_mortality_rates()` | dimensions, `union_type`, `previous_union` | `sex` |
| `calculate_fertility_rates()` | dimensions, `union_type`, `previous_union` | women only; `parity_breaks` optional |
| `calculate_union_rates()` | dimensions, `previous_union`, `union_type` (type entered) | `sex` |
| `calculate_marriage_rates()` | dimensions, `previous_union` | `sex` |
| `calculate_transition_rates()` | `sex`, `union_type`, `previous_union`, other dimensions | `origin`, `destination` (the dimension chosen by `group`) |
| `calculate_population_size()` | `sex`, dimensions, `union_type`, `previous_union` | — |
| `calculate_dissolution_rates()` | `<dim>_w`, `<dim>_h`, `union_type` | — |
| `calculate_cohab_convert_rate()` | `<dim>_w`, `<dim>_h` | — |

"Dimensions" are declared group dimensions by name, such as `race` and
`region` (D134); the examples below use a single dimension named `group`.
Unknown or repeated names are errors that list the accepted values, naming
the declared dimensions.

**Several dimensions (D134 stage 4).** Any number of dimensions may be listed,
each giving its own column; for union rates any combination of partners and
dimensions, such as `by = c("race_w", "region_h")`. Each dimension's
categories appear in order of first appearance (in persons, then its
transition history), as with one dimension, and the output is the full grid
of their categories, including combinations without exposure. Exposure is
split at every move in any listed dimension; a person's moves in several
dimensions at one instant are applied in declaration order. Internally the
listed dimensions form one combined category with a merged history; results
are then split back into one column per dimension, and union rates sum over
partner-dimension combinations that were not requested. Rows follow the
usual rule: the first `by` variable varies fastest after the fixed columns.
In `calculate_transition_rates()`, other dimensions in `by` split the origin
exposure by the person's category there; only moves in the chosen dimension
count, and duration counts from entering its current category, unaffected by
moves in other dimensions. The chosen dimension itself cannot be in `by`.

**Variables.** Each is the person's (or union's) state at the time, changing
only at recorded events:

- `sex`: `m` or `f`.
- A dimension's name (such as `race`): the current category, following that
  dimension's transition history.
- `union_type`: current union state, `unpartnered`, `cohab` or `marriage` for
  typed histories and `unpartnered` or `partnered` for untyped ones, always in
  that order. In `calculate_union_rates()` it is instead the type of union
  entered, `cohab` or `marriage`: it splits the formations, and both rows share
  the same unpartnered exposure.
- `previous_union`: how the most recent ended union ended, `none`,
  `separated` or `widowed`. Conversion is not an ending, and cohabitation
  endings follow the population's recorded `cohabitation_counts_as_previous_union`
  policy (D085). It applies to partnered time too: a woman in her first union
  is `none`. `pool_dissolved = TRUE` merges `separated` and `widowed` into
  `dissolved` in any calculator when `previous_union` is in `by`.
- `<dim>_w`, `<dim>_h` (such as `race_w`): the wife's and husband's current
  categories, for union-level rates. Directional: A-wife/B-husband and B-wife/A-husband are
  separate rows; sum them for an unordered pair.

Dimension names require declared groups. `union_type` and `previous_union` require
a `hammel_population` with union histories, and `previous_union` requires
nonoverlapping histories for each person. Union rates by entry type and
conversion rates require typed histories. Bare union tables passed to
`calculate_dissolution_rates()` need `group_w`, `group_h` or `union_type`
columns for those values.

**Exposure and events.** Exposure is split at every change in the requested
variables, and each event counts under the state immediately before it in the
simulation's event order (death, birth, dissolution, conversion, search, group
transition). Changes at the same instant from earlier events apply first:

- A death counts under the state just before it, with no same-time changes, so
  a death that ends a union counts as partnered.
- A birth counts in a union that starts at the same time (on-demand formation
  happens at the birth) and in a union ending at the same time by dissolution
  or conversion (birth precedes them), but not in one ending at the same time by
  the partner's death. It uses the mother's group before any same-time
  transition. Simultaneous births advance parity one at a time, in row order.
- A union formation counts under the previous-union category and group before
  it.
- A group transition counts under the union state after same-time union
  changes, which precede transitions.
- `calculate_population_size()` counts states after all events at each time.

**Row order.** Time varies fastest, then the age or duration band, then `sex`
(where it is always included), then parity, then the transition destination
and origin (destinations vary faster, so each origin's rows are adjacent),
then the `by` variables in the order given, earlier ones varying faster. In `calculate_union_rates()` with `dur_breaks`,
`previous_union` is added to `by` (last, unless listed), `none` rows use age
bands and `separated`/`widowed` rows use duration since the union ended.

**Example.** A woman born at 0 starts cohabiting at 20, converts to marriage at
22, and dissolves the marriage at 30, giving birth at 30. Her exposure by
`union_type` is unpartnered to 20, cohab 20–22 and marriage 22–30, then
unpartnered; by `previous_union` it is `none` to 30 and `separated` after.
The birth at 30 counts as `marriage` and `none`, because it precedes the
dissolution.

## Union-years and dissolution rates

Use `calculate_dissolution_rates()` to observe union episodes rather than people:

```r
union_counts <- result |>
  calculate_dissolution_rates(
    time_breaks = seq(2025, 2100, 5),
    dur_breaks = c(0, 5, 10, 15, 20, Inf)
  )
```

The tibble contains `time_start`, `time_end`, `dur_start`, `dur_end`,
`union_years`, `dissolutions`, and `rate_dissolution`. The rate is dissolutions
per union-year, directly comparable with duration-specific input hazards.
Each recorded episode contributes once, even if both partners are present;
repeat episodes are separate observations. Death stops exposure but does not
count as dissolution. Missing ends are censored at the observation endpoint.

Exposure is split exactly at calendar and duration boundaries. Events at an
internal boundary belong to the period/band beginning there; the last time
boundary is excluded. Finite final duration breaks automatically gain `Inf`.
Zero exposure gives an NA rate; zero dissolutions with positive exposure gives
zero. This is a cause-specific rate, not a cumulative probability of dissolution
when death can also end unions.

`by = c("group_w", "group_h")` splits by the wife's and husband's current
groups (D133). Rows are directional: sum A-wife/B-husband and B-wife/A-husband
rows for the unordered pair that dissolution inputs use (D079). A full
period/duration grid is retained for every combination of population groups.
The function also accepts a bare union data frame with `union_id`, `start_time`,
`end_time`, and `end_reason`; caller-provided follow-up must cover the window.
Bare tables need `group_w` and `group_h` columns for those variables, with
labels sorted in radix order. Summarize independent replications separately before binding
results, as for person-level Lexis counts. See
[union dissolution](union-dissolution.md#union-year-counts-and-dissolution-rates)
for worked accounting and censoring rules.

## Total fertility rate

```r
input_tfr <- model$processes$fertility |> calculate_tfr()
output_tfr <- fertility_summary |> calculate_tfr()
```

`calculate_tfr()` selects female rows when sex is present, omits open age bands,
and sums ASFR times each finite band's actual width by period. There is no
`age_range` argument; subset the summary to the desired age range first.
Remaining bands must be contiguous. The output reports the finite range used.
Any column that is not a band or measure, such as a group dimension (`race`,
`region`), automatically separates TFR calculations within each period and
union status, preserving labels and their type. This applies to process rates
as well as Lexis summaries; groups are not averaged together. Life tables and
PATFRs separate strata the same way. (Before stage 4 only a column named
`group` did, so tables split by a dimension with another name were pooled.)

If any ASFR in the included finite age bands is missing, TFR is `NA` for
that period/group/status, including when all ASFRs are missing. Known zero
rates contribute zero. For example, two five-year bands with rates 0.1 and
`NA` yield `NA`, not 0.5; rates 0.1 and zero yield 0.5. Missing rates in
excluded open age bands do not affect TFR. Inspect exposure and missing rates
alongside TFR; when summarizing with `na.rm = TRUE`, report the number of
nonmissing TFRs because the summary includes only complete estimates.
TFR warns if it omits a known nonzero upper-open rate, since the finite-range
result underestimates full-age fertility. Zero or missing omitted rates do not
trigger this warning. HFD preparation now maps 12- to [12,13) and 55+ to
[55,Inf), preserving their rates; TFR drops only the upper-open row and need
not equal HFD's full published TFR.

## Life tables

```r
input_life_table <- model$processes$mortality |> calculate_life_table()
output_life_table <- summary |> calculate_life_table()

output_life_table |>
  dplyr::filter(age_start == 0) |>
  dplyr::select(time_start, time_end, sex, ex)
```

Each period/sex/group/status combination present needs contiguous age bands from zero through an explicit open
interval. Lexis counts add that interval automatically; manually supplied rate
tables must define it themselves. The function treats mortality as piecewise
constant hazards, integrates exponential survival/person-years, and calculates
`ex = Tx/lx`. It does not apply conventional HMD infant adjustments.

The output includes input identifiers and `mx`, `n`, `ax`, `qx`, `lx`, `lx_end`,
`dx`, `Lx`, `Tx`, and `ex`. Missing mortality rates propagate into life expectancy.
A positive-exposure terminal cell with zero deaths has estimated rate zero;
under the open-interval model this gives infinite remaining life expectancy
if life-table survivorship there is positive. `na.rm = TRUE` does not remove
infinity. Report finite, missing, and infinite estimates separately rather than
silently excluding difficult replications.

Sparse older ages can make a 110+ life table unestimable. Pooling at 100+ may
make the summary usable, but a reference retaining detailed ages through 110
has different bands. Label that life-expectancy comparison descriptive unless
the terminal treatment is aligned. A pooled target ASMR requires exposure
weights within the pooled band, not an unweighted average of input rates.

TFR and life-table functions also accept their corresponding `$rates` tables.
TFR uses `rate_fertility` or fertility `rate`; life tables use `mx`, `rate_mortality`, or mortality
`rate`. Typed processes prevent passing mortality to TFR or fertility to a life
table. A bare `rate` column has no such type protection.

Input-rate summaries, including those passed a process object, summarize its
stored rate table, including strata added by `stratify_rates()` (D132). They do
not expand a shared table to population groups.

## Population size and observation boundaries

```r
population_sizes <- result |>
  calculate_population_size(times = seq(2000, 2100, by = 10))
```

This returns `time` and `population`. Counts include births and exclude deaths
at the requested instant, preserving input order and repeated times.
`calculate_population_size(result, times, by = "group")` adds `group` and
counts each group separately; `by` also accepts `sex`, `union_type` and
`previous_union`. Group requires group data, retains the group's input
type, and lists groups in first-occurrence order; times vary fastest.
Groups present in the history get zero counts where no members are alive.
There is no extra pooled row. Empty time vectors and empty group histories
return zero-row tibbles with the appropriate column types.

Lexis exposure splits `[birth_time, death_time)` directly at period/age
boundaries; fractional years and unequal bands are supported. Population
objects must have a known endpoint and cannot be queried beyond it. With bare
persons tables, the caller must ensure follow-up through the requested endpoint.
Missing death times are censored there.

## Replications and interpretation

Analyze each independent history before combining summaries. Stratification
is available through `by` (D133); arbitrary custom variables are not supported.
Stacking person tables
with repeated IDs would make parent references ambiguous and is rejected.
Retain a replication ID when binding derived tables:

```r
population_sizes <- purrr::map(
  results, calculate_population_size, times = seq(2000, 2100, 10)
) |>
  dplyr::bind_rows(.id = "replication")
```

For unnamed lists, `.id` uses character positions; named lists use their names.
Quantiles across runs describe replication variability, not confidence intervals
for the mean. For sparse cells, show the number of contributing replications and
exposure/events. Zero ASMRs cannot appear on a log scale. Longer periods can
improve precision for calendar-constant inputs, but do not establish absence
of age-specific bias.

Use aligned ASFR/ASMR alongside aggregate TFR/life expectancy. Distinguish
sampling variability from persistent differences due to unrealized births in
sparse cells (D101), initialization, or grouped ages. See [testing and validation](testing.md) and
the [accepted milestone 1 report](../validation/milestone1.qmd) and
[accepted milestone 2 report](../validation/milestone2.qmd).


When input fertility rates include `union_type`, `calculate_tfr()` returns a
separate row per status and period, preserving the status column. Each sum is
`rate * age-band width` over that status's finite contiguous age range. These
conditional TFRs describe continuous exposure to a status, not additive
components of the overall TFR. Tables without status retain their existing
behavior; default Lexis summaries use all female exposure, while
`by = "union_type"` also returns status-specific exposure as described below.

## Parent information

`get_parent_info(population)` enriches the full persons table with
`age_mother` and `age_father` at each person's birth, in years, and
`group_mother` and `group_father` when the population records groups.
Parents are matched by ID, including deceased parents. Unknown parents yield
missing values; founder rows remain present. For example, birth in 2000.5
and maternal birth in 1975.25 give `age_mother = 25.25`.

```r
parent_info <- get_parent_info(result)
cross_group_births <- parent_info |>
  dplyr::filter(!is.na(group_mother), !is.na(group_father),
                group_mother != group_father)
```

All original person rows and columns are preserved, as are group types.
Conflicting names cause an error. Empty populations return the enriched
zero-row tibble. No endpoint is required and the population is not modified.
Parent groups are recorded identities, not reconstructed past transitions.

## Union exposure and formation counts

### Detailed union records

`get_detailed_unions(population)` accepts a `hammel_population` and returns a
separate tibble containing all union rows and columns, in their original order,
with partner attributes matched through `wife_id` and `husband_id`. It adds
`birth_time_w`, `death_time_w`, `birth_time_h`, and `death_time_h`, plus
`group_w` and `group_h` when persons has a group column. Groups retain their
input type. Partner groups are reconstructed at each union's start from
transition histories; later moves do not relabel earlier unions. An exact-time
transition follows formation, so its origin applies. Without transition
history, use recorded person groups. A surviving partner's
death time remains NA. Dates remain in years and can be combined with union
`start_time` to derive ages at formation and the age gap.

An empty union table returns the same enriched zero-row schema. The input is
unchanged, including its union table and custom columns. Existing columns with
the names to be added cause an error rather than being overwritten. No known
population endpoint is required merely to attach recorded partner information.

For example, a union whose wife ID is 10 and husband ID is 30 receives birth
dates from those IDs regardless of their row positions in persons; repeated
unions reuse those attributes while retaining distinct episode records.

```r
detailed_unions <- get_detailed_unions(result)
```

### Current union type exposure and events

`by = "union_type"` returns rows labeled by `union_type`: `unpartnered`,
`cohab` and `marriage` for typed histories, or `unpartnered` and `partnered`
for legacy untyped ones, always all of them in that order.
There is no overall row. Omit the option to obtain pooled demographic rates
without a type column. Type-specific counts and exposure sum to pooled totals;
rates must be recomputed after pooling counts and exposure.

Exposure follows half-open episode intervals. At an exact conversion or
dissolution tie, births retain the ending type to match birth-first execution.
A birth at formation is assigned to the new episode; births at other endings
are unpartnered unless another episode covers that time. Deaths use the type
immediately before death, so union-ending deaths retain that type. Typed
histories must not overlap; untyped overlapping coverage is merged.

```r
counts <- calculate_fertility_rates(result, seq(2020, 2100, 20), 0:100,
                                by = "union_type")
tfr <- calculate_tfr(counts)
```

`calculate_tfr()` and `calculate_life_table()` automatically keep union types
separate. Conditional TFRs describe continuous exposure to one type and are
not additive components of pooled TFR. Life tables similarly describe
conditional mortality without modeling transitions between types. Empty cells
are retained, and zero exposure gives NA rates. Unknown-mother births remain
in the single period-level `unassigned_births` attribute, without a type guess.

For an untyped history where a woman enters at 0.5, loses her partner at 1.5,
and enters again at 2.5, [0,3) contributes 1.5 partnered and 1.5 unpartnered
years. Her partner's death at 1.5 is partnered; a birth at that death ending is
unpartnered. A birth at an ending caused by dissolution is partnered.

## Union entry rates by previous union

Use the separate helper for formations per unpartnered person-year:

```r
marriages <- calculate_union_rates(result, seq(2025, 2100, 5),
                                      0:100, by = c("group", "previous_union"))
```

By default, union counts pool previous categories and omit `previous_union`.
Add `"previous_union"` to `by` as above to stratify independently of group.
The stratified output returns period, age, sex, `previous_union`, optional individual `group`,
`person_years`, `unions`, and `rate_union`. Categories are `none`,
`separated`, and `widowed`, determined by the most recent recorded episode.
Each spouse contributes a start in their category immediately before formation.
Partnered time contributes no exposure. Men's rates describe matching outcomes;
female search inputs can exceed realized rates when searches fail.

For the widowhood example above, the woman has 0.5 unpartnered years in `none`
and 1 in `widowed`, with one formation in each category. The `separated` row has
zero exposure and starts, with an NA rate. Pool counts and exposure before
calculating a combined rate: `(1 + 1) / (0.5 + 1)`, rather than averaging the
category rates. Summing counts/exposure also gives pooled results across groups
or replications; Monte Carlo means and SEs across independent simulations should
be calculated from their separate replication estimates.

Histories before the window remain necessary to identify prior category.
No recorded union means `none`, conditional on supplied history. Overlapping
episodes are rejected for this helper because previous status is ambiguous;
the existing Lexis exposure helper still merges overlaps. Touching episodes are
allowed, and zero-duration episodes can contribute starts without exposure.
All time and age bins are half-open. See [remarriage](remarriage.md) for worked
examples, tie conventions, and schedule requirements.

This replaces the former `unions` and `rate_union` columns from
`calculate_fertility_rates()` and `calculate_mortality_rates()`, which provide
current-status exposure and their respective event counts and rates.


Optional `dur_breaks` on `calculate_union_rates()` measures separated and
widowed rates since the latest union ending (D084). `none` still uses
`age_breaks`. Supplying duration breaks automatically enables previous-category
stratification, even if its flag is FALSE; unlike clocks are not pooled.
The output has both boundary pairs, with NA in the unused pair. Age breaks
restrict only none rows; duration rows include all attained ages. Duration
origins use complete pre-window histories and reset after each ending.

```r
marriages <- calculate_union_rates(result, c(2025, 2100), seq(15, 70, 5),
                                      dur_breaks = c(0, 1, 3, 5, 10, Inf))
```

Time/clock bands are half-open; finite final duration breaks gain Inf. With
`dur_breaks = NULL`, the age-only schema and pooling behavior are unchanged.


With typed histories, `by = "union_type"` adds separate `cohab` and `marriage`
rows. The numerator counts entry from an unpartnered state; conversion is
excluded. Each type uses the same unpartnered exposure, so rates sum to the
pooled `rate_union`, but repeated exposure must not be summed across types.
Without this option the helper pools both entry types. Untyped histories support
pooled rates only. See [cohabitation](cohabitation.md).


## Marriage entry including conversion

Use `calculate_marriage_rates(result, time_breaks, age_breaks, by = NULL)`
to measure all marriage starts, including conversions from cohabitation. It
returns `marriages`, `person_years`, and `rate_marriage`. Person-years include
both unpartnered and cohabiting time, and exclude marriage. This helper supports
age/period/sex and optional individual group; it optionally stratifies by `previous_union`, but not duration. Cohabiting
exposure and conversion retain the category held before cohabitation began. Legacy untyped unions are interpreted as marriages.

| Helper | Events | Exposure |
| --- | --- | --- |
| `calculate_union_rates()` | Direct marriage or cohabitation entry | Unpartnered years |
| `calculate_marriage_rates()` | Direct marriage or conversion | Unpartnered plus cohabiting years |

For example, unpartnered from 0–1, cohabiting from 1–3, then married from 3–5
gives one partnership entry per unpartnered year and one marriage per three
nonmarried years. Conversions do not add events to `calculate_union_rates()`.
These are realized demographic rates, not exact copies of female search or
conversion hazard inputs. Age/group composition and failed partner searches
can affect comparisons with input schedules.


Dissolution summaries also accept `"union_type"` in `by`, optionally combined
with `"group_w"` and `"group_h"`. Separate cohabitation and marriage rows use their own
union-years and dissolution counts; summing counts and exposure across types
recovers pooled totals. Conversion ends cohabitation exposure without a
dissolution event. Marriage duration starts at conversion. Typed union histories
are required for splitting; the default pools types.


## Cohabitation-to-marriage conversion rates

```r
calculate_cohab_convert_rate(result,
  time_breaks = c(2025, 2100),
  dur_breaks = c(0, 1, 2, 5, Inf)
)
calculate_cohab_convert_rate(result,
  time_breaks = c(2025, 2100),
  age_breaks = seq(15, 70, 5),
  by = "group_w"
)
```

Supply exactly one axis. Duration is time since cohabitation entry; age is
woman's age, matching conversion schedules. `by` accepts `"group_w"` (the
woman's group, matching conversion inputs) and `"group_h"`.
Typed population histories are required. Output contains period and selected
axis boundaries, `union_years`, `conversions`, and `rate_cohab_convert`.
Exposure is cohabitation only, once per union episode. Only conversions count
as events; dissolution and death censor exposure, and marriage episodes are
excluded. This is a cause-specific rate per cohabitation union-year, not the
probability of eventual conversion. Empty cells remain; zero exposure gives NA.
Finite last axis breaks gain Inf, as with the other rate helpers.


## Births with fathers selected outside partnerships

`get_unpartnered_births(population)` returns parent-enriched person rows for
births with known mother and father where the mother had no active partnership.
It excludes both marriage and cohabitation coverage, without deduplicating
parental pairs: each birth is a father-selection event. All input person columns
and the parent ages/groups from `get_parent_info()` are retained in person-row
order. Founders and births missing either parent are omitted. Complete histories
are needed; no recorded partnership is interpreted as unpartnered, not proof of
never partnering. Formation coverage is half-open; exact dissolution/conversion
ending births remain partnered (birth-first), while death-ending births are
unpartnered unless another episode covers that instant. No observation window
is imposed; users can filter birth_time afterward. The input is not modified.


```r
father_crosstabs <- results |>
  purrr::map(get_unpartnered_births) |>
  dplyr::bind_rows(.id = "replication") |>
  dplyr::count(replication, group_mother, group_father)
```


## Parity-adjusted fertility summaries

`calculate_patfr(rates)` accepts conditional age/current-parity fertility rates.
`calculate_tfr(rates)` automatically delegates when parity bounds are present;
its output then has `patfr`, rather than `tfr`. Period, group, union_type and
legacy union_status grouping are the same as for ordinary TFR, and only female
rows are used when sex is present. Parity is progressed through, not grouped
into separate raw hazard integrals. A conditional index for one group/union type
holds that stratum fixed throughout the synthetic cohort's reproductive life.

```r
rates <- get_fertility_hfd("USA", 2000, by_parity = TRUE)
calculate_patfr(rates)
calculate_tfr(rates) # Same result, with patfr column.
```

The calculation uses HFD's period fertility-table midpoint convention, not
microsimulation. Start at parity zero at the reported lower age, count births
and move exposure through parity states; the terminal category admits repeated
births. Return NA for a stratum with any missing included rate. Exclude upper-
open age rows, warning when known nonzero fertility is omitted. Report the finite
age range used. The user's HFD zero-fill import choice and endpoint convention
remain distinct from the summary's missing-data/upper-age rules. See
[parity fertility](parity-fertility.md) for formulae, coverage and limitations.


## Fertility rates by observed parity

```r
asfr <- calculate_fertility_rates(
  result, time_breaks = c(2025, 2100), age_breaks = seq(15, 50, 5),
  parity_breaks = c(0, 1, 2, 3, 4, Inf)
)
```

`parity_breaks = NULL` pools parities. Otherwise, half-open integer bands start
at zero; finite final breaks gain Inf. `parity_start`/`parity_end` label current
parity before a birth. Exposure changes parity immediately after each birth,
including births before the observation window and births of children who later
died. Custom person parity columns are ignored. Complete maternal histories are
required; absent births cannot be reconstructed. Minimum-spacing time remains
in the denominator. Group and union-type options combine with parity bands.
Sum births and exposure, then divide, to pool bands; do not average their rates.
See [parity fertility](parity-fertility.md#observed-fertility-by-parity).

## Unrealized births under cell accounting

For populations returned by `simulate()` (D101),
`calculate_fertility_rates(..., unrealized = TRUE)` adds `unrealized_births`
and `rate_fertility_adjusted = (births + unrealized_births) / person_years`.
Unrealized births are counts of the times the fertility clock drew a cell that
had members but no eligible woman (`simulation$unrealized_births`; D128; before
November 4, 2026 they were expected values), so the adjusted
rate shows what the input schedule would have produced without that loss. The
observed `rate_fertility` remains the simulated outcome.

Each record covers one one-year period from the run start and one fertility
cell (an age year, actual parity when schedules are parity-specific, group and
union type). A record is added to the observation row that contains it, so
actual parities pool into broader parity bands and age years into wider age
bands. Values are `NA`, rather than approximated, for periods before the most
recent simulation's start (a restarted population keeps only its latest run's
record), for time or age breaks that cut through a record, and for union-type or
parity strata the fertility schedule did not use. Populations without an
unrealized-births record (never simulated, or simulated before D101) are
rejected.

```r
asfr <- calculate_fertility_rates(result, time_breaks = c(2050, 2150),
  age_breaks = 12:56, parity_breaks = 0:4, unrealized = TRUE)
```

## Directed group-transition Lexis rates

Use `calculate_transition_rates()` with exactly one age or duration axis:

```r
calculate_transition_rates(result, time_breaks = seq(2000, 2100, 10),
  age_breaks = seq(0, 100, 10), by = c("sex", "union_type"))
calculate_transition_rates(result, time_breaks = seq(2000, 2100, 10),
  dur_breaks = c(0, 1, 5, 10, Inf))
```

Output retains every directed pair among identities in the population and
transition history, with `person_years`, `transitions`, and `rate_transition`.
Age uses time since birth; duration uses time since latest group entry, or
birth without earlier moves. Exposure splits at moves and, when requested,
union-type boundaries. Directional event counts use the origin being exited;
all destinations from one origin share its denominator. Pool sex/type counts
and exposure before computing rates, and do not sum denominator copies across
destinations. Zero exposure gives NA; a positive-exposure cell without moves
has rate zero. See the [full contract](group-transitions.md#observed-transition-rates)
for endpoint, empty-output and history-completeness rules.

## Optional Poisson rate intervals

All seven event-rate calculators (mortality, fertility, union formation,
marriage, dissolution, cohabitation conversion and group transitions) accept
`conf.level = NULL` by default, retaining their existing output. A finite numeric
scalar strictly between zero and one requests two-sided exact Poisson intervals
at that confidence level, appended as numeric `conf_low` and `conf_high` columns.
For example:

```r
calculate_fertility_rates(
  result, time_breaks = c(2000, 2010), age_breaks = c(15, 20, 25, 30, 35, 40, 45, 50),
  conf.level = 0.95
)
```

This requests 95% intervals. Each positive-exposure
row uses `stats::poisson.test(x = events, T = exposure, conf.level = conf.level)`.
Exposure is person-years or union-years according to the existing measure.
Zero exposure produces two NA bounds, even when an endpoint event is recorded;
empty results retain numeric zero-length confidence columns. Zero events with
positive exposure have a zero lower bound and a positive upper bound: zero
events in 100 years of exposure gives approximately [0, 0.036889] at 95%.

Intervals use the final stratum's counts and exposure, after pooling where
requested. They preserve all existing rate definitions, attributes and ordering.
Fertility intervals describe observed `births / person_years`, including when
`unrealized = TRUE`; unrealized births are would-be events rather than
observed births and do not enter the interval. These are marginal Poisson-model
intervals, not simultaneous intervals, Monte Carlo SEs, or a guarantee of
coverage under demographic feedback, shared partners or endogenous exposure.

Validation compares every calculator with `stats::poisson.test()` and checks
zero events, zero exposure, empty results, pooling/strata, fertility parity and
unrealized outputs, invalid confidence levels and unchanged default outputs.

## Detailed and pooled previous unions (D126)

Union and marriage rate calculators use `none`, `separated`, and `widowed`
when `by = "previous_union"`. Set `pool_dissolved = TRUE` to return
`none` and `dissolved` instead, together with `by = "previous_union"`.
Duration breaks also enable previous-union stratification.
Dissolved pools separation and widowhood. Events and exposure are summed
within each period, sex, group, age/duration and union-type cell before computing
rates and optional exact Poisson intervals. Confidence bounds are never averaged.
Duration remains elapsed time since the most recent ending. Actual histories
and the recorded cohabitation policy are unchanged. Rate and probability inputs
accept the same detailed or pooled category sets.
