# Demographic analysis

Updated October 19, 2026. All outputs are tibbles; base data frames are valid
inputs. Dates and age boundaries use years, and exposure uses person-years.
See the [scientific contract](scientific-contract.md) for model assumptions.

## Exposure, counts, and rates

For a population simulated from 2000 through 2100:

```r
summary <- result |>
  calculate_demog_rates(
    time_breaks = seq(2020, 2100, by = 20),
    age_breaks = 0:100
  )
```

The function appends `Inf` to finite final age breaks. Thus `0:100` ends with
`[100, Inf)` and `0:110` ends with `[110, Inf)`; those ages are pooled, not
removed. Existing `Inf` breaks are retained.

| Column | Meaning |
| --- | --- |
| `time_start`, `time_end` | Half-open observation period |
| `sex`, `age_start`, `age_end` | Sex and half-open age band |
| `person_years` | Exact recorded-lifetime exposure in the cell |
| `deaths`, `rate_mortality` | Deaths and deaths/person-years |
| `births`, `rate_fertility` | Maternal births and births/all female person-years |

Birth counts and ASFR are `NA` for males. Female exposure includes the time
excluded by birth spacing. Zero events with positive exposure give a zero rate;
zero exposure gives `NA`. Keep those situations distinct.

`by_group = TRUE` returns a separate grid for every group in the person
history, preserving character labels or numeric codes. The default FALSE pools
groups. TRUE requires a `group` column; no pooled group row is added. Exposure
and deaths follow each person's group; births follow the mother's group, even
when the child's group differs. It can be combined with `by_union_type = TRUE`:
each spouse's formations and exposure follow their own group. Complete histories
are retained during extraction, so cross-group family links remain available.

```r
group_counts <- calculate_demog_rates(
  result, seq(2020, 2100, 20), 0:100,
  by_group = TRUE, by_union_type = TRUE
)
group_tfr <- calculate_tfr(group_counts)
group_life_tables <- calculate_life_table(group_counts)
```

Groups appear in first-occurrence order, with all period/age/sex cells retained,
including zero exposure. Empty grouped histories produce typed empty tables.
See [groups](groups.md) for exact accounting examples and observation rules.

Births with unknown mothers cannot be assigned maternal age or group. Their period
counts are in `attr(summary, "unassigned_births")`, also a tibble. Initial
persons born within the observation window can appear there if parents are
unknown; a baseline is not complete historical birth follow-up. The attribute
remains one pooled count per period when `by_group = TRUE`; the child's group
does not establish an unknown mother's group.

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

`by_group = TRUE` returns symmetric `group_pair` strata such as `A:A`, `A:B`,
and `B:B` (D079), pooling A-wife/B-husband with B-wife/A-husband. A full
period/duration grid is retained for every unordered pair of population groups.
The function also accepts a bare union data frame with `union_id`, `start_time`,
`end_time`, and `end_reason`; caller-provided follow-up must cover the window.
Bare tables need a `group_pair` column when grouping is requested; reversed
labels are normalized. Pair labels are canonical character strings in radix
order, including when the underlying group identities are numeric. Summarize independent replications separately before binding
results, as for person-level Lexis counts. See
[union dissolution](union-dissolution.md#union-year-counts-and-dissolution-rates)
for worked accounting and censoring rules.

## Total fertility rate

```r
input_tfr <- model$processes$fertility |> calculate_tfr()
output_tfr <- summary |> calculate_tfr()
```

`calculate_tfr()` selects female rows when sex is present, omits open age bands,
and sums ASFR times each finite band's actual width by period. There is no
`age_range` argument; subset the summary to the desired age range first.
Remaining bands must be contiguous. The output reports the finite range used.
A `group` column automatically separates TFR calculations within each period
and union status, preserving group labels and their type. This applies to
process rates as well as Lexis summaries; groups are not averaged together.

If any ASFR in the included finite age bands is missing, TFR is `NA` for
that period/group/status, including when all ASFRs are missing. Known zero
rates contribute zero. For example, two five-year bands with rates 0.1 and
`NA` yield `NA`, not 0.5; rates 0.1 and zero yield 0.5. Missing rates in
excluded open age bands do not affect TFR. Inspect exposure and missing rates
alongside TFR; when summarizing with `na.rm = TRUE`, report the number of
nonmissing TFRs because the summary includes only complete estimates.
HFD input TFR excludes the omitted open endpoint rows, so do not compare it
uncritically with the full published HFD TFR.

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
stored rate table. They do not apply `group_prop_hazard` settings or expand a
shared table to population groups. Account for those factors when preparing
effective input-rate comparisons; summaries of simulated histories already
reflect the configured multipliers.

## Population size and observation boundaries

```r
population_sizes <- result |>
  calculate_population_size(times = seq(2000, 2100, by = 10))
```

This returns `time` and `population`. Counts include births and exclude deaths
at the requested instant, preserving input order and repeated times.
`calculate_population_size(result, times, by_group = TRUE)` adds `group` and
counts each group separately. It requires group data, retains the group's input
type, and lists groups in first-occurrence order within each requested time.
Groups present in the history get zero counts where no members are alive.
There is no extra pooled row. Empty time vectors and empty group histories
return zero-row tibbles with the appropriate column types.

Lexis exposure splits `[birth_time, death_time)` directly at period/age
boundaries; fractional years and unequal bands are supported. Population
objects must have a known endpoint and cannot be queried beyond it. With bare
persons tables, the caller must ensure follow-up through the requested endpoint.
Missing death times are censored there.

## Replications and interpretation

Analyze each independent history before combining summaries. Group stratification
is available through `by_group`; arbitrary custom stratification is not supported.
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
sampling variability from persistent differences due to midpoint compensation,
initialization, or grouped ages. See [testing and validation](testing.md) and
the [accepted milestone 1 report](../validation/milestone1.qmd) and
[accepted milestone 2 report](../validation/milestone2.qmd).


When input fertility rates include `union_type`, `calculate_tfr()` returns a
separate row per status and period, preserving the status column. Each sum is
`rate * age-band width` over that status's finite contiguous age range. These
conditional TFRs describe continuous exposure to a status, not additive
components of the overall TFR. Tables without status retain their existing
behavior; default Lexis summaries use all female exposure, while
`by_union_type = TRUE` also returns status-specific exposure as described below.

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
input type. These are recorded person attributes, not inferred historical group
transitions; current group identity is fixed over life. A surviving partner's
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

`by_union_type = TRUE` returns adjacent rows labeled by `union_type`:
`unpartnered` plus each type discovered in the union histories. Typed histories
can contain `cohab` and `marriage`; legacy untyped histories use `partnered`.
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
counts <- calculate_demog_rates(result, seq(2020, 2100, 20), 0:100,
                                by_union_type = TRUE)
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
                                      0:100, by_group = TRUE,
                                      by_previous_union = TRUE)
```

By default, union counts pool previous categories and omit `previous_union`.
Set `by_previous_union = TRUE` as above to stratify independently of `by_group`.
The stratified output returns period, age, sex, `previous_union`, optional individual `group`,
`person_years`, `unions`, and `rate_union`. Categories are `none`,
`dissolved`, and `widowed`, determined by the most recent recorded episode.
Each spouse contributes a start in their category immediately before formation.
Partnered time contributes no exposure. Men's rates describe matching outcomes;
female search inputs can exceed realized rates when searches fail.

For the widowhood example above, the woman has 0.5 unpartnered years in `none`
and 1 in `widowed`, with one formation in each category. The `dissolved` row has
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
`calculate_demog_rates(by_union_type = TRUE)`, which retains current-status
exposure, births, deaths, fertility and mortality summaries.


Optional `dur_breaks` on `calculate_union_rates()` measures dissolved and
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


With typed histories, `by_union_type = TRUE` adds separate `cohab` and `marriage`
rows. The numerator counts entry from an unpartnered state; conversion is
excluded. Each type uses the same unpartnered exposure, so rates sum to the
pooled `rate_union`, but repeated exposure must not be summed across types.
Without this option the helper pools both entry types. Untyped histories support
pooled rates only. See [cohabitation](cohabitation.md).


## Marriage entry including conversion

Use `calculate_marriage_rates(result, time_breaks, age_breaks, by_group = FALSE, by_previous_union = FALSE)`
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


Dissolution summaries also accept `by_union_type = TRUE`, optionally combined
with `by_group = TRUE`. Separate cohabitation and marriage rows use their own
union-years and dissolution counts; summing counts and exposure across types
recovers pooled totals. Conversion ends cohabitation exposure without a
dissolution event. Marriage duration starts at conversion. Typed union histories
are required for splitting; the default pools types.
