# Demographic analysis

Updated October 19, 2026. All outputs are tibbles; base data frames are valid
inputs. Dates and age boundaries use years, and exposure uses person-years.
See the [scientific contract](scientific-contract.md) for model assumptions.

## Exposure, counts, and rates

For a population simulated from 2000 through 2100:

```r
summary <- result |>
  calculate_lexis_counts(
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

Births with unknown mothers cannot be assigned maternal age. Their period
counts are in `attr(summary, "unassigned_births")`, also a tibble. Initial
persons born within the observation window can appear there if parents are
unknown; a baseline is not complete historical birth follow-up.

## Total fertility rate

```r
input_tfr <- model$processes$fertility |> calculate_tfr()
output_tfr <- summary |> calculate_tfr()
```

`calculate_tfr()` selects female rows when sex is present, omits open age bands,
and sums ASFR times each finite band's actual width by period. There is no
`age_range` argument; subset the summary to the desired age range first.
Remaining bands must be contiguous. The output reports the finite range used.

Missing ASFR contributions are omitted with `na.rm = TRUE`; an all-missing set
sums to zero. That is the accepted calculation rule, not evidence of zero
fertility in unobserved ages. Inspect exposure and missing rates alongside TFR.
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

Each period/sex needs contiguous age bands from zero through an explicit open
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

## Population size and observation boundaries

```r
population_sizes <- result |>
  calculate_population_size(times = seq(2000, 2100, by = 10))
```

This returns `time` and `population`. Counts include births and exclude deaths
at the requested instant, preserving input order and repeated times.

Lexis exposure splits `[birth_time, death_time)` directly at period/age
boundaries; fractional years and unequal bands are supported. Population
objects must have a known endpoint and cannot be queried beyond it. With bare
persons tables, the caller must ensure follow-up through the requested endpoint.
Missing death times are censored there.

## Replications and interpretation

Analyze each independent history before combining summaries. There is no
`groups` argument or automatic custom stratification. Stacking person tables
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


When input fertility rates include `union_status`, `calculate_tfr()` returns a
separate row per status and period, preserving the status column. Each sum is
`rate * age-band width` over that status's finite contiguous age range. These
conditional TFRs describe continuous exposure to a status, not additive
components of the overall TFR. Tables without status retain their existing
behavior; default Lexis summaries use all female exposure, while
`include_unions = TRUE` also returns status-specific exposure as described below.

## Union exposure and formation counts

`include_unions = TRUE` returns three adjacent rows per period/age/sex cell,
identified by `union_status`: `overall`, `partnered`, and `unpartnered`. All use
`person_years`, `deaths`, `births`, `rate_mortality`, and `rate_fertility`. Births and ASFR are NA
for men. The option also adds `unions` and `rate_marriage`, populated only for
`unpartnered`; both are NA on `overall` and `partnered` rows. The default FALSE
preserves the original output and computation. A population object with union
histories is required; bare person tables cannot establish status.

Exposure is the overlap of recorded lifetimes, observation cells, and union
episodes, with ongoing episodes censored at the observation endpoint. Partnered
means covered by at least one recorded episode; overlapping episodes do not
double-count exposure. Unpartnered exposure is the remaining lifetime exposure,
including widowhood and time between unions. These are current partnership
statuses, not never-married/ever-married categories. Missing earlier histories
cannot establish true historical status; select windows with complete union
follow-up. Union records do not distinguish legal marriage.

Each recorded union start counts once for each spouse at their own attained
age and the formation time, using half-open age and period cells. Repeated
unions count separately. Starts before the window contribute exposure but no
formation; starts at its upper endpoint are excluded. All formation counts are
reported on the unpartnered row, representing entry into a union.

Birth status uses `[start_time, end_time)`: births at formation are partnered,
including those triggering on-demand formation. Births at ending are unpartnered
unless another union covers that time. Unknown mothers remain in the
`unassigned_births` attribute; their status is not inferred. That attribute
remains one count per period, not repeated by status.

Deaths use status immediately before death: a union covers the death when
`start_time < death_time <= end_time`. Thus a death ending a union is partnered;
a later death following widowhood is unpartnered. A union starting exactly at
death or a zero-duration union contributes no preceding partnered time. This
observation convention uses the left limit even for coincident recorded events;
it does not reconstruct their within-time execution order. Event period and age
classification still use the death time and usual half-open cells.

Worked example: a woman aged 20 at time 0 enters a union at 0.5, becomes
widowed at 1.5, and forms another union at 2.5. Across [0,3), her overall row
has 3 person-years; partnered and unpartnered rows each have 1.5. Two union
formations on the unpartnered row give `rate_marriage = 2 / 1.5`. Births at
0.25 and 1.5 are unpartnered; births at 0.5 and 2.5 are partnered, giving
conditional ASFR 2 / 1.5 in each status and overall ASFR 4 / 3. Her first
husband's death at 1.5 is partnered. Each husband has his own age-specific
exposure and formation count.

All rates use each row's `person_years`: `rate_mortality = deaths / person_years`,
`rate_fertility = births / person_years`, and, on unpartnered rows only,
`rate_marriage = unions / person_years`. Rates are per person-year, with ASFR
per woman-year, and are NA when exposure is zero. Positive exposure with zero
events gives rate zero. Exposure includes time excluded by birth spacing.

```r
counts <- calculate_lexis_counts(result, seq(2020, 2100, 20), 0:100,
                                include_unions = TRUE)
# Existing TFR grouping reads union_status without a new option:
tfr <- calculate_tfr(counts)
```

`calculate_tfr()` already groups by `union_status`, yielding overall, partnered,
and unpartnered TFRs per period. Conditional TFRs describe continuous exposure
to the selected status, and do not add to the overall TFR. Existing missing-ASFR
summation rules still apply. Life tables also keep statuses separate; conditional
mortality life tables do not model transitions between partnership statuses.

The overall row duplicates the combined status rows. Do not sum all three:
partnered plus unpartnered exposure, deaths, and maternal births equal overall.
Validation covers exact overlap, boundaries, widowhood/remarriage, empty and
overlapping histories, zero-duration episodes, count conservation, rate
formulas, downstream grouping, and unchanged default output.
