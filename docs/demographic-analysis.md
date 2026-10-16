# Demographic analysis

Updated October 13, 2026. All outputs are tibbles; base data frames are valid
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
| `deaths`, `asmr` | Deaths and deaths/person-years |
| `births`, `asfr` | Maternal births and births/all female person-years |

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
TFR uses `asfr` or fertility `rate`; life tables use `mx`, `asmr`, or mortality
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
[the report template](../validation/milestone1.qmd).


When input fertility rates include `union_status`, `calculate_tfr()` returns a
separate row per status and period, preserving the status column. Each sum is
`rate * age-band width` over that status's finite contiguous age range. These
conditional TFRs describe continuous exposure to a status, not additive
components of the overall TFR. Tables without status retain their existing
behavior; Lexis summaries currently use all female exposure.

## Union exposure and formation counts

The optional extension retains the existing period/age/sex rows and
adds `person_years_unpartnered`, `person_years_partnered`, and `unions` when
`include_unions = TRUE`. It requires a population object containing union
histories; bare person tables cannot establish union status. The default is
FALSE and preserves existing output and computation.

Exposure is the overlap of recorded lifetimes, observation cells, and union
episodes, with ongoing episodes censored at the observation endpoint. Partnered
means covered by at least one recorded episode; overlapping episodes must not
double-count exposure. Unpartnered exposure is the remaining lifetime exposure,
including widowhood and time between unions. These are current partnership
statuses, not never-married/ever-married categories. Missing earlier histories
cannot establish true historical status; users must select windows with complete
union follow-up. Union records do not distinguish legal marriage.

Each recorded union start counts once for each spouse, at that spouse's attained
age and the formation time, using the usual half-open age and period cells.
Repeated unions count separately. Starts before the observation window contribute
exposure but no formation; starts at its upper endpoint are excluded.

Worked exposure example: a woman aged 20 at time 0 enters a union at 0.5,
becomes widowed at 1.5, and forms another union at 2.5. Across [0, 3), she
contributes 1.5 partnered and 1.5 unpartnered person-years and two formations.
Splitting at time 1 gives partnered exposure 0.5 and 1.0 in the two periods.
Each husband contributes his own age-specific exposure and formation count.

`births_partnered` and `births_unpartnered` split maternal births by recorded
status, with NA for men. A union covers `[start_time, end_time)`: births exactly
at formation are partnered, including births triggering on-demand formation.
Births exactly at ending are unpartnered unless another union covers that time.
Unknown mothers remain in `unassigned_births`; their status is not inferred.

Validation checks exact interval overlap, age/period boundaries, widowhood,
remarriage, ongoing and zero-duration episodes, empty histories, conservation
of total exposure and births, spouse-specific formation ages, and unchanged
default output. Rates require positive denominators: status-specific fertility
uses births in that status divided by female exposure in the same status;
formation rates among unpartnered people use `unions / person_years_unpartnered`.

```r
counts <- calculate_lexis_counts(result, seq(2020, 2100, 20), 0:100,
                                include_unions = TRUE)
counts <- counts |> dplyr::mutate(
  asfr_partnered = dplyr::if_else(person_years_partnered > 0,
    births_partnered / person_years_partnered, NA_real_),
  asfr_unpartnered = dplyr::if_else(person_years_unpartnered > 0,
    births_unpartnered / person_years_unpartnered, NA_real_),
  union_formation_rate = dplyr::if_else(person_years_unpartnered > 0,
    unions / person_years_unpartnered, NA_real_)
)
```

`calculate_tfr(counts)` continues to use overall `asfr`. To calculate a
status-specific TFR, select female rows and replace `asfr` with the desired
conditional rate before calling `calculate_tfr()`. The optional columns do not
change the overall mortality or fertility measures.
