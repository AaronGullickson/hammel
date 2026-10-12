# Extracting demographic measures

Load the development package with `devtools::load_all()`. `calculate_lexis_counts()` accepts a `hammel_population` such as `result`, or a persons data frame representing one independent population history. Dates and age boundaries are in years; exposure is in person-years.

All table outputs are tibbles, including population sizes, Lexis counts, TFR,
life tables, and auxiliary unassigned-birth counts. Population histories and
process rate tables are also tibbles; base data frames are accepted as input.

## One summary for fertility and mortality

For the USA/2000 example simulated from 2000 to 2100:

```r
summary <- result |>
  calculate_lexis_counts(
    time_breaks = seq(2025, 2100, by = 25),
    age_breaks = c(0, 1, 5, 10, 13:55, seq(60, 100, by = 5))
  )

tfr <- summary |> calculate_tfr()
life_tables <- summary |> calculate_life_table()

life_tables[life_tables$age_start == 0,
            c("time_start", "time_end", "sex", "ex")]
```

The summary contains `time_start`, `time_end`, `sex`, `age_start`, `age_end`, `person_years`, `deaths`, `births`, `asmr`, and `asfr`. ASMR is deaths divided by person-years. Births are classified by maternal age and ASFR uses all female exposure, including time temporarily excluded by birth spacing. Births and ASFR are NA for men.

TFR automatically selects female rows, omits open-ended age bands, and sums ASFR times each remaining finite age-band width separately for each period. No age-range argument is needed. The finite age range reported in the output tells you which bands contributed. Missing ASFRs are omitted from the TFR sum using `na.rm = TRUE`. You can subset the summary first if you want a particular finite age range. HFD open endpoint rows remain omitted by our rate adapter.

Life tables are constructed separately for each period and sex. The age schedule must begin at zero, be contiguous, and end at Inf. `calculate_lexis_counts()` automatically appends Inf if the supplied age breaks end at a finite value: breaks ending at 100 create a final `[100, Inf)` band. Manually supplied rate tables must still specify their actual open-ended interval. Columns include `n`, `ax`, `mx`, `qx`, `lx`, `lx_end`, `dx`, `Lx`, `Tx`, and `ex`. Life expectancy is Tx/lx. The method uses piecewise constant hazards and corresponding integrated life-table person-years; it does not apply conventional HMD infant adjustments. See the [HMD Methods Protocol](https://www.mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf) for that database's methodology.

Input-derived reference tables are also supported:

```r
input_tfr <- model$processes$fertility |> calculate_tfr()
input_life_tables <- model$processes$mortality |>
  calculate_life_table()
```

Both functions also accept their process's `$rates` table. TFR uses `asfr` from
Lexis counts or `rate` from a fertility input table; life tables use `asmr`,
`mx`, or mortality `rate`. Typed mortality processes are rejected by TFR,
and typed fertility processes are rejected by the life-table function. A bare
table's `rate` column must represent the corresponding demographic process.

## Observation rules

Calculate total living population at selected time points:

```r
population_sizes <- result |>
  calculate_population_size(times = seq(2000, 2100, by = 10))
```

The result has `time` and `population` columns. Births at a requested time are
included and deaths at that time are excluded. Requested order and repeated
times are preserved. The function also accepts a persons table for one history;
as with Lexis counts, the caller must ensure complete follow-up when supplying
a bare table.

Intervals are half-open. Exposure is recorded lifespan overlap, split directly at age and time boundaries, supporting fractional periods and unequal age bands. A population's known endpoint limits observation. For a bare persons table, ensure complete follow-up through the requested endpoint yourself.

Births with unknown mothers cannot be classified by maternal age; inspect `attr(summary, "unassigned_births")` for period counts. Initial persons with unknown parents may appear there if their birth dates lie inside the observation window. Use periods with complete birth histories for fertility estimates.

Empty cells remain in the summary. Zero exposure gives NA rates. TFR omits missing ASFRs from its sum; missing mortality rates still propagate into life tables. A zero terminal mortality rate implies infinite remaining life expectancy under the hazard model. Wider bands or longer periods can help with sparse data, but missing estimates should not silently become zero.

Analyze each independent simulation separately. There is no `groups` parameter at this stage. Stacked histories with duplicated IDs must be split before calling the summary function so parent references remain unambiguous. Extra person attributes are not automatically treated as analysis strata.

These functions support formal validation rather than establish it. Step-start event recording, birth-spacing compensation, initialization, and sampling uncertainty still affect comparisons with inputs.
