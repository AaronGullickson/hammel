# United States: fertility 2000, marriage transitions 2001

- `usa_2000_fertility_overall.csv`: six all-races, all-marital-status rates
  from Table 1, ages 15-19 through 40-44, with `age_start`, `age_end`, and
  births per woman-year (`rate`), ready for `set_fertility()`. Omits published
  10-14 and 45-49 groups to match the conditional input's age boundaries.
- `usa_2000_fertility_by_marital_status.csv`: 12 rows, all races combined,
  original five-year bands, legal marital status, published rates per 1,000
  women and converted births per woman-year (`rate`).
- `usa_2000_fertility_by_union_status.csv`: the same fertility values in the
  England/Wales input format for `set_fertility()`, with married mapped to
  `partnered` and unmarried to `unpartnered`.
- `usa_2001_marriage_transitions.csv`: 24 rows, both sexes, original age bands,
  starting unmarried population (thousands), transitions per 1,000, 90% margin
  of error per 1,000, and annual transition proportion. Full spreadsheet
  precision is retained. `m` means men; `f` means women. The last band is 70+.

- `usa_2001_marriage_rates_female.csv`: 12 rows with `age_start`, `age_end`,
  and `rate`, ready for female rate-based initiation in `set_union_formation()`.
- `usa_2001_marriage_rates_male.csv`: the same three columns for men, retained
  for comparison with simulated male marriage rates. The current engine only
  schedules female searches; do not pass the male file as a female schedule.
  Both files preserve five-year bands from 15-19 through 65-69 and the 70+
  open band (`age_end = Inf`). They use the direct-rate approximation below.

## Sources

Fertility: NCHS, *Revised Birth and Fertility Rates for the 1990s and New Rates
for Hispanic Populations, 2000 and 2001: United States*, NVSR 51(12), August 4,
2003, revised June 2004. Table 1, printed page 13 (PDF page 13), supplies
overall rates. Tables 7 and 8, printed pages 54 and 56 (PDF pages
54 and 56), supply marital-status rates. Archived as `usa_fertility_revised_nchs_2003.pdf`.
[Official PDF](https://www.cdc.gov/nchs/data/nvsr/nvsr51/nvsr51_12.pdf).
These are revised rates based on population estimates consistent with the
2000 Census, rather than the original final-2000 report's denominators.

Marriage: Census Bureau, *Number, Timing, and Duration of Marriages and
Divorces: 2001*, P70-97, February 2005, supplemental Table 2 (Internet release
February 10, 2005), SIPP 2001 Panel, Wave 2 Topical Module. Archived as
`usa_marriage_transitions_2001_census.xls`.
[Report page](https://www.census.gov/library/publications/2005/demo/p70-097.html)
and [official spreadsheet](https://www2.census.gov/programs-surveys/demo/tables/marriage-and-divorce/2001/tab02.xls).
Events occurred during the year before the survey, rather than precisely
calendar 2001. Marriage entry includes first marriages and remarriages.

## Interpretation and limits

Overall fertility denominators are women in each age group; Table 1 specifies
April 1 enumeration for 2000. Marital-status rates use the corresponding
age/status denominators (the table headers specify July 1 estimates). Unmarried includes never married, widowed and divorced women. Legal
marital status does not distinguish cohabitation, so its mapping to simulation
union status is an illustrative legal-marriage scenario. Rates count live
births, not deliveries. The source's 40-44 rates include births to women
**40 and older**, divided by women aged 40-44. The finite output band [40,45)
uses that published rate; no separate 45+ schedule is inferred. No rates are
inferred below 15. Age ends are exclusive.

The overall Table 1 rates per 1,000 are 47.7, 109.7, 113.5, 91.2, 39.7,
and 8.0, respectively. Its 40-44 numerator covers births at ages 40-44;
Tables 7 and 8 instead use births at 40+ divided by women 40-44. Matching
output age bands therefore does not make that final cell exactly comparable.
No alternative overall 40+ rate is inferred from rounded rates or from
unavailable marital-status population weights. The omitted endpoint rates
are 0.9 (10-14) and 0.5 (45-49) per 1,000. The retained schedule integrates
to 2.049 births per woman across [15,45), rather than the report's full-range
TFR of 2.056. Omitting endpoint groups is a truncation, not an estimate of
zero fertility at those ages.

Marriage denominators are people unmarried at the beginning of the period,
not person-years of unmarried exposure. `annual_transition_proportion` is
the published estimate divided by 1,000; it is **not a continuous-time hazard**.
For the two new model-format files, the user selected the approximation on
October 16, 2026: `rate = transitions_per_1000 / 1000`, interpreted directly as
per-person-year rates. For example, the female 25-29 estimate 137.760548112966
per 1,000 becomes approximately 0.137760548112966 per person-year. No
`-log(1-p)` transformation is applied. This is an explicit modeling approximation,
not a claim that the source measured person-time hazards. At constant hazard h,
the one-year search probability is `1-exp(-h)`, rather than h itself.
Observed marriage entries also differ from search attempts that may fail, so
supplying these female rates does not guarantee matching realized marriage rates.
Sampling uncertainty is substantial in some groups; the original combined CSV
retains the source's 90% margins of error without clipping intervals. All source
values and the combined transitions CSV are retained unchanged. No interpolation,
smoothing or fertility compensation is applied.

```r
marriage_rates <- read.csv("validation/data/usa_2001_marriage_rates_female.csv")
model <- set_union_formation(model, "rate", matching, rates = marriage_rates)
```

Rebuild offline from the package root:

```sh
Rscript validation/data/prepare_usa_2000_rates.R
```

Requires `readxl` and Poppler's `pdftotext`. The script checks source labels,
row counts, unique marriage cells, finite values and the extracted 2000
fertility values. Source URLs, retrieval time and SHA-256 hashes are recorded
in `usa_2000_rates.source.json`. These are U.S. federal statistical sources;
the NCHS report explicitly permits reproduction with source citation.
