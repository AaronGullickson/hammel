# US marital-rate endpoints: 1980 and approximately 2020

Prepared November 1, 2026 for the demographic-transition demonstration. Female
rates are supplied because Hammel carries union searches and dissolution on the
woman. These are published empirical endpoints, not new demographic defaults.
No interpolation, rebinning, missing-age imputation or simulation is performed.

| Measure | Early endpoint | Later endpoint | Denominator within age band |
| --- | --- | --- | --- |
| First marriage | 1980, NCHS Table 4 | 2017, NCFMR Figure 1 | Never-married women |
| Remarriage | 1980, NCHS Table 4 | 2018, NCFMR Figure 3 | Divorced or widowed women, pooled |
| Divorce | 1980, NCHS Table 3 | 2019, NCFMR Figure 2, women | Married women |

## Sources and provenance

- NCHS, *Advance Report of Final Marriage Statistics, 1981*, Monthly Vital
  Statistics Report 32(11), supplement, February 29, 1984. Table 4, printed page 8
  (PDF page 8), contains 1980 and 1981 side by side. Extract the **women's 1980
  rates**, excluding totals and the overlapping 25–44 subtotal. The report is
  archived as `usa_marriage_1980_1981_nchs.pdf`.
  [Official PDF](https://www.cdc.gov/nchs/data/mvsr/supp/mv32_11s.pdf).
- Wilson, B. F., NCHS, *Teenage Marriage and Divorce: United States, 1970–81*,
  Vital and Health Statistics 21(43), September 1985. Despite its title, Table 3,
  printed page 19 (PDF page 24), supplies estimated divorce rates for **all adult
  ages**. Extract the wife/1980 column. The report is archived as
  `usa_divorce_1970_1981_nchs.pdf`.
  [Official PDF](https://stacks.cdc.gov/view/cdc/12848/cdc_12848_DS1.pdf).
- Eickmeyer, K. J. (2019), *Age Variation in the First Marriage Rate, 1990 &
  2017*, NCFMR FP-19-05, Figure 1: 2017 series, excluding Total.
  [DOI](https://doi.org/10.25035/ncfmr/fp-19-05).
- Schweizer, V. (2020), *Age Variation in the Remarriage Rate, 1990 & 2018*,
  NCFMR FP-20-13, Figure 3: women, 2018 series.
  [DOI](https://doi.org/10.25035/ncfmr/fp-20-13).
- Carlson, L. (2021), *Age Variation in the Divorce Rate, 1990 & 2019*, NCFMR
  FP-21-16, Figure 2: women, 2019 series, not the combined-sex Figure 1.
  [DOI](https://doi.org/10.25035/ncfmr/fp-21-16).

Later estimates are NCFMR tabulations of ACS one-year data. Printed labels were
visually verified; no values were estimated from bar heights. Only numerical
facts are transcribed from NCFMR, with attribution. Its figures and prose are
not redistributed and no open licence is claimed. Federal NCHS reports are
public-domain government publications; the divorce report explicitly states
permission to reproduce with attribution appreciated. No upstream code is used.
Source URLs, downloaded-source hashes and output hashes are recorded in
`usa_1980_2019_marital_rates.source.json`.

## Files and units

Six individual CSVs have `age_start`, `age_end`, `rate`:

- `usa_1980_first_marriage_rates_female.csv`
- `usa_1980_remarriage_rates_female.csv`
- `usa_1980_divorce_rates_female.csv`
- `usa_2017_first_marriage_rates_female.csv`
- `usa_2018_remarriage_rates_female.csv`
- `usa_2019_divorce_rates_female.csv`

`usa_1980_2019_marital_rates.csv` combines these with measure, actual source year,
sex, original age labels, denominator, source ID and original rate per 1,000.
Do not pass the combined table directly as one process schedule: select a
measure/year first. Source periods are not shifted to nominal 1980/2020 dates.
All tables can be read with `readr::read_csv()` or
`tibble::as_tibble(read.csv(...))`.

Each rate is the published annual events per 1,000 women in the relevant risk
population, divided by 1,000. Age ends are exclusive, including Inf for the
source's open final band. For example, the 1980 first-marriage rate at 20–24 is
119.8 per 1,000 never-married women; its CSV rate is 0.1198. This is a central
population-based rate approximation, not an exactly measured person-time hazard.
No `-log(1-p)` conversion is applied.

## Coverage and interpretation

1980 marriage rates refer to the marriage-registration area (MRA), with all
rates excluding Iowa. They are not full-country registration rates. The
1980 remarriage rate for 15–24 is a single shared source rate (231.0 per 1,000),
although the source shows separate marriage counts for 15–19 and 20–24. It must
not be read as a rate available separately for each band. This pooled remarriage
series includes widowed and divorced women; separate status rows are not used.
The youngest marriage numerator includes under-15 marriages, but its denominator
excludes under-15s, as stated in Table 4.

1980 divorce rates are national **estimates**: NCHS applies the DRA's age
composition to the national divorce total. Its appendix explains this procedure
and the population denominators. The underlying registration area does not
cover all states. Table 3's overall female rate is 22.5 per 1,000; do not mix it
with the differently estimated 22.6 in the earlier final-divorce report.

The later first-marriage schedule ends with a single 50+ band. Later remarriage
starts at age 20; ages below 20 are unavailable here, not zero. Native source
bands differ among measures and endpoints. No unsupported age detail is added.
The later ACS estimates have survey uncertainty not supplied by these charts;
no SEs or confidence bounds are invented. Sources and denominators differ,
so the endpoints are suitable as illustrative anchors, not a fully harmonized
historical series. They pool marriage orders for divorce and include remarriages
of any order for remarriage.

With cohabitation enabled, legal first marriage is not the same as first union;
legal remarriage is not the same as repartnering. Source marriage events include
marriages following cohabitation. The tables do not separately identify direct
marriage, cohabitation entry, or conversion hazards. In particular, do not treat
legal never-married status as automatically equivalent to Hammel's previous-union
`none` category. A mapping would be an explicit demonstration assumption.

The divorce tables use attained **age**, while current Hammel dissolution
schedules use **union duration**. These divorce CSVs are reference data, not
immediately usable dissolution-process tables. No age-based dissolution feature
or age-to-duration mapping is introduced by this preparation.

The observed divorce trend differs by age: younger-age rates decline while
older-age rates rise. Using these endpoints need not implement a uniform
increase in dissolution. Cohabitation estimates remain separate future work.

## Rebuild and checks

Run `Rscript validation/data/usa_1980_2019/prepare_usa_1980_2019_marital_rates.R` from the
repository root. It rebuilds the CSVs offline from the attributed printed-value
transcription in the script, checks nonoverlap/contiguity over each supplied age
range, and verifies units and endpoint cells. Source PDFs are retained for
inspection, not parsed by fragile OCR. The preparation script downloads nothing.
Neither the engine nor the user's milestone reports are modified.
