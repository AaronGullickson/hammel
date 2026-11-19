# United States demographic rates around 2000

## Rebuild the bundled data

Run [usa_2000.R](usa_2000.R) from the package root with
`readxl`, `tibble`, `usethis`, and Poppler's `pdftotext` installed.

```sh
Rscript data-raw/usa_2000/usa_2000.R
```

## Mortality, 1999–2001

`usa_2000_mortality` contains 44 rows of death rates per person-year by sex
(`f`, `m`) in bands 0-1, 1-5, five-year bands from 5 to 100, and an open 100+
band (`age_end = Inf`). Columns `sex`, `age_start`, `age_end` and `rate`,
ready for `mortality()`.

### Source

Arias E, Curtin LR, Wei R, Anderson RN. *U.S. Decennial Life Tables for
1999-2001, United States Life Tables*. National Vital Statistics Reports
57(1), National Center for Health Statistics, August 5, 2008 (updated to
correct calculation errors in Tables 1-13). Tables 2 (males) and 3 (females).

- Published report:
  <https://www.cdc.gov/nchs/data/nvsr/nvsr57/nvsr57_01.pdf>
- `usa_life_table_1999_2001_males.xls`:
  <https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/57_01/Table02.xls>
- `usa_life_table_1999_2001_females.xls`:
  <https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/57_01/Table03.xls>

Retrieved 2026-11-09. SHA-256:

| File | SHA-256 |
| --- | --- |
| `usa_life_tables_1999_2001_nchs_2008.pdf` | `40f6c144feee4c9748c202e8d36b7605c48ca0bf09765c9cc08b93a966f96a57` |
| `usa_life_table_1999_2001_males.xls` | `fb3e42521fcab50fc5aaf6b8cb9fdd5c95c47236c28acd05a8cd0205e292681d` |
| `usa_life_table_1999_2001_females.xls` | `51c75fee618ef6dacd4c1ca243a58209a62b7604e9c97bd559ad4c734f6418f6` |

The tables use deaths to U.S. residents in 1999-2001 and population
estimates based on the 2000 census, consistent with the revised 2000
fertility denominators in this folder. Ages over 85 draw on Medicare data.
The spreadsheets carry more precision than the printed tables and are the
inputs. As a U.S. federal publication the report is in the public domain;
cite the source.

### Method

From the single-year rows (ages 0-109) of each table, a closed band's rate
is its life-table deaths divided by its person-years lived,
`sum(ndx) / sum(nLx)`. The open band's rate is `l100 / T100`, which keeps the
source's remaining life expectancy at 100. The script checks table titles,
column labels, ages and `lx` - `dx` consistency.

### Caveats

- **Constant hazards within bands.** Hammel treats each band's rate as a
  constant hazard, so its survival within a band differs slightly from the
  source's single-year table. `calculate_life_table()` gives e0 of 74.05
  (men) and 79.38 (women) against published 74.13 and 79.47; e65 is 16.03
  and 19.03 against 16.11 and 19.12.
- **Infancy.** The 0-1 rate is a constant hazard over the first year, so
  infant deaths are not concentrated in the first weeks. The source also has
  day-level infant rows (0-1, 1-7, 7-28 and 28-365 days) if a finer schedule
  is needed.
- **Period.** The rates pool three years centered on 2000, rather than
  calendar 2000 alone.

## Fertility in 2000 and marriage transitions in 2001

- `usa_2000_fertility` contains six all-races, all-marital-status rates from
  Table 1, in births per woman-year.
- `usa_2000_fertility_union` contains 12 marital-status rates from Tables 7
  and 8, with married mapped to `partnered` and unmarried to `unpartnered`.
- `usa_2000_marriage` contains 12 female marriage search rates with five-year
  age bands from 15–19 through 65–69 and an open 70+ band.

### Sources

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

Source URLs, retrieval time and SHA-256 hashes are recorded in
[sources.json](sources.json). The NCHS report
permits reproduction with source citation.

### Interpretation and limits

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
The user selected the approximation on
October 16, 2026: `rate = transitions_per_1000 / 1000`, interpreted directly as
per-person-year rates. For example, the female 25-29 estimate 137.760548112966
per 1,000 becomes approximately 0.137760548112966 per person-year. No
`-log(1-p)` transformation is applied. This is an explicit modeling approximation,
not a claim that the source measured person-time hazards. At constant hazard h,
the one-year search probability is `1-exp(-h)`, rather than h itself.
Observed marriage entries also differ from search attempts that may fail, so
supplying these female rates does not guarantee matching realized marriage rates.
Sampling uncertainty is substantial in some groups. The original spreadsheet
contains the 90% margins of error. No interpolation, smoothing or fertility
compensation is applied.

## Migration between census regions

`usa_2000_migration` contains ballpark annual migration rates between the
four census regions (`NE` Northeast, `MW` Midwest, `S` South, `W` West) by
age. Columns `origin`, `destination`, `age_start`, `age_end` and `rate` (moves
per person-year), ready for
`set_group_transition(model, rates, group = "region")` with a region
dimension using these labels. 12 directed pairs x 11 age bands; the last band
is open (`age_end = Inf`).

### Source

U.S. Census Bureau, *Geographical Mobility: March 1999 to March 2000*, Current
Population Reports P20-538 (Jason Schachter, 2001), from the March 2000
Current Population Survey. The report used for transcription is archived
in commit `c104bc1` as `usa_mobility_1999_2000_census_p20-538.pdf`
(SHA-256 `26635b340a99b26d3283f9c26dd97d314778e449a598048d7c346ae3ffbd9d62`).

- Figure 3: migrants between regions, March 1999 to March 2000 (thousands).
  The figure's labels do not survive text extraction, so its values were
  assigned to flows using Table D's regional in- and out-migrant totals; the
  assignment is the only one consistent with all eight totals and reproduces
  them exactly (the script checks this).
- Table B: percent moving to a different state, by age.
- Denominators: Census 2000 regional populations (NE 53.594, MW 64.393,
  S 100.237, W 63.198 million).

### Method

rate(origin, destination, age) = flow(origin, destination) / population(origin)
x m(age), where m(age) is the age group's interstate moving rate divided by
the all-ages interstate rate (3.1 percent). About 11.5 per 1,000 people moved
between regions that year. Annual proportions are used as rates per
person-year, which is close at these levels.

### Caveats

- **Age pattern.** Moves between regions by age are not published here, so
  the national interstate age pattern is applied to every pair of regions.
  This understates older-age moves to the South (retirement migration) and is
  only a ballpark.
- **First band.** The source's youngest group is ages 1-4; its multiplier is
  applied from age 0.
- **Sampling error.** One year of CPS data; individual flows carry sampling
  error. Only the South's net gain (+227,000) was reported as statistically
  significant.
- **Individuals, not households.** Hammel moves people individually (household
  moves are deferred, D114), so applying these rates to everyone, children
  included, separates children from parents and spouses from each other.
- Five-year census flows (Census 2000 special report CENSR-7) or microdata
  with age (the 2000 census sample or the March 2000 CPS) would give better
  pair-specific age patterns.
