# United States: mortality, 1999-2001

`usa_1999_2001_mortality.csv`: 44 rows of death rates per person-year by sex
(`f`, `m`) in bands 0-1, 1-5, five-year bands from 5 to 100, and an open 100+
band (`age_end = Inf`). Columns `sex`, `age_start`, `age_end` and `rate`,
ready for `mortality()`. Bundled in the package as `usa_2000_mortality`.
Rebuild with `prepare_usa_1999_2001_mortality.R`; requires `readxl`.

## Source

Arias E, Curtin LR, Wei R, Anderson RN. *U.S. Decennial Life Tables for
1999-2001, United States Life Tables*. National Vital Statistics Reports
57(1), National Center for Health Statistics, August 5, 2008 (updated to
correct calculation errors in Tables 1-13). Tables 2 (males) and 3 (females).

- `usa_life_tables_1999_2001_nchs_2008.pdf`:
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

## Method

From the single-year rows (ages 0-109) of each table, a closed band's rate
is its life-table deaths divided by its person-years lived,
`sum(ndx) / sum(nLx)`. The open band's rate is `l100 / T100`, which keeps the
source's remaining life expectancy at 100. The script checks table titles,
column labels, ages and `lx` - `dx` consistency.

## Caveats

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
