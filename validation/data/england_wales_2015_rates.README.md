# England and Wales, 2015

Source: Office for National Statistics. Contains public sector information
licensed under the Open Government Licence v3.0.

- `england_wales_2015_fertility_by_union_status.csv`: 14 rows, five-year fertility
  schedules for `partnered` and `unpartnered`, ready for `set_fertility()`.
- `england_wales_2015_overall_fertility.csv`: 7 rows, overall fertility alternative.
- `england_wales_2015_female_union_rates.csv`: 11 rows, female marriage schedule,
  original bands (16-19, five-year bands 20-64, and 65+), ready for rate initiation.
- `england_wales_fertility_ons_table6.csv`: tidy archival Table 6, 1999-2015,
  including original birth counts, source age labels, per-1,000 rates and
  converted per-woman-year rates. Includes `total` rows: select a year and either
  total or both conditional schedules before passing data to Hammel.

All model-ready CSV rates are per woman-year in the corresponding status.
Original fertility rates are **per 1,000**, not per 100,000. No interpolation
is used. Partnered means married/civil partnered; unpartnered means legally
unmarried, including cohabiting women. This is a legal-marriage example.

The source's under-20 and 45+ fertility rates use women aged 15-19 and 45-49
as denominators. The CSV uses those finite exposure bands, but the source
birth counts include births below 15 and above 49 respectively. Marriage risk
starts at 16 in its source. Small numbers of births to same-sex couples are
included in fertility, while marriage rates here refer to opposite-sex couples.
Without divorce, the simulation need not reproduce observed status proportions.

```r
fertility_rates <- read.csv(
  "validation/data/england_wales_2015_fertility_by_union_status.csv"
)
union_rates <- read.csv(
  "validation/data/england_wales_2015_female_union_rates.csv"
)
model <- model |>
  set_fertility(fertility_rates) |>
  set_union_formation("rate", matching, rates = union_rates)
```

Rebuild offline from the package root with
`Rscript validation/data/prepare_england_wales_2015_rates.R` (requires `readxl`).
The archived XLS and marriage chart CSV are inputs. Source links, hashes,
transformation details and checks are in `england_wales_2015_rates.source.json`.
