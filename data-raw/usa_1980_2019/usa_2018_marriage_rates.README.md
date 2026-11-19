# US female marriage rates, 2018

Prepared November 1, 2026 for the second demographic-transition demonstration.
No interpolation is included; the existing 2001 inputs are unchanged.

Source: Schweizer, V. (2019), *Marriage to Divorce Ratio in the U.S.: Demographic
Variation, 2018*, NCFMR Family Profiles FP-19-27, Figure 3.
[Report](https://www.bgsu.edu/ncfmr/resources/data/family-profiles/schweizer-marriage-divorce-ratio-demo-variation-fp-19-27),
[DOI](https://doi.org/10.25035/ncfmr/fp-19-27). The underlying estimates are
NCFMR analyses of the 2018 ACS one-year data.

| Age | Marriages per 1,000 unmarried women | CSV rate |
| --- | ---: | ---: |
| 15–24 | 23.0 | 0.0230 |
| 25–34 | 75.3 | 0.0753 |
| 35–44 | 45.9 | 0.0459 |
| 45–54 | 27.9 | 0.0279 |
| 55–64 | 12.1 | 0.0121 |
| 65+ | 3.3 | 0.0033 |

`usa_2018_marriage_rates_female.csv` contains `age_start`, `age_end`, and `rate`.
The six values were transcribed from printed marriage-rate labels (dark green),
not the orange divorce-rate labels or the grey marriage-to-divorce ratios.
The conversion divides by 1,000 and uses exclusive age ends; 65+ ends at Inf.
The CSV can be read as a tibble with `readr::read_csv()` or
`tibble::as_tibble(read.csv(...))`. No schedule is inferred below age 15.

## Comparison with the earlier schedule

The [existing source notes](usa_2000_rates.README.md) describe the approximate
2000 baseline, actually labelled 2001 because that is its SIPP survey year.
These 2018 rates pool first marriages and remarriages, as does that baseline.
They do not separate previously divorced and widowed women.

The 2018 bands are broader. The earlier five-year rates can be pooled using
their `unmarried_at_start_thousands` weights in
`usa_2001_marriage_transitions.csv`; a simple unweighted average would not
preserve the pooled source transition proportion. Neither pooling nor
interpolation has been applied here.

The denominators differ: the old table uses people unmarried at the start of
the period, whereas the later refined rate uses unmarried women at interview.
The estimates therefore provide illustrative endpoints, not a fully harmonized
historical hazard series. Neither denominator is exact unmarried person-years.
Interpreting the divided values directly as simulation hazards would use the
same kind of approximation previously chosen for the old inputs; this file
does not change the engine or establish a new demographic default.

With cohabitation enabled, observed marriages include conversion as well as
direct entry. These rates are empirical marriage targets, not independently
estimated union-search or conversion hazards. Applying them to total union
search would be an additional demonstration assumption.

Only factual numerical values are reproduced with attribution; the source
figure is not redistributed and no open licence is claimed. Source URLs,
retrieval timestamp, chart hash, output hash and original rounded values are
recorded in `usa_2018_marriage_rates.source.json`.
