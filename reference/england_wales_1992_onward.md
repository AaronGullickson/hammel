# England and Wales inputs for a second demographic transition example

A named list of six tibbles used by the second demographic transition
article vignette. Fertility and marriage tables have Hammel's schedule
column names; cohabitation and marital divorce/death observations retain
source values for the article's transformation examples. Cohabitation
cohort outcomes cover Great Britain; the other tables cover England and
Wales.

## Usage

``` r
england_wales_1992_onward
```

## Format

A named list containing:

- fertility:

  18 rows: `union_type`, `age_start`, `age_end`, `rate`.

- marriage:

  672 rows: `time_start`, `time_end`, `age_start`, `age_end`,
  `previous_union`, `rate`.

- premarital_cohabitation:

  373 rows: `year`, `age_groups`, `prop_cohabit`.

- cohabitation_outcomes:

  3 rows: `cohort_start`, `cohort_end`, `marriage_5`, `separation_5`,
  `continuing_5`, `marriage_10`, `separation_10`, `continuing_10`.

- marital_divorce_death_2005:

  52 rows: `year`, `anniversary`, `cumulative_divorce`,
  `cumulative_death`.

- marital_divorce_death_2010:

  60 rows, with the same columns as 2005.

Ages, durations and calendar times are in years. Rates are annual;
proportions are on the zero-to-one scale. Source cohorts are half-open
intervals, not calendar periods of a rate schedule.

## Source

Office for National Statistics (ONS). Contains public sector information
licensed under the Open Government Licence v3.0:
<https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/>

Fertility: O'Leary L, Natamba E, Jefferies J, Wilson B (2010). Fertility
and partnership status in the last two decades. Population Trends 140,
Table B1. [doi:10.1057/pt.2010.10](https://doi.org/10.1057/pt.2010.10) .

Marriage and premarital cohabitation: ONS Marriages in England and
Wales, 2019 workbook; female marital population estimates, 1971–2001 and
2002–2019. Exact source URLs, retrieval dates and checksums are recorded
in the repository's source-folder `sources.csv`.

Cohabitation outcomes: Beaujouan E, Ni Bhrolchain M (2011). Cohabitation
and marriage in Britain since the 1970s. Population Trends 145, Table 3.
[doi:10.1057/pt.2011.16](https://doi.org/10.1057/pt.2011.16) .

Marital divorce/death, 2005: Wilson B, Smallwood S (2008). The
proportion of marriages ending in divorce. Population Trends 131, Table
3; cumulative proportions reconstructed from life-table event counts.
Marital divorce/death, 2010: ONS Divorces in England and Wales, 2012
bulletin, Figure 5. Publication tables and PDFs are not bundled with
this dataset.

## Details

Fertility rates are births per woman-year in 1992 by `union_type`:
`"marriage"`, `"cohab"` and `"unpartnered"`. Age bounds use the source's
denominator ages, 16–45.

Marriage rates are annual central rates of first marriages and
remarriages for 1992–2019. `previous_union` maps these to `"none"` and
`"dissolved"`; this describes legal marital history, not prior
cohabitation. Calendar intervals are `[time_start, time_end)`. The 48
unavailable rates remain `NA`. These marriage observations are not
direct measurements of union-entry hazards including cohabitation.

Premarital cohabitation proportions condition on marrying; Cohabitation
outcomes are cumulative proportions at five and ten years for women
under 45 at cohabitation start. Source rounding and missing ten-year
values for the latest cohort are retained.

Marital divorce/death tables are cumulative competing-risk proportions,
not annual hazards. No hazard conversion, normalization, extrapolation
or interpolation is applied to these tables.

Original numerical sources and their provenance remain in
`validation/data/england_wales_1992_onward/` in the repository. The
offline preparation script is `data-raw/england_wales_1992_onward.R`.

## Examples

``` r
england_wales_1992_onward$fertility
#> # A tibble: 18 × 4
#>    union_type  age_start age_end   rate
#>    <chr>           <int>   <int>  <dbl>
#>  1 marriage           16      20 0.287 
#>  2 marriage           20      25 0.196 
#>  3 marriage           25      30 0.170 
#>  4 marriage           30      35 0.103 
#>  5 marriage           35      40 0.0355
#>  6 marriage           40      45 0.0057
#>  7 cohab              16      20 0.192 
#>  8 cohab              20      25 0.128 
#>  9 cohab              25      30 0.108 
#> 10 cohab              30      35 0.101 
#> 11 cohab              35      40 0.0623
#> 12 cohab              40      45 0.0166
#> 13 unpartnered        16      20 0.022 
#> 14 unpartnered        20      25 0.031 
#> 15 unpartnered        25      30 0.0327
#> 16 unpartnered        30      35 0.0263
#> 17 unpartnered        35      40 0.0133
#> 18 unpartnered        40      45 0.0029
fertility(england_wales_1992_onward$fertility)
#> $name
#> [1] "fertility"
#> 
#> $rates
#> # A tibble: 18 × 4
#>    union_type  age_start age_end   rate
#>    <chr>           <int>   <int>  <dbl>
#>  1 marriage           16      20 0.287 
#>  2 marriage           20      25 0.196 
#>  3 marriage           25      30 0.170 
#>  4 marriage           30      35 0.103 
#>  5 marriage           35      40 0.0355
#>  6 marriage           40      45 0.0057
#>  7 cohab              16      20 0.192 
#>  8 cohab              20      25 0.128 
#>  9 cohab              25      30 0.108 
#> 10 cohab              30      35 0.101 
#> 11 cohab              35      40 0.0623
#> 12 cohab              40      45 0.0166
#> 13 unpartnered        16      20 0.022 
#> 14 unpartnered        20      25 0.031 
#> 15 unpartnered        25      30 0.0327
#> 16 unpartnered        30      35 0.0263
#> 17 unpartnered        35      40 0.0133
#> 18 unpartnered        40      45 0.0029
#> 
#> $birth_interval
#> [1] 0.75
#> 
#> $ffab
#> [1] 0.4886
#> 
#> attr(,"class")
#> [1] "hammel_fertility" "hammel_process"  
england_wales_1992_onward$cohabitation_outcomes
#> # A tibble: 3 × 8
#>   cohort_start cohort_end marriage_5 separation_5 continuing_5 marriage_10
#>          <int>      <int>      <dbl>        <dbl>        <dbl>       <dbl>
#> 1         1990       1995      0.463        0.232        0.305       0.561
#> 2         1995       2000      0.399        0.305        0.296       0.494
#> 3         2000       2005      0.395        0.354        0.252      NA    
#> # ℹ 2 more variables: separation_10 <dbl>, continuing_10 <dbl>
```
