# Document the example demographic datasets bundled with the package.

# DATASETS -----------------------------------------------------------------

#' United States example rates, around 2000
#'
#' Published U.S. demographic rates from around 2000, prepared as input tables
#' for Hammel's process constructors. They give realistic examples without
#' downloading rates. All rates are per person-year, and ages are in years,
#' with half-open bands `[age_start, age_end)`. Each source is a U.S. federal
#' statistical publication. Source files, preparation scripts and further
#' notes are in the package repository under
#' `validation/data/usa_2000_2001/`; `data-raw/usa_2000.R` builds the
#' datasets.
#'
#' @section Mortality:
#' `usa_2000_mortality` gives death rates by sex for ages 0-1, 1-5, five-year
#' bands to 100, and an open 100+ band, for use with [mortality()]. Rates are
#' abridged from the single-year U.S. decennial life tables for 1999-2001 as
#' life-table deaths divided by person-years lived in each band (`l100 /
#' T100` for 100+). Under Hammel's constant hazards within bands,
#' [calculate_life_table()] gives life expectancy at birth of 74.05 years for
#' men and 79.38 for women, against the source's 74.13 and 79.47.
#'
#' @section Fertility:
#' `usa_2000_fertility` gives all-women birth rates for 2000 in five-year
#' bands from 15 to 45, for use with [fertility()]. The published 10-14 and
#' 45-49 rates are omitted, so the schedule sums to 2.049 births per woman
#' rather than the published total fertility rate of 2.056.
#'
#' `usa_2000_fertility_union` gives the same bands by `union_type`, mapping
#' legally married women to `"partnered"` and unmarried women to
#' `"unpartnered"`. Legal marital status does not distinguish cohabitation.
#' The 40-45 rates divide births to women 40 and older by women 40-44.
#'
#' @section Marriage:
#' `usa_2000_marriage` gives female marriage rates in five-year bands from 15
#' to 70 and an open 70+ band, for rate-based initiation in
#' [union_formation()]. They are the published annual proportions of
#' unmarried women who married in the year before a 2001 survey, used
#' directly as search rates per person-year without conversion to hazards.
#' Marriages include first marriages and remarriages. Search attempts can fail
#' to find a partner, so simulated marriage rates need not match these.
#'
#' @section Migration:
#' `usa_2000_migration` gives approximate migration rates between the four
#' census regions (`"NE"` Northeast, `"MW"` Midwest, `"S"` South, `"W"` West)
#' by age, March 1999 to March 2000, for [group_transition()] on a region
#' dimension with these categories. Each origin-destination rate is the
#' all-ages flow divided by the origin's Census 2000 population, scaled by
#' the national age pattern of interstate moves. The same age pattern is
#' applied to every pair, which understates older-age moves to the South.
#' The rates apply to individuals; see [household()] for moving family
#' members together.
#'
#' @format Tibbles with numeric `age_start`, `age_end` and `rate` columns,
#'   plus:
#' \describe{
#'   \item{`usa_2000_mortality`}{44 rows; `sex` (`"f"` or `"m"`).}
#'   \item{`usa_2000_fertility`}{6 rows; no other columns.}
#'   \item{`usa_2000_fertility_union`}{12 rows; `union_type`
#'     (`"partnered"` or `"unpartnered"`).}
#'   \item{`usa_2000_marriage`}{12 rows; no other columns.}
#'   \item{`usa_2000_migration`}{132 rows; `origin` and `destination`
#'     regions, 12 directed pairs by 11 age bands.}
#' }
#' @source
#' Mortality: Arias E, Curtin LR, Wei R, Anderson RN. U.S. decennial life
#' tables for 1999-2001, United States life tables. *National Vital
#' Statistics Reports* 57(1). National Center for Health Statistics; 2008.
#' Tables 2 and 3. <https://www.cdc.gov/nchs/data/nvsr/nvsr57/nvsr57_01.pdf>
#'
#' Fertility: Hamilton BE, Sutton PD, Ventura SJ. Revised birth and fertility
#' rates for the 1990s and new rates for Hispanic populations, 2000 and 2001:
#' United States. *National Vital Statistics Reports* 51(12). National Center
#' for Health Statistics; 2003. Tables 1, 7 and 8.
#' <https://www.cdc.gov/nchs/data/nvsr/nvsr51/nvsr51_12.pdf>
#'
#' Marriage: Kreider RM. Number, timing, and duration of marriages and
#' divorces: 2001. *Current Population Reports* P70-97. U.S. Census Bureau;
#' 2005. Supplemental Table 2.
#' <https://www.census.gov/library/publications/2005/demo/p70-097.html>
#'
#' Migration: Schachter J. Geographical mobility: March 1999 to March 2000.
#' *Current Population Reports* P20-538. U.S. Census Bureau; 2001. Figure 3
#' and Tables B and D.
#' @name usa_2000
#' @keywords datasets
#' @examples
#' usa_2000_mortality
#' calculate_life_table(usa_2000_mortality)
#' calculate_tfr(usa_2000_fertility)
#' process <- mortality(usa_2000_mortality)
"usa_2000_mortality"

#' @rdname usa_2000
"usa_2000_fertility"

#' @rdname usa_2000
"usa_2000_fertility_union"

#' @rdname usa_2000
"usa_2000_marriage"

#' @rdname usa_2000
"usa_2000_migration"

#' England and Wales inputs for a second demographic transition example
#'
#' A named list of six tibbles used by the second demographic transition
#' article. Fertility and marriage tables have Hammel's schedule column names;
#' cohabitation and marital divorce/death observations retain their source
#' columns for the article's transformation examples. Cohabitation cohort
#' outcomes cover Great Britain; the other tables cover England and Wales.
#'
#' @details
#' Fertility rates are births per woman-year in 1992 by `union_type`:
#' `"marriage"`, `"cohab"` and `"unpartnered"`. Age bounds use the source's
#' denominator ages, 16--45, as in the article. The youngest and oldest birth
#' numerators are not strictly limited to these denominator ages; the final
#' band is finite, not open-ended. Source labels and duplicate rate units
#' have been removed without changing the rates.
#'
#' Marriage rates are annual central rates of first marriages and remarriages
#' for 1992--2019. `previous_union` maps these to `"none"` and `"dissolved"`;
#' this describes legal marital history, not prior cohabitation. Calendar
#' intervals are `[time_start, time_end)`. The 48 unavailable rates remain
#' `NA`; the article's filling and multiplicative tuning are not applied.
#' These marriage observations are not direct measurements of union-entry
#' hazards including cohabitation. Counts, denominators and source-method
#' columns have been removed.
#'
#' Premarital cohabitation proportions condition on marrying; they are not
#' probabilities that a newly formed union is a cohabitation. Cohabitation
#' outcomes are cumulative proportions at five and ten years for women under
#' 45 at cohabitation start. Source rounding and missing ten-year values for
#' the latest cohort are retained. Marital divorce/death tables are cumulative
#' competing-risk proportions, not annual hazards. No hazard conversion,
#' normalization, extrapolation or interpolation is applied to these tables.
#' HMD mortality is retrieved separately and is not included in this list.
#'
#' Original numerical sources and their provenance remain in
#' `validation/data/england_wales_1992_onward/` in the repository. The offline
#' preparation script is `data-raw/england_wales_1992_onward.R`.
#'
#' @format A named list containing:
#' \describe{
#'   \item{fertility}{18 rows: `union_type`, `age_start`, `age_end`, `rate`.}
#'   \item{marriage}{672 rows: `time_start`, `time_end`, `age_start`,
#'     `age_end`, `previous_union`, `rate`.}
#'   \item{premarital_cohabitation}{373 rows: `year`, `source_age_label`,
#'     `proportion_cohabiting_before_marriage`.}
#'   \item{cohabitation_outcomes}{3 rows: `cohort_start`, `cohort_end`,
#'     `marriage_5`, `separation_5`, `continuing_5`, `marriage_10`,
#'     `separation_10`, `continuing_10`.}
#'   \item{marital_divorce_death_2005}{52 rows: `year`, `anniversary`,
#'     `cumulative_divorce`, `cumulative_death`.}
#'   \item{marital_divorce_death_2010}{60 rows, with the same columns as 2005.}
#' }
#' Ages, durations and calendar times are in years. Rates are annual;
#' proportions are on the zero-to-one scale. Source cohorts are half-open
#' intervals, not calendar periods of a rate schedule.
#'
#' @source
#' Office for National Statistics (ONS). Contains public sector information
#' licensed under the Open Government Licence v3.0:
#' <https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/>
#'
#' Fertility: O'Leary L, Natamba E, Jefferies J, Wilson B (2010). Fertility and
#' partnership status in the last two decades. Population Trends 140, Table
#' B1. \doi{10.1057/pt.2010.10}.
#'
#' Marriage and premarital cohabitation: ONS Marriages in England and Wales,
#' 2019 workbook; female marital population estimates, 1971--2001 and
#' 2002--2019. Exact source URLs, retrieval dates and checksums are recorded in
#' the repository's source-folder `sources.csv`.
#'
#' Cohabitation outcomes: Beaujouan E, Ni Bhrolchain M (2011). Cohabitation and
#' marriage in Britain since the 1970s. Population Trends 145, Table 3.
#' \doi{10.1057/pt.2011.16}.
#'
#' Marital divorce/death, 2005: Wilson B, Smallwood S (2008). The proportion
#' of marriages ending in divorce. Population Trends 131, Table 3; cumulative
#' proportions reconstructed from life-table event counts. Marital
#' divorce/death, 2010: ONS Divorces in England and Wales, 2012 bulletin,
#' Figure 5. Publication tables and PDFs are not bundled with this dataset.
#' @name england_wales_1992_onward
#' @keywords datasets
#' @examples
#' england_wales_1992_onward$fertility
#' fertility(england_wales_1992_onward$fertility)
#' england_wales_1992_onward$cohabitation_outcomes
"england_wales_1992_onward"
