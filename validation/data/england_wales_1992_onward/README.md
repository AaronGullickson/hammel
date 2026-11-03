# England and Wales demonstration evidence, 1992 onward

Prepared November 2, 2026. The user selected England and Wales as the main
empirical basis, starting in 1992, with outside evidence permitted for gaps.
These files collect observations and calibration targets. They are **not a
selected or calibrated simulation scenario**. No temporal interpolation,
extrapolation, search calibration, or foreign-country substitution is applied.
The user's milestone 7 report is unchanged.

## Available files

| File | Content and interpretation |
| --- | --- |
| `fertility_by_union_type_1992.csv` | 18 observations: six age groups by marriage, cohabitation and other/lone (mapped to unpartnered); candidate constant fertility baseline |
| `fertility_by_union_type_1992_2001_2007.csv` | All 54 published Table B1 estimates, retaining original age labels and denominator age limits |
| `female_union_type_population_1992_2001_2007.csv` | Table A4 female population stocks, useful as age-specific composition targets |
| `female_marital_population_1992_2019.csv` | Annual age-specific never-married, divorced and widowed female population estimates, including available survey precision flags |
| `female_marriage_rates_1992_2019.csv` | Annual first-marriage, pooled remarriage and all-marriage counts and derived central rates; 1,008 rows, 912 estimable rates and 96 missing rates |
| `female_published_pooled_marriage_rates_1992_2019.csv` | 308 published ONS pooled age-specific marriage rates in their native bands, including where the separate derivation has missing population cells |
| `female_living_arrangements_2002_2019.csv` | Annual female population stocks by age and detailed living arrangement; not entry probabilities |
| `female_premarital_cohabitation_1994_2019.csv` | ONS female shared-address proportions at marriage by age; conditional on marrying, all marriage orders |
| `marital_divorce_probabilities_by_duration_2010.csv` | 51 annual conditional divorce probabilities at exact durations 0 through 50 from ONS Figure 6; probabilities, not central rates or instantaneous hazards |
| `marital_divorce_death_cumulative_2010.csv` | Companion Figure 5 cumulative incidence of divorce and partner death, anniversaries 1 through 60 |
| `marital_life_table_2005.csv` | Full numerical transcription of the ONS 2005 life table: conditional probabilities, survivors, events and published cumulative percentages, with explicit interval boundaries |
| `marital_divorce_death_probabilities_2005.csv` | Separate conditional divorce and death probabilities in each native duration interval |
| `marital_divorce_death_cumulative_2005.csv` | Cumulative incidence reconstructed from the more precise life-table event counts, with `anniversary` referring to interval ends |
| `marital_divorce_hazard_anchors_2005_2010.csv` | Annual divorce hazard estimates in five-year duration bands for both verified period anchors; 45–50 hazard extended to 45+ |
| `marital_divorce_benchmarks_1993_1994.csv` | Selected early-period cumulative divorce benchmarks reproduced in the 2008 paper; insufficient for a complete competing-risk hazard schedule |
| `marriage_cohort_cumulative_divorce_1992_2018.csv` | Cumulative proportions divorced by anniversary for marriage cohorts beginning in 1992 or later, observed through the 2019 release |
| `great_britain_cohabitation_cohort_outcomes.csv` | Published female competing-risk outcomes for cohabitations beginning in 1990–94, 1995–99 and 2000–04; marriage, separation, or continuing cohabitation at five/ten years |
| `great_britain_premarital_cohabitation.csv` | Female first-marriage premarital cohabitation measures in 1990–94, 2000–03 and 2004–07; the first interval spans 1992 |

The `year` columns identify observations, not automatically scheduled inputs.
CSV files can be read with `readr::read_csv()` or
`tibble::as_tibble(read.csv(...))`. The original ONS Excel workbooks are retained
locally under the Open Government Licence v3.0. Published paper tables are
transcribed as numerical observations with attribution; their PDFs are not
redistributed here. `sources.csv` records exact download URLs, retrieval date,
licence and MD5 checksums for the six retained ONS source files.

## Fertility ages, open bands and TFR

Source: O'Leary, Natamba, Jefferies and Wilson (2010),
[Fertility and partnership status in the last two decades](https://link.springer.com/content/pdf/10.1057/pt.2010.10.pdf),
Population Trends 140, Tables B1 and A4, printed pp. 31 and 29; methods pp. 21–23.

The source uses female population denominators aged **16–44**. Accordingly,
`denominator_age_start`/`denominator_age_end` give `[16,20)`, `[20,25)`,
`[25,30)`, `[30,35)`, `[35,40)`, `[40,45)`. These deliberately are not named
`age_start`/`age_end`: mapping them into a simulation schedule requires an
explicit choice about end bands. Original labels are retained as `Under 20`
and `40 and over`. Birth numerators are not necessarily restricted to these
denominator ages; the source explicitly includes births below age 16. Treating
the last row as a conventional 40–44 ASFR also requires caution about births
above 44. Do not extend the pooled rate to infinite ages.

For finite age intervals the TFR contribution is rate times age width.
A genuinely open-ended ASFR has no identifiable contribution without an
upper age limit or a model for the age pattern. A pooled 40+ rate cannot reveal
the separate 40–44 and 45–49 ASFRs. Candidate approaches are to use a documented
finite reproductive-age window or borrow an overall age pattern to split the
tail while preserving the observed birth total with appropriate exposures.
Neither has been silently selected here.

Hammel's current `calculate_tfr()` warns and omits nonzero open age bands.
That is not a complete TFR when such a band contains births. Simply replacing
`Inf` in a pooled output table with 45 does not undo pooling of exposures at
older ages: calculate rates with the intended finite boundaries first.
Conditional status-specific sums describe hypothetical continuous exposure to
each status, not additive components of an overall TFR.

Rates are births per woman-year after dividing the published values by 1,000.
Cohabitation is inferred from joint birth registration at a shared address,
using cohabiting population estimates (Census denominators in 2001).
These are estimates, not perfect observed partnership histories. Married
denominators use legal status, so separation/repartnering can cause mismatch
with the birth classification. Population counts are rounded and need not sum
exactly. No constant married/cohabiting fertility ratio is imposed.

## Marriage rates and their denominators

Sources:

- ONS [Marriages in England and Wales, 2019 edition](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/marriagecohabitationandcivilpartnerships/datasets/marriagesinenglandandwales2013), Tables 3a and 11a.
- ONS [historical marital-status population backseries](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/backseriespopulationestimatesbymaritalstatusselectedyearsfrommid1971tomid2001singleyearofageandsexforenglandandwales), annual sheets 1992–2001; source units thousands, multiplied by 1,000.
- ONS [population by marital status and living arrangements](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesbymaritalstatusandlivingarrangements/current), archived July 2020 workbook (file `2020ewpublicationtablefinal.xlsx`, covering 2002–2019), female sections of Tables 1 and 2. Archive version v7.

First marriage = count of previously single women marrying / never-married
female midyear population. Remarriage = divorced plus widowed women marrying /
divorced plus widowed female midyear population. All-marriage counts here sum
these known previous statuses. Numerators cover opposite-sex marriages.
Later source definitions include dissolved civil partnerships among divorced
and surviving civil partners among widowed. No same-sex marriage counts are
added. Counts refer to marriages occurring in England/Wales, whereas the
denominator is resident population, so marriages abroad and nonresident
marriages can produce a mismatch.

Native count bands are 16–19, five-year bands 20–59, 60–69, 70–79, 80+.
Single-age or five-year population cells are summed to match these exactly.
Unstated ages are excluded after 1993; ONS already redistributed unstated
ages for 1961–93 (workbook note 10). Unknown previous status is excluded.
The rebuild checks that each year's band counts reconcile with its
status-specific age-known totals.

There is a **population-method break at 2002**: the later estimates use survey
marital-status distributions. These calculations are our derived central rates,
not a claim to reproduce a particular vintage of published ONS rates. Missing
or suppressed denominator cells remain missing, including when only one
component of the pooled remarriage denominator is missing. No suppressed
widow count is replaced by zero. Source CV codes and confidence half-widths
are retained in the population file; rate uncertainty is not estimated here.

These legal marriage rates are **targets**, not direct Hammel search inputs.
They include marriages between cohabiting partners. The simulator separates
new unions from same-partner conversion, and its search population excludes
current cohabiters. Never legally married also differs from never previously
partnered if cohabitation counts as previous union. Search/type/conversion
inputs therefore need a joint mapping, rather than copying legal marriage
rates into both search and conversion schedules.

## Dissolution and cohabitation evidence

ONS [2012 divorce bulletin](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/bulletins/divorcesinenglandandwales/2014-02-06),
Figures 5 and 6, supplies the retained **2010 period** divorce evidence.
Figure 6's conditional annual divorce probabilities arise from the ONS
divorce/death life-table construction. They must not be passed as hazards
without stating a conversion approximation; companion mortality figures are
rounded and do not uniquely recover unrounded instantaneous hazards.

ONS [2019 divorce workbook](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/datasets/divorcesinenglandandwales),
Table 6, supplies the later marriage-cohort curves. ONS
[methods](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/methodologies/userguidetodivorcestatistics)
explain that the cohort denominator is the original marriage count and does
not remove marriages ended by death. These curves are not death-censored
survival curves and are not period hazards. Legal divorce also occurs after
physical separation, which matters for partnership status and fertility.

Supplementary source: Beaujouan and Ni Bhrolchain (2011),
[Cohabitation and marriage in Britain since the 1970s](https://link.springer.com/content/pdf/10.1057/pt.2011.16.pdf),
Population Trends 145, Tables 2 and 3, printed pp. 42 and 47. This covers
**Great Britain**, not exclusively England and Wales, and uses retrospective
GHS histories. The retained Table 3 measures refer to women under 45 at union
start. They are competing-risk outcomes: marriage ends the cohabitation spell,
even if a subsequent divorce occurs. Rounded outcomes can total 100.1%.
They can support an explicit approximate hazard fit, but are not annual rates.
The ten-year observation is unavailable for the latest cohort.

ONS Table 11a's shared address at marriage and the GHS premarital histories
are different measures. Both condition on marrying. Neither equals Hammel's
`cohab_prob`, the probability a **new union** starts as cohabitation. Population
prevalence also cannot be substituted for this entry probability.

## Remaining scenario choices and gaps

- Hold the 1992 fertility estimates fixed, after selecting finite age support
  or an explicitly borrowed tail pattern.
- Retrieve the chosen fixed mortality baseline separately; no HMD login,
  credential access or mortality download was performed in this preparation.
- Map legal marriage targets to union searches, entry type and conversion.
  Some small remarriage denominator cells need an explicit approximation or
  a different/coarser source; missingness is retained for review.
- Choose cohabitation entry probabilities using composition/premarital targets
  as evidence. Direct entry probabilities are not identified by these tables.
- Use British cohabitation outcomes to inform conversion and breakup patterns;
  evidence beyond the early 2000s remains missing here.
- The 2005 and 2010 marital duration patterns are available. The requested
  early anchor uses 1993–94 as a stand-in for 1992, but its full life table has
  not been retrieved. Interpolation to 2005 cannot yet be completed from the
  verified data. Hold the last anchor constant after 2010 as requested.
  Danish SKI5
  duration-specific observations are a possible external supplement if a
  historical period series is needed, but no Danish values have been mixed in.
- Reproducing historical British TFR is not the target: fertility within union
  types will deliberately be held fixed while composition changes.

## Additional marital divorce anchors: retrieval status

The user requested the 1993–94 and 2005 period estimates, assigning the early
estimate to 1992, interpolating between anchors, and holding rates fixed after
2010. The verified 2005 extraction is from Wilson and Smallwood (2008),
[The proportion of marriages ending in divorce](https://sociologytwynham.com/wp-content/uploads/2008/07/proportion-of-marriages-ending-in-divorce.pdf),
Table 3, printed p.33 (PDF page 6). Attribution, download checksum, extraction
method and requested scenario conventions are in
`marital_divorce_additional_sources.json`. The paper is downloaded only to a
temporary location; the numerical observations are retained here.

Source duration `x` labels the **start** of an interval. Columns g/h report
cumulative incidence **before the next stated anniversary**, at the end of that
interval. For example, row x=4 ends at anniversary 5. The last native intervals
are 50–55 and 55–60, not single years. The cumulative output reconstructs
incidence from the rounded integer event counts per 100,000 initial marriages,
which preserve substantially more precision than the one-decimal published
percentages. Life-table survivor accounting reconciles within one marriage
between adjacent rows. The reconstruction matches each published cumulative
percentage to its stated rounding precision.

`prepare_marital_divorce_anchors.R` produces both 2005 cumulative/probability
files and five-year divorce hazard anchors for 2005 and 2010. It uses the
constant competing-hazards formula within each five-year band and verifies
recovery of both conditional outcomes. The death hazard is used in this
calculation only; the output `rate` is divorce, not divorce plus death.
2005 uses event-count precision; 2010 uses the already retained rounded
cumulative curves. The 45–50 estimate is extended to 45+ following the user's
tail convention; `observed_dur_end` preserves the actual estimation endpoint.
No interpolated calendar schedule has been fabricated.

For 1993–94, Tables 4 and 5 reproduce cumulative divorce figures labelled at
durations 10, 20, 30, 40 and 50. These are preserved with their **source labels**,
not silently converted into matching five-year interval endpoints. Table 4's
2005 comparisons also differ in anniversary alignment from Table 3's interval-end
cumulative columns, so original-table definitions need checking before pooling.
Corresponding early mortality values and the full duration profile have not
been recovered. These few divorce benchmarks do not identify both hazards.
The original source is Haskey (1996),
[The proportion of married couples who divorce: past patterns and current prospects](https://pubmed.ncbi.nlm.nih.gov/8729793/),
Population Trends 83, pp.25–36. Public source and author-profile searches located
the citation and abstract, but not a retrievable full table. The 2008 authors
also describe a change in the treatment of mortality and migration compared
with earlier analyses; even with the early table, interpolation would be a
demonstration approximation across estimation methods.

A complete 1992–2010 interpolation still needs that original early table, or
an explicit alternative assumption for the early mortality/duration pattern.
Neither missing death probabilities nor intermediate five-year values are
filled by default.

## Rebuild and checks

From the package root run:

```r
source("validation/data/england_wales_1992_onward/prepare.R")
```

This rebuild is offline. It validates year coverage, nonnegative rates,
probability bounds and rounded competing-risk totals, and reconciles marriage
band counts against the workbook's age-known totals. Existing package tests
are not changed or run: no package behavior is modified. A separate empirical
or simulation validation has not been performed.
