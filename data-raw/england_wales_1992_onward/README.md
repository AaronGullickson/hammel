# England and Wales demonstration evidence, 1992 onward

Run `Rscript data-raw/england_wales_1992_onward/england_wales_1992_onward.R` from the
package root with `here`, `readr`, `readxl`, and `tibble` installed. The script
reads the original sources below
and saves the six-table bundle to `data/england_wales_1992_onward.rda`.

| Input | Use |
| --- | --- |
| `ons_marriages_2019.xlsx` | Annual female marriage counts and premarital cohabitation |
| `ons_marital_population_1971_2001.xls` | Female marriage rate denominators through 2001 |
| `ons_marital_population_2002_2019.xlsx` | Female marriage rate denominators from 2002 |
| `marital_life_table_2005.csv` | Transcribed life-table observations for cumulative divorce and death |
| `ons_divorce_death_cumulative_2010.csv` | Original ONS cumulative divorce and death estimates |

The script includes the transcribed fertility and cohabitation cohort values.
`sources.json` records the exact ONS download URLs, retrieval dates, licences,
and checksums. Contains public sector information licensed under the
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Fertility ages, open bands and TFR

Source: O'Leary, Natamba, Jefferies and Wilson (2010),
[Fertility and partnership status in the last two decades](https://link.springer.com/content/pdf/10.1057/pt.2010.10.pdf),
Population Trends 140, Table B1, printed p. 31; methods pp. 21–23.

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
- ONS [population by marital status and living arrangements](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesbymaritalstatusandlivingarrangements/current), archived July 2020 workbook (file `2020ewpublicationtablefinal.xlsx`, covering 2002–2019), female section of Table 1. Archive version v7.

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
widow count is replaced by zero. The bundled rates carry no uncertainty estimates.

These legal marriage rates are **targets**, not direct Hammel search inputs.
They include marriages between cohabiting partners. The simulator separates
new unions from same-partner conversion, and its search population excludes
current cohabiters. Never legally married also differs from never previously
partnered if cohabitation counts as previous union. Search/type/conversion
inputs therefore need a joint mapping, rather than copying legal marriage
rates into both search and conversion schedules.

## Cohabitation outcomes

Source: Beaujouan and Ni Bhrolchain (2011),
[Cohabitation and marriage in Britain since the 1970s](https://link.springer.com/content/pdf/10.1057/pt.2011.16.pdf),
Population Trends 145, Table 3, printed p. 47. This covers
**Great Britain**, not exclusively England and Wales, and uses retrospective
GHS histories. The retained Table 3 measures refer to women under 45 at union
start. They are competing-risk outcomes: marriage ends the cohabitation spell,
even if a subsequent divorce occurs. Rounded outcomes can total 100.1%.
They can support an explicit approximate hazard fit, but are not annual rates.
The ten-year observation is unavailable for the latest cohort.

## Cumulative divorce and death

The 2010 period estimates come from Figure 5 of the ONS
[2012 divorce bulletin](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/divorce/bulletins/divorcesinenglandandwales/2014-02-06).
The published cumulative proportions are rounded.

The verified 2005 extraction is from Wilson and Smallwood (2008),
[The proportion of marriages ending in divorce](https://sociologytwynham.com/wp-content/uploads/2008/07/proportion-of-marriages-ending-in-divorce.pdf),
Table 3, printed p.33 (PDF page 6). Attribution, download checksum, and extraction
method are in
`sources.json`. The paper is downloaded only to a
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
