# Historical U.S. immigration counts

`us_immigration_c89_c119.csv` transcribes all 31 series for 1820–1970 from
the U.S. Bureau of the Census, *Historical Statistics of the United States,
Colonial Times to 1970*, part 1 (1975), pp. 105–109. See pp. 97–99 for
definitions and comparability notes. The transcription used the PDF text
layer and page images, with regional totals checked against component counts.
[Source PDF](https://fraser.stlouisfed.org/files/docs/publications/histstatus/hstat1970_cen_1975_v1.pdf)

Each row gives the printed reporting year and counts of immigrants.
During 1880–1960 the reporting year ends June 30. Earlier exceptions
are listed in the source table. The year printed as 1922 on p. 108
between 1823 and 1821 is aligned to 1822 using the other panels.

Zero represents a source dash. Empty cells indicate unavailable or
combined reporting, including Poland in 1899–1919, Mexico in 1886–1893,
and the source's exceptions for Asia, America, Africa and the Pacific.
Published totals and component counts are retained where they disagree.
The Europe components exceed C90 by 10 in both 1855 and 1870.
The Asia components fall short of C102 by 3 in 1872.

| Column | Source series |
| --- | --- |
| c89 | All countries |
| c90 | Europe, total |
| c91 | Great Britain |
| c92 | Ireland |
| c93 | Scandinavia |
| c94 | Other Northwestern Europe |
| c95 | Germany |
| c96 | Poland |
| c97 | Other Central Europe |
| c98 | U.S.S.R. and Baltic States |
| c99 | Other Eastern Europe |
| c100 | Italy |
| c101 | Other Southern Europe |
| c102 | Asia, total |
| c103 | Turkey in Asia |
| c104 | China |
| c105 | India |
| c106 | Japan |
| c107 | Korea |
| c108 | Philippines |
| c109 | Other Asia |
| c110 | America, total |
| c111 | Canada and Newfoundland |
| c112 | Mexico |
| c113 | West Indies |
| c114 | Other America |
| c115 | Africa, total |
| c116 | Australasia, total |
| c117 | Australia and New Zealand |
| c118 | Other Pacific Islands |
| c119 | All other countries |

Before 1906 the country series refer to the country from which immigrants
came. From 1906 onward they refer to last permanent residence. The source
notes changes in borders and tabulation, including the reporting of Northern
Ireland, Poland, and Austria. For crude U.S. immigration rates, divide the
counts by estimated U.S. person-years over the corresponding intervals.

## Annual U.S. resident population

`us_resident_population_a7.csv` transcribes Series A7, total resident
population, from the same volume, p. 8. It has one row per year from 1790
through 1970. The estimates are for July 1. The source reports in thousands,
and `resident_population` gives the equivalent number of people. Estimates
before 1900 are interpolated between decennial censuses. The report marks
1960–1970 as preliminary.

The source prints two 1959 values. The CSV keeps 176,289,000 on the same
geographic basis as the preceding years. The alternate estimate including
Alaska and Hawaii is 177,135,000. Estimates from 1960 onward include both
states. To approximate person-years for the July 1959 to June 1960
immigration reporting interval on a consistent geographic basis, use the
alternate 1959 value with the 1960 estimate.
