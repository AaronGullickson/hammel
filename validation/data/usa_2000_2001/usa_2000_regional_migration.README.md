# United States: migration between census regions, around 2000

`usa_2000_regional_migration.csv`: ballpark annual migration rates between the
four census regions (`NE` Northeast, `MW` Midwest, `S` South, `W` West) by
age. Columns `origin`, `destination`, `age_start`, `age_end` and `rate` (moves
per person-year), ready for
`set_group_transition(model, rates, group = "region")` with a region
dimension using these labels. 12 directed pairs x 11 age bands; the last band
is open (`age_end = Inf`). Rebuild with `prepare_usa_2000_regional_migration.R`.

## Source

U.S. Census Bureau, *Geographical Mobility: March 1999 to March 2000*, Current
Population Reports P20-538 (Jason Schachter, 2001), from the March 2000
Current Population Survey; archived here as
`usa_mobility_1999_2000_census_p20-538.pdf`.

- Figure 3: migrants between regions, March 1999 to March 2000 (thousands).
  The figure's labels do not survive text extraction, so its values were
  assigned to flows using Table D's regional in- and out-migrant totals; the
  assignment is the only one consistent with all eight totals and reproduces
  them exactly (the script checks this).
- Table B: percent moving to a different state, by age.
- Denominators: Census 2000 regional populations (NE 53.594, MW 64.393,
  S 100.237, W 63.198 million).

## Method

rate(origin, destination, age) = flow(origin, destination) / population(origin)
x m(age), where m(age) is the age group's interstate moving rate divided by
the all-ages interstate rate (3.1 percent). About 11.5 per 1,000 people moved
between regions that year. Annual proportions are used as rates per
person-year, which is close at these levels.

## Caveats

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
