# Preparing HMD and HFD rate schedules

Specification for the user-requested multi-year getters and native HMD shapes,
October 30, 2026 (D119), extended by user-selected open-age pooling (D120).
Preparation is separate from simulation; no new demographic event mechanism
or interpolation is introduced.

## Selecting years

Both `get_mortality_hmd()` and `get_fertility_hfd()` accept `year` as one
whole-number calendar year or a nonempty vector of distinct whole-number years.
Requests are sorted chronologically; duplicates and unavailable requested
years are errors. A length-one request retains the existing static rate-table
schema without `time_start` or `time_end`. A request with more than one year
returns those columns, using the source periods as half-open intervals.
Annual year 2000 becomes `[2000, 2001)`.

Each requested product is downloaded once per call, then all requested years
are selected locally. Nonconsecutive requests retain their calendar gaps;
the simulator rejects gaps in a run rather than interpolating or extending a
schedule. Extra downloaded years are not returned. Rates retain their original
per-person-year units; they are not multiplied or divided by interval widths.

## Native HMD shapes

`get_mortality_hmd(..., shape = "1x1")` selects `Mx_<shape>`. Supported native
shapes are `1x1`, `1x5`, `1x10`, `5x1`, `5x5`, and `5x10`, where the first
number is age width and the second is period width. The database's five-year
age layout has infant age `[0,1)`, childhood `[1,5)`, then `[5,10)`, etc.,
and an upper-open age band. Raw inclusive labels such as `5-9` become `[5,10)`;
parsed numeric ages are interpreted using that documented layout and the
reader's `OpenInterval` flag. Raw labels inconsistent with the requested shape
are errors. Missing rates are rejected as before.

For grouped periods, `year` selects the source period's **starting year**.
For example, `year = c(2000, 2005), shape = "5x5"` returns source periods
`[2000,2005)` and `[2005,2010)` when those periods are available. Requesting
2001 does not silently select the period starting in 2000. A singleton request
for 2000 returns a static schedule estimated over its source period, without
calendar columns; source metadata still identifies that period.

Source period labels take precedence over nominal shape widths. A shortened
first or last source period is retained without inventing additional years.
HMDHFDplus 2.0.8's web reader parses labels even with `fixup = FALSE`, discarding
upper bounds. For grouped-time products with parsed labels, the getter also
downloads `Mx_1x1` once to establish annual coverage. Native period starts,
nominal width, the following period start and annual coverage jointly recover
the source interval; inconsistent or missing annual coverage is an error.
This auxiliary request is recorded in provenance. No credentials are retained.

These are the source's pooled mortality rates, interpreted as constant hazards
within the supplied bands. They are not local averages of annual rates or an
attempt to reconstruct the original single-age, single-year schedules.

## Optional pooled final age band

`get_mortality_hmd(..., open_age = 95)` retains source rates below age 95 and
replaces all ages 95 and above with one `[95, Inf)` band. The default `NULL`
retains the original age layout and missing-rate rejection. This is an open
age category, not a maximum lifespan or a forced-death age.

When requested, download `Deaths_<shape>` and `Exposures_<shape>` once each,
in addition to the products described above. For each selected source period
and sex, compute `sum(deaths) / sum(person_years)` over the entire tail.
Do not average age-specific rates, extend the cutoff-age rate, or treat missing
counts/exposures as zeros. Individual source Mx values in the replaced tail
may be missing; retained rates below the cutoff must still be finite.

The cutoff must be a finite nonnegative whole age aligned with a native age
band start: any whole age in single-age products, or 0, 1, 5, 10, etc. in
five-year-age products. It must exist in each selected period. Do not split
count/exposure bands or automatically move an unaligned cutoff. Counts and
exposures must match the rates' age/period bands, with a contiguous tail from
the cutoff through the source upper-open band. Selected tail values must be
finite and nonnegative, and pooled exposure must be positive. Zero-death,
zero-exposure individual cells contribute zero to both sums; a wholly
zero-exposure tail is an error. Unselected years are not pooled or validated
for numeric count/exposure values.

Worked example: ages 95, 96, and 97+ have deaths 2, 3, 0 and exposures
100, 20, 0 person-years. The 97+ source rate can be undefined, but the pooled
95+ rate is `5 / 120 = 0.0416667` per person-year. Zero pooled deaths with
positive exposure produce rate zero under the existing mortality contract.

Provenance records the cutoff, method, auxiliary product names and a tibble of
pooled deaths, person-years and rates by sex/source period. Singleton outputs
remain static; multiple-period outputs retain calendar bands. No source counts
or exposures are imputed. This is an explicit coarsening of old-age mortality;
the pooled rate is constant throughout the resulting open interval.

## HFD products and parity

HFD's existing `asfrRR` and conditional `mi` products remain on their native
annual, single-age layout. They do not offer HMD's six shape products, so no
HFD `shape` argument or local aggregation is introduced. The accepted endpoint
mapping remains: `12-` is `[12,13)` and `55+` is `[55,Inf)` (D095).

For `by_parity = TRUE`, each selected year retains all five current-parity
bands and the accepted zero-fill policy for unavailable conditional cells
(D094/D096). Source metadata counts replacements only across selected years,
including retained endpoint ages. Missing overall fertility rates remain errors.

## Provenance, worked examples and validation

Returned tables are tibbles. Source metadata records the requested years,
database, population, product, represented source periods, retrieval timestamp,
reader version, HMD shape and any auxiliary coverage product. HFD retains
revision, endpoint and conditional missing-rate metadata. RDS save/load and
process attachment retain this provenance.

Examples (downloads require the user's database credentials):

```r
# Existing static output, default source shape:
mortality_2000 <- get_mortality_hmd("SWE", 2000)

# Annual historical schedules; each product is downloaded once:
mortality_history <- get_mortality_hmd("SWE", 1891:2000, open_age = 95)
fertility_history <- get_fertility_hfd("SWE", 1891:2000)

# Native five-year ages and calendar periods:
mortality_grouped <- get_mortality_hmd("SWE", seq(1900, 1995, 5),
  shape = "5x5", open_age = 95)

# Retain ages below 95 and pool source deaths/exposure into a 95+ band:
mortality_truncated <- get_mortality_hmd("SWE", 1800:1950, open_age = 95)
```

Offline tests use synthetic raw/parsed tables and mocked web readers. Check
one versus multiple years, chronological ordering, missing/duplicate years,
download counts, all native HMD shapes, infant/child/open age widths, inclusive
period labels, shortened source periods, parsed-label coverage recovery,
unchanged rates, all HFD parity bands, selected missing-cell counts, invalid
requests before downloads, provenance and direct simulation with prepared tables.
For pooling, check the ratio of sums, sex/period separation, unchanged younger
rates, missing tail Mx, zero-exposure individual cells, invalid aggregate
exposure, missing/negative counts, source alignment, cutoff validation, download
counts, pooled provenance and direct simulation with a complete age schedule.
Live database tests remain explicitly opt-in.

Sources: [HMD explanatory notes](https://www.mortality.org/Data/ExplanatoryNotes),
[HMD product layout](https://www.mortality.org/Country/Country?cntr=SWE),
[HFD product layout](https://www.humanfertility.org/Country/Country?cntr=SWE),
and the installed HMDHFDplus 2.0.8 reader implementation. No upstream code is
copied into Hammel.
