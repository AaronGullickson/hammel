# Parity-specific fertility and PATFR

Milestone 5 accepted complete October 29, 2026 (D105). The user selected half-open parity
bands and HFD preparation, then authorized history-derived parity and parity
execution (D098–D099). Births use cell-accounted fertility with one cell per
actual parity (D101), which replaced D099's SOCSIM-style spacing correction.
Persistent lognormal fertility weights and maternal inheritance are implemented
under D102/D104; see [heterogeneity](fertility-heterogeneity-notes.md). The
[accepted report](../validation/milestone5.qmd) retains sparse/unrealized-cell
and reproductive-selection limitations; group transitions are next.

## Selected representation

Schedules add `parity_start` and `parity_end` alongside `age_start`, `age_end`
and `rate`. Parity bands are half-open, with nonnegative integer finite bounds
and `Inf` for the open final category. HFD conditional birth-order rates map
onto current maternal parity as follows:

| HFD rate | Current parity | parity_start | parity_end |
| --- | --- | --- | --- |
| First birth | 0 | 0 | 1 |
| Second birth | 1 | 1 | 2 |
| Third birth | 2 | 2 | 3 |
| Fourth birth | 3 | 3 | 4 |
| Fifth or higher birth | 4+ | 4 | Inf |

Rates are births per woman-year at the indicated age and current parity.
Unconditional birth-order ASFRs have the wrong exposure denominator for this
purpose. See [HFD definitions](https://www.humanfertility.org/Data/ExplanatoryNotes)
and [file formats](https://www.humanfertility.org/File/GetDocument/Docs/formats.pdf).

## HFD preparation

`get_fertility_hfd("USA", 2000, by_parity = TRUE)` requests the conditional
period product `mi` (`m1x`, `m2x`, `m3x`, `m4x`, `m5px`), selects one year, and returns a tibble of age/parity bounds
and unadjusted annual rates. The default FALSE retains the existing total-rate
`asfrRR` import. An optional `update` pins the HFD revision; source metadata
records the product and the birth-order/current-parity mapping without retaining
credentials. Missing products or years are errors, with no fallback to totals.
Under D095 both products retain endpoint rows: 12- maps to [12,13) and
55+ maps to [55,Inf), preserving supplied rates. This supersedes their earlier
omission; source metadata records the convention.

Under D096, unavailable HFD conditional age/parity cells are replaced with
zero at import, as selected by the user. Source metadata records
`missing_rates = "zero"` and the number of cells replaced in the selected year.
This is an explicit preparation assumption, not a claim that HFD observed zero
risk. Explicit zeros and known values remain unchanged; no interpolation is used.
Conditional rates are not transformed into probabilities or
adjusted for birth spacing during import.

The prepared parity table can be supplied to `set_fertility()` or `fertility()`.
D099 defines schedule coverage and actual parity tracking; D101 defines the birth mechanism.

## Selected simulation parity initialization (D098)

When parity-specific fertility schedules are detected, initialize an internal
integer parity vector from recorded maternal birth histories. Count every
recorded child linked by mother_id, including children who have died and births
with an unknown father. Women without recorded births start at zero, including
all women in a random baseline without birth histories. Use stable engine person
indices for storage, with the existing ID/index map for initialization. Do not
use an optional persons parity column to override the observed history.

Parity is the actual recorded count, not capped at the lower bound of an open
schedule band. Update it only on successful births; cancelled attempts do not
change it. Retained complete person histories reconstruct parity at restart.
Incomplete earlier histories consequently understate parity; no hidden births
or externally supplied parity values are inferred. This is the implemented initialization contract under D099.

## Observation and remaining scope

Parity-stratified observations use `calculate_fertility_rates(parity_breaks = ...)`
(D100), described below.
D102 defines persistent lognormal weights and maternal inheritance. Other
distributions, cohort normalization, group transitions and calendar-time inputs
remain future scope. The HFD-only zero-fill rule does not set
missing-data defaults for arbitrary rate tables. Remaining young high-parity
shortfalls are births in cells without an eligible woman; see below.

## PATFR calculation (D097)

`calculate_patfr(rates)` constructs a deterministic synthetic fertility cohort;
`calculate_tfr()` delegates to it when either parity-bound column is present
(and requires both). Output uses `patfr`, not `tfr`. Both accept `rate` or
`rate_fertility`, select female rows when sex is present, and separate the same
period, group, union_type and legacy union_status strata. They preserve group
types and the unassigned_births attribute. Parity is integrated over, not an
output grouping variable. Each group/union-type stratum has an independent
cohort held in that stratum; transitions between groups or union types are not
modeled by this conditional index. No engine initialization or mortality is involved.

Use the HFD period fertility-table midpoint convention (Methods Protocol
2015, sections 8.3–8.4), not a stochastic microsimulation or exact continuous-time
matrix exponential. Start with unit mass at parity zero at the first included
age. On a step of width d, use q = d*m/(1+d*m/2) for progression from each
nonterminal parity. For starting masses l[k], interval exposure is:

- k=0: d*l[0]*(1-q[0]/2);
- interior k: d*l[k]*(1-q[k]/2) + d*l[k-1]*q[k-1]/2;
- pooled upper parity P: d*l[P] + d*l[P-1]*q[P-1]/2.

Births at each parity are its exposure times its supplied rate. Move those
birth counts to the next parity, retaining women in the pooled final parity
while counting its repeated births. Sum births per starting woman to get PATFR.
Split at integer ages and every supplied age boundary, giving steps at most
one year and unchanged results when constant five-year bands are expanded to
single years. Reject progression probabilities above one rather than clipping
invalid table calculations. This generalizes HFD's annual midpoint convention
to fractional boundaries. Exact HFD agreement additionally requires identical
age endpoints, source revision and missing-rate assumptions.

Finite parity intervals have integer nonnegative boundaries; broader intervals
apply the same rate to every contained parity. At each included age, coverage
must run from zero through one final open parity interval, without overlap or
gaps. All parity states must cover the same contiguous finite age range. Birth
rates must be numeric/nonnegative; any missing included cell returns NA for its
whole stratum. Import-time zero filling (D096) remains distinct from this rule.

As selected by the user, upper-open age rows are excluded, with a warning if a
known omitted female rate is nonzero, because the finite-range PATFR
underestimates full-age fertility. Missing omitted rows do not invalidate the
finite result. Output includes the finite age_start/age_end used; absent finite
rows return a typed empty tibble. The resulting index depends on the selected
age range and the assumed parity-zero cohort at its beginning. No arbitrary
finite upper age or replacement for missing included rates is inferred.

## Parity-dependent execution

Schedule rules come from D099; the birth mechanism is cell-accounted fertility
(D101). Both parity bounds must be present or both absent. Within every
group/union-type stratum, distinct parity bands must partition [0,Inf) without
overlap or gaps. Finite bounds are nonnegative whole numbers within the engine's
integer-count capacity. Age bands within each parity band follow existing rules:
no overlaps, with uncovered ages having zero fertility. Parity partitions may
differ between strata. Overall and two-/three-type fertility and group
multipliers compose with parity; multipliers apply before cell accounting.
On-demand still requires overall union-type fertility but may distinguish
parity.

Track actual counts internally by person index under D098. Infer from all
recorded maternal births, ignoring any custom persons parity column. Each
successful birth increments the mother once; newborn girls start at zero;
cancellations, union changes and deaths do not increment parity. Retained
histories reconstruct the state at restart. No new authoritative public parity
column is added and no missing prehistory is inferred.

Each actual parity has its own fertility cells, so a pooled band such as
[4,Inf) applies its rate to actual parities 4, 5, 6, ... separately; cells are
created when a woman first reaches a parity. A birth moves the mother to the
next parity's cell as an ineligible member until her birth interval ends.
Father-choice calibration uses the cell-derived hazards of unpartnered eligible
women. See [cell-accounted fertility](fertility-cell-accounting.md) for the
mechanism, worked examples and validation.

For parity-specific schedules with positive spacing, `fertility()` warns when
positive rates fall in one-year age cells that no eligible woman can occupy;
parity k is first eligible after k singleton births at least `birth_interval`
apart, beginning at the schedule's first positive rates.

## Recovery with HFD inputs

With the user's HFD USA 2000 conditional rates (50,000-person baselines, 100
replications, years 2050–2150), conditional rates are recovered at most ages
and parities. Young high-parity cells fall short (for example parity 4+ at age
21 at .55 of input) because those cells hold very few women, nearly all of whom
have just given birth, so they often have no eligible member. As an accounting diagnostic, adding the
reported unrealized expected births with `calculate_fertility_rates(unrealized = TRUE)`
gives .98 of input for that cell; these are not additional realized births. Unrealized births were about 0.1% of
all births and fell to about 0.04% in a 200,000-person check, so the shortfall
shrinks with population size. The earlier preceding-parity correction (D099)
gave .53 for that cell with nothing to add back; see
[the investigation](parity-spacing-investigation.md) and Git history at commit
`89d451b` for its controlled recovery check.

## Observed fertility by parity

`calculate_fertility_rates()` measures maternal births per woman-year. Optional
`parity_breaks` defines half-open current-parity bands; a finite final boundary
receives `Inf`. Bounds are nonnegative integers starting at zero. NULL pools
parities. Reconstruct parity from every recorded maternal birth, including
pre-window births and deceased children; do not use custom parity columns.
Each birth counts at parity immediately before that birth and advances parity
for subsequent exposure. Simultaneous births advance in person-row order.
Exposure includes birth-spacing exclusions. Group and union-type options retain
the existing maternal attribution and birth-first endpoint rules. Fertility
outputs contain female rows only; mortality retains both sexes. Separate
`calculate_mortality_rates()` reports deaths per person-year. Neither public
function returns the other process's counts or rates. Unknown-mother births
remain in the fertility table's period-level `unassigned_births` attribute.
Summing parity counts and exposure recovers pooled female counts and exposure;
zero exposure gives NA rates. These are conditional on complete supplied history.


## Spacing-calibration investigation

[The investigation](parity-spacing-investigation.md) records the failures of
the preceding-parity correction, the independent SOCSIM reference experiment,
the composition diagnosis and the evidence for cell accounting (D101).
