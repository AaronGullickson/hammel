# Parity-specific fertility and PATFR

Milestone 5 design, October 29, 2026 (D094). The user selected half-open parity
bands and HFD preparation, then authorized history-derived parity and the
SOCSIM-style birth-spacing correction (D098–D099). Parity-specific execution is
implemented; fertility heterogeneity remains unselected.

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
D099 defines schedule coverage, actual parity tracking and spacing compensation.

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

## Remaining scope

Parity-stratified observations use `calculate_fertility_rates(parity_breaks = ...)`
(D100), described below.
Heterogeneity distribution, persistence, association with other attributes and
calibration remain unselected. The HFD-only zero-fill rule does not set
missing-data defaults for arbitrary rate tables. The SOCSIM-style spacing
correction is approximate; see the validation evidence below.

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

## Parity-dependent execution (D099)

The user authorized the SOCSIM-style correction on October 29, 2026. This
implements parity-specific fertility, not fertility heterogeneity. Both parity
bounds must be present or both absent. Within every group/union-type stratum,
distinct parity bands must partition [0,Inf) without overlap or gaps. Finite
bounds are nonnegative whole numbers within the engine's integer-count capacity.
Age bands within each parity band follow existing rules: no overlaps, with
uncovered ages having zero fertility. Parity partitions may differ between
strata. Overall and two-/three-type fertility and group multipliers compose
with parity; multipliers apply before compensation. On-demand still requires
overall union-type fertility but may distinguish parity.

Track actual counts internally by person index under D098. Infer from all
recorded maternal births, ignoring any custom persons parity column. Each
successful birth increments the mother once; newborn girls start at zero;
cancellations, union changes and deaths do not increment parity. Retained
histories reconstruct the state at restart. No new authoritative public parity
column is added and no missing prehistory is inferred.

For actual parity k=0, use the unadjusted age schedule. For k>0, let f_k be the
current band's input rate after group multiplication. Evaluate the integral of
the preceding actual parity's unadjusted schedule over the birth interval at
the current age band's midpoint; adjusted hazard is
f_k / (1 - integral(previous_parity_rate)). At an upper-open age band, use the
preceding parity schedule's interior tail rate times the interval instead.
A missing upper tail contributes zero. Zero input fertility remains zero;
zero spacing leaves hazards unchanged. Reject nonpositive eligibility for a
positive rate and nonrepresentable adjusted hazards. Do not recursively use
already-adjusted rates. This is an approximation, not the exact current-parity
eligible fraction and not a guarantee of arbitrary conditional-rate recovery.
It retains the same group/union-type stratum for the preceding schedule; prior
status transitions are not reconstructed in compensation (D045).

Lookup uses actual k and k-1 even when input parity bands pool several counts.
For [4,Inf), actual parity 4 uses the parity-3 schedule for correction; parity
5 and above use the pooled schedule itself. A broad [0,Inf) parity input also
leaves first births unadjusted; omission of parity columns retains legacy
unstratified compensation. These are deliberately different interpretations.

After birth, select the new schedule and retain the actual-birth spacing clock.
Refresh any unpartnered woman's calibration contribution and future boundary
entry using the new adjusted hazard, without adding parity dimensions to the
calibration bins. Union-search intensities remain governed by search schedules;
father-choice intensities use current parity fertility. Returned histories and
ordinary schemas stay unchanged when parity columns are absent.

Worked software checks: parity-0 rate .2 is unchanged with spacing .75;
parity-1 rate .3 with prior .2 in a constant region becomes .3/(1-.75*.2).
For a final [2,Inf) band at .4, actual parity 2 corrects from prior .3 while
actual parity 3 corrects from .4. History-derived parity counts deceased
children and ignores a conflicting custom column. Tests must verify age-band
crossings, upper-open ages, spacing, cancellation, union types, groups, restarts,
calibration updates and unstratified fixed-seed controls. Statistical checks
report the approximation's recovery or bias rather than setting scientific
acceptance thresholds.


## Controlled recovery check and limitations

The offline diagnostic at
[validation/milestone5/check-parity-spacing.R](../validation/milestone5/check-parity-spacing.R)
uses ten independent seeds, each with 3,000 women aged 20, all initially
nulliparous, followed for ten years without mortality. Inputs are constant
.2 at parity 0, .3 at parity 1 and .4 at parity 2+. Cohort parity exposure is
reconstructed independently from recorded birth times. This is an initialized
cohort check, not a stable-population or empirical HFD validation.

| Birth interval | Parity | Mean observed rate | MCSE | Input |
| --- | --- | --- | --- | --- |
| 0 | 0 | .199760 | .001479 | .2 |
| 0 | 1 | .299198 | .002280 | .3 |
| 0 | 2+ | .403456 | .001528 | .4 |
| .75 | 0 | .199760 | .001479 | .2 |
| .75 | 1 | .266465 | .001925 | .3 |
| .75 | 2+ | .333460 | .000637 | .4 |

The zero-spacing control supports schedule switching and exposure accounting.
The selected preceding-parity approximation does not recover the higher-parity
inputs in this initialized-cohort case: actual ineligible fractions depend on
parity composition, which this correction does not estimate. Successful software
checks are therefore not scientific acceptance of its rate recovery. Do not
silently add population-ratio corrections; a different compensation model needs
user agreement. Replication estimates and the summary CSV are retained beside
the diagnostic script. The main validation report remains under user development.

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


## Current spacing-calibration investigation

The user's later-window empirical-schedule checks also show large high-parity
shortfalls with positive spacing and close recovery with spacing disabled. See
[the investigation and collaboration handoff](parity-spacing-investigation.md)
for those reported numbers, the independent SOCSIM reference experiment, a
composition diagnosis of the preceding-parity correction with a candidate
cohort-based alternative, and unresolved questions. A general eligibility-calibration framework is proposed
for investigation only; no replacement algorithm has been selected. D099
remains the implemented approximation with a known rate-recovery limitation.
