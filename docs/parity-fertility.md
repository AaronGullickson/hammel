# Parity-specific fertility: selected format and data preparation

Milestone 5 design, October 29, 2026 (D094). The user selected half-open parity
bands and requested an HFD importer. This does not specify or implement parity-
dependent simulation or fertility heterogeneity.

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

The parity table is data preparation for milestone 5. `set_fertility()` and
`fertility()` reject parity columns until the engine implements their meaning,
rather than silently treating them as ordinary unstratified age rates.

## Unresolved implementation choices

Define parity state/initialization for incomplete birth histories, grouping and
coverage validation, transitions after birth, interactions with birth spacing
and compensation, calibration contributions, and parity-stratified observation
before implementing execution. Heterogeneity distribution, persistence,
association with other attributes and calibration remain unselected. The HFD-only zero-fill rule does not set missing-data defaults for arbitrary
rate tables or select a parity maximum or heterogeneity distribution.
