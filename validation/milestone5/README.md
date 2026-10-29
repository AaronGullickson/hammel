# Milestone 5 supporting diagnostics

The user-authored report remains at [../milestone5.qmd](../milestone5.qmd).
This folder contains supplementary software/method checks, not an accepted
milestone assessment.

Run `Rscript validation/milestone5/check-parity-spacing.R` from the repository
root to regenerate the two CSVs. It uses offline constant inputs and ten
independent seeds (9501–9510), each with 3,000 women aged 20 and ten years of
follow-up, with no mortality and no births among descendants before observation
ends. Birth intervals are 0 and .75 years. All women start at parity zero;
this is an initialized cohort, not a stable-population experiment. Rates and
exposure are reconstructed from event histories independently of engine parity.

`parity-spacing-replications.csv` contains each rate/event/exposure estimate;
`parity-spacing-summary.csv` gives means and MCSEs. The zero-spacing comparison
supports schedule switching. With spacing, the selected SOCSIM-style correction
does not recover the higher-parity target rates in this case. See
[the parity specification](../../docs/parity-fertility.md) and
[testing](../../docs/testing.md) for the interpretation and limits. No empirical
rates, credentials or live downloads are needed; no SOCSIM run is implied.


The [SOCSIM reference experiment](socsim/README.md) adds independent executable
evidence: initialized cohorts, later generational windows, and age-18 rates
with/without spacing, including MCSEs and explicit zero-spacing limitations.


`check-composition-spacing.R` is an independent toy simulation (not engine
code) comparing no spacing, the D099/SOCSIM preceding-parity correction, and a
candidate composition-based correction derived from the target synthetic
cohort. Run `Rscript validation/milestone5/check-composition-spacing.R` from the
repository root to regenerate `composition-spacing-summary.csv` (dt = .0025,
seed 9601); pass a step size, e.g. `0.001`, to print a sensitivity rerun
without writing. The candidate is not an accepted model; see
[the investigation](../../docs/parity-spacing-investigation.md).

`check-runtime-spacing.R` extends that toy to a schedule change at t = 5,
comparing no spacing, the composition-based correction (table from the original
schedule), and candidate runtime cell accounting (hazard m * N / N_eligible in
each input parity cell). Run `Rscript validation/milestone5/check-runtime-spacing.R`
to regenerate `runtime-spacing-summary.csv` (seed 7), which includes expected
births lost in cells with no eligible women. Neither candidate is an accepted
model.

Engine checks of cell-accounted fertility (D101, `spacing_method = "cell"`):

- `check-cell-spacing.R [cores]` compares cell and midpoint methods over ten
  seeds in five synthetic scenarios (initialized cohort, generational young
  ages, restart to a high-fertility schedule, union types, groups) and writes
  `cell-spacing-replications.csv` and `cell-spacing-summary.csv`.
- `check-cell-resolution.R` reruns the generational scenario with the same
  rates supplied as one 25-year band, 1-year bands, 1-year actual-parity bands
  and 0.25-year bands, writing `cell-resolution-summary.csv`. Cells are split
  to at most one year and by actual parity, so the first three are identical;
  supplied 0.25-year bands give 0.25-year cells.
- `benchmark-cell-spacing.R LIB` times both methods in fresh processes using
  an optimized build installed in LIB (`R CMD INSTALL --library=LIB .`) and
  writes `cell-spacing-benchmark.csv`.

See [cell-accounted fertility](../../docs/fertility-cell-accounting.md#implementation-and-validation).

Midpoint compensation was removed after these checks (D101 is the only
method). `check-parity-spacing.R`, `check-cell-spacing.R` and
`benchmark-cell-spacing.R` compare against midpoint and run only at commit
`89d451b`; their CSVs are the recorded results. `check-cell-resolution.R`,
`check-composition-spacing.R` and `check-runtime-spacing.R` still run.
