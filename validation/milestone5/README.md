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
