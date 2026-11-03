# SOCSIM parity/spacing reference experiment

Executed October 29, 2026 on the locally available MPIDR/rsocsim 1.9.18 source
at commit `125450cff968bb3ff762fa99d4c926a8136062a2`, clean checkout.
Root license GPLv3. Engine source was not modified, copied into Hammel, or
adapted. A scratch copy was compiled and installed under `/tmp`; no rates were
downloaded. These are narrow synthetic experiments, not milestone acceptance.

## Design and observation

Annual input hazards are .2 at current parity zero, .3 at parity one, and .4 at
parity 2+, constant from age 15 to 40 and zero elsewhere. SOCSIM fertility
`add_rate_block()` uses `lambda = -input`, so supplied fertility values are
annual hazard/12, **not** `1-exp(-annual/12)`. Mortality and marriage use a
different input conversion in that source. SOCSIM remains a monthly engine;
exact recovery of continuous annual hazards is not assumed.

Compare `bint 0` and `bint 9` (nine months = .75 years), ten seeds 9801–9810
per scenario. Run each in a fresh R process, avoiding process-global reference
state between runs. Initial persons are women aged about 20 at parity zero;
a single male newborn anchors the last observed month at 1200. Simulation runs
from month 1201. No marriages, random fathers, group differences or status
transitions; births outside partnerships are retained. Mortality is zero before
60 and .1 monthly probability thereafter. Birth spacing affects births only.

Set `hetfert 0` and, independently, fix all fertility multipliers to one with
founder fmult=1, `alpha 1`, `betaT 1`. Saved histories confirm unity for all women
and single marital status. A preliminary generational pilot using only the
toggle had variable recorded descendant multipliers and poor zero-spacing
recovery; it is excluded. No claim about the toggle's implementation is made
without a separate investigation. All CSVs here use explicitly fixed multipliers.

Two experiments:

- **Initialized cohort:** 3,000 founder women, ten years, founder exposure only.
  Descendants cannot reproduce in that interval. Exposure sums to 30,000
  woman-years in each replication.
- **Generational:** 300 founder women, 100 years, including descendant women.
  Rates measured at ages [20,35) in successive 20-year windows. These are later
  windows, not an assertion that the population is mathematically stationary.
  A separate measurement uses age [18,19) during years [80,100), with actual
  parities 0,1,2,3,4+ (the last three share input .4).

Exposure and parity are reconstructed independently from raw maternal birth
histories, including deceased children. A monthly risk census attributes the
birth month to pre-birth parity and advances exposure the following month;
multiple same-month births advance successively. Counts and exposure use
half-open observation/age windows. These monthly observations are not Hammel's
continuous Lexis kernel. Small discrepancies in the zero-spacing control must
be retained when interpreting the extra spacing effect.

## Results

Means across ten separate replication rates; MCSE = their sample SD/sqrt(10).
All cells below have ten contributing replications.

| Experiment / window | Parity | Input | No spacing, mean (MCSE) | Nine months, mean (MCSE) |
| --- | --- | --- | --- | --- |
| Initialized, years 0–10 | 0 | .20 | .197759 (.001074) | .197759 (.001074) |
| Initialized, years 0–10 | 1 | .30 | .308838 (.002105) | .270157 (.001591) |
| Initialized, years 0–10 | 2+ | .40 | .414309 (.001831) | .344593 (.001353) |
| Generational, years 80–100, ages 20–35 | 0 | .20 | .197649 (.000905) | .199626 (.000759) |
| Generational, years 80–100, ages 20–35 | 1 | .30 | .299209 (.000650) | .304430 (.000762) |
| Generational, years 80–100, ages 20–35 | 2+ | .40 | .408200 (.000469) | .382117 (.000389) |

The initialized-cohort reference run reproduces the substantial spacing-related
higher-parity shortfall observed in Hammel's corresponding test. Later-window
results show that its size depends on composition: parity one is close to input,
while the pooled top category remains below input. Zero-spacing controls have
small systematic offsets as well, so this does not isolate every monthly engine
or observation effect. The evidence supports a limitation of SOCSIM's correction
for these inputs, not identical numerical behavior to Hammel, universal failure
at every parity, or a diagnosed new source-code defect.

The later-window age-[18,19) comparison is more pronounced, with all ten
replications contributing in every cell:

| Actual current parity | Input | No spacing, mean (MCSE) | Nine months, mean (MCSE) |
| --- | --- | --- | --- |
| 0 | .20 | .198435 (.000926) | .199657 (.001340) |
| 1 | .30 | .302024 (.002211) | .272932 (.001755) |
| 2 | .40 | .413259 (.003382) | .284820 (.003407) |
| 3 | .40 | .428620 (.006431) | .186698 (.005988) |
| 4+ | .40 | .435636 (.008090) | .054263 (.009823) |

This reproduces the user's qualitative pattern after long burn-in: first births
are close, while higher-parity fertility at young ages is strongly suppressed
with spacing. The zero-spacing offsets, particularly at the highest young-age
parities, mean the reference experiment is not a clean validation of all other
SOCSIM behavior. They are much smaller than the additional positive-spacing
shortfalls here; do not claim that no-spacing rates match exactly.

See `age18-summary.csv` for the later-window young-age comparison and
`age18-replications.csv` for its exposure/event denominators. The separate
[collaboration note](../../../design/parity-spacing-investigation.md) distinguishes
these reference results from the user's empirical-schedule results.

## Files and reproduction

- `run-toy.R`: generates inputs, runs one fresh-process job, independently
  reconstructs parity and exposure, writes raw persons and `rates.csv` to scratch.
- `summarize.R`: combines the 40 jobs into the saved `replications.csv` (360 rows)
  and `summary.csv` (36 rows), including intermediate observation windows.
- `check-young-age.R`: remeasures the 20 saved generational histories at age 18,
  verifies female multipliers are one, saves 100 replication cells and 10 summaries.

From the Hammel repository root, with the pinned local checkout at `../rsocsim`:

```sh
cp -a ../rsocsim /tmp/hammel-rsocsim-reference
mkdir -p /tmp/hammel-reference-library
R CMD INSTALL --no-multiarch --library=/tmp/hammel-reference-library /tmp/hammel-rsocsim-reference
```

Run the following R orchestration; it invokes a new Rscript process per job.
Use an empty scratch directory to prevent mixing prior pilot outputs.

```r
scratch <- "/tmp/hammel-socsim-fixed-runs"
dir.create(scratch, recursive = TRUE)
for (experiment in c("cohort", "long")) {
  months <- if (experiment == "cohort") 120 else 1200
  women <- if (experiment == "cohort") 3000 else 300
  for (bint in c(0, 9)) for (seed in 9801:9810) {
    folder <- file.path(scratch, sprintf("%s-b%d-s%d", experiment, bint, seed))
    status <- system2("Rscript", c("validation/milestone5/socsim/run-toy.R",
      "/tmp/hammel-reference-library", folder, bint, seed, months, women),
      stdout = paste0(folder, ".log"), stderr = paste0(folder, ".log"))
    stopifnot(status == 0)
  }
}
system2("Rscript", c("validation/milestone5/socsim/summarize.R", scratch,
                    "validation/milestone5/socsim"))
system2("Rscript", c("validation/milestone5/socsim/check-young-age.R", scratch,
                    "validation/milestone5/socsim"))
```

Raw output directories and compiled reference library remain scratch artifacts;
committed inputs/runner and aggregate evidence are sufficient to reproduce them.
No source-code licensing claim is made beyond recording the reference provenance.
