# SOCSIM rate conversion: initial source trace

Inspected October 12, 2026. Reference: https://github.com/MPIDR/rsocsim at commit `125450cff968bb3ff762fa99d4c926a8136062a2`, full repository cloned into `/home/aarong/Projects/rsocsim`. This is a narrow static source trace, not an executable validation or a complete scheduling audit. No source code was copied into Hammel.

## Ordinary rate-file inputs

In `src/src/load.cpp`, `add_rate_block()` (lines 1643–1698) treats births and other events differently:

- Birth input `r` is stored as `lambda = -r`.
- Death, marriage, divorce, and transition input `p` is stored as `lambda = log(1 - p)`.
- Zero inputs are replaced with a small positive value; births also have a floor on the magnitude of their negative lambda. These numerical conventions require separate evaluation before any adoption in Hammel.

In `src/src/events.cpp`, `modify_rates()` applies individual and other multipliers to these values. `datev()` (lines 1570–1678) draws a uniform variate, takes its logarithm, divides by the modified negative lambda, and truncates the positive waiting time to integer months, walking across age blocks as needed. Thus this path schedules a future event rather than simply testing a Bernoulli draw each month.

For a constant block without modifiers, boundaries, or competing events, the implied probability of a candidate event in the first month is `1 - exp(lambda)`. Consequently:

- Fertility input `r` implies `1 - exp(-r)`.
- Ordinary mortality input `p` recovers `p`, since `1 - exp(log(1 - p)) = p`.

These statements concern candidate-event timing. Eligibility, birth spacing, rate adjustments, competing events, and rescheduling can change realized demographic outcomes; they remain to be audited.

## Separate Lee–Carter mortality path

`lc_datev()` in `src/src/events.cpp` (around lines 1747–1752) calculates `exp(ax + bx * kt + ext)`, then uses `log(1 - temp_rate) / 12` (with a numerical floor). Algebraically, away from that floor, this corresponds to a monthly probability `1 - (1 - temp_rate)^(1/12)`. This differs from treating the calculated value as a hazard and dividing it by twelve. The intended interpretation of the Lee–Carter inputs needs further investigation.

## Implication for Hammel

SOCSIM does not apply one universal rate-to-probability conversion. Hammel's accepted person-year rate interface should distinguish exposure-unit conversion (`r / 12`) from event-generation semantics. Under a constant hazard interpretation, the monthly event probability is `1 - exp(-r / 12)`. Central mortality rates need a specified relationship to hazards; annual probabilities require different treatment. This source finding does not settle Hammel's algorithm.
