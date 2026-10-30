# SOCSIM parity-specific birth-spacing correction

Static inspection, October 29, 2026. Reference: local `../rsocsim` at commit
`125450cff968bb3ff762fa99d4c926a8136062a2`, the same pinned MPIDR/rsocsim source
as the earlier audits. Its root LICENSE.md contains GPLv3. No source was copied
or adapted. The initial audit was static; an executable follow-up is recorded
in [the reference experiment](../validation/milestone5/socsim/README.md). This note covers parity and
spacing only, not fertility heterogeneity or all variants.

## Verified source behavior

- `src/src/utils.c`, get_parity(), lines 1416 onward: female parity is counted
  from the maternal linked child history, not from surviving children only.
- `src/src/events.cpp`, lines 1853–1885: choose the woman's current-parity rate
  schedule. Earliest birth scheduling is current_month for no previous child,
  otherwise max(current_month, last_birth_month + bint). Time since birth is an
  eligibility restriction; the selected fertility schedule remains age-specific.
- `src/src/load.cpp`, adjust_birth_for_bint(), lines 2659–2819: when schedules
  distinguish parities, do not adjust parity-zero rates. For higher parity p,
  sum the previous parity's monthly event probabilities over the preceding
  birth interval and divide the current-parity log-rate parameter by one minus
  that sum. Work from high parity downward so previous inputs are unadjusted.
  This does not estimate parity-population ratios or actual eligible fractions.
- The same routine (lines 2674–2720) creates an additional copy of the final
  supplied parity schedule. Its purpose is to distinguish first entry to the
  pooled highest category (correction from the preceding parity) from later
  births within that category (correction from its own supplied rate).
- If only parity-zero rates were supplied as an unstratified schedule, they
  still receive spacing correction; that case does not mean first births only.

## Implication for Hammel

SOCSIM uses a previous-parity probability-sum approximation, not an exact
population-exposure correction. An exact current-parity eligible fraction would
also depend on parity composition and earlier transitions. These are different
methods; the source trace does not establish that the SOCSIM approximation
recovers arbitrary supplied conditional rates without bias.

A proposed Hammel analogue would retain annual input rates, continuous event
timing and existing finite-band midpoint/upper-open interior-limit evaluation,
using the preceding actual parity's unadjusted age schedule in the correction.
The internal actual-parity count can distinguish entry to the highest open band
from repeated births within it without changing the public input format. The user selected this approximation under D099; the Hammel implementation
is independent code. Its synthetic recovery check reports higher-parity bias
in the initialized-cohort case. This static audit itself remains separate from
executable SOCSIM comparison. D101 subsequently replaced this correction with
cell accounting; milestone 5 acceptance is recorded in D105 and does not
validate the superseded correction.


## Executable follow-up

The user subsequently requested toy SOCSIM runs to test whether the shortfall
also occurs in the reference engine. See [reference results and reproduction](../validation/milestone5/socsim/README.md).
`load.cpp::add_rate_block()` stores fertility input as negative monthly hazard
(`lambda = -input`), whereas several other event types use `log(1 - input)`.
Use annual fertility hazard / 12 for this experiment; do not apply the
probability conversion used for mortality. Source inspection alone did not
establish recovery; the new reference evidence is separate from that audit.
