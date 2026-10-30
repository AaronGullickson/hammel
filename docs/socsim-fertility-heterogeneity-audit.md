# SOCSIM fertility heterogeneity and inheritance audit

October 29, 2026. Reference: MPIDR/rsocsim 1.9.18 at commit
`125450cff968bb3ff762fa99d4c926a8136062a2`, inspected in the clean local
`../rsocsim` checkout. Root license: GPLv3. No reference engine source was copied
or adapted into Hammel. Static findings below concern the ordinary non-ENHANCED
engine; four controlled executable checks supplement the trace. This does not
select a Hammel heterogeneity model.

## Persistent individual hazard multipliers

`src/src/utils.c::usefertmult()` (lines 165–170) returns the schedule's lambda
unchanged when `hetfert` is zero, otherwise multiplied by the woman's stored
`fmult`. `src/src/events.cpp::modify_rates()` (lines 1544–1555) applies this
alongside group birth-adjustment factors. Schedules have already undergone the
reference birth-spacing correction; the actual minimum birth interval remains
an additional eligibility restriction.

The multiplier is persistent, rather than redrawn at each event or age change.
For the same eligible schedule, a multiplier of two doubles the woman's hazard.
It does not redistribute a fixed cell intensity among women and does not
normalize observed age/parity rates. Inspecting this code establishes mechanism,
not demographic recovery under arbitrary schedules or composition.

## Distribution of fresh multipliers

`src/src/random.c::fertmult()` (lines 97–103) takes a uniform draw U and returns

    F = 2.178 U - 1.764 U^2 + 1.995 U^3.

This is the guidebook's cubic approximation to a beta distribution, not a direct
beta variate with configurable shape parameters. Under an ideal continuous
uniform draw, analytical polynomial moments give:

| Quantity | Value |
| --- | --- |
| Mean | 0.99975 |
| Variance | 0.4166301375 |
| Range endpoints of the cubic on [0,1] | 0 and 2.409 |

These agree with the manual's rounded mean 1, variance .416 and range 0–2.4.
The implementation uses its package-owned finite-state RNG rather than an ideal
continuous uniform, so these are the polynomial's idealized moments, not an
empirical certification of its RNG. Alpha and betaT do not parameterize this
fresh-draw distribution.

## Maternal inheritance

`src/src/xevents.c::birth()` (lines 114–121) assigns a daughter

    x = alpha * F_mother + (1 - alpha) * gamma
    F_daughter = 2.5 * exp(betaT * log(x / 2.5)),

where gamma is a fresh `fertmult()` draw. The father is not used. Boys receive
zero. For positive x, the second equation is equivalently

    F_daughter = 2.5 * (x / 2.5)^betaT.

| Parameters | Consequence |
| --- | --- |
| alpha=0, betaT=1 | Daughter receives a fresh independent multiplier |
| alpha=1, betaT=1 | Daughter copies her mother's multiplier |
| 0<alpha<1, betaT=1 | Arithmetic blend of maternal multiplier and fresh draw |
| betaT=0, positive x | Every daughter receives 2.5, not one |

For 0<x<2.5, betaT between zero and one raises values toward 2.5, and betaT
above one lowers them toward zero. On the log scale relative to 2.5, betaT is a
scale factor. This changes both level and dispersion; there is no mean-one
renormalization afterward. It applies to newborn daughters even when alpha=0,
but it does not transform multipliers loaded/generated for initial women.

The arithmetic blend also changes dispersion: if mother and innovation are
independent with variance v, Var(x) = [alpha^2 + (1-alpha)^2] v. For alpha=.5,
that is v/2. Maternal selection and unequal input variances complicate this
identity in an actual population. Alpha is not generally the observed
mother–daughter correlation coefficient.

## Defaults, keywords and validation

`src/src/events.cpp` (lines 218–220) initializes:

    hetfert = 1
    alpha = 0
    betaT = 1

The default is heterogeneity without maternal inheritance. The user's first
quoted guidebook table gave alpha=1, beta=0; that does not match this pinned
implementation. The local detailed guide describes the same two-stage formula
as the user's later excerpt.

`src/src/load.cpp` recognizes the literal keyword `betaT`, not `beta`. Alpha is
checked against [0,1] at lines 1039–1048; betaT is parsed without a corresponding
range check at lines 1050–1054. `hetfert` is read as an integer and used as a
truth value, rather than explicitly restricted to zero or one. These findings
should not be generalized to all SOCSIM versions or enhanced builds.

## Initialization, disabling and restart

`src/src/io.c::read_pop()` (lines 275–309):

- A 14-column initial female record with nonzero fmult retains it.
- A female fmult of zero requests a fresh draw; it cannot encode permanent
  infertility through a fixed zero multiplier.
- A 13-column record without fmult generates a fresh draw for each woman.
- These assignments are not conditional on hetfert. Births also generate stored
  daughter multipliers even when their effect on hazards is disabled.

`write_person()` (lines 944–959) writes the multiplier to the population file.
It is read again on restart rather than systematically redrawn, subject to
output rounding and the special treatment of zero. Disabling hetfert bypasses
its hazard effect, not storage, inheritance or all associated RNG consumption.

## Executable switch check and loader finding

An unmodified scratch copy of the pinned package was installed under
`/tmp/hammel-reference-library`. Each check used a fresh R process, seed 9901,
2,000 founder women aged about 20, ten years, fertility input .1/12 per month
on ages 15–50, no birth spacing, negligible reference zero mortality, no unions,
and births retained without fathers. Descendants could not reproduce during
follow-up. Initial multipliers were fixed at either .5 or 2. These are
mechanism checks at one seed, not replicated rate-recovery validation.

| hetfert | Founder multiplier | Founder births | Births / founder woman-year |
| --- | --- | --- | --- |
| 0 | .5 | 2069 | .10345 |
| 0 | 2 | 2069 | .10345 |
| 1 | .5 | 1025 | .05125 |
| 1 | 2 | 4103 | .20515 |

The switch and proportional hazard effect operate as the source describes.
Scratch runner: `/tmp/hammel-hetfert-audit.R`; controlled output directories:
`/tmp/hammel-hetfert-safe-{0,1}-{0.5,2}`. These are local scratch artifacts and
are not portable project dependencies.

**Important loader quirk:** at lines 337–346, `load.cpp::load()` consumes the
first line while positioning the file and starts processing at the next line.
With `hetfert 0` on the first line, the directive was skipped: preliminary
checks returned 1025 births at multiplier .5 and 4103 at multiplier 2 for both
switch settings. Placing a harmless `*` comment first produced the controlled
table above. This also explains the earlier variable-multiplier spacing pilot's
apparently ineffective switch. The retained parity/spacing reference runs fixed
all female multipliers explicitly to one, so that specific omission did not
introduce heterogeneity into their birth hazards. Do not report the switch as
broken independently of this loader issue.

## Consequences for Hammel design

A near-mean-one fresh multiplier does not guarantee mean-one weights among
women at a given parity or among mothers. Fertility selects higher-multiplier
women into births and higher parities; inheritance can transmit that selection
across generations. SOCSIM's mean-one draw and arithmetic inheritance do not
by themselves preserve conditional observed rates or a common marginal
weight distribution across generations.

Hammel's accepted cell-intensity semantics differ: future individual weights
redistribute a cell's fixed m*N birth intensity among eligible women. The
[distribution/inheritance notes](fertility-heterogeneity-notes.md) separate
possible dispersion and inheritance controls from this audited legacy method.
No distribution, parameter defaults or inheritance implementation is accepted
by this audit.

## Primary references

- [Pinned source tree](https://github.com/MPIDR/rsocsim/tree/125450cff968bb3ff762fa99d4c926a8136062a2/src/src).
- [Fresh multiplier generator](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/random.c#L97).
- [Hazard multiplier switch](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/utils.c#L165).
- [Newborn inheritance](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/xevents.c#L114).
- [Pinned guidebook, multiplier and inheritance sections](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/socsim_oversimplified.md#fertility-multiplier-fmult).
