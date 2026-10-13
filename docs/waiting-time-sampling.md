# How waiting-time sampling works

The engine draws a random amount of accumulated risk needed before an event occurs. Time consumes that amount faster when the event hazard is high. The explanation below uses mortality; the same sampling principle applies to the fertility hazard while a woman is eligible.

## Draw a risk threshold

In R, `rexp(1)` draws one exponential random value with rate 1. Call it `E`. It has the same distribution as `-log(U)`, where `U` is uniform between zero and one:

\[
E = -\log(U), \qquad U \sim \operatorname{Uniform}(0,1).
\]

These are equivalent distributions, not a promise that `rexp()` and `runif()` produce matching values under the same seed. The C++ engine uses R's exponential generator to draw this threshold.

**The threshold is accumulated risk, not a waiting time in years.** If the mortality hazard is `mx` per person-year, spending `duration` years at that hazard accumulates `mx * duration` units of risk. This product, called cumulative hazard, is dimensionless.

If the age band contributes less risk than the remaining threshold, subtract `mx * duration` and carry the remainder into the next band. In the band where the threshold is reached, divide the remaining threshold by that band's hazard to obtain the remaining waiting time in years. A zero-hazard band consumes time without reducing the threshold.

## Worked example

Suppose the threshold is `E = 0.05`. There are two months left in the current age band, with hazard 0.12 per year. The next age band has hazard 0.24 per year.

| Age band | Hazard per year | Time available | Accumulated risk |
| --- | ---: | ---: | ---: |
| Current band | 0.12 | Two months, or `2/12` years | 0.02 |
| Next band | 0.24 | Until the event | Remaining threshold: 0.03 |

The first band does not reach the threshold:

\[
0.05 - 0.12(2/12) = 0.03.
\]

In the next band, the remaining waiting time is

\[
0.03 / 0.24 = 0.125\text{ years} = 1.5\text{ months}.
\]

The event therefore occurs **3.5 months from now**: two months in the first band and 1.5 months in the next. Under Hammel's continuous-time execution, that fractional time is preserved. There is no need to draw a new threshold at the age-band boundary.

## Why this gives the survival distribution

Let `H(t)` be the accumulated hazard by elapsed time `t`. A person survives until `t` when the sampled threshold exceeds that accumulated hazard. An exponential threshold satisfies `P(E > z) = exp(-z)`, so

\[
P(T > t) = P(E > H(t)) = \exp(-H(t)).
\]

For a constant hazard `mx`, `H(t) = mx * t`, giving survival `exp(-mx * t)`. For successive age bands, accumulated hazard is the sum of each band's hazard times the time spent in it. Inverting that sum yields the waiting time under the specified piecewise-constant hazard model.

## Connection to SOCSIM

The ordinary `datev()` path in the audited SOCSIM source uses the same underlying log-survival calculation, expressed with negative numbers. It starts with `log(U)` and stores a negative quantity called `lambda`. If the corresponding positive hazard is `mx`, then

\[
\lambda = -m_x, \qquad \log(U) = -E.
\]

SOCSIM's division therefore yields the same positive within-band wait:

\[
\frac{\log(U)}{\lambda} = \frac{-E}{-m_x} = \frac{E}{m_x}.
\]

When crossing an age band, its update is `logu = logu - duration * lambda`. Substitution gives

\[
-E - (-m_x d) = -(E - m_x d).
\]

This is Hammel's threshold subtraction expressed with the opposite sign. SOCSIM's `datev()` converts its within-band wait to an integer and returns an integer-month date; Hammel now preserves continuous event times. This mathematical connection does not assert identical boundary behavior, random-number consumption, or realized events.

There is also an important input distinction. For ordinary mortality rate files, SOCSIM transforms a monthly death probability `p` to `lambda = log(1 - p)`. Its equivalent positive hazard per month is therefore `-log(1 - p)`. Hammel takes supplied annual `mx` directly as a hazard and uses durations in years. Treating an empirical central mortality rate as a hazard is Hammel's documented modeling assumption; the two input conventions must not be confused.

Reference: [rsocsim source at commit 125450cff968bb3ff762fa99d4c926a8136062a2](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp), functions `datev()` and `date_and_event()`. See also the [rate conversion audit](socsim-rate-audit.md) for input transformations and the separate Lee–Carter path. These findings are based on static source inspection, not executable equivalence testing.
