# How event timing works

Updated November 5, 2026. This explains how the engine in `src/simulate.cpp`
times events since milestone 8 (D128): one clock per rate process, spent at the
process's total intensity, with cells choosing who experiences the event. The
[cell scheduling specification](binned-group-hazards.md#cell-scheduling-d128)
gives the data structures; this note gives the reasoning.

## Draw a risk threshold

In R, `rexp(1)` draws one exponential random value with rate 1. Call it `E`. It has the same distribution as `-log(U)`, where `U` is uniform between zero and one:

\[
E = -\log(U), \qquad U \sim \operatorname{Uniform}(0,1).
\]

These are equivalent distributions, not a promise that `rexp()` and `runif()` produce matching values under the same seed. The C++ engine uses R's exponential generator to draw this threshold.

**The threshold is accumulated risk, not a waiting time in years.** Spending `duration` years at an intensity of `Lambda` events per year accumulates `Lambda * duration` units of risk. This product, called cumulative hazard, is dimensionless. The event occurs when accumulated risk reaches the threshold. If the intensity changes, the elapsed interval is charged at the old intensity and the remainder is spent at the new one; no new threshold is needed.

## Why this gives the right distribution

Let `H(t)` be the accumulated intensity by elapsed time `t`. No event has occurred by `t` when the threshold exceeds it. An exponential threshold satisfies `P(E > z) = exp(-z)`, so

\[
P(T > t) = P(E > H(t)) = \exp(-H(t)),
\]

the waiting-time distribution of a Poisson process with that piecewise-constant intensity.

## One clock per process

A process's cells are independent Poisson streams: cell `c` produces events at
`m_c * N_c`, its input rate times its members. Their sum is one Poisson stream
at the total intensity

\[
\Lambda = \sum_c m_c N_c ,
\]

in which each event belongs to cell `c` with probability `m_c N_c / Lambda`.
So each process keeps one threshold, spent at `Lambda`. When it is reached, one
uniform draw picks the cell in proportion to its intensity and another picks
the member, then a new threshold is drawn. Within a cell every member has
hazard `m_c` (with weights, `m_c N_c w_i / sum(w)`), the same as giving each
person their own clock.

Whenever `Lambda` changes (someone dies, is born, moves to another age band or
changes union status, or a rate changes at a calendar boundary), the clock
charges the elapsed time at the old total and re-times the rest at the new one.
Because exponential waiting times are memoryless, this is exact.

## Worked example

A mortality process has two cells: rate 0.01 with 300 members and rate 0.05
with 100 members, so `Lambda = 3 + 5 = 8` deaths per year. Suppose the
threshold is `E = 0.4`.

| Interval | Total intensity | Time | Accumulated risk |
| --- | ---: | ---: | ---: |
| Until a member of the second cell leaves | 8 | 0.02 years | 0.16 |
| Afterwards | 7.95 | until the event | remaining 0.24 |

The remaining wait is `0.24 / 7.95 = 0.0302` years, so the next death occurs
0.0502 years from the start. A uniform draw then picks the second cell with
probability `4.95 / 7.95`, and another picks one of its 99 members.

## Eligibility, boundaries, and restarts

Members move between cells at band edges through bookkeeping boundaries, which
use no random draws. In fertility, women inside the birth interval remain cell
members and contribute intensity, but births go to eligible women (D101); a
cell with no eligible woman counts an unrealized birth. Zero intensity means no
event. Restarts recover spacing from histories and draw new thresholds rather
than restore a checkpoint. See the [scientific contract](scientific-contract.md)
for event competition and RNG isolation.

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

This is the threshold subtraction expressed with the opposite sign. SOCSIM draws one threshold per person and event type and spends it across that person's age bands, as Hammel's engine did before milestone 8. Hammel now spends one threshold per process at the process's total intensity (D128), and SOCSIM's `datev()` returns integer-month dates where Hammel keeps continuous times. Both give the same piecewise-constant hazard model in distribution; this does not assert identical boundary behavior, random-number consumption, or realized events.

There is also an important input distinction. For ordinary mortality rate files, SOCSIM transforms a monthly death probability `p` to `lambda = log(1 - p)`. Its equivalent positive hazard per month is therefore `-log(1 - p)`. Hammel takes supplied annual `mx` directly as a hazard and uses durations in years. Treating an empirical central mortality rate as a hazard is Hammel's documented modeling assumption; the two input conventions must not be confused.

Reference: [rsocsim source at commit 125450cff968bb3ff762fa99d4c926a8136062a2](https://github.com/MPIDR/rsocsim/blob/125450cff968bb3ff762fa99d4c926a8136062a2/src/src/events.cpp), functions `datev()` and `date_and_event()`. See also the [rate conversion audit](socsim-rate-audit.md) for input transformations and the separate Lee–Carter path. These findings are based on static source inspection, not executable equivalence testing.
